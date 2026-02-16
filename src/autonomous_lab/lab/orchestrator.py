"""
Autonomous Lab - Multi-Agent Orchestrator

Coordinates multiple LLM agents (each with its own context window) to
execute PI, Trainee, and Reviewer roles. Opt-in via `orchestration: multi`
in .autolab/config.yaml. Falls through to single-agent mode on any failure.

Supports:
- Single trainee (backward compatible `trainee` key)
- Multiple trainees in parallel (`trainees` key, runs via asyncio.gather)
- Mixed providers (PI on claude-cli, trainees on codex-cli / cursor-cli / API)

CLI providers (claude-cli, codex-cli, cursor-cli) are the primary path:
zero API key configuration, leverages existing AI subscriptions.
API providers (anthropic, openai) are the fallback for per-token billing.
"""

import asyncio
import logging
import os
from dataclasses import dataclass, field

from .event_log import log_event
from .providers.base import LLMResponse
from .providers.cli_provider import CLIProvider, is_cli_available
from .providers.factory import CLI_PROVIDER_MAP, create_provider
from .state import (
    COMPRESS_EVERY,
    append_meeting_log,
    compress_old_meetings,
    drain_feedback_queue,
    get_domain_config,
    get_meeting_summaries,
    get_paper_progress,
    get_recent_meetings,
    load_config,
    load_idea,
    load_profile,
    load_state,
    save_state,
    scan_project_files,
)

logger = logging.getLogger(__name__)


@dataclass
class AgentConfig:
    """Configuration for a single agent role."""

    role: str
    provider_name: str
    model: str
    api_key: str = ""


@dataclass
class TraineeConfig:
    """Configuration for one trainee in a multi-trainee setup."""

    name: str
    provider_name: str
    model: str
    focus: str = ""
    api_key: str = ""


class Orchestrator:
    """Multi-agent orchestrator for Autonomous Lab.

    Loads agent configuration from .autolab/config.yaml, creates
    providers, and runs scoped turns for each role. Supports both
    single-trainee and multi-trainee configurations.
    """

    def __init__(self, project_dir: str):
        self.project_dir = os.path.abspath(project_dir)
        self.pi_config: AgentConfig | None = None
        self.trainees: list[TraineeConfig] = []
        self._load_config()

    def _load_config(self) -> None:
        """Load agent configs from config.yaml, resolve API keys."""
        config = load_config(self.project_dir)
        agents_cfg = config.get("agents", self._default_agents())
        api_keys_env = config.get("api_keys_env", {})
        orch_cfg = config.get("orchestrator", {})
        self.trainee_timeout = orch_cfg.get("trainee_timeout", 600)
        self.max_tool_calls = orch_cfg.get("max_tool_calls_per_turn", 20)
        self.tool_timeout = orch_cfg.get("tool_timeout", 120)

        # PI config
        pi_cfg = agents_cfg.get("pi", {"provider": "claude-cli", "model": "sonnet"})
        self.pi_config = self._make_agent_config("pi", pi_cfg, api_keys_env)

        # Trainee config — supports both singular and plural
        if "trainees" in agents_cfg:
            # Multiple trainees
            for t in agents_cfg["trainees"]:
                self.trainees.append(TraineeConfig(
                    name=t.get("name", "Trainee"),
                    provider_name=t.get("provider", "claude-cli"),
                    model=t.get("model", "sonnet"),
                    focus=t.get("focus", ""),
                    api_key=self._resolve_api_key(
                        t.get("provider", "claude-cli"), api_keys_env
                    ),
                ))
        else:
            # Single trainee (backward compatible)
            t_cfg = agents_cfg.get(
                "trainee", {"provider": "claude-cli", "model": "sonnet"}
            )
            self.trainees.append(TraineeConfig(
                name="Trainee",
                provider_name=t_cfg.get("provider", "claude-cli"),
                model=t_cfg.get("model", "sonnet"),
                focus="",
                api_key=self._resolve_api_key(
                    t_cfg.get("provider", "claude-cli"), api_keys_env
                ),
            ))

    def _make_agent_config(
        self, role: str, cfg: dict, api_keys_env: dict
    ) -> AgentConfig:
        """Create an AgentConfig from a config dict."""
        provider_name = cfg.get("provider", "claude-cli")
        return AgentConfig(
            role=role,
            provider_name=provider_name,
            model=cfg.get("model", "sonnet"),
            api_key=self._resolve_api_key(provider_name, api_keys_env),
        )

    @staticmethod
    def _resolve_api_key(provider_name: str, api_keys_env: dict) -> str:
        """Resolve API key from environment for API providers."""
        if provider_name in ("anthropic", "openai"):
            env_var = api_keys_env.get(
                provider_name, f"{provider_name.upper()}_API_KEY"
            )
            return os.environ.get(env_var, "")
        return ""

    @staticmethod
    def _default_agents() -> dict:
        """Default: both roles use claude-cli with sonnet."""
        return {
            "pi": {"provider": "claude-cli", "model": "sonnet"},
            "trainee": {"provider": "claude-cli", "model": "sonnet"},
        }

    # ------------------------------------------------------------------
    # Availability checks
    # ------------------------------------------------------------------

    def is_available(self, role: str) -> bool:
        """Check if the agent(s) for a role are usable.

        For PI: checks if the PI provider is available.
        For trainee: checks if at least one trainee provider is available.
        """
        if role == "pi":
            return self._provider_available(self.pi_config)
        elif role == "trainee":
            return any(
                self._trainee_provider_available(t) for t in self.trainees
            )
        return False

    def _provider_available(self, agent: AgentConfig) -> bool:
        """Check if a single agent's provider is usable."""
        if agent.provider_name in CLI_PROVIDER_MAP:
            cli_key = CLI_PROVIDER_MAP[agent.provider_name]
            return is_cli_available(cli_key)
        return bool(agent.api_key)

    def _trainee_provider_available(self, trainee: TraineeConfig) -> bool:
        """Check if a trainee's provider is usable."""
        if trainee.provider_name in CLI_PROVIDER_MAP:
            cli_key = CLI_PROVIDER_MAP[trainee.provider_name]
            return is_cli_available(cli_key)
        return bool(trainee.api_key)

    # ------------------------------------------------------------------
    # Turn execution
    # ------------------------------------------------------------------

    async def run_turn(self, role: str) -> str:
        """Execute a full turn for the given role (backward-compatible API).

        Routes to run_pi_turn() or run_trainee_turn() based on role.
        """
        if role == "pi":
            return await self.run_pi_turn()
        else:
            return await self.run_trainee_turn()

    async def run_pi_turn(self) -> str:
        """Run PI turn: spawns one agent session, returns summary."""
        agent = self.pi_config
        provider = create_provider(
            agent.provider_name, agent.model,
            self.project_dir, agent.api_key,
        )

        prompt = self._build_scoped_prompt("pi")

        if isinstance(provider, CLIProvider):
            output = await provider.run_agent(prompt)
            summary = self._extract_summary(output, "pi")
            self._record_turn("pi", summary, output)
            log_event(
                self.project_dir, "multi_agent_turn",
                role="pi", provider=agent.provider_name,
                model=agent.model, output_length=len(output),
            )
            return summary
        else:
            return await self._run_api_turn("pi", provider, prompt)

    async def run_trainee_turn(self) -> str:
        """Run trainee turn(s). Multiple trainees run in PARALLEL."""
        if len(self.trainees) == 1:
            summary, content = await self._run_single_trainee(self.trainees[0])
            self._record_turn("trainee", summary, content)
            return summary

        # Multiple trainees — run in parallel via asyncio.gather
        tasks = []
        for trainee in self.trainees:
            tasks.append(self._run_single_trainee(trainee))

        results = await asyncio.gather(*tasks, return_exceptions=True)

        # Combine results
        combined_summary = []
        combined_content = []
        for trainee, result in zip(self.trainees, results):
            if isinstance(result, BaseException):
                combined_summary.append(f"{trainee.name}: FAILED ({result})")
                log_event(
                    self.project_dir, "trainee_error",
                    trainee=trainee.name, error=str(result),
                )
            else:
                summary, content = result
                combined_summary.append(f"{trainee.name}: {summary}")
                combined_content.append(f"## {trainee.name}\n\n{content}")

        full_summary = " | ".join(combined_summary)
        full_content = "\n\n---\n\n".join(combined_content)
        self._record_turn("trainee", full_summary, full_content)
        return full_summary

    async def _run_single_trainee(
        self, trainee: TraineeConfig
    ) -> tuple[str, str]:
        """Run one trainee, return (summary, full_content)."""
        provider = create_provider(
            trainee.provider_name, trainee.model,
            self.project_dir, trainee.api_key,
            timeout=self.trainee_timeout,
        )

        prompt = self._build_trainee_prompt(trainee)

        if isinstance(provider, CLIProvider):
            output = await provider.run_agent(prompt)
        else:
            output = await self._run_api_turn_raw(
                "trainee", provider, prompt
            )

        summary = self._extract_summary(output, "trainee")
        log_event(
            self.project_dir, "multi_agent_turn",
            role="trainee", trainee_name=trainee.name,
            provider=trainee.provider_name, model=trainee.model,
            output_length=len(output),
        )
        return summary, output

    # ------------------------------------------------------------------
    # API provider agentic loop
    # ------------------------------------------------------------------

    async def _run_api_turn(self, role: str, provider, prompt: str) -> str:
        """Agentic loop for API providers — records turn and returns summary."""
        output = await self._run_api_turn_raw(role, provider, prompt)
        summary = self._extract_summary(output, role)
        self._record_turn(role, summary, output)
        log_event(
            self.project_dir, "multi_agent_turn",
            role=role, provider=provider.provider_name,
            model=provider.model, output_length=len(output),
        )
        return summary

    async def _run_api_turn_raw(
        self, role: str, provider, prompt: str
    ) -> str:
        """Agentic loop for API providers — returns raw output text."""
        from .agent_tools import ToolExecutor

        tool_executor = ToolExecutor(
            self.project_dir, shell_timeout=self.tool_timeout
        )
        tools = tool_executor.get_tools_for_role(role)
        formatted_tools = provider.format_tools(tools)

        system_prompt = (
            "You are an AI research agent working on an Autonomous Lab project. "
            "Follow the instructions in the user message precisely. "
            "Use the available tools to read/write files and run scripts."
        )
        messages = [{"role": "user", "content": prompt}]
        full_text = ""

        for _i in range(self.max_tool_calls):
            response: LLMResponse = await provider.complete(
                system_prompt=system_prompt,
                messages=messages,
                tools=formatted_tools,
                max_tokens=8192,
            )
            full_text += response.content or ""

            if not response.tool_calls:
                break

            # Execute each tool call
            tool_results = []
            for tc in response.tool_calls:
                try:
                    result = await tool_executor.execute(
                        tc.name, tc.arguments, role
                    )
                except Exception as e:
                    result = f"Tool error: {e}"
                    logger.warning("Tool %s failed: %s", tc.name, e)
                tool_results.append(
                    {"tool_call_id": tc.id, "content": result}
                )

            # Append messages in provider-specific format
            if provider.provider_name == "anthropic":
                self._append_anthropic_messages(
                    messages, response, tool_results
                )
            elif provider.provider_name == "openai":
                self._append_openai_messages(
                    messages, response, tool_results
                )

        return full_text

    @staticmethod
    def _append_anthropic_messages(
        messages: list[dict], response: LLMResponse, tool_results: list[dict]
    ) -> None:
        """Append Anthropic-format tool use/result messages."""
        assistant_content = []
        if response.content:
            assistant_content.append(
                {"type": "text", "text": response.content}
            )
        for tc in response.tool_calls:
            assistant_content.append({
                "type": "tool_use",
                "id": tc.id,
                "name": tc.name,
                "input": tc.arguments,
            })
        messages.append({"role": "assistant", "content": assistant_content})
        tool_content = [
            {
                "type": "tool_result",
                "tool_use_id": tr["tool_call_id"],
                "content": tr["content"],
            }
            for tr in tool_results
        ]
        messages.append({"role": "user", "content": tool_content})

    @staticmethod
    def _append_openai_messages(
        messages: list[dict], response: LLMResponse, tool_results: list[dict]
    ) -> None:
        """Append OpenAI-format tool call/result messages."""
        import json

        messages.append({
            "role": "assistant",
            "content": response.content or None,
            "tool_calls": [
                {
                    "id": tc.id,
                    "type": "function",
                    "function": {
                        "name": tc.name,
                        "arguments": json.dumps(tc.arguments),
                    },
                }
                for tc in response.tool_calls
            ],
        })
        for tr in tool_results:
            messages.append({
                "role": "tool",
                "tool_call_id": tr["tool_call_id"],
                "content": tr["content"],
            })

    # ------------------------------------------------------------------
    # Prompt building
    # ------------------------------------------------------------------

    def _build_scoped_prompt(self, role: str) -> str:
        """Build a context-scoped prompt for PI or single-trainee role."""
        from .prompts import build_pi_prompt, build_trainee_prompt

        state = load_state(self.project_dir)
        idea = load_idea(self.project_dir)
        iteration = state["iteration"]
        dcfg = get_domain_config(self.project_dir)
        user_feedback = drain_feedback_queue(self.project_dir)
        meeting_history = get_recent_meetings(self.project_dir, n=3)
        summaries = get_meeting_summaries(self.project_dir)
        file_listings = scan_project_files(self.project_dir)
        profile = load_profile(self.project_dir, role)

        if role == "pi":
            prompt = build_pi_prompt(
                idea=idea,
                profile=profile,
                meeting_history=meeting_history,
                summaries=summaries,
                file_listings=file_listings,
                user_feedback=user_feedback,
                iteration=iteration,
                domain_config=dcfg,
            )
        else:
            prompt = build_trainee_prompt(
                idea=idea,
                profile=profile,
                meeting_history=meeting_history,
                summaries=summaries,
                file_listings=file_listings,
                user_feedback=user_feedback,
                iteration=iteration,
                domain_config=dcfg,
            )

        dcfg_label = (
            dcfg.get("senior_label", "PI")
            if role == "pi"
            else dcfg.get("junior_label", "Trainee")
        )

        header = (
            f"# Autonomous Lab -- {dcfg_label} Turn (Iteration {iteration})\n\n"
            f"You are the {dcfg_label} in an autonomous research lab.\n"
            f"Project directory: {self.project_dir}\n\n"
            f"Execute your role actions below. Write code to scripts/, "
            f"figures to figures/, results to results/, "
            f"and paper sections to paper/sections/.\n\n"
            f"---\n\n"
        )

        footer = (
            "\n\n---\n\n"
            "IMPORTANT: Complete all actions described above. "
            "Your output will be recorded as this turn's contribution."
        )

        return header + prompt + footer

    def _build_trainee_prompt(self, trainee: TraineeConfig) -> str:
        """Build a scoped prompt for a specific trainee, including focus area.

        When multiple trainees exist, each gets a tailored prompt with
        their focus area prepended and abbreviated context.
        """
        from .prompts import build_trainee_prompt

        state = load_state(self.project_dir)
        idea = load_idea(self.project_dir)
        iteration = state["iteration"]
        dcfg = get_domain_config(self.project_dir)
        meeting_history = get_recent_meetings(self.project_dir, n=3)
        summaries = get_meeting_summaries(self.project_dir)
        file_listings = scan_project_files(self.project_dir)
        profile = load_profile(self.project_dir, "trainee")

        # If only one trainee with no focus, use the standard prompt
        if len(self.trainees) == 1 and not trainee.focus:
            return self._build_scoped_prompt("trainee")

        # For multi-trainee or focused trainee, build a tailored prompt
        base_prompt = build_trainee_prompt(
            idea=idea,
            profile=profile,
            meeting_history=meeting_history,
            summaries=summaries,
            file_listings=file_listings,
            user_feedback="",  # Feedback already drained for PI
            iteration=iteration,
            domain_config=dcfg,
        )

        focus_note = ""
        if trainee.focus:
            focus_note = (
                f"YOUR FOCUS AREA: {trainee.focus}\n"
                f"Concentrate your work on tasks related to your focus area. "
                f"Other trainees are handling other aspects of the project.\n\n"
            )

        other_trainees = [
            t for t in self.trainees if t.name != trainee.name
        ]
        coordination_note = ""
        if other_trainees:
            others = ", ".join(
                f"{t.name} ({t.focus})" if t.focus else t.name
                for t in other_trainees
            )
            coordination_note = (
                f"NOTE: Other trainees working in parallel: {others}. "
                f"Avoid duplicating their work.\n\n"
            )

        header = (
            f"# Autonomous Lab -- {trainee.name} Turn (Iteration {iteration})\n\n"
            f"You are '{trainee.name}', a trainee in an autonomous research lab.\n"
            f"Project directory: {self.project_dir}\n\n"
            f"{focus_note}"
            f"{coordination_note}"
            f"Execute your role actions below. Write code to scripts/, "
            f"figures to figures/, results to results/, "
            f"and paper sections to paper/sections/.\n\n"
            f"---\n\n"
        )

        footer = (
            "\n\n---\n\n"
            "IMPORTANT: Complete all actions described above. "
            "Your output will be recorded as this turn's contribution."
        )

        return header + base_prompt + footer

    # ------------------------------------------------------------------
    # Summary extraction
    # ------------------------------------------------------------------

    def _extract_summary(self, output: str, role: str) -> str:
        """Extract a brief summary from agent output.

        Looks for common summary markers, falls back to first lines.
        """
        if not output:
            return f"{role} turn completed (no output captured)"

        for marker in ["### Summary", "## Summary", "**Summary**", "Summary:"]:
            idx = output.find(marker)
            if idx != -1:
                after = output[idx + len(marker):].strip()
                for end_marker in ["\n## ", "\n### ", "\n---", "\n**"]:
                    end_idx = after.find(end_marker)
                    if end_idx != -1:
                        after = after[:end_idx]
                        break
                if len(after) > 500:
                    after = after[:500] + "..."
                return after.strip()

        lines = output.strip().split("\n")
        preview = "\n".join(lines[:5])
        if len(preview) > 300:
            preview = preview[:300] + "..."
        return preview

    # ------------------------------------------------------------------
    # Turn recording
    # ------------------------------------------------------------------

    def _record_turn(self, role: str, summary: str, content: str) -> None:
        """Record the completed turn via existing state management."""
        state = load_state(self.project_dir)
        iteration = state["iteration"]

        append_meeting_log(self.project_dir, role, iteration, summary, content)

        # Flip role and advance iteration
        next_role = "trainee" if role == "pi" else "pi"
        new_iteration = iteration + (1 if role == "trainee" else 0)

        state["iteration"] = new_iteration
        state["next_role"] = next_role
        save_state(self.project_dir, state)

        # Compress old meetings periodically
        if new_iteration > 0 and new_iteration % COMPRESS_EVERY == 0:
            try:
                compress_old_meetings(self.project_dir)
            except Exception as e:
                logger.warning("Meeting compression failed: %s", e)
