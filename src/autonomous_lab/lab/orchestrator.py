"""
Autonomous Lab - Multi-Agent Orchestrator

Coordinates multiple LLM agents (each with its own context window) to
execute PI, Trainee, and Reviewer roles. Opt-in via `orchestration: multi`
in .autolab/config.yaml. Falls through to single-agent mode on any failure.

CLI providers (claude-cli, codex-cli, cursor-cli) are the primary path:
zero API key configuration, leverages existing AI subscriptions.
API providers (anthropic, openai) are the fallback for per-token billing.
"""

import logging
import os
import shutil
from dataclasses import dataclass

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


class Orchestrator:
    """Multi-agent orchestrator for Autonomous Lab.

    Loads agent configuration from .autolab/config.yaml, creates
    providers, and runs scoped turns for each role.
    """

    def __init__(self, project_dir: str):
        self.project_dir = os.path.abspath(project_dir)
        self.agents: dict[str, AgentConfig] = {}
        self._load_agents()

    def _load_agents(self) -> None:
        """Load agent configs from config.yaml, resolve API keys."""
        config = load_config(self.project_dir)
        agents_cfg = config.get("agents", self._default_agents())
        api_keys_env = config.get("api_keys_env", {})

        for role, acfg in agents_cfg.items():
            provider_name = acfg.get("provider", "claude-cli")
            model = acfg.get("model", "sonnet")
            api_key = ""

            if provider_name in ("anthropic", "openai"):
                env_var = api_keys_env.get(
                    provider_name, f"{provider_name.upper()}_API_KEY"
                )
                api_key = os.environ.get(env_var, "")

            self.agents[role] = AgentConfig(
                role=role,
                provider_name=provider_name,
                model=model,
                api_key=api_key,
            )

    @staticmethod
    def _default_agents() -> dict:
        """Default: both roles use claude-cli with sonnet."""
        return {
            "pi": {"provider": "claude-cli", "model": "sonnet"},
            "trainee": {"provider": "claude-cli", "model": "sonnet"},
        }

    def is_available(self, role: str) -> bool:
        """Check if the agent for a role is usable.

        CLI providers: checks if the binary is on PATH.
        API providers: checks if the API key env var is set.
        """
        agent = self.agents.get(role)
        if not agent:
            return False

        if agent.provider_name in CLI_PROVIDER_MAP:
            cli_key = CLI_PROVIDER_MAP[agent.provider_name]
            return is_cli_available(cli_key)

        # API providers need an API key
        return bool(agent.api_key)

    async def run_turn(self, role: str) -> str:
        """Execute a full turn for the given role.

        For CLI providers: builds a scoped prompt, spawns the CLI,
        captures output, records the turn, and advances state.

        For API providers: runs an agentic tool loop with the
        ToolExecutor sandbox.

        Args:
            role: 'pi' or 'trainee'.

        Returns:
            Summary string of what the agent accomplished.
        """
        agent = self.agents[role]
        provider = create_provider(
            agent.provider_name,
            agent.model,
            self.project_dir,
            agent.api_key,
        )

        prompt = self._build_scoped_prompt(role)

        if isinstance(provider, CLIProvider):
            output = await provider.run_agent(prompt)
            summary = self._extract_summary(output, role)
            self._record_turn(role, summary, output)
            log_event(
                self.project_dir,
                "multi_agent_turn",
                role=role,
                provider=agent.provider_name,
                model=agent.model,
                output_length=len(output),
            )
            return summary
        else:
            return await self._run_api_turn(role, provider, prompt)

    async def _run_api_turn(self, role: str, provider, prompt: str) -> str:
        """Agentic loop for API providers with tool use."""
        from .agent_tools import ToolExecutor

        agent = self.agents[role]
        config = load_config(self.project_dir)
        max_iterations = config.get("orchestrator", {}).get(
            "max_tool_calls_per_turn", 20
        )

        tool_executor = ToolExecutor(self.project_dir)
        tools = tool_executor.get_tools_for_role(role)
        formatted_tools = provider.format_tools(tools)

        system_prompt = (
            "You are an AI research agent working on an Autonomous Lab project. "
            "Follow the instructions in the user message precisely. "
            "Use the available tools to read/write files and run scripts."
        )
        messages = [{"role": "user", "content": prompt}]
        full_text = ""

        for _i in range(max_iterations):
            response: LLMResponse = await provider.complete(
                system_prompt=system_prompt,
                messages=messages,
                tools=formatted_tools,
                max_tokens=8192,
            )
            full_text += response.content or ""

            if not response.tool_calls:
                break

            # Execute each tool call and collect results
            tool_results = []
            for tc in response.tool_calls:
                try:
                    result = await tool_executor.execute(
                        tc.name, tc.arguments, role
                    )
                except Exception as e:
                    result = f"Tool error: {e}"
                    logger.warning("Tool %s failed: %s", tc.name, e)
                tool_results.append({"tool_call_id": tc.id, "content": result})

            # Append assistant message with tool calls and tool results
            # Anthropic format
            if agent.provider_name == "anthropic":
                assistant_content = []
                if response.content:
                    assistant_content.append(
                        {"type": "text", "text": response.content}
                    )
                for tc in response.tool_calls:
                    assistant_content.append(
                        {
                            "type": "tool_use",
                            "id": tc.id,
                            "name": tc.name,
                            "input": tc.arguments,
                        }
                    )
                messages.append({"role": "assistant", "content": assistant_content})
                # Tool results as user message
                tool_content = [
                    {
                        "type": "tool_result",
                        "tool_use_id": tr["tool_call_id"],
                        "content": tr["content"],
                    }
                    for tr in tool_results
                ]
                messages.append({"role": "user", "content": tool_content})

            # OpenAI format
            elif agent.provider_name == "openai":
                import json

                messages.append(
                    {
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
                    }
                )
                for tr in tool_results:
                    messages.append(
                        {
                            "role": "tool",
                            "tool_call_id": tr["tool_call_id"],
                            "content": tr["content"],
                        }
                    )

        summary = self._extract_summary(full_text, role)
        self._record_turn(role, summary, full_text)
        log_event(
            self.project_dir,
            "multi_agent_turn",
            role=role,
            provider=agent.provider_name,
            model=agent.model,
            output_length=len(full_text),
            tool_iterations=_i + 1 if '_i' in dir() else 0,
        )
        return summary

    def _build_scoped_prompt(self, role: str) -> str:
        """Build a context-scoped prompt for a role.

        Each role sees only what it needs:
        - PI: full idea + summaries + recent meetings + file listings + feedback
        - Trainee: full idea + last PI directive + file listings + feedback
        """
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

        # Wrap with execution context for the agent
        dcfg_label = (
            dcfg.get("senior_label", "PI")
            if role == "pi"
            else dcfg.get("junior_label", "Trainee")
        )

        header = (
            f"# Autonomous Lab — {dcfg_label} Turn (Iteration {iteration})\n\n"
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

    def _extract_summary(self, output: str, role: str) -> str:
        """Extract a brief summary from agent output.

        Looks for common summary markers, falls back to first lines.
        """
        if not output:
            return f"{role} turn completed (no output captured)"

        # Look for explicit summary markers
        for marker in ["### Summary", "## Summary", "**Summary**", "Summary:"]:
            idx = output.find(marker)
            if idx != -1:
                # Take text after marker up to next heading or 500 chars
                after = output[idx + len(marker) :].strip()
                # Find end of summary section
                for end_marker in ["\n## ", "\n### ", "\n---", "\n**"]:
                    end_idx = after.find(end_marker)
                    if end_idx != -1:
                        after = after[:end_idx]
                        break
                if len(after) > 500:
                    after = after[:500] + "..."
                return after.strip()

        # Fallback: first 300 chars
        lines = output.strip().split("\n")
        preview = "\n".join(lines[:5])
        if len(preview) > 300:
            preview = preview[:300] + "..."
        return preview

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
