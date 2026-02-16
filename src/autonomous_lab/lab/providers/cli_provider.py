"""
Autonomous Lab - CLI Agent Provider

Spawns CLI agents (claude, codex, cursor-agent) as subprocesses.
Each CLI handles its own tool calls internally — no ToolExecutor needed.
This is the primary provider path: zero API key configuration,
leverages existing AI subscriptions (Claude Pro, ChatGPT Plus, etc.).
"""

import asyncio
import logging
import os
import shutil
from typing import Any

from .base import LLMProvider, LLMResponse

logger = logging.getLogger(__name__)

# Default subprocess timeout (seconds)
_DEFAULT_TIMEOUT = 600

# CLI binary configs: how to build the command for each CLI
_CLI_CONFIGS: dict[str, dict[str, Any]] = {
    "claude": {
        "binary": "claude",
        "build_cmd": lambda model, prompt, project_dir: [
            "claude",
            "-p",
            prompt,
            "--model",
            model,
            "--output-format",
            "text",
            "--dangerously-skip-permissions",
            "--add-dir",
            project_dir,
        ],
    },
    "codex": {
        "binary": "codex",
        "build_cmd": lambda model, prompt, project_dir: [
            "codex",
            "exec",
            prompt,
            "--model",
            model,
            "-s",
            "workspace-write",
        ],
    },
    "cursor-agent": {
        "binary": "cursor-agent",
        "build_cmd": lambda model, prompt, project_dir: [
            "cursor-agent",
            "-p",
            prompt,
            "--model",
            model,
            "--output-format",
            "text",
            "--workspace",
            project_dir,
            "--force",
            "--trust",
        ],
    },
}


def is_cli_available(cli: str) -> bool:
    """Check if a CLI binary is available on PATH."""
    config = _CLI_CONFIGS.get(cli)
    if not config:
        return False
    return shutil.which(config["binary"]) is not None


class CLIProvider(LLMProvider):
    """Runs a CLI agent (claude, codex, or cursor-agent) as a subprocess.

    The CLI handles tool execution internally — we just pass a prompt and
    capture the output. This makes CLI providers zero-config: users only
    need their existing AI subscription CLI installed.
    """

    def __init__(
        self,
        cli: str,
        model: str,
        project_dir: str,
        timeout: int = _DEFAULT_TIMEOUT,
    ):
        if cli not in _CLI_CONFIGS:
            raise ValueError(
                f"Unknown CLI: {cli!r}. Supported: {list(_CLI_CONFIGS.keys())}"
            )
        self.cli = cli
        self.model = model
        self.project_dir = os.path.abspath(project_dir)
        self.timeout = timeout

    @property
    def provider_name(self) -> str:
        return f"{self.cli}-cli"

    async def complete(
        self,
        system_prompt: str,
        messages: list[dict],
        tools: list[dict] | None = None,
        max_tokens: int = 8192,
    ) -> LLMResponse:
        """Run the CLI agent. Ignores tools param (CLI handles its own).

        Combines system_prompt and messages into a single prompt string,
        then delegates to run_agent().
        """
        # Build a flat prompt from system + messages
        parts = []
        if system_prompt:
            parts.append(system_prompt)
        for msg in messages:
            content = msg.get("content", "")
            if isinstance(content, list):
                # Handle structured content blocks
                content = "\n".join(
                    block.get("text", "") for block in content if isinstance(block, dict)
                )
            if content:
                parts.append(content)

        prompt = "\n\n".join(parts)
        output = await self.run_agent(prompt)

        return LLMResponse(
            content=output,
            tool_calls=[],  # CLI handles tools internally
            stop_reason="end_turn",
            model=self.model,
        )

    async def run_agent(self, prompt: str) -> str:
        """Spawn CLI agent with prompt, return its full output.

        Args:
            prompt: The complete prompt to send to the agent.

        Returns:
            The agent's stdout output as a string.

        Raises:
            TimeoutError: If the agent exceeds the timeout.
            RuntimeError: If the agent process fails to start or exits with error.
        """
        config = _CLI_CONFIGS[self.cli]
        cmd = config["build_cmd"](self.model, prompt, self.project_dir)

        # Build clean environment
        env = {**os.environ}
        # Allow nested Claude sessions by removing the marker
        env.pop("CLAUDECODE", None)
        # Ensure the project dir is the working directory for codex
        cwd = self.project_dir if self.cli == "codex" else None

        logger.info(
            "Spawning %s agent: model=%s, project=%s",
            self.cli,
            self.model,
            self.project_dir,
        )

        try:
            proc = await asyncio.create_subprocess_exec(
                *cmd,
                cwd=cwd or self.project_dir,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                env=env,
            )
        except FileNotFoundError:
            raise RuntimeError(
                f"CLI binary '{config['binary']}' not found. "
                f"Install it or switch to a different provider."
            )

        try:
            stdout, stderr = await asyncio.wait_for(
                proc.communicate(), timeout=self.timeout
            )
        except asyncio.TimeoutError:
            proc.kill()
            await proc.communicate()  # drain pipes
            raise TimeoutError(
                f"{self.cli} agent timed out after {self.timeout}s"
            )

        output = stdout.decode("utf-8", errors="replace").strip()
        err_output = stderr.decode("utf-8", errors="replace").strip()

        if proc.returncode != 0 and not output:
            # Only raise if we got no output — some CLIs write warnings to stderr
            raise RuntimeError(
                f"{self.cli} agent exited with code {proc.returncode}: {err_output[:500]}"
            )

        if err_output:
            logger.debug("%s stderr: %s", self.cli, err_output[:200])

        return output
