"""
Autonomous Lab - Agent Tool Sandbox

Provides tool definitions and a sandboxed executor for API providers.
CLI providers handle their own tools — this module is only used when
the orchestrator manages the agentic tool loop for API-based agents.
"""

import asyncio
import json
import logging
import os
from pathlib import Path

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Role-based tool permissions
# ---------------------------------------------------------------------------
TOOL_PERMISSIONS: dict[str, set[str]] = {
    "pi": {"autolab_cite", "file_read"},
    "trainee": {"autolab_cite", "file_read", "file_write", "shell_execute"},
    "reviewer": {"file_read", "autolab_cite"},
}

# ---------------------------------------------------------------------------
# Canonical tool specifications (JSON Schema)
# ---------------------------------------------------------------------------
TOOL_SPECS: list[dict] = [
    {
        "name": "file_read",
        "description": "Read the contents of a file within the project directory.",
        "parameters": {
            "type": "object",
            "properties": {
                "path": {
                    "type": "string",
                    "description": "Relative path to the file from the project root.",
                },
            },
            "required": ["path"],
        },
    },
    {
        "name": "file_write",
        "description": "Write content to a file within the project directory. Creates parent directories as needed.",
        "parameters": {
            "type": "object",
            "properties": {
                "path": {
                    "type": "string",
                    "description": "Relative path to the file from the project root.",
                },
                "content": {
                    "type": "string",
                    "description": "The content to write to the file.",
                },
            },
            "required": ["path", "content"],
        },
    },
    {
        "name": "shell_execute",
        "description": "Execute a shell command in the project directory. Use for running Python scripts, installing packages, etc.",
        "parameters": {
            "type": "object",
            "properties": {
                "command": {
                    "type": "string",
                    "description": "The shell command to execute.",
                },
                "timeout": {
                    "type": "integer",
                    "description": "Timeout in seconds (default 120).",
                    "default": 120,
                },
            },
            "required": ["command"],
        },
    },
    {
        "name": "autolab_cite",
        "description": "Search for papers, look up a DOI, or validate references.bib. Action: 'search' (find papers), 'doi' (get BibTeX from DOI), 'validate' (check .bib file).",
        "parameters": {
            "type": "object",
            "properties": {
                "action": {
                    "type": "string",
                    "enum": ["search", "doi", "validate"],
                    "description": "Citation action to perform.",
                },
                "query": {
                    "type": "string",
                    "description": "Search query, DOI string, or .bib file path.",
                },
                "count": {
                    "type": "integer",
                    "description": "Number of search results (max 10).",
                    "default": 5,
                },
            },
            "required": ["action", "query"],
        },
    },
]


class ToolExecutor:
    """Sandboxed tool executor for API-based agent turns.

    All file operations are confined to the project directory.
    Shell commands run with the project directory as cwd.
    """

    def __init__(self, project_dir: str, shell_timeout: int = 120):
        self.project_dir = os.path.abspath(project_dir)
        self.shell_timeout = shell_timeout

    def get_tools_for_role(self, role: str) -> list[dict]:
        """Return tool specs filtered by role permissions."""
        allowed = TOOL_PERMISSIONS.get(role, set())
        return [spec for spec in TOOL_SPECS if spec["name"] in allowed]

    async def execute(self, tool_name: str, arguments: dict, role: str) -> str:
        """Execute a tool call and return the result string.

        Args:
            tool_name: Name of the tool to execute.
            arguments: Tool arguments dict.
            role: The role executing this tool (for permission checking).

        Returns:
            Tool output as a string.

        Raises:
            PermissionError: If the role doesn't have access to this tool.
        """
        allowed = TOOL_PERMISSIONS.get(role, set())
        if tool_name not in allowed:
            raise PermissionError(
                f"Role '{role}' does not have permission to use tool '{tool_name}'"
            )

        if tool_name == "file_read":
            return self._file_read(arguments.get("path", ""))
        elif tool_name == "file_write":
            return self._file_write(
                arguments.get("path", ""), arguments.get("content", "")
            )
        elif tool_name == "shell_execute":
            return await self._shell_execute(
                arguments.get("command", ""),
                arguments.get("timeout", self.shell_timeout),
            )
        elif tool_name == "autolab_cite":
            return self._autolab_cite(
                arguments.get("action", "search"),
                arguments.get("query", ""),
                arguments.get("count", 5),
            )
        else:
            return f"Unknown tool: {tool_name}"

    def _resolve_path(self, relative_path: str) -> Path:
        """Resolve a relative path within the project dir, preventing traversal."""
        resolved = (Path(self.project_dir) / relative_path).resolve()
        project_resolved = Path(self.project_dir).resolve()
        if not str(resolved).startswith(str(project_resolved)):
            raise PermissionError(
                f"Path traversal blocked: {relative_path} resolves outside project"
            )
        return resolved

    def _file_read(self, path: str) -> str:
        """Read a file within the project directory."""
        try:
            resolved = self._resolve_path(path)
            if not resolved.exists():
                return f"Error: File not found: {path}"
            content = resolved.read_text(encoding="utf-8", errors="replace")
            # Cap at 50k chars to avoid blowing up context
            if len(content) > 50_000:
                return content[:50_000] + f"\n\n... (truncated, {len(content)} chars total)"
            return content
        except PermissionError as e:
            return f"Error: {e}"
        except Exception as e:
            return f"Error reading {path}: {e}"

    def _file_write(self, path: str, content: str) -> str:
        """Write a file within the project directory."""
        try:
            resolved = self._resolve_path(path)
            resolved.parent.mkdir(parents=True, exist_ok=True)
            resolved.write_text(content, encoding="utf-8")
            return f"Written {len(content)} chars to {path}"
        except PermissionError as e:
            return f"Error: {e}"
        except Exception as e:
            return f"Error writing {path}: {e}"

    async def _shell_execute(self, command: str, timeout: int) -> str:
        """Execute a shell command in the project directory."""
        timeout = min(timeout, 300)  # Hard cap at 5 minutes

        try:
            proc = await asyncio.create_subprocess_shell(
                command,
                cwd=self.project_dir,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            stdout, stderr = await asyncio.wait_for(
                proc.communicate(), timeout=timeout
            )
        except asyncio.TimeoutError:
            return f"Error: Command timed out after {timeout}s"
        except Exception as e:
            return f"Error executing command: {e}"

        output = stdout.decode("utf-8", errors="replace")
        err = stderr.decode("utf-8", errors="replace")
        result = ""
        if output:
            result += output
        if err:
            result += f"\nSTDERR:\n{err}"
        if proc.returncode != 0:
            result += f"\nExit code: {proc.returncode}"

        # Cap output length
        if len(result) > 30_000:
            result = result[:30_000] + "\n... (truncated)"
        return result or "(no output)"

    def _autolab_cite(self, action: str, query: str, count: int) -> str:
        """Execute a citation operation via the citations module."""
        try:
            from .citations import search_papers, validate_bibtex_file
        except ImportError:
            return "Error: citations module not available"

        if action == "search":
            results = search_papers(query, rows=min(count, 10))
            if not results:
                return "No papers found."
            parts = []
            for r in results:
                if "error" in r:
                    parts.append(f"Error: {r['error']}")
                else:
                    parts.append(
                        f"**{r['title']}** ({r['year']})\n"
                        f"DOI: {r['doi']}\n"
                        f"```bibtex\n{r['bibtex']}\n```"
                    )
            return "\n\n".join(parts)

        elif action == "doi":
            results = search_papers(f"doi:{query}", rows=1)
            if results and "bibtex" in results[0]:
                return f"```bibtex\n{results[0]['bibtex']}\n```"
            return f"Could not resolve DOI: {query}"

        elif action == "validate":
            bib_path = query or os.path.join(
                self.project_dir, "paper", "references.bib"
            )
            report = validate_bibtex_file(bib_path)
            return json.dumps(report, indent=2)

        return f"Unknown citation action: {action}"
