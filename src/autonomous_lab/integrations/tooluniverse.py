"""
ToolUniverse integration for Autonomous Lab.

ToolUniverse (https://aiscientist.tools) provides 1000+ scientific tools
(ML models, datasets, APIs, analysis packages) with a standardized
AI-Tool Interaction Protocol.

ARCHITECTURE
  This module wraps ToolUniverse in two ways:
    1. Local SDK (preferred): pip install tooluniverse
    2. HTTP API fallback: https://aiscientist.tools

  During RECRUITMENT, tools are discovered and auto-converted to SKILL.md
  files that integrate with the character/skill system.

  During the ACTIVE PHASE, trainees execute tools via:
    from tooluniverse import ToolUniverse
    tu = ToolUniverse()
    result = tu.tools.<tool_name>(**kwargs)

REFERENCE
  Gao et al., "Democratizing AI scientists using ToolUniverse" (2025)
  https://github.com/mims-harvard/ToolUniverse
"""

from __future__ import annotations

import importlib
import json
import logging
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)

# ── singleton state ─────────────────────────────────────────────
_sdk_available: Optional[bool] = None
_tu_instance: Any = None

# ── config ──────────────────────────────────────────────────────
TOOLS_CFG_NAME = "tools_config.json"
_default_config: dict = {
    "enabled": True,
    "prefer_local": True,
    "api_url": "https://aiscientist.tools",
    "api_key": None,
}

API_BASE = "https://aiscientist.tools"


# ==================================================================
# Detection
# ==================================================================
def is_sdk_available() -> bool:
    """Return True if ``tooluniverse`` is importable."""
    global _sdk_available
    if _sdk_available is None:
        try:
            importlib.import_module("tooluniverse")
            _sdk_available = True
        except ImportError:
            _sdk_available = False
    return _sdk_available


def get_sdk_version() -> Optional[str]:
    """Return the installed ToolUniverse version string, or None."""
    if not is_sdk_available():
        return None
    try:
        mod = importlib.import_module("tooluniverse")
        return getattr(mod, "__version__", "unknown")
    except Exception:
        return None


def _get_tu_instance() -> Any:
    """Get or create a singleton ToolUniverse instance."""
    global _tu_instance
    if _tu_instance is None and is_sdk_available():
        try:
            from tooluniverse import ToolUniverse
            _tu_instance = ToolUniverse()
        except Exception as exc:
            logger.warning("Failed to create ToolUniverse instance: %s", exc)
    return _tu_instance


# ==================================================================
# Installation
# ==================================================================
def install_sdk(from_source: bool = True) -> dict:
    """Install tooluniverse at runtime. Users must explicitly opt in."""
    global _sdk_available
    try:
        target = (
            "git+https://github.com/mims-harvard/ToolUniverse.git@main"
            if from_source
            else "tooluniverse --upgrade"
        )
        cmd = [sys.executable, "-m", "pip", "install"] + target.split()
        logger.info("Installing ToolUniverse: %s", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode == 0:
            _sdk_available = None  # force re-check
            if is_sdk_available():
                return {"success": True, "message": f"ToolUniverse installed (v{get_sdk_version()})."}
            return {"success": False, "message": "pip succeeded but import still fails."}
        snippet = (result.stderr or result.stdout)[-500:]
        return {"success": False, "message": f"pip failed:\n{snippet}"}
    except subprocess.TimeoutExpired:
        return {"success": False, "message": "Installation timed out (300 s)."}
    except Exception as exc:
        return {"success": False, "message": str(exc)}


# ==================================================================
# Configuration (per-project, stored in .autolab/)
# ==================================================================
def load_config(project_directory: str = ".") -> dict:
    cfg = dict(_default_config)
    cfg_path = Path(project_directory) / ".autolab" / TOOLS_CFG_NAME
    if cfg_path.exists():
        try:
            with open(cfg_path) as f:
                cfg.update(json.load(f))
        except Exception:
            pass
    return cfg


def save_config(project_directory: str = ".", **overrides: Any) -> dict:
    cfg = load_config(project_directory)
    cfg.update(overrides)
    cfg_path = Path(project_directory) / ".autolab" / TOOLS_CFG_NAME
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cfg_path, "w") as f:
        json.dump(cfg, f, indent=2)
    return cfg


# ==================================================================
# Tool Discovery
# ==================================================================
def search_tools(query: str, limit: int = 10) -> list[dict]:
    """Search ToolUniverse for tools matching a natural language query.

    Uses local SDK if available, falls back to HTTP API.
    Returns list of dicts: [{name, description, category}, ...]
    """
    if is_sdk_available():
        try:
            return _sdk_search_tools(query, limit)
        except Exception as exc:
            logger.warning("SDK search failed, trying HTTP: %s", exc)

    try:
        return _http_search_tools(query, limit)
    except Exception as exc:
        logger.warning("HTTP search failed: %s", exc)
        return []


def _sdk_search_tools(query: str, limit: int) -> list[dict]:
    """Search using local ToolUniverse SDK."""
    tu = _get_tu_instance()
    if tu is None:
        return []

    results = []
    try:
        from tooluniverse.tool_registry import get_tool_registry
        registry = get_tool_registry()
        query_lower = query.lower()
        for name, tool_cls in registry.items():
            desc = getattr(tool_cls, "__doc__", "") or ""
            if query_lower in name.lower() or query_lower in desc.lower():
                results.append({
                    "name": name,
                    "description": desc.strip().split("\n")[0],
                })
            if len(results) >= limit:
                break
    except Exception:
        pass
    return results


def _http_search_tools(query: str, limit: int) -> list[dict]:
    """Search using ToolUniverse HTTP API."""
    encoded = urllib.parse.quote(query)
    url = f"{API_BASE}/api/search?q={encoded}&limit={limit}"
    req = urllib.request.Request(url, headers={
        "User-Agent": "AutonomousLab/0.5",
        "Accept": "application/json",
    })
    with urllib.request.urlopen(req, timeout=15) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    tools = data if isinstance(data, list) else data.get("tools", data.get("results", []))
    return [{"name": t.get("name", ""), "description": t.get("description", "")} for t in tools[:limit]]


def get_tool_spec(tool_name: str) -> Optional[dict]:
    """Get full specification for a named tool.

    Returns dict with: name, description, parameters, etc.
    """
    if is_sdk_available():
        try:
            return _sdk_get_tool_spec(tool_name)
        except Exception:
            pass

    try:
        return _http_get_tool_spec(tool_name)
    except Exception:
        return None


def _sdk_get_tool_spec(tool_name: str) -> Optional[dict]:
    """Get tool spec via local SDK."""
    tu = _get_tu_instance()
    if tu is None:
        return None
    try:
        spec = tu.tool_specification(tool_name)
        return spec if spec else None
    except Exception:
        return None


def _http_get_tool_spec(tool_name: str) -> Optional[dict]:
    """Get tool spec via HTTP API."""
    encoded = urllib.parse.quote(tool_name)
    url = f"{API_BASE}/api/tools/{encoded}"
    req = urllib.request.Request(url, headers={
        "User-Agent": "AutonomousLab/0.5",
        "Accept": "application/json",
    })
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except Exception:
        return None


# ==================================================================
# SKILL.md Generation (for recruitment pipeline)
# ==================================================================
def generate_skill_md(spec: dict) -> str:
    """Generate a SKILL.md from a ToolUniverse tool specification.

    The generated file is compatible with the Autonomous Lab skill system
    and uses content-only validation (no validation.py needed).
    """
    name = spec.get("name", "unknown_tool")
    description = spec.get("description", "No description available.")
    parameters = spec.get("parameters", [])

    lines = [
        f"# {name}",
        "",
        f"Source: [ToolUniverse](https://aiscientist.tools) — auto-generated from tool catalog.",
        "",
        "## Description",
        "",
        description,
        "",
        "## When to Use",
        "",
        f"Use this tool when you need to {description.lower().rstrip('.')}.",
        "This tool is provided by the ToolUniverse scientific tool catalog.",
        "",
        "## Usage",
        "",
        "```python",
        "from tooluniverse import ToolUniverse",
        "",
        "tu = ToolUniverse()",
        f"result = tu.tools.{name}(",
    ]

    if parameters:
        param_args = []
        for p in parameters:
            pname = p.get("name", "param")
            ptype = p.get("type", "str")
            required = p.get("required", False)
            if required:
                param_args.append(f"    {pname}=<{ptype}>,  # required")
            else:
                param_args.append(f"    {pname}=<{ptype}>,  # optional")
        lines.extend(param_args)
    lines.extend([
        ")",
        "```",
        "",
    ])

    if parameters:
        lines.extend([
            "## Parameters",
            "",
            "| Name | Type | Required | Description |",
            "|------|------|----------|-------------|",
        ])
        for p in parameters:
            pname = p.get("name", "param")
            ptype = p.get("type", "string")
            required = "Yes" if p.get("required", False) else "No"
            pdesc = p.get("description", "")
            lines.append(f"| `{pname}` | {ptype} | {required} | {pdesc} |")
        lines.append("")

    lines.extend([
        "## Execution",
        "",
        "If the ToolUniverse SDK is installed locally (`pip install tooluniverse`),",
        "the tool runs locally. Otherwise, it executes via the ToolUniverse HTTP API",
        "at https://aiscientist.tools.",
        "",
        "## Standard Workflow",
        "",
        "1. Import ToolUniverse and create an instance",
        "2. Call the tool with the required parameters",
        "3. Parse the result dictionary for your analysis",
        "4. Save results to the `results/` directory",
        "5. Document findings in the paper",
        "",
        "## Key Decisions",
        "",
        "- Always validate tool output before using in downstream analysis",
        "- If the tool fails or times out, document the failure and proceed with alternatives",
        "- Prefer local SDK execution for reproducibility; use HTTP API for heavy compute",
    ])

    return "\n".join(lines) + "\n"


# ==================================================================
# Status summary
# ==================================================================
def get_status(project_directory: str = ".") -> dict:
    """Return a dict summarising the ToolUniverse integration state."""
    cfg = load_config(project_directory)
    installed = is_sdk_available()

    return {
        "installed": installed,
        "version": get_sdk_version(),
        "enabled": cfg.get("enabled", True),
        "prefer_local": cfg.get("prefer_local", True),
        "api_url": cfg.get("api_url", API_BASE),
    }
