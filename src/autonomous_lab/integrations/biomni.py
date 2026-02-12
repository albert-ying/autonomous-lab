"""
Optional Biomni integration for Autonomous Lab.

Biomni (https://github.com/snap-stanford/Biomni) provides a curated
collection of 100+ biomedical tools, databases, and a data lake for
computational biology research.

ARCHITECTURE
  This module is a *thin wrapper* — it NEVER ships or bundles Biomni.
  Instead it:
    1. Detects whether Biomni is importable in the current env.
    2. Offers to pip-install Biomni from GitHub on-the-fly.
    3. Exposes Biomni's *tools and databases* for direct use by
       the Autonomous Lab Trainee in their Python scripts.

  IMPORTANT: We do NOT instantiate the Biomni A1 agent or use their
  LLM pipeline.  We treat Biomni purely as a tool/database library
  with a curated Python environment.  The Trainee imports and calls
  individual functions directly.

LICENSING NOTE
  Biomni itself is Apache 2.0, but some integrated tools, databases,
  and third-party software may carry more restrictive licenses.
  Review each component before commercial use.
"""

from __future__ import annotations

import importlib
import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)

# ── singleton state ─────────────────────────────────────────────
_biomni_available: Optional[bool] = None

# ── config ──────────────────────────────────────────────────────
BIOMNI_CFG_NAME = "biomni_config.json"
_default_config: dict = {
    "enabled": False,
    "data_path": "./data",
    "skip_datalake": False,
}


# ==================================================================
# Detection
# ==================================================================
def is_biomni_available() -> bool:
    """Return True if ``biomni`` is importable."""
    global _biomni_available
    if _biomni_available is None:
        try:
            importlib.import_module("biomni")
            _biomni_available = True
        except ImportError:
            _biomni_available = False
    return _biomni_available


def get_biomni_version() -> Optional[str]:
    """Return the installed Biomni version string, or None."""
    if not is_biomni_available():
        return None
    try:
        mod = importlib.import_module("biomni")
        return getattr(mod, "__version__", "unknown")
    except Exception:
        return None


def check_environment() -> dict:
    """
    Check what Biomni components are available in the current env.

    Returns a dict describing:
      - whether biomni is importable
      - whether key sub-packages exist (tools, data loaders, etc.)
      - whether the datalake directory exists and has files
      - which conda env is active (if any)
    """
    import os
    info: dict[str, Any] = {
        "installed": is_biomni_available(),
        "version": get_biomni_version(),
        "conda_env": os.environ.get("CONDA_DEFAULT_ENV", None),
        "python": sys.executable,
    }

    if is_biomni_available():
        # Check key sub-packages
        for pkg in ("biomni.tools", "biomni.data", "biomni.utils"):
            try:
                importlib.import_module(pkg)
                info[pkg] = True
            except ImportError:
                info[pkg] = False
    return info


# ==================================================================
# Installation
# ==================================================================
def install_biomni(from_source: bool = True) -> dict:
    """
    Install Biomni at runtime.  Users must explicitly opt in.

    Parameters
    ----------
    from_source : bool
        True (default) → pip install from GitHub main.
        False → install latest PyPI release.

    Returns
    -------
    dict  {"success": bool, "message": str}
    """
    global _biomni_available
    try:
        target = (
            "git+https://github.com/snap-stanford/Biomni.git@main"
            if from_source
            else "biomni --upgrade"
        )
        cmd = [sys.executable, "-m", "pip", "install"] + target.split()
        logger.info("Installing Biomni: %s", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode == 0:
            _biomni_available = None  # force re-check
            if is_biomni_available():
                return {"success": True, "message": f"Biomni installed (v{get_biomni_version()})."}
            return {"success": False, "message": "pip succeeded but import still fails."}
        snippet = (result.stderr or result.stdout)[-500:]
        return {"success": False, "message": f"pip failed:\n{snippet}"}
    except subprocess.TimeoutExpired:
        return {"success": False, "message": "Installation timed out (600 s)."}
    except Exception as exc:
        return {"success": False, "message": str(exc)}


# ==================================================================
# Configuration (per-project, stored in .autolab/)
# ==================================================================
def load_biomni_config(project_directory: str = ".") -> dict:
    cfg = dict(_default_config)
    cfg_path = Path(project_directory) / ".autolab" / BIOMNI_CFG_NAME
    if cfg_path.exists():
        try:
            with open(cfg_path) as f:
                cfg.update(json.load(f))
        except Exception:
            pass
    return cfg


def save_biomni_config(project_directory: str = ".", **overrides: Any) -> dict:
    cfg = load_biomni_config(project_directory)
    cfg.update(overrides)
    cfg_path = Path(project_directory) / ".autolab" / BIOMNI_CFG_NAME
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cfg_path, "w") as f:
        json.dump(cfg, f, indent=2)
    return cfg


# ==================================================================
# Tool & database introspection
# ==================================================================
def list_available_tools() -> list[dict]:
    """
    List Biomni's available tools with descriptions.

    Each tool is returned as:
      {"name": "...", "module": "biomni.tools.xxx", "description": "..."}
    """
    if not is_biomni_available():
        return []
    tools: list[dict] = []
    try:
        # Attempt to access tools package
        tools_pkg = importlib.import_module("biomni.tools")
        if hasattr(tools_pkg, "__all__"):
            for name in tools_pkg.__all__:
                mod_path = f"biomni.tools.{name}"
                desc = ""
                try:
                    m = importlib.import_module(mod_path)
                    desc = getattr(m, "__doc__", "") or ""
                    desc = desc.strip().split("\n")[0]  # first line
                except Exception:
                    pass
                tools.append({"name": name, "module": mod_path, "description": desc})
        else:
            # Fallback: inspect package directory
            import pkgutil
            pkg_path = getattr(tools_pkg, "__path__", [])
            for importer, modname, ispkg in pkgutil.iter_modules(pkg_path):
                if modname.startswith("_"):
                    continue
                mod_path = f"biomni.tools.{modname}"
                desc = ""
                try:
                    m = importlib.import_module(mod_path)
                    desc = (getattr(m, "__doc__", "") or "").strip().split("\n")[0]
                except Exception:
                    pass
                tools.append({"name": modname, "module": mod_path, "description": desc})
    except Exception as exc:
        tools.append({"name": "(error)", "module": "", "description": str(exc)})
    return tools


def list_available_databases() -> list[dict]:
    """
    List Biomni's available database connectors/data loaders.
    """
    if not is_biomni_available():
        return []
    databases: list[dict] = []
    try:
        data_pkg = importlib.import_module("biomni.data")
        import pkgutil
        pkg_path = getattr(data_pkg, "__path__", [])
        for importer, modname, ispkg in pkgutil.iter_modules(pkg_path):
            if modname.startswith("_"):
                continue
            mod_path = f"biomni.data.{modname}"
            desc = ""
            try:
                m = importlib.import_module(mod_path)
                desc = (getattr(m, "__doc__", "") or "").strip().split("\n")[0]
            except Exception:
                pass
            databases.append({"name": modname, "module": mod_path, "description": desc})
    except Exception:
        pass
    return databases


def get_datalake_path(project_directory: str = ".") -> Optional[str]:
    """Return the path to Biomni's data lake, or None."""
    cfg = load_biomni_config(project_directory)
    dp = Path(cfg.get("data_path", "./data"))
    if dp.exists() and any(dp.iterdir()):
        return str(dp.resolve())
    return None


def get_import_snippet(tool_name: str) -> str:
    """
    Return a Python code snippet the Trainee can use to import
    and call a Biomni tool in their scripts.
    """
    return (
        f"# Biomni tool: {tool_name}\n"
        f"from biomni.tools.{tool_name} import *\n"
        f"# See: https://github.com/snap-stanford/Biomni\n"
    )


# ==================================================================
# Status summary
# ==================================================================
def get_status(project_directory: str = ".") -> dict:
    """Return a dict summarising the Biomni integration state."""
    cfg = load_biomni_config(project_directory)
    installed = is_biomni_available()
    data_path = cfg.get("data_path", "./data")

    dl_exists = False
    if installed:
        dp = Path(data_path)
        dl_exists = dp.exists() and any(dp.iterdir()) if dp.exists() else False

    tool_count = len(list_available_tools()) if installed else 0
    db_count = len(list_available_databases()) if installed else 0

    return {
        "installed": installed,
        "version": get_biomni_version(),
        "enabled": cfg.get("enabled", False),
        "data_path": data_path,
        "skip_datalake": cfg.get("skip_datalake", False),
        "datalake_exists": dl_exists,
        "tool_count": tool_count,
        "database_count": db_count,
    }
