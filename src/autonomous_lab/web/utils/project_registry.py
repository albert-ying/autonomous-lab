"""
Cross-instance project registry.

All MCP server instances share a single registry file at
~/.autolab/registry.json. Each instance writes its own entry
(keyed by host:port) listing which projects it serves.

The Editorial Center reads this file to show ALL projects
across ALL running instances, linking to the correct port.
"""

import json
import os
import time
from pathlib import Path
from typing import Any

from ...debug import web_debug_log as debug_log

_REGISTRY_DIR = Path.home() / ".autolab"
_REGISTRY_FILE = _REGISTRY_DIR / "registry.json"
_STALE_SECONDS = 120  # entries older than 2 min without heartbeat = stale


def _read_registry() -> dict:
    """Read the shared registry, returning {} on any error."""
    try:
        if _REGISTRY_FILE.exists():
            data = json.loads(_REGISTRY_FILE.read_text(encoding="utf-8"))
            if isinstance(data, dict):
                return data
    except Exception as e:
        debug_log(f"Registry read error: {e}")
    return {}


def _write_registry(data: dict) -> None:
    """Write the registry atomically."""
    try:
        _REGISTRY_DIR.mkdir(parents=True, exist_ok=True)
        tmp = _REGISTRY_FILE.with_suffix(".tmp")
        tmp.write_text(json.dumps(data, indent=2), encoding="utf-8")
        tmp.replace(_REGISTRY_FILE)
    except Exception as e:
        debug_log(f"Registry write error: {e}")


def register_instance(host: str, port: int) -> None:
    """Register this MCP server instance in the shared registry."""
    reg = _read_registry()
    key = f"{host}:{port}"
    if key not in reg:
        reg[key] = {
            "host": host,
            "port": port,
            "projects": {},
            "started_at": time.time(),
        }
    reg[key]["last_seen"] = time.time()
    _write_registry(reg)
    debug_log(f"Registered instance {key} in shared registry")


def register_project(
    host: str,
    port: int,
    project_dir: str,
    project_name: str,
    state: dict | None = None,
    config: dict | None = None,
) -> None:
    """Register or update a project under this instance."""
    reg = _read_registry()
    key = f"{host}:{port}"

    if key not in reg:
        reg[key] = {
            "host": host,
            "port": port,
            "projects": {},
            "started_at": time.time(),
        }

    state = state or {}
    config = config or {}
    editorial = state.get("editorial", {})

    reg[key]["last_seen"] = time.time()
    reg[key]["projects"][project_dir] = {
        "name": project_name,
        "iteration": state.get("iteration", 0),
        "next_role": state.get("next_role", "pi"),
        "status": state.get("status", "active"),
        "progress": state.get("progress", 0),
        "domain": config.get("domain", "research"),
        "created_at": state.get("created_at", ""),
        "editorial_phase": editorial.get("phase", "none"),
        "editorial_decision": editorial.get("decision", ""),
        "updated_at": time.time(),
    }

    _write_registry(reg)


def heartbeat(host: str, port: int) -> None:
    """Update last_seen timestamp for this instance."""
    reg = _read_registry()
    key = f"{host}:{port}"
    if key in reg:
        reg[key]["last_seen"] = time.time()
        _write_registry(reg)


def unregister_instance(host: str, port: int) -> None:
    """Remove this instance from the registry (on shutdown)."""
    reg = _read_registry()
    key = f"{host}:{port}"
    if key in reg:
        del reg[key]
        _write_registry(reg)
        debug_log(f"Unregistered instance {key} from shared registry")


def get_all_projects() -> list[dict[str, Any]]:
    """Get all projects across all instances, cleaning stale entries.

    Returns a flat list of project dicts, each with host/port/url fields
    so the Editorial Center can link to the correct instance.
    """
    reg = _read_registry()
    now = time.time()
    result = []
    stale_keys = []

    for key, instance in reg.items():
        last_seen = instance.get("last_seen", 0)
        if now - last_seen > _STALE_SECONDS:
            # Check if the port is actually still listening before removing
            host = instance.get("host", "127.0.0.1")
            port = instance.get("port", 0)
            if not _is_port_alive(host, port):
                stale_keys.append(key)
                continue

        display_host = "127.0.0.1" if instance.get("host") in ("0.0.0.0", "::") else instance.get("host", "127.0.0.1")
        port = instance.get("port", 8765)
        base_url = f"http://{display_host}:{port}"

        for pdir, pinfo in instance.get("projects", {}).items():
            result.append({
                **pinfo,
                "project_dir": pdir,
                "host": display_host,
                "port": port,
                "base_url": base_url,
                "lab_url": f"{base_url}/lab?project={pdir}",
                "instance_key": key,
            })

    # Clean stale entries
    if stale_keys:
        for k in stale_keys:
            del reg[k]
        _write_registry(reg)
        debug_log(f"Cleaned {len(stale_keys)} stale instance(s) from registry")

    return result


def _is_port_alive(host: str, port: int) -> bool:
    """Quick check if a port is still listening."""
    import socket
    check_host = "127.0.0.1" if host in ("0.0.0.0", "::") else host
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.settimeout(1)
            s.connect((check_host, port))
            return True
    except Exception:
        return False
