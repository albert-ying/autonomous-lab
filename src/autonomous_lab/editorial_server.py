"""
Standalone Editorial Center server.

A lightweight FastAPI app that serves the Editorial Center dashboard
on a fixed, bookmarkable port. It reads from the shared registry at
~/.autolab/registry.json to show all projects across all MCP instances.

Usage:
    autonomous-lab editorial              # http://127.0.0.1:8700
    autonomous-lab editorial --port 9000  # http://127.0.0.1:9000
"""

import json
import sys
from pathlib import Path

import uvicorn
from fastapi import FastAPI
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi import Request

from . import __version__


def run_editorial_server(host: str = "127.0.0.1", port: int = 8700):
    """Start the standalone editorial center."""
    app = FastAPI(title="Autonomous Lab — Editorial Center")

    # Reuse the existing templates and static files
    web_dir = Path(__file__).parent / "web"
    templates_dir = web_dir / "templates"
    static_dir = web_dir / "static"

    if static_dir.exists():
        app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")

    templates = Jinja2Templates(directory=str(templates_dir))

    @app.get("/", response_class=HTMLResponse)
    @app.get("/editorial-center", response_class=HTMLResponse)
    async def editorial_center(request: Request):
        return templates.TemplateResponse(
            request,
            "editorial_center.html",
            {"title": "Editorial Center", "version": __version__},
        )

    @app.get("/api/projects")
    async def list_projects():
        """Read the shared registry and return all projects."""
        from .web.utils.project_registry import get_all_projects
        from .lab.state import load_config, load_state

        all_projects = get_all_projects()

        result = []
        for p in all_projects:
            pdir = p.get("project_dir", "")
            try:
                state = load_state(pdir)
                config = load_config(pdir)
                editorial = state.get("editorial", {})
                result.append({
                    "project_dir": pdir,
                    "name": (
                        config.get("title")
                        or state.get("title")
                        or Path(pdir).name
                    ),
                    "iteration": state.get("iteration", 0),
                    "next_role": state.get("next_role", "pi"),
                    "status": state.get("status", "active"),
                    "progress": state.get("progress", 0),
                    "domain": config.get("domain", "research"),
                    "created_at": state.get("created_at", ""),
                    "editorial_phase": editorial.get("phase", "none"),
                    "editorial_decision": editorial.get("decision", ""),
                    "host": p.get("host", "127.0.0.1"),
                    "port": p.get("port", 8765),
                    "base_url": p.get("base_url", ""),
                    "lab_url": p.get("lab_url", ""),
                })
            except Exception:
                result.append({
                    **p,
                    "status": p.get("status", "unknown"),
                })

        return JSONResponse(content={"projects": result})

    print(f"\n  Autonomous Lab — Editorial Center v{__version__}")
    print(f"  http://{host}:{port}")
    print(f"  Bookmark this URL. It shows all projects across all instances.\n")

    try:
        uvicorn.run(app, host=host, port=port, log_level="warning")
    except KeyboardInterrupt:
        print("\nEditorial Center stopped.")
        sys.exit(0)
