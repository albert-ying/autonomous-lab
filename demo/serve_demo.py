#!/usr/bin/env python3
"""Direct FastAPI server for demo — bypasses WebUIManager complexity."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import uvicorn
from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse, JSONResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
import json, os

DEMO_DIR = "/tmp/autolab_demo_full"
WEB_DIR = Path(__file__).parent.parent / "src" / "autonomous_lab" / "web"

app = FastAPI()

# Static files
app.mount("/static", StaticFiles(directory=str(WEB_DIR / "static")), name="static")
templates = Jinja2Templates(directory=str(WEB_DIR / "templates"))

def read_state():
    p = Path(DEMO_DIR) / ".autolab" / "state.json"
    if p.exists():
        return json.loads(p.read_text())
    return {}

def read_editorial():
    s = read_state()
    return s.get("editorial", {"phase": "none"})

@app.get("/lab", response_class=HTMLResponse)
async def lab_page(request: Request):
    return templates.TemplateResponse("lab.html", {"request": request, "version": "0.5.0"})

@app.get("/api/autolab/state")
async def get_state():
    try:
        state = read_state()
        base = Path(DEMO_DIR)

        # Build file listings
        files = {"figures": [], "scripts": [], "results": [], "data": [], "paper": []}
        for cat in files:
            d = base / cat
            if d.exists():
                for f in sorted(d.rglob("*")):
                    if f.is_file():
                        files[cat].append(str(f.relative_to(base)))

        # Parse meeting log
        idea = ""
        idea_path = base / ".autolab" / "idea.md"
        if idea_path.exists():
            idea = idea_path.read_text()

        return JSONResponse(content={
            "active": True,
            "project_dir": DEMO_DIR,
            "iteration": state.get("iteration", 0),
            "next_role": state.get("next_role", "pi"),
            "status": state.get("status", "active"),
            "progress": state.get("progress", 0),
            "experts": state.get("experts", []),
            "editorial": state.get("editorial", {"phase": "none"}),
            "files": files,
            "idea": idea,
            "user_feedback": state.get("user_feedback", ""),
        })
    except Exception as e:
        return JSONResponse(content={"error": str(e)}, status_code=500)

@app.get("/api/autolab/meeting-log")
async def get_meeting_log():
    try:
        log_path = Path(DEMO_DIR) / ".autolab" / "meeting_log.md"
        if not log_path.exists():
            return JSONResponse(content={"turns": []})
        text = log_path.read_text()

        # Parse turns
        turns = []
        import re
        parts = re.split(r"---\s*\n", text)
        for part in parts:
            part = part.strip()
            if not part:
                continue
            m = re.search(r"Iteration\s+(\d+)\s+.*?(PI|TRAINEE|REVIEWER\S*)\s+Turn", part, re.I)
            if m:
                role = m.group(2).lower()
                iteration = int(m.group(1))
                # Extract summary
                sm = re.search(r"### Summary\s*\n(.+?)(?:\n###|\Z)", part, re.S)
                summary = sm.group(1).strip() if sm else ""
                # Extract details
                dm = re.search(r"### Details\s*\n(.+?)(?:\Z)", part, re.S)
                details = dm.group(1).strip() if dm else part
                turns.append({
                    "role": role,
                    "iteration": iteration,
                    "summary": summary,
                    "content": details,
                })
        return JSONResponse(content={"turns": turns})
    except Exception as e:
        return JSONResponse(content={"turns": [], "error": str(e)})

@app.get("/api/autolab/file")
async def get_file(path: str = ""):
    try:
        full = os.path.join(DEMO_DIR, path)
        full = os.path.normpath(full)
        if not full.startswith(DEMO_DIR):
            return JSONResponse(status_code=403, content={"error": "forbidden"})
        if os.path.isfile(full):
            return FileResponse(full)
        return JSONResponse(status_code=404, content={"error": "not found"})
    except Exception as e:
        return JSONResponse(status_code=500, content={"error": str(e)})

@app.get("/api/autolab/editorial")
async def get_editorial():
    return JSONResponse(content=read_editorial())

@app.post("/api/autolab/editorial/invite-reviewers")
async def invite_reviewers(request: Request):
    body = await request.json()
    reviewers = body.get("reviewers", [])
    state = read_state()
    ed = state.get("editorial", {})
    ed["phase"] = "reviewers_invited"
    ed["reviewers"] = [{"name": r["name"], "role": r["role"], "id": f"reviewer_{i}"}
                       for i, r in enumerate(reviewers)]
    state["editorial"] = ed
    state["status"] = "under_review"
    p = Path(DEMO_DIR) / ".autolab" / "state.json"
    p.write_text(json.dumps(state, indent=2))
    return JSONResponse(content={"ok": True, "editorial": ed})

@app.post("/api/autolab/editorial/decision")
async def editorial_decision(request: Request):
    body = await request.json()
    decision = body.get("decision", "")
    feedback = body.get("feedback", "")
    print(f"  [DECISION] decision={decision}, feedback={feedback[:60] if feedback else '(none)'}...")
    state = read_state()
    ed = state.get("editorial", {})
    ed["phase"] = "decision_made"
    ed["decision"] = decision
    ed["decision_feedback"] = feedback
    state["editorial"] = ed
    p = Path(DEMO_DIR) / ".autolab" / "state.json"
    p.write_text(json.dumps(state, indent=2))
    return JSONResponse(content={"ok": True})

@app.post("/api/autolab/editorial/desk-reject")
async def desk_reject(request: Request):
    state = read_state()
    ed = state.get("editorial", {})
    ed["phase"] = "decision_made"
    ed["decision"] = "reject"
    ed["decision_feedback"] = "Desk rejected by editor."
    state["editorial"] = ed
    state["status"] = "rejected"
    p = Path(DEMO_DIR) / ".autolab" / "state.json"
    p.write_text(json.dumps(state, indent=2))
    return JSONResponse(content={"ok": True})

if __name__ == "__main__":
    # Kill any old process on 8765
    import subprocess
    subprocess.run("lsof -ti:8765 | xargs kill -9 2>/dev/null", shell=True, capture_output=True)

    print("\n" + "="*50)
    print("  Autonomous Lab Demo")
    print(f"  Project: {DEMO_DIR}")
    print(f"  URL: http://127.0.0.1:8765/lab")
    print("="*50 + "\n")

    uvicorn.run(app, host="127.0.0.1", port=8765, log_level="info")
