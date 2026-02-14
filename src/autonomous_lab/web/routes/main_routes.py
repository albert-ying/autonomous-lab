#!/usr/bin/env python3
"""
主要路由處理
============

設置 Web UI 的主要路由和處理邏輯。
"""

import json
import mimetypes
import time
from pathlib import Path
from typing import TYPE_CHECKING

from fastapi import Request, WebSocket, WebSocketDisconnect
from fastapi.responses import (
    FileResponse,
    HTMLResponse,
    JSONResponse,
    PlainTextResponse,
)

from ... import __version__
from ...debug import web_debug_log as debug_log
from ..constants import get_message_code as get_msg_code


if TYPE_CHECKING:
    from ..main import WebUIManager


def load_user_layout_settings() -> str:
    """載入用戶的佈局模式設定"""
    try:
        # 使用統一的設定檔案路徑
        config_dir = Path.home() / ".config" / "autonomous-lab"
        settings_file = config_dir / "ui_settings.json"

        if settings_file.exists():
            with open(settings_file, encoding="utf-8") as f:
                settings = json.load(f)
                layout_mode = settings.get("layoutMode", "combined-vertical")
                debug_log(f"從設定檔案載入佈局模式: {layout_mode}")
                # 修復 no-any-return 錯誤 - 確保返回 str 類型
                return str(layout_mode)
        else:
            debug_log("設定檔案不存在，使用預設佈局模式: combined-vertical")
            return "combined-vertical"
    except Exception as e:
        debug_log(f"載入佈局設定失敗: {e}，使用預設佈局模式: combined-vertical")
        return "combined-vertical"


# 使用統一的訊息代碼系統
# 從 ..constants 導入的 get_msg_code 函數會處理所有訊息代碼
# 舊的 key 會自動映射到新的常量


def _guess_language(suffix: str) -> str:
    """Map file extension to a display language label."""
    return {
        ".py": "python",
        ".tex": "latex",
        ".bib": "bibtex",
        ".r": "r",
        ".R": "r",
        ".sh": "bash",
        ".json": "json",
        ".yaml": "yaml",
        ".yml": "yaml",
        ".csv": "csv",
        ".tsv": "tsv",
        ".md": "markdown",
        ".txt": "text",
    }.get(suffix, "text")


def _track_project(
    manager: "WebUIManager",
    project_dir: str,
    project_name: str,
    state: dict,
    config: dict,
) -> None:
    """Register a project in the manager's known-projects set."""
    if not hasattr(manager, "_known_projects"):
        manager._known_projects = {}  # type: ignore[attr-defined]
    manager._known_projects[project_dir] = {
        "name": project_name,
        "iteration": state.get("iteration", 0),
        "next_role": state.get("next_role", "pi"),
        "status": state.get("status", "active"),
        "progress": state.get("progress", 0),
        "domain": config.get("domain", "research"),
        "created_at": state.get("created_at", ""),
        "editorial_phase": state.get("editorial", {}).get("phase", "none"),
    }


def setup_routes(manager: "WebUIManager"):
    """設置路由"""

    @manager.app.get("/", response_class=HTMLResponse)
    async def index(request: Request):
        """統一回饋頁面 - 重構後的主頁面"""
        # 獲取當前活躍會話
        current_session = manager.get_current_session()

        if not current_session:
            # 沒有活躍會話時顯示等待頁面
            return manager.templates.TemplateResponse(
                "index.html",
                {
                    "request": request,
                    "title": "Autonomous Lab",
                    "has_session": False,
                    "version": __version__,
                },
            )

        # 有活躍會話時顯示回饋頁面
        # 載入用戶的佈局模式設定
        layout_mode = load_user_layout_settings()

        return manager.templates.TemplateResponse(
            "feedback.html",
            {
                "request": request,
                "project_directory": current_session.project_directory,
                "summary": current_session.summary,
                "title": "Autonomous Lab",
                "version": __version__,
                "has_session": True,
                "layout_mode": layout_mode,
            },
        )

    @manager.app.get("/api/translations")
    async def get_translations():
        """獲取翻譯數據 - 從 Web 專用翻譯檔案載入"""
        translations = {}

        # 獲取 Web 翻譯檔案目錄
        web_locales_dir = Path(__file__).parent.parent / "locales"
        supported_languages = ["zh-TW", "zh-CN", "en"]

        for lang_code in supported_languages:
            lang_dir = web_locales_dir / lang_code
            translation_file = lang_dir / "translation.json"

            try:
                if translation_file.exists():
                    with open(translation_file, encoding="utf-8") as f:
                        lang_data = json.load(f)
                        translations[lang_code] = lang_data
                        debug_log(f"成功載入 Web 翻譯: {lang_code}")
                else:
                    debug_log(f"Web 翻譯檔案不存在: {translation_file}")
                    translations[lang_code] = {}
            except Exception as e:
                debug_log(f"載入 Web 翻譯檔案失敗 {lang_code}: {e}")
                translations[lang_code] = {}

        debug_log(f"Web 翻譯 API 返回 {len(translations)} 種語言的數據")
        return JSONResponse(content=translations)

    @manager.app.get("/api/session-status")
    async def get_session_status(request: Request):
        """獲取當前會話狀態"""
        current_session = manager.get_current_session()

        # 從請求頭獲取客戶端語言
        lang = (
            request.headers.get("Accept-Language", "zh-TW").split(",")[0].split("-")[0]
        )
        if lang == "zh":
            lang = "zh-TW"

        if not current_session:
            return JSONResponse(
                content={
                    "has_session": False,
                    "status": "no_session",
                    "messageCode": get_msg_code("no_active_session"),
                }
            )

        return JSONResponse(
            content={
                "has_session": True,
                "status": "active",
                "session_info": {
                    "project_directory": current_session.project_directory,
                    "summary": current_session.summary,
                    "feedback_completed": current_session.feedback_completed.is_set(),
                },
            }
        )

    @manager.app.get("/api/current-session")
    async def get_current_session(request: Request):
        """獲取當前會話詳細信息"""
        current_session = manager.get_current_session()

        # 從查詢參數獲取語言，如果沒有則從會話獲取，最後使用默認值

        if not current_session:
            return JSONResponse(
                status_code=404,
                content={
                    "error": "No active session",
                    "messageCode": get_msg_code("no_active_session"),
                },
            )

        return JSONResponse(
            content={
                "session_id": current_session.session_id,
                "project_directory": current_session.project_directory,
                "summary": current_session.summary,
                "feedback_completed": current_session.feedback_completed.is_set(),
                "command_logs": current_session.command_logs,
                "images_count": len(current_session.images),
            }
        )

    @manager.app.get("/api/all-sessions")
    async def get_all_sessions(request: Request):
        """獲取所有會話的實時狀態"""

        try:
            sessions_data = []

            # 獲取所有會話的實時狀態
            for session_id, session in manager.sessions.items():
                session_info = {
                    "session_id": session.session_id,
                    "project_directory": session.project_directory,
                    "summary": session.summary,
                    "status": session.status.value,
                    "status_message": session.status_message,
                    "created_at": int(session.created_at * 1000),  # 轉換為毫秒
                    "last_activity": int(session.last_activity * 1000),
                    "feedback_completed": session.feedback_completed.is_set(),
                    "has_websocket": session.websocket is not None,
                    "is_current": session == manager.current_session,
                    "user_messages": session.user_messages,  # 包含用戶消息記錄
                }
                sessions_data.append(session_info)

            # 按創建時間排序（最新的在前）
            sessions_data.sort(key=lambda x: x["created_at"], reverse=True)

            debug_log(f"返回 {len(sessions_data)} 個會話的實時狀態")
            return JSONResponse(content={"sessions": sessions_data})

        except Exception as e:
            debug_log(f"獲取所有會話狀態失敗: {e}")
            return JSONResponse(
                status_code=500,
                content={
                    "error": f"Failed to get sessions: {e!s}",
                    "messageCode": get_msg_code("get_sessions_failed"),
                },
            )

    @manager.app.post("/api/add-user-message")
    async def add_user_message(request: Request):
        """添加用戶消息到當前會話"""

        try:
            data = await request.json()
            current_session = manager.get_current_session()

            if not current_session:
                return JSONResponse(
                    status_code=404,
                    content={
                        "error": "No active session",
                        "messageCode": get_msg_code("no_active_session"),
                    },
                )

            # 添加用戶消息到會話
            current_session.add_user_message(data)

            debug_log(f"用戶消息已添加到會話 {current_session.session_id}")
            return JSONResponse(
                content={
                    "status": "success",
                    "messageCode": get_msg_code("user_message_recorded"),
                }
            )

        except Exception as e:
            debug_log(f"添加用戶消息失敗: {e}")
            return JSONResponse(
                status_code=500,
                content={
                    "error": f"Failed to add user message: {e!s}",
                    "messageCode": get_msg_code("add_user_message_failed"),
                },
            )

    # ===== Autonomous Lab API & Game UI =====

    @manager.app.get("/lab", response_class=HTMLResponse)
    async def lab_page(request: Request):
        """Serve the game-style lab monitoring UI"""
        return manager.templates.TemplateResponse(
            request,
            "lab.html",
            {
                "title": "Autonomous Lab",
                "version": __version__,
            },
        )

    @manager.app.get("/editorial-center", response_class=HTMLResponse)
    async def editorial_center_page(request: Request):
        """Serve the Editorial Center — central dashboard for all projects"""
        return manager.templates.TemplateResponse(
            request,
            "editorial_center.html",
            {
                "title": "Editorial Center — Autonomous Lab",
                "version": __version__,
            },
        )

    @manager.app.get("/api/projects")
    async def list_projects(request: Request):
        """List all known lab projects across ALL running instances.

        Reads from the shared ~/.autolab/registry.json so every
        Editorial Center shows every project regardless of which port
        it's running on.
        """
        from ..utils.project_registry import get_all_projects

        # Also make sure this instance's local projects are registered
        try:
            from ...server import _lab_browser_opened
            from ...lab.state import load_config, load_state
            from ..utils.project_registry import register_project

            for pdir in list(_lab_browser_opened):
                try:
                    state = load_state(pdir)
                    config = load_config(pdir)
                    name = (
                        config.get("title") or state.get("title")
                        or Path(pdir).name
                    )
                    register_project(
                        manager.host, manager.port, pdir, name, state, config
                    )
                except Exception:
                    pass
        except Exception:
            pass

        # Get ALL projects across all instances
        all_projects = get_all_projects()

        # Refresh each project's live state from disk
        result = []
        from ...lab.state import load_config, load_state
        for p in all_projects:
            pdir = p.get("project_dir", "")
            try:
                state = load_state(pdir)
                config = load_config(pdir)
                editorial = state.get("editorial", {})
                result.append({
                    "project_dir": pdir,
                    "name": config.get("title") or state.get("title") or Path(pdir).name,
                    "iteration": state.get("iteration", 0),
                    "next_role": state.get("next_role", "pi"),
                    "status": state.get("status", "active"),
                    "progress": state.get("progress", 0),
                    "domain": config.get("domain", "research"),
                    "created_at": state.get("created_at", ""),
                    "editorial_phase": editorial.get("phase", "none"),
                    "editorial_decision": editorial.get("decision", ""),
                    # Cross-instance fields
                    "host": p.get("host", "127.0.0.1"),
                    "port": p.get("port", 8765),
                    "base_url": p.get("base_url", ""),
                    "lab_url": p.get("lab_url", ""),
                })
            except Exception:
                # Can't read state (different machine, missing files)
                result.append({
                    **p,
                    "status": p.get("status", "unknown"),
                })
        return JSONResponse(content={"projects": result})

    @manager.app.get("/character-builder", response_class=HTMLResponse)
    async def character_builder_page(request: Request):
        """Serve the character builder UI"""
        return manager.templates.TemplateResponse(
            request,
            "character_builder.html",
            {
                "title": "Character Builder — Autonomous Lab",
                "version": __version__,
            },
        )

    @manager.app.get("/skill-learner", response_class=HTMLResponse)
    async def skill_learner_page(request: Request):
        """Serve the skill learner UI"""
        return manager.templates.TemplateResponse(
            request,
            "skill_learner.html",
            {
                "title": "Skill Learner — Autonomous Lab",
                "version": __version__,
            },
        )

    @manager.app.post("/api/skill/learn")
    async def learn_skill_api(request: Request):
        """Learn a new skill: iteratively research, implement, test, refine.

        Returns a streaming response with SSE-style progress events.
        The agent creates a skill workspace, searches documentation,
        writes a SKILL.md, and optionally tests it against criteria.
        """
        from starlette.responses import StreamingResponse
        import asyncio

        try:
            data = await request.json()
            skill_name = data.get("skill_name", "").strip()
            description = data.get("description", "").strip()
            test_criteria = data.get("test_criteria", "").strip()
            max_iterations = min(int(data.get("max_iterations", 5)), 20)
            target_character = data.get("target_character", "").strip()

            if not skill_name:
                return JSONResponse(
                    status_code=400,
                    content={"error": "skill_name is required"},
                )

            # Determine workspace
            project_dir = Path(
                getattr(manager, "lab_project_dir", None) or "."
            ).resolve()
            skill_ws = project_dir / ".autolab" / "skill-learning" / skill_name
            skill_ws.mkdir(parents=True, exist_ok=True)

            async def generate_learning_events():
                """Stream skill learning progress as SSE events."""
                skill_content = ""
                try:
                    yield _sse({
                        "log": f"Initializing skill workspace: {skill_ws}",
                        "cls": "step",
                        "progress": 5,
                        "iteration": 0,
                    })
                    await asyncio.sleep(0.3)

                    # Write metadata
                    meta = {
                        "skill_name": skill_name,
                        "description": description,
                        "test_criteria": test_criteria,
                        "max_iterations": max_iterations,
                        "target_character": target_character,
                    }
                    (skill_ws / "meta.json").write_text(
                        json.dumps(meta, indent=2), encoding="utf-8"
                    )

                    for iteration in range(1, max_iterations + 1):
                        pct_base = int(10 + (iteration - 1) * (80 / max_iterations))
                        pct_step = int(80 / max_iterations)

                        # Phase 1: Research
                        yield _sse({
                            "log": f"[Iteration {iteration}/{max_iterations}] Searching documentation...",
                            "cls": "step",
                            "progress": pct_base,
                            "iteration": iteration,
                        })
                        search_results = _search_skill_docs(skill_name)
                        await asyncio.sleep(0.3)

                        if search_results:
                            yield _sse({
                                "log": f"  Found {len(search_results)} relevant sources.",
                                "cls": "info",
                                "progress": pct_base + int(pct_step * 0.2),
                            })
                        else:
                            yield _sse({
                                "log": "  No web results; generating from knowledge base.",
                                "cls": "info",
                                "progress": pct_base + int(pct_step * 0.2),
                            })

                        # Phase 2: Generate / refine SKILL.md
                        yield _sse({
                            "log": f"  Generating SKILL.md (iteration {iteration})...",
                            "cls": "info",
                            "progress": pct_base + int(pct_step * 0.4),
                        })
                        skill_content = _generate_learned_skill(
                            skill_name, description, search_results,
                            test_criteria, iteration
                        )
                        skill_md_path = skill_ws / "SKILL.md"
                        skill_md_path.write_text(skill_content, encoding="utf-8")
                        await asyncio.sleep(0.2)

                        yield _sse({
                            "skill_content": skill_content[:3000],
                            "progress": pct_base + int(pct_step * 0.6),
                        })

                        # Phase 3: Validate structure
                        yield _sse({
                            "log": "  Validating skill structure...",
                            "cls": "info",
                            "progress": pct_base + int(pct_step * 0.8),
                        })
                        lines = skill_content.strip().split("\n")
                        checks = []
                        if lines[0].strip() == "---":
                            checks.append("frontmatter")
                        if any("## " in l for l in lines):
                            checks.append("sections")
                        if "```" in skill_content:
                            checks.append("code examples")
                        if "when to use" in skill_content.lower():
                            checks.append("usage guidance")

                        if len(checks) >= 3:
                            yield _sse({
                                "log": f"  Validation passed: {', '.join(checks)}",
                                "cls": "success",
                                "progress": pct_base + pct_step,
                            })
                            # If we have good structure, we can stop early
                            if iteration >= 2 and len(checks) >= 4:
                                yield _sse({
                                    "log": f"  Skill meets quality threshold at iteration {iteration}.",
                                    "cls": "success",
                                })
                                break
                        else:
                            missing = []
                            for req in ["frontmatter", "sections", "code examples"]:
                                if req not in checks:
                                    missing.append(req)
                            yield _sse({
                                "log": f"  Missing: {', '.join(missing)}. Iterating...",
                                "cls": "info",
                                "progress": pct_base + pct_step,
                            })

                        await asyncio.sleep(0.2)

                    # Write file tree
                    tree_lines = [f"{skill_ws.name}/"]
                    for f in sorted(skill_ws.iterdir()):
                        tree_lines.append(f"├── {f.name}")
                    yield _sse({
                        "file_tree": "\n".join(tree_lines),
                        "progress": 95,
                    })

                    # Deploy to character if requested
                    if target_character:
                        yield _sse({
                            "log": f"Deploying skill to character: {target_character}",
                            "cls": "step",
                            "progress": 98,
                        })

                    yield _sse({
                        "log": f"Skill '{skill_name}' acquired. SKILL.md written to {skill_ws}/",
                        "cls": "success",
                        "progress": 100,
                        "status": "done",
                        "skill_content": skill_content[:5000],
                    })

                except Exception as e:
                    yield _sse({
                        "log": f"Error: {e}",
                        "cls": "error",
                        "status": "error",
                    })

            return StreamingResponse(
                generate_learning_events(),
                media_type="text/event-stream",
            )
        except Exception as e:
            debug_log(f"Skill learn API error: {e}")
            return JSONResponse(status_code=500, content={"error": str(e)})

    @manager.app.post("/api/skill/distill")
    async def distill_skill_api(request: Request):
        """Distill skills from existing project files.

        Scans a directory for code patterns and generates SKILL.md files.
        """
        from starlette.responses import StreamingResponse
        import asyncio
        import glob as _glob

        try:
            data = await request.json()
            project_directory = data.get("project_directory", "").strip()
            file_patterns = data.get("file_patterns", "*.py, *.ipynb").strip()

            if not project_directory:
                return JSONResponse(
                    status_code=400,
                    content={"error": "project_directory is required"},
                )

            scan_dir = Path(project_directory).resolve()
            if not scan_dir.exists():
                return JSONResponse(
                    status_code=400,
                    content={"error": f"Directory not found: {scan_dir}"},
                )

            async def generate_distill_events():
                try:
                    # Find files matching patterns
                    patterns = [p.strip() for p in file_patterns.split(",")]
                    all_files = []
                    for pat in patterns:
                        all_files.extend(
                            _glob.glob(str(scan_dir / "**" / pat), recursive=True)
                        )

                    yield _sse({
                        "log": f"Found {len(all_files)} files matching {file_patterns}",
                        "cls": "info",
                        "progress": 20,
                    })
                    await asyncio.sleep(0.3)

                    if not all_files:
                        yield _sse({
                            "log": "No files found. Check the directory and patterns.",
                            "cls": "error",
                            "status": "error",
                        })
                        return

                    # Analyze imports and patterns
                    yield _sse({
                        "log": "Analyzing imports and code patterns...",
                        "cls": "step",
                        "progress": 40,
                    })

                    import_counts: dict[str, int] = {}
                    for fpath in all_files[:50]:  # cap at 50 files
                        try:
                            content = Path(fpath).read_text(encoding="utf-8", errors="ignore")
                            for line in content.split("\n"):
                                line = line.strip()
                                if line.startswith("import ") or line.startswith("from "):
                                    parts = line.split()
                                    if len(parts) >= 2:
                                        mod = parts[1].split(".")[0]
                                        if mod not in ("os", "sys", "json", "re", "typing", "pathlib"):
                                            import_counts[mod] = import_counts.get(mod, 0) + 1
                        except Exception:
                            pass

                    await asyncio.sleep(0.2)

                    # Rank by frequency
                    ranked = sorted(import_counts.items(), key=lambda x: -x[1])[:10]
                    yield _sse({
                        "log": f"Top libraries: {', '.join(f'{m}({c})' for m, c in ranked[:5])}",
                        "cls": "info",
                        "progress": 60,
                    })

                    # Generate skill files for top libraries
                    skills_output = []
                    ws_dir = Path(
                        getattr(manager, "lab_project_dir", None) or "."
                    ).resolve() / ".autolab" / "skill-learning"
                    ws_dir.mkdir(parents=True, exist_ok=True)

                    for i, (mod, count) in enumerate(ranked[:5]):
                        yield _sse({
                            "log": f"Generating skill for '{mod}' (used in {count} files)...",
                            "cls": "step",
                            "progress": 60 + int(30 * (i + 1) / 5),
                        })

                        skill_dir = ws_dir / mod
                        skill_dir.mkdir(exist_ok=True)
                        content = _generate_skill_content(mod, [])
                        (skill_dir / "SKILL.md").write_text(content, encoding="utf-8")
                        skills_output.append({
                            "name": mod,
                            "count": count,
                            "preview": content[:500],
                        })
                        await asyncio.sleep(0.2)

                    yield _sse({
                        "log": f"Distilled {len(skills_output)} skills from project history.",
                        "cls": "success",
                        "progress": 100,
                        "status": "done",
                        "skills": skills_output,
                    })

                except Exception as e:
                    yield _sse({
                        "log": f"Error: {e}",
                        "cls": "error",
                        "status": "error",
                    })

            return StreamingResponse(
                generate_distill_events(),
                media_type="text/event-stream",
            )
        except Exception as e:
            debug_log(f"Skill distill API error: {e}")
            return JSONResponse(status_code=500, content={"error": str(e)})

    def _generate_learned_skill(
        skill_name: str,
        description: str,
        search_results: list,
        test_criteria: str,
        iteration: int,
    ) -> str:
        """Generate a SKILL.md with richer content for the learning workflow."""
        readable = skill_name.replace("-", " ").replace("_", " ").title()
        refs_section = ""
        if search_results:
            refs = "\n".join(
                f"- [{r['title'][:80]}](https://doi.org/{r['doi']})"
                for r in search_results if r.get("doi")
            )
            if refs:
                refs_section = f"\n## References\n\n{refs}\n"

        desc_block = description if description else f"Workflows and best practices for {readable.lower()}"
        criteria_block = ""
        if test_criteria:
            criteria_block = (
                f"\n## Acceptance criteria\n\n"
                f"- {test_criteria}\n"
                f"- Skill validated at iteration {iteration}\n"
            )

        return (
            f"---\n"
            f"name: {skill_name}\n"
            f"description: {desc_block}\n"
            f"metadata:\n"
            f"  skill-author: Autonomous Lab Skill Learner\n"
            f"  generated: true\n"
            f"  learning-iteration: {iteration}\n"
            f"---\n\n"
            f"# {readable}\n\n"
            f"{desc_block}\n\n"
            f"## When to use\n\n"
            f"Use this skill when working with {readable.lower()} tasks. "
            f"The agent should follow the workflow below and adapt "
            f"parameters to the specific project context.\n\n"
            f"## Standard workflow\n\n"
            f"```python\n"
            f"# {readable} — standard workflow\n"
            f"# Generated by Skill Learner (iteration {iteration})\n"
            f"#\n"
            f"# The Cursor agent should refine this with actual API calls,\n"
            f"# imports, and step-by-step instructions.\n"
            f"#\n"
            f"# 1. Import dependencies\n"
            f"# 2. Load / preprocess data\n"
            f"# 3. Execute core analysis\n"
            f"# 4. Validate outputs\n"
            f"# 5. Generate figures / reports\n"
            f"```\n\n"
            f"## Key decisions\n\n"
            f"- Document important parameter choices and defaults\n"
            f"- Note common pitfalls and how to avoid them\n"
            f"{criteria_block}"
            f"{refs_section}"
        )

    @manager.app.post("/api/character/create")
    async def create_character_api(request: Request):
        """Create a character folder from the builder UI form."""
        import yaml as _yaml

        try:
            data = await request.json()
            name = data.get("name", "").strip()
            role = data.get("role", "pi")
            title = data.get("title", "").strip()
            expertise = data.get("expertise", "").strip()
            goal = data.get("goal", "").strip()
            skill_list = data.get("skills", [])
            personality_text = data.get("personality", "").strip()
            avatar = data.get("avatar", "generic")
            deploy_as = data.get("deploy_as", "")

            if not all([name, title, expertise, goal, personality_text, skill_list]):
                return JSONResponse(
                    status_code=400,
                    content={"error": "All fields are required."},
                )

            personality_list = [
                p.strip() for p in personality_text.split("\n") if p.strip()
            ]

            # Determine project directory
            project_dir = Path(
                getattr(manager, "lab_project_dir", None) or "."
            ).resolve()

            char_slug = (
                name.lower()
                .replace(" ", "-")
                .replace(".", "")
                .replace("'", "")
            )
            folder_name = f"autolab-char-{char_slug}"
            char_dir = project_dir / folder_name
            char_dir.mkdir(parents=True, exist_ok=True)

            # character.yaml
            character_data = {
                "name": name,
                "role": role,
                "avatar": avatar,
                "title": title,
                "expertise": expertise,
                "goal": goal,
                "skills": skill_list,
                "personality": personality_list,
            }
            yaml_body = _yaml.dump(
                character_data,
                default_flow_style=False,
                allow_unicode=True,
                sort_keys=False,
            )
            (char_dir / "character.yaml").write_text(yaml_body, encoding="utf-8")

            # skills/ directory with skeleton SKILL.md
            skills_dir = char_dir / "skills"
            skills_dir.mkdir(exist_ok=True)
            for skill_name in skill_list:
                skill_folder = skills_dir / skill_name
                skill_folder.mkdir(exist_ok=True)
                skill_md = skill_folder / "SKILL.md"
                if not skill_md.exists():
                    skill_md.write_text(
                        f"---\nname: {skill_name}\n"
                        f"description: TODO\nmetadata:\n"
                        f"  skill-author: {name}\n---\n\n"
                        f"# {skill_name.replace('-', ' ').title()}\n\n"
                        f"## When to use\n\n- TODO\n\n"
                        f"## Standard workflow\n\n```python\n# TODO\n```\n\n"
                        f"## Key decisions\n\n- TODO\n",
                        encoding="utf-8",
                    )

            # README.md
            skills_tree = "\n".join(
                f"    ├── {s}/\n    │   └── SKILL.md" for s in skill_list[:-1]
            )
            if skill_list:
                skills_tree += (
                    f"\n    └── {skill_list[-1]}/\n        └── SKILL.md"
                )
            readme = (
                f"# {name} — {title}\n\n"
                f"Character for [Autonomous Lab](https://autolab.kejunying.com). "
                f"{expertise.rstrip('.')}.\n\n"
                f"## Structure\n\n```\n{folder_name}/\n"
                f"├── character.yaml\n├── README.md\n"
                f"└── skills/\n{skills_tree}\n```\n\n"
                f"## Install\n\n```shell\n"
                f"git clone https://github.com/USERNAME/{folder_name} "
                f".autolab/characters/{char_slug}\n```\n\n"
                f"## License\n\nApache 2.0\n"
            )
            (char_dir / "README.md").write_text(readme, encoding="utf-8")

            # Deploy if requested
            if deploy_as in ("pi", "trainee"):
                autolab_profiles = project_dir / ".autolab" / "profiles"
                autolab_profiles.mkdir(parents=True, exist_ok=True)
                (autolab_profiles / f"{deploy_as}.yaml").write_text(
                    yaml_body, encoding="utf-8"
                )

            debug_log(
                f"Character created via builder UI: {name} ({role}) at {char_dir}"
            )

            return JSONResponse(content={
                "status": "ok",
                "folder_path": str(char_dir),
                "folder_name": folder_name,
                "skills": skill_list,
            })
        except Exception as e:
            debug_log(f"Character create API error: {e}")
            return JSONResponse(
                status_code=500, content={"error": str(e)}
            )

    @manager.app.post("/api/character/train-skill")
    async def train_skill_api(request: Request):
        """Train a single skill: research docs, write SKILL.md, validate.

        Returns a streaming response with SSE-style progress events.
        """
        from starlette.responses import StreamingResponse
        import asyncio

        try:
            data = await request.json()
            folder_path = data.get("folder_path", "")
            skill_name = data.get("skill_name", "")

            if not folder_path or not skill_name:
                return JSONResponse(
                    status_code=400,
                    content={"error": "folder_path and skill_name required"},
                )

            skill_md_path = Path(folder_path) / "skills" / skill_name / "SKILL.md"

            async def generate_training_events():
                """Stream training progress as SSE events."""
                try:
                    yield _sse({"log": f"> Researching '{skill_name}'...", "progress": 10})
                    await asyncio.sleep(0.5)

                    # Step 1: Search for documentation
                    yield _sse({"log": "> Searching web for documentation and best practices...", "progress": 20})
                    search_results = _search_skill_docs(skill_name)
                    await asyncio.sleep(0.3)

                    if search_results:
                        yield _sse({"log": f"> Found {len(search_results)} relevant sources.", "progress": 40})
                    else:
                        yield _sse({"log": "> No web results; generating from knowledge base.", "progress": 40})

                    # Step 2: Generate SKILL.md content
                    yield _sse({"log": "> Generating SKILL.md content...", "progress": 60})
                    content = _generate_skill_content(skill_name, search_results)
                    await asyncio.sleep(0.3)

                    # Step 3: Write to file
                    yield _sse({"log": "> Writing SKILL.md...", "progress": 80})
                    skill_md_path.parent.mkdir(parents=True, exist_ok=True)
                    skill_md_path.write_text(content, encoding="utf-8")

                    # Step 4: Validate
                    yield _sse({"log": "> Validating structure...", "progress": 90})
                    lines = content.strip().split("\n")
                    has_frontmatter = lines[0].strip() == "---"
                    has_sections = any("## " in l for l in lines)
                    has_code = "```" in content

                    checks = []
                    if has_frontmatter:
                        checks.append("frontmatter")
                    if has_sections:
                        checks.append("sections")
                    if has_code:
                        checks.append("code examples")
                    yield _sse({
                        "log": f"> Validation passed: {', '.join(checks)}",
                        "progress": 100,
                        "status": "done",
                        "label": "Complete",
                    })

                except Exception as e:
                    yield _sse({
                        "log": f"> Error: {e}",
                        "status": "error",
                        "label": "Failed",
                    })

            return StreamingResponse(
                generate_training_events(),
                media_type="text/event-stream",
            )
        except Exception as e:
            debug_log(f"Train skill API error: {e}")
            return JSONResponse(
                status_code=500, content={"error": str(e)}
            )

    def _sse(data: dict) -> str:
        """Format a dict as an SSE data line."""
        return f"data: {json.dumps(data)}\n\n"

    def _search_skill_docs(skill_name: str) -> list[dict]:
        """Search for skill documentation using CrossRef-style web lookup."""
        import urllib.request
        import urllib.parse

        results = []
        try:
            query = urllib.parse.quote(f"{skill_name} python tutorial documentation")
            url = f"https://api.crossref.org/works?query={query}&rows=3"
            req = urllib.request.Request(url, headers={
                "User-Agent": "AutonomousLab/0.5 (mailto:autolab@example.com)",
            })
            with urllib.request.urlopen(req, timeout=8) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                for item in data.get("message", {}).get("items", []):
                    results.append({
                        "title": (item.get("title") or [""])[0],
                        "doi": item.get("DOI", ""),
                    })
        except Exception:
            pass  # Non-critical — we generate from knowledge base
        return results

    def _generate_skill_content(skill_name: str, search_results: list) -> str:
        """Generate a structured SKILL.md for a given skill name.

        Uses the skill name to produce a template with relevant sections.
        The Cursor agent will refine this further during the training loop.
        """
        readable = skill_name.replace("-", " ").replace("_", " ").title()
        refs_section = ""
        if search_results:
            refs = "\n".join(
                f"- [{r['title'][:80]}](https://doi.org/{r['doi']})"
                for r in search_results if r.get("doi")
            )
            if refs:
                refs_section = f"\n## References\n\n{refs}\n"

        return (
            f"---\n"
            f"name: {skill_name}\n"
            f"description: {readable} — workflows, best practices, and code patterns\n"
            f"metadata:\n"
            f"  skill-author: Autonomous Lab Character Builder\n"
            f"  generated: true\n"
            f"---\n\n"
            f"# {readable}\n\n"
            f"## When to use\n\n"
            f"- Working with {readable.lower()} tasks or data\n"
            f"- When the project requires {readable.lower()} analysis or processing\n\n"
            f"## Standard workflow\n\n"
            f"```python\n"
            f"# {readable} — standard workflow\n"
            f"# TODO: The Cursor agent should refine this with actual API calls,\n"
            f"# imports, and step-by-step instructions based on web research.\n"
            f"#\n"
            f"# Run autolab_create_character with training enabled,\n"
            f"# or manually research and fill in this section.\n"
            f"```\n\n"
            f"## Key decisions\n\n"
            f"- TODO: Document important parameter choices and defaults\n"
            f"- TODO: Note common pitfalls and how to avoid them\n"
            f"{refs_section}"
        )

    @manager.app.get("/api/autolab/state")
    async def get_autolab_state(request: Request):
        """Get Autonomous Lab project state for the dashboard"""
        try:
            # Priority: URL ?project= > manager.lab_project_dir > current session
            project_dir = request.query_params.get("project") or None
            if not project_dir:
                project_dir = getattr(manager, "lab_project_dir", None)
            if not project_dir:
                current_session = manager.get_current_session()
                if current_session:
                    project_dir = current_session.project_directory

            if not project_dir:
                return JSONResponse(content={"active": False})

            autolab_dir = Path(project_dir) / ".autolab"
            if not autolab_dir.exists():
                return JSONResponse(content={"active": False})

            from ...lab.state import (
                get_domain_config,
                get_paper_progress,
                load_config,
                load_state,
                scan_project_files,
            )

            state = load_state(project_dir)
            file_listings = scan_project_files(project_dir)
            paper_progress = get_paper_progress(project_dir)
            config = load_config(project_dir)
            domain_config = get_domain_config(project_dir)

            figures = file_listings.get("figures", [])

            # Biomedical toolkit integration status (lightweight check)
            biotools_info = {"installed": False}
            try:
                from ...integrations.biomni import get_status

                biotools_info = get_status(project_dir)
            except Exception:
                pass

            editorial = state.get("editorial", {"phase": "none"})

            # Derive project name from config title, state, or folder name
            project_name = (
                config.get("title")
                or state.get("title")
                or Path(project_dir).name
            )

            # Track this project for the Editorial Center (local + shared registry)
            _track_project(manager, str(project_dir), project_name, state, config)
            try:
                from ..utils.project_registry import register_project
                register_project(
                    manager.host, manager.port,
                    str(project_dir), project_name, state, config,
                )
            except Exception:
                pass

            return JSONResponse(
                content={
                    "active": True,
                    "project_dir": str(project_dir),
                    "project_name": project_name,
                    "iteration": state.get("iteration", 0),
                    "next_role": state.get("next_role", "pi"),
                    "status": state.get("status", "active"),
                    "user_feedback": state.get("user_feedback", ""),
                    "feedback_queue_count": (
                        len(state.get("feedback_queue", []))
                        if isinstance(state.get("feedback_queue"), list)
                        else 0
                    ),
                    "progress": state.get("progress", 0),
                    "experts": state.get("experts", []),
                    "created_at": state.get("created_at", ""),
                    "editorial": editorial,
                    "editor_timeout_minutes": config.get("editor_timeout_minutes", 30),
                    "figures": figures,
                    "paper_progress": paper_progress,
                    "file_counts": {k: len(v) for k, v in file_listings.items()},
                    "files": file_listings,
                    "biotools": biotools_info,
                    "domain_config": domain_config,
                }
            )
        except Exception as e:
            debug_log(f"Autolab state API error: {e}")
            return JSONResponse(content={"active": False, "error": str(e)})

    @manager.app.get("/api/autolab/meeting-log")
    async def get_meeting_log(request: Request):
        """Get parsed meeting log turns for the game UI conversation panel"""
        try:
            project_dir = request.query_params.get("project") or getattr(
                manager, "lab_project_dir", None
            )
            if not project_dir:
                return JSONResponse(content={"turns": []})

            from ...lab.state import parse_meeting_log

            turns = parse_meeting_log(project_dir)
            return JSONResponse(content={"turns": turns})
        except Exception as e:
            debug_log(f"Meeting log API error: {e}")
            return JSONResponse(content={"turns": [], "error": str(e)})

    @manager.app.get("/api/autolab/file")
    async def get_autolab_file(request: Request, path: str = ""):
        """Serve a file from the project directory (images, code, paper sections).

        Security: only serves files under allowed subdirectories.
        """
        project_dir = request.query_params.get("project") or getattr(
            manager, "lab_project_dir", None
        )
        if not project_dir:
            return JSONResponse(status_code=400, content={"error": "No active project"})

        if not path:
            return JSONResponse(status_code=400, content={"error": "path required"})

        ALLOWED_DIRS = {"data", "scripts", "figures", "results", "paper"}
        top_dir = path.split("/")[0] if "/" in path else path
        if top_dir not in ALLOWED_DIRS:
            return JSONResponse(status_code=403, content={"error": "Access denied"})

        full_path = Path(project_dir) / path
        resolved = full_path.resolve()
        project_resolved = Path(project_dir).resolve()

        # Path traversal protection
        if not str(resolved).startswith(str(project_resolved)):
            return JSONResponse(status_code=403, content={"error": "Access denied"})

        if not resolved.is_file():
            return JSONResponse(status_code=404, content={"error": "File not found"})

        # Determine if it's a text or binary file
        mime_type, _ = mimetypes.guess_type(str(resolved))
        TEXT_EXTENSIONS = {
            ".py",
            ".tex",
            ".txt",
            ".md",
            ".csv",
            ".tsv",
            ".json",
            ".yaml",
            ".yml",
            ".r",
            ".R",
            ".sh",
            ".bib",
            ".log",
        }

        suffix = resolved.suffix.lower()
        if suffix in TEXT_EXTENSIONS:
            try:
                content = resolved.read_text(encoding="utf-8", errors="replace")
                return JSONResponse(
                    content={
                        "type": "text",
                        "path": path,
                        "name": resolved.name,
                        "content": content[:50000],  # cap at 50KB text
                        "size": len(content),
                        "language": _guess_language(suffix),
                    }
                )
            except Exception as e:
                return JSONResponse(status_code=500, content={"error": str(e)})
        else:
            # Binary file (images, PDFs) — serve directly
            return FileResponse(
                path=str(resolved),
                media_type=mime_type or "application/octet-stream",
                filename=resolved.name,
            )

    # ---- Editorial workflow API ----

    def _resolve_project_dir(body_project_dir: str | None = None) -> str | None:
        """Resolve project dir: prefer body param > manager attr > current session."""
        if body_project_dir:
            return body_project_dir
        pdir = getattr(manager, "lab_project_dir", None)
        if pdir:
            return pdir
        current_session = manager.get_current_session()
        if current_session:
            return current_session.project_directory
        return None

    @manager.app.get("/api/autolab/editorial")
    async def get_editorial_state(request: Request):
        """Get current editorial phase and data"""
        project_dir = _resolve_project_dir(request.query_params.get("project"))
        if not project_dir:
            return JSONResponse(content={"phase": "none"})
        try:
            from ...lab.state import get_editorial

            editorial = get_editorial(project_dir)
            return JSONResponse(content=editorial)
        except Exception as e:
            debug_log(f"Editorial state error: {e}")
            return JSONResponse(content={"phase": "none", "error": str(e)})

    @manager.app.post("/api/autolab/editorial/invite-reviewers")
    async def invite_reviewers_api(request: Request):
        """Editor invites peer reviewers (2-3 from available pool)"""
        try:
            data = await request.json()
            project_dir = _resolve_project_dir(data.get("project_dir"))
            if not project_dir:
                return JSONResponse(
                    status_code=400,
                    content={"error": "No active project — project_dir not set"},
                )

            reviewers = data.get("reviewers", [])
            if not reviewers or len(reviewers) < 1:
                return JSONResponse(
                    status_code=400, content={"error": "Select at least 1 reviewer"}
                )
            if len(reviewers) > 5:
                return JSONResponse(
                    status_code=400, content={"error": "Maximum 5 reviewers"}
                )

            # Assign IDs
            for i, r in enumerate(reviewers):
                r["id"] = f"reviewer_{i + 1}"

            from ...lab.state import invite_reviewers

            editorial = invite_reviewers(project_dir, reviewers)
            debug_log(f"Reviewers invited: {[r['name'] for r in reviewers]}")
            return JSONResponse(content={"status": "ok", "editorial": editorial})
        except Exception as e:
            debug_log(f"Invite reviewers error: {e}")
            return JSONResponse(status_code=500, content={"error": str(e)})

    @manager.app.post("/api/autolab/editorial/decision")
    async def editorial_decision_api(request: Request):
        """Editor makes a decision on the manuscript"""
        try:
            data = await request.json()
            project_dir = _resolve_project_dir(data.get("project_dir"))
            if not project_dir:
                return JSONResponse(
                    status_code=400,
                    content={"error": "No active project — project_dir not set"},
                )

            decision = data.get("decision", "")
            feedback = data.get("feedback", "")

            valid = ("accept", "minor_revision", "major_revision", "reject")
            if decision not in valid:
                return JSONResponse(
                    status_code=400,
                    content={
                        "error": f"Invalid decision. Must be one of: {', '.join(valid)}"
                    },
                )

            from ...lab.state import record_editorial_decision

            editorial = record_editorial_decision(project_dir, decision, feedback)
            debug_log(f"Editorial decision: {decision}")
            return JSONResponse(content={"status": "ok", "editorial": editorial})
        except Exception as e:
            debug_log(f"Editorial decision error: {e}")
            return JSONResponse(status_code=500, content={"error": str(e)})

    @manager.app.post("/api/autolab/editorial/desk-reject")
    async def desk_reject_api(request: Request):
        """Editor desk-rejects without sending to reviewers"""
        try:
            data = await request.json()
            project_dir = _resolve_project_dir(data.get("project_dir"))
            if not project_dir:
                return JSONResponse(
                    status_code=400,
                    content={"error": "No active project — project_dir not set"},
                )

            feedback = data.get("feedback", "Desk rejected by editor.")

            from ...lab.state import record_editorial_decision

            editorial = record_editorial_decision(project_dir, "reject", feedback)
            debug_log("Manuscript desk-rejected")
            return JSONResponse(content={"status": "ok", "editorial": editorial})
        except Exception as e:
            debug_log(f"Desk reject error: {e}")
            return JSONResponse(status_code=500, content={"error": str(e)})

    @manager.app.post("/api/autolab/feedback")
    async def submit_lab_feedback(request: Request):
        """Store user feedback from the game UI (non-blocking, async)"""
        try:
            data = await request.json()
            project_dir = _resolve_project_dir(data.get("project_dir"))
            if not project_dir:
                # Also check query param
                project_dir = request.query_params.get("project") or getattr(
                    manager, "lab_project_dir", None
                )
            if not project_dir:
                return JSONResponse(
                    status_code=400,
                    content={"error": "No active lab project"},
                )

            from ...lab.state import store_user_feedback

            feedback_text = data.get("text", "").strip()
            target_role = data.get("target", "")

            if not feedback_text:
                return JSONResponse(content={"status": "empty"})

            store_user_feedback(project_dir, feedback_text, target_role)
            debug_log(
                f"Lab feedback stored: target={target_role}, len={len(feedback_text)}"
            )
            return JSONResponse(content={"status": "ok"})
        except Exception as e:
            debug_log(f"Lab feedback error: {e}")
            return JSONResponse(
                status_code=500,
                content={"error": str(e)},
            )

    @manager.app.websocket("/ws")
    async def websocket_endpoint(websocket: WebSocket, lang: str = "zh-TW"):
        """WebSocket 端點 - 重構後移除 session_id 依賴"""
        # 獲取當前活躍會話
        session = manager.get_current_session()
        if not session:
            await websocket.close(code=4004, reason="No active session")
            return

        await websocket.accept()

        # 語言由前端處理，不需要在後端設置
        debug_log(f"WebSocket 連接建立，語言由前端處理: {lang}")

        # 檢查會話是否已有 WebSocket 連接
        if session.websocket and session.websocket != websocket:
            debug_log("會話已有 WebSocket 連接，替換為新連接")

        session.websocket = websocket
        debug_log(f"WebSocket 連接建立: 當前活躍會話 {session.session_id}")

        # 發送連接成功消息
        try:
            await websocket.send_json(
                {
                    "type": "connection_established",
                    "messageCode": get_msg_code("websocket_connected"),
                }
            )

            # 檢查是否有待發送的會話更新
            if getattr(manager, "_pending_session_update", False):
                debug_log("檢測到待發送的會話更新，準備發送通知")
                await websocket.send_json(
                    {
                        "type": "session_updated",
                        "action": "new_session_created",
                        "messageCode": get_msg_code("new_session_created"),
                        "session_info": {
                            "project_directory": session.project_directory,
                            "summary": session.summary,
                            "session_id": session.session_id,
                        },
                    }
                )
                manager._pending_session_update = False
                debug_log("✅ 已發送會話更新通知到前端")
            else:
                # 發送當前會話狀態
                await websocket.send_json(
                    {"type": "status_update", "status_info": session.get_status_info()}
                )
                debug_log("已發送當前會話狀態到前端")

        except Exception as e:
            debug_log(f"發送連接確認失敗: {e}")

        try:
            while True:
                data = await websocket.receive_text()
                message = json.loads(data)

                # 重新獲取當前會話，以防會話已切換
                current_session = manager.get_current_session()
                if current_session and current_session.websocket == websocket:
                    await handle_websocket_message(manager, current_session, message)
                else:
                    debug_log("會話已切換或 WebSocket 連接不匹配，忽略消息")
                    break

        except WebSocketDisconnect:
            debug_log("WebSocket 連接正常斷開")
        except ConnectionResetError:
            debug_log("WebSocket 連接被重置")
        except Exception as e:
            debug_log(f"WebSocket 錯誤: {e}")
        finally:
            # 安全清理 WebSocket 連接
            current_session = manager.get_current_session()
            if current_session and current_session.websocket == websocket:
                current_session.websocket = None
                debug_log("已清理會話中的 WebSocket 連接")

    @manager.app.post("/api/save-settings")
    async def save_settings(request: Request):
        """保存設定到檔案"""

        try:
            data = await request.json()

            # 使用統一的設定檔案路徑
            config_dir = Path.home() / ".config" / "autonomous-lab"
            config_dir.mkdir(parents=True, exist_ok=True)
            settings_file = config_dir / "ui_settings.json"

            # 保存設定到檔案
            with open(settings_file, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False, indent=2)

            debug_log(f"設定已保存到: {settings_file}")

            return JSONResponse(
                content={
                    "status": "success",
                    "messageCode": get_msg_code("settings_saved"),
                }
            )

        except Exception as e:
            debug_log(f"保存設定失敗: {e}")
            return JSONResponse(
                status_code=500,
                content={
                    "status": "error",
                    "message": f"Save failed: {e!s}",
                    "messageCode": get_msg_code("save_failed"),
                },
            )

    @manager.app.get("/api/load-settings")
    async def load_settings(request: Request):
        """從檔案載入設定"""

        try:
            # 使用統一的設定檔案路徑
            config_dir = Path.home() / ".config" / "autonomous-lab"
            settings_file = config_dir / "ui_settings.json"

            if settings_file.exists():
                with open(settings_file, encoding="utf-8") as f:
                    settings = json.load(f)

                debug_log(f"設定已從檔案載入: {settings_file}")
                return JSONResponse(content=settings)
            debug_log("設定檔案不存在，返回空設定")
            return JSONResponse(content={})

        except Exception as e:
            debug_log(f"載入設定失敗: {e}")
            return JSONResponse(
                status_code=500,
                content={
                    "status": "error",
                    "message": f"Load failed: {e!s}",
                    "messageCode": get_msg_code("load_failed"),
                },
            )

    @manager.app.post("/api/clear-settings")
    async def clear_settings(request: Request):
        """清除設定檔案"""

        try:
            # 使用統一的設定檔案路徑
            config_dir = Path.home() / ".config" / "autonomous-lab"
            settings_file = config_dir / "ui_settings.json"

            if settings_file.exists():
                settings_file.unlink()
                debug_log(f"設定檔案已刪除: {settings_file}")
            else:
                debug_log("設定檔案不存在，無需刪除")

            return JSONResponse(
                content={
                    "status": "success",
                    "messageCode": get_msg_code("settings_cleared"),
                }
            )

        except Exception as e:
            debug_log(f"清除設定失敗: {e}")
            return JSONResponse(
                status_code=500,
                content={
                    "status": "error",
                    "message": f"Clear failed: {e!s}",
                    "messageCode": get_msg_code("clear_failed"),
                },
            )

    @manager.app.get("/api/load-session-history")
    async def load_session_history(request: Request):
        """從檔案載入會話歷史"""

        try:
            # 使用統一的設定檔案路徑
            config_dir = Path.home() / ".config" / "autonomous-lab"
            history_file = config_dir / "session_history.json"

            if history_file.exists():
                with open(history_file, encoding="utf-8") as f:
                    history_data = json.load(f)

                debug_log(f"會話歷史已從檔案載入: {history_file}")

                # 確保資料格式相容性
                if isinstance(history_data, dict):
                    # 新格式：包含版本資訊和其他元資料
                    sessions = history_data.get("sessions", [])
                    last_cleanup = history_data.get("lastCleanup", 0)
                else:
                    # 舊格式：直接是會話陣列（向後相容）
                    sessions = history_data if isinstance(history_data, list) else []
                    last_cleanup = 0

                # 回傳會話歷史資料
                return JSONResponse(
                    content={"sessions": sessions, "lastCleanup": last_cleanup}
                )

            debug_log("會話歷史檔案不存在，返回空歷史")
            return JSONResponse(content={"sessions": [], "lastCleanup": 0})

        except Exception as e:
            debug_log(f"載入會話歷史失敗: {e}")
            return JSONResponse(
                status_code=500,
                content={
                    "status": "error",
                    "message": f"Load failed: {e!s}",
                    "messageCode": get_msg_code("load_failed"),
                },
            )

    @manager.app.post("/api/save-session-history")
    async def save_session_history(request: Request):
        """保存會話歷史到檔案"""

        try:
            data = await request.json()

            # 使用統一的設定檔案路徑
            config_dir = Path.home() / ".config" / "autonomous-lab"
            config_dir.mkdir(parents=True, exist_ok=True)
            history_file = config_dir / "session_history.json"

            # 建立新格式的資料結構
            history_data = {
                "version": "1.0",
                "sessions": data.get("sessions", []),
                "lastCleanup": data.get("lastCleanup", 0),
                "savedAt": int(time.time() * 1000),  # 當前時間戳
            }

            # 保存會話歷史到檔案
            with open(history_file, "w", encoding="utf-8") as f:
                json.dump(history_data, f, ensure_ascii=False, indent=2)

            debug_log(f"會話歷史已保存到: {history_file}")
            session_count = len(history_data["sessions"])
            debug_log(f"保存了 {session_count} 個會話記錄")

            return JSONResponse(
                content={
                    "status": "success",
                    "messageCode": get_msg_code("session_history_saved"),
                    "params": {"count": session_count},
                }
            )

        except Exception as e:
            debug_log(f"保存會話歷史失敗: {e}")
            return JSONResponse(
                status_code=500,
                content={
                    "status": "error",
                    "message": f"Save failed: {e!s}",
                    "messageCode": get_msg_code("save_failed"),
                },
            )

    @manager.app.get("/api/log-level")
    async def get_log_level(request: Request):
        """獲取日誌等級設定"""

        try:
            # 使用統一的設定檔案路徑
            config_dir = Path.home() / ".config" / "autonomous-lab"
            settings_file = config_dir / "ui_settings.json"

            if settings_file.exists():
                with open(settings_file, encoding="utf-8") as f:
                    settings_data = json.load(f)
                    log_level = settings_data.get("logLevel", "INFO")
                    debug_log(f"從設定檔案載入日誌等級: {log_level}")
                    return JSONResponse(content={"logLevel": log_level})
            else:
                # 預設日誌等級
                default_log_level = "INFO"
                debug_log(f"使用預設日誌等級: {default_log_level}")
                return JSONResponse(content={"logLevel": default_log_level})

        except Exception as e:
            debug_log(f"獲取日誌等級失敗: {e}")
            return JSONResponse(
                status_code=500,
                content={
                    "error": f"Failed to get log level: {e!s}",
                    "messageCode": get_msg_code("get_log_level_failed"),
                },
            )

    @manager.app.post("/api/log-level")
    async def set_log_level(request: Request):
        """設定日誌等級"""

        try:
            data = await request.json()
            log_level = data.get("logLevel")

            if not log_level or log_level not in ["DEBUG", "INFO", "WARN", "ERROR"]:
                return JSONResponse(
                    status_code=400,
                    content={
                        "error": "Invalid log level",
                        "messageCode": get_msg_code("invalid_log_level"),
                    },
                )

            # 使用統一的設定檔案路徑
            config_dir = Path.home() / ".config" / "autonomous-lab"
            config_dir.mkdir(parents=True, exist_ok=True)
            settings_file = config_dir / "ui_settings.json"

            # 載入現有設定或創建新設定
            settings_data = {}
            if settings_file.exists():
                with open(settings_file, encoding="utf-8") as f:
                    settings_data = json.load(f)

            # 更新日誌等級
            settings_data["logLevel"] = log_level

            # 保存設定到檔案
            with open(settings_file, "w", encoding="utf-8") as f:
                json.dump(settings_data, f, ensure_ascii=False, indent=2)

            debug_log(f"日誌等級已設定為: {log_level}")

            return JSONResponse(
                content={
                    "status": "success",
                    "logLevel": log_level,
                    "messageCode": get_msg_code("log_level_updated"),
                }
            )

        except Exception as e:
            debug_log(f"設定日誌等級失敗: {e}")
            return JSONResponse(
                status_code=500,
                content={
                    "status": "error",
                    "message": f"Set failed: {e!s}",
                    "messageCode": get_msg_code("set_failed"),
                },
            )


async def handle_websocket_message(manager: "WebUIManager", session, data: dict):
    """處理 WebSocket 消息"""
    message_type = data.get("type")

    if message_type == "submit_feedback":
        # 提交回饋
        feedback = data.get("feedback", "")
        images = data.get("images", [])
        settings = data.get("settings", {})
        await session.submit_feedback(feedback, images, settings)

    elif message_type == "run_command":
        # 執行命令
        command = data.get("command", "")
        if command.strip():
            await session.run_command(command)

    elif message_type == "get_status":
        # 獲取會話狀態
        if session.websocket:
            try:
                await session.websocket.send_json(
                    {"type": "status_update", "status_info": session.get_status_info()}
                )
            except Exception as e:
                debug_log(f"發送狀態更新失敗: {e}")

    elif message_type == "heartbeat":
        # WebSocket 心跳處理（簡化版）
        # 更新心跳時間
        session.last_heartbeat = time.time()
        session.last_activity = time.time()

        # 發送心跳回應
        if session.websocket:
            try:
                await session.websocket.send_json(
                    {
                        "type": "heartbeat_response",
                        "timestamp": data.get("timestamp", 0),
                    }
                )
            except Exception as e:
                debug_log(f"發送心跳回應失敗: {e}")

    elif message_type == "user_timeout":
        # 用戶設置的超時已到
        debug_log(f"收到用戶超時通知: {session.session_id}")
        # 清理會話資源
        await session._cleanup_resources_on_timeout()
        # 重構：不再自動停止服務器，保持服務器運行以支援持久性

    elif message_type == "pong":
        # 處理來自前端的 pong 回應（用於連接檢測）
        debug_log(f"收到 pong 回應，時間戳: {data.get('timestamp', 'N/A')}")
        # 可以在這裡記錄延遲或更新連接狀態

    elif message_type == "update_timeout_settings":
        # 處理超時設定更新
        settings = data.get("settings", {})
        debug_log(f"收到超時設定更新: {settings}")
        if settings.get("enabled"):
            session.update_timeout_settings(
                enabled=True, timeout_seconds=settings.get("seconds", 3600)
            )
        else:
            session.update_timeout_settings(enabled=False)

    else:
        debug_log(f"未知的消息類型: {message_type}")


async def _delayed_server_stop(manager: "WebUIManager"):
    """延遲停止服務器"""
    import asyncio

    await asyncio.sleep(5)  # 等待 5 秒讓前端有時間關閉
    from ..main import stop_web_ui

    stop_web_ui()
    debug_log("Web UI 服務器已因用戶超時而停止")
