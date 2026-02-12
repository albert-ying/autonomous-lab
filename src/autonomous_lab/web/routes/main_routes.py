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
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, PlainTextResponse

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
        ".py": "python", ".tex": "latex", ".bib": "bibtex",
        ".r": "r", ".R": "r", ".sh": "bash", ".json": "json",
        ".yaml": "yaml", ".yml": "yaml", ".csv": "csv",
        ".tsv": "tsv", ".md": "markdown", ".txt": "text",
    }.get(suffix, "text")


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

    @manager.app.get("/api/autolab/state")
    async def get_autolab_state():
        """Get Autonomous Lab project state for the dashboard"""
        try:
            # Use lab_project_dir if set, else fall back to current session
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
                get_paper_progress,
                load_state,
                scan_project_files,
            )

            state = load_state(project_dir)
            file_listings = scan_project_files(project_dir)
            paper_progress = get_paper_progress(project_dir)

            figures = file_listings.get("figures", [])

            # Biomni integration status (lightweight check)
            biomni_info = {"installed": False}
            try:
                from ...integrations.biomni import get_status
                biomni_info = get_status(project_dir)
            except Exception:
                pass

            editorial = state.get("editorial", {"phase": "none"})

            return JSONResponse(content={
                "active": True,
                "iteration": state.get("iteration", 0),
                "next_role": state.get("next_role", "pi"),
                "status": state.get("status", "active"),
                "user_feedback": state.get("user_feedback", ""),
                "progress": state.get("progress", 0),
                "experts": state.get("experts", []),
                "editorial": editorial,
                "figures": figures,
                "paper_progress": paper_progress,
                "file_counts": {k: len(v) for k, v in file_listings.items()},
                "files": file_listings,
                "biomni": biomni_info,
            })
        except Exception as e:
            debug_log(f"Autolab state API error: {e}")
            return JSONResponse(content={"active": False, "error": str(e)})

    @manager.app.get("/api/autolab/meeting-log")
    async def get_meeting_log():
        """Get parsed meeting log turns for the game UI conversation panel"""
        try:
            project_dir = getattr(manager, "lab_project_dir", None)
            if not project_dir:
                return JSONResponse(content={"turns": []})

            from ...lab.state import parse_meeting_log

            turns = parse_meeting_log(project_dir)
            return JSONResponse(content={"turns": turns})
        except Exception as e:
            debug_log(f"Meeting log API error: {e}")
            return JSONResponse(content={"turns": [], "error": str(e)})

    @manager.app.get("/api/autolab/file")
    async def get_autolab_file(path: str = ""):
        """Serve a file from the project directory (images, code, paper sections).

        Security: only serves files under allowed subdirectories.
        """
        project_dir = getattr(manager, "lab_project_dir", None)
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
            ".py", ".tex", ".txt", ".md", ".csv", ".tsv", ".json",
            ".yaml", ".yml", ".r", ".R", ".sh", ".bib", ".log",
        }

        suffix = resolved.suffix.lower()
        if suffix in TEXT_EXTENSIONS:
            try:
                content = resolved.read_text(encoding="utf-8", errors="replace")
                return JSONResponse(content={
                    "type": "text",
                    "path": path,
                    "name": resolved.name,
                    "content": content[:50000],  # cap at 50KB text
                    "size": len(content),
                    "language": _guess_language(suffix),
                })
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

    @manager.app.get("/api/autolab/editorial")
    async def get_editorial_state():
        """Get current editorial phase and data"""
        project_dir = getattr(manager, "lab_project_dir", None)
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
        project_dir = getattr(manager, "lab_project_dir", None)
        if not project_dir:
            return JSONResponse(status_code=400, content={"error": "No active project"})
        try:
            data = await request.json()
            reviewers = data.get("reviewers", [])
            if not reviewers or len(reviewers) < 1:
                return JSONResponse(status_code=400, content={"error": "Select at least 1 reviewer"})
            if len(reviewers) > 5:
                return JSONResponse(status_code=400, content={"error": "Maximum 5 reviewers"})

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
        project_dir = getattr(manager, "lab_project_dir", None)
        if not project_dir:
            return JSONResponse(status_code=400, content={"error": "No active project"})
        try:
            data = await request.json()
            decision = data.get("decision", "")
            feedback = data.get("feedback", "")

            valid = ("accept", "minor_revision", "major_revision", "reject")
            if decision not in valid:
                return JSONResponse(status_code=400, content={
                    "error": f"Invalid decision. Must be one of: {', '.join(valid)}"
                })

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
        project_dir = getattr(manager, "lab_project_dir", None)
        if not project_dir:
            return JSONResponse(status_code=400, content={"error": "No active project"})
        try:
            data = await request.json()
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
            project_dir = getattr(manager, "lab_project_dir", None)
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
            debug_log(f"Lab feedback stored: target={target_role}, len={len(feedback_text)}")
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
