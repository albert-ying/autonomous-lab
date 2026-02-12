#!/usr/bin/env python3
"""
MCP Feedback Enhanced 伺服器主要模組

此模組提供 MCP (Model Context Protocol) 的增強回饋收集功能，
支援智能環境檢測，自動使用 Web UI 介面。

主要功能：
- MCP 工具實現
- 介面選擇（Web UI）
- 環境檢測 (SSH Remote, WSL, Local)
- 國際化支援
- 圖片處理與上傳
- 命令執行與結果展示
- 專案目錄管理

主要 MCP 工具：
- interactive_feedback: 收集用戶互動回饋
- get_system_info: 獲取系統環境資訊

作者: Fábio Ferreira (原作者)
增強: Minidoracat (Web UI, 圖片支援, 環境檢測)
重構: 模塊化設計
"""

import base64
import io
import json
import os
import sys
from typing import Annotated, Any

from fastmcp import FastMCP
from fastmcp.utilities.types import Image as MCPImage
from mcp.types import TextContent
from pydantic import Field

# 導入統一的調試功能
from .debug import server_debug_log as debug_log

# 導入多語系支援
# 導入錯誤處理框架
from .utils.error_handler import ErrorHandler, ErrorType

# 導入資源管理器
from .utils.resource_manager import create_temp_file


# ===== 編碼初始化 =====
def init_encoding():
    """初始化編碼設置，確保正確處理中文字符"""
    try:
        # Windows 特殊處理
        if sys.platform == "win32":
            import msvcrt

            # 設置為二進制模式
            msvcrt.setmode(sys.stdin.fileno(), os.O_BINARY)
            msvcrt.setmode(sys.stdout.fileno(), os.O_BINARY)

            # 重新包裝為 UTF-8 文本流，並禁用緩衝
            # 修復 union-attr 錯誤 - 安全獲取 buffer 或 detach
            stdin_buffer = getattr(sys.stdin, "buffer", None)
            if stdin_buffer is None and hasattr(sys.stdin, "detach"):
                stdin_buffer = sys.stdin.detach()

            stdout_buffer = getattr(sys.stdout, "buffer", None)
            if stdout_buffer is None and hasattr(sys.stdout, "detach"):
                stdout_buffer = sys.stdout.detach()

            sys.stdin = io.TextIOWrapper(
                stdin_buffer, encoding="utf-8", errors="replace", newline=None
            )
            sys.stdout = io.TextIOWrapper(
                stdout_buffer,
                encoding="utf-8",
                errors="replace",
                newline="",
                write_through=True,  # 關鍵：禁用寫入緩衝
            )
        else:
            # 非 Windows 系統的標準設置
            if hasattr(sys.stdout, "reconfigure"):
                sys.stdout.reconfigure(encoding="utf-8", errors="replace")
            if hasattr(sys.stdin, "reconfigure"):
                sys.stdin.reconfigure(encoding="utf-8", errors="replace")

        # 設置 stderr 編碼（用於調試訊息）
        if hasattr(sys.stderr, "reconfigure"):
            sys.stderr.reconfigure(encoding="utf-8", errors="replace")

        return True
    except Exception:
        # 如果編碼設置失敗，嘗試基本設置
        try:
            if hasattr(sys.stdout, "reconfigure"):
                sys.stdout.reconfigure(encoding="utf-8", errors="replace")
            if hasattr(sys.stdin, "reconfigure"):
                sys.stdin.reconfigure(encoding="utf-8", errors="replace")
            if hasattr(sys.stderr, "reconfigure"):
                sys.stderr.reconfigure(encoding="utf-8", errors="replace")
        except:
            pass
        return False


# 初始化編碼（在導入時就執行）
_encoding_initialized = init_encoding()

# ===== 常數定義 =====
SERVER_NAME = "Autonomous Lab MCP"
SSH_ENV_VARS = ["SSH_CONNECTION", "SSH_CLIENT", "SSH_TTY"]
REMOTE_ENV_VARS = ["REMOTE_CONTAINERS", "CODESPACES"]


# 初始化 MCP 服務器
from . import __version__


# 確保 log_level 設定為正確的大寫格式
fastmcp_settings = {}

# 檢查環境變數並設定正確的 log_level
env_log_level = os.getenv("FASTMCP_LOG_LEVEL", "").upper()
if env_log_level in ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"):
    fastmcp_settings["log_level"] = env_log_level
else:
    # 預設使用 INFO 等級
    fastmcp_settings["log_level"] = "INFO"

mcp: Any = FastMCP(SERVER_NAME)


# ===== 工具函數 =====
def is_wsl_environment() -> bool:
    """
    檢測是否在 WSL (Windows Subsystem for Linux) 環境中運行

    Returns:
        bool: True 表示 WSL 環境，False 表示其他環境
    """
    try:
        # 檢查 /proc/version 文件是否包含 WSL 標識
        if os.path.exists("/proc/version"):
            with open("/proc/version") as f:
                version_info = f.read().lower()
                if "microsoft" in version_info or "wsl" in version_info:
                    debug_log("偵測到 WSL 環境（通過 /proc/version）")
                    return True

        # 檢查 WSL 相關環境變數
        wsl_env_vars = ["WSL_DISTRO_NAME", "WSL_INTEROP", "WSLENV"]
        for env_var in wsl_env_vars:
            if os.getenv(env_var):
                debug_log(f"偵測到 WSL 環境變數: {env_var}")
                return True

        # 檢查是否存在 WSL 特有的路徑
        wsl_paths = ["/mnt/c", "/mnt/d", "/proc/sys/fs/binfmt_misc/WSLInterop"]
        for path in wsl_paths:
            if os.path.exists(path):
                debug_log(f"偵測到 WSL 特有路徑: {path}")
                return True

    except Exception as e:
        debug_log(f"WSL 檢測過程中發生錯誤: {e}")

    return False


def is_remote_environment() -> bool:
    """
    檢測是否在遠端環境中運行

    Returns:
        bool: True 表示遠端環境，False 表示本地環境
    """
    # WSL 不應被視為遠端環境，因為它可以訪問 Windows 瀏覽器
    if is_wsl_environment():
        debug_log("WSL 環境不被視為遠端環境")
        return False

    # 檢查 SSH 連線指標
    for env_var in SSH_ENV_VARS:
        if os.getenv(env_var):
            debug_log(f"偵測到 SSH 環境變數: {env_var}")
            return True

    # 檢查遠端開發環境
    for env_var in REMOTE_ENV_VARS:
        if os.getenv(env_var):
            debug_log(f"偵測到遠端開發環境: {env_var}")
            return True

    # 檢查 Docker 容器
    if os.path.exists("/.dockerenv"):
        debug_log("偵測到 Docker 容器環境")
        return True

    # Windows 遠端桌面檢查
    if sys.platform == "win32":
        session_name = os.getenv("SESSIONNAME", "")
        if session_name and "RDP" in session_name:
            debug_log(f"偵測到 Windows 遠端桌面: {session_name}")
            return True

    # Linux 無顯示環境檢查（但排除 WSL）
    if (
        sys.platform.startswith("linux")
        and not os.getenv("DISPLAY")
        and not is_wsl_environment()
    ):
        debug_log("偵測到 Linux 無顯示環境")
        return True

    return False


def save_feedback_to_file(feedback_data: dict, file_path: str | None = None) -> str:
    """
    將回饋資料儲存到 JSON 文件

    Args:
        feedback_data: 回饋資料字典
        file_path: 儲存路徑，若為 None 則自動產生臨時文件

    Returns:
        str: 儲存的文件路徑
    """
    if file_path is None:
        # 使用資源管理器創建臨時文件
        file_path = create_temp_file(suffix=".json", prefix="feedback_")

    # 確保目錄存在
    directory = os.path.dirname(file_path)
    if directory and not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)

    # 複製數據以避免修改原始數據
    json_data = feedback_data.copy()

    # 處理圖片數據：將 bytes 轉換為 base64 字符串以便 JSON 序列化
    if "images" in json_data and isinstance(json_data["images"], list):
        processed_images = []
        for img in json_data["images"]:
            if isinstance(img, dict) and "data" in img:
                processed_img = img.copy()
                # 如果 data 是 bytes，轉換為 base64 字符串
                if isinstance(img["data"], bytes):
                    processed_img["data"] = base64.b64encode(img["data"]).decode(
                        "utf-8"
                    )
                    processed_img["data_type"] = "base64"
                processed_images.append(processed_img)
            else:
                processed_images.append(img)
        json_data["images"] = processed_images

    # 儲存資料
    with open(file_path, "w", encoding="utf-8") as f:
        json.dump(json_data, f, ensure_ascii=False, indent=2)

    debug_log(f"回饋資料已儲存至: {file_path}")
    return file_path


def create_feedback_text(feedback_data: dict) -> str:
    """
    建立格式化的回饋文字

    Args:
        feedback_data: 回饋資料字典

    Returns:
        str: 格式化後的回饋文字
    """
    text_parts = []

    # 基本回饋內容
    if feedback_data.get("interactive_feedback"):
        text_parts.append(f"=== 用戶回饋 ===\n{feedback_data['interactive_feedback']}")

    # 命令執行日誌
    if feedback_data.get("command_logs"):
        text_parts.append(f"=== 命令執行日誌 ===\n{feedback_data['command_logs']}")

    # 圖片附件概要
    if feedback_data.get("images"):
        images = feedback_data["images"]
        text_parts.append(f"=== 圖片附件概要 ===\n用戶提供了 {len(images)} 張圖片：")

        for i, img in enumerate(images, 1):
            size = img.get("size", 0)
            name = img.get("name", "unknown")

            # 智能單位顯示
            if size < 1024:
                size_str = f"{size} B"
            elif size < 1024 * 1024:
                size_kb = size / 1024
                size_str = f"{size_kb:.1f} KB"
            else:
                size_mb = size / (1024 * 1024)
                size_str = f"{size_mb:.1f} MB"

            img_info = f"  {i}. {name} ({size_str})"

            # 為提高兼容性，添加 base64 預覽信息
            if img.get("data"):
                try:
                    if isinstance(img["data"], bytes):
                        img_base64 = base64.b64encode(img["data"]).decode("utf-8")
                    elif isinstance(img["data"], str):
                        img_base64 = img["data"]
                    else:
                        img_base64 = None

                    if img_base64:
                        # 只顯示前50個字符的預覽
                        preview = (
                            img_base64[:50] + "..."
                            if len(img_base64) > 50
                            else img_base64
                        )
                        img_info += f"\n     Base64 預覽: {preview}"
                        img_info += f"\n     完整 Base64 長度: {len(img_base64)} 字符"

                        # 如果 AI 助手不支援 MCP 圖片，可以提供完整 base64
                        debug_log(f"圖片 {i} Base64 已準備，長度: {len(img_base64)}")

                        # 檢查是否啟用 Base64 詳細模式（從 UI 設定中獲取）
                        include_full_base64 = feedback_data.get("settings", {}).get(
                            "enable_base64_detail", False
                        )

                        if include_full_base64:
                            # 根據檔案名推斷 MIME 類型
                            file_name = img.get("name", "image.png")
                            if file_name.lower().endswith((".jpg", ".jpeg")):
                                mime_type = "image/jpeg"
                            elif file_name.lower().endswith(".gif"):
                                mime_type = "image/gif"
                            elif file_name.lower().endswith(".webp"):
                                mime_type = "image/webp"
                            else:
                                mime_type = "image/png"

                            img_info += f"\n     完整 Base64: data:{mime_type};base64,{img_base64}"

                except Exception as e:
                    debug_log(f"圖片 {i} Base64 處理失敗: {e}")

            text_parts.append(img_info)

        # 添加兼容性說明
        text_parts.append(
            "\n💡 注意：如果 AI 助手無法顯示圖片，圖片數據已包含在上述 Base64 信息中。"
        )

    return "\n\n".join(text_parts) if text_parts else "用戶未提供任何回饋內容。"


def process_images(images_data: list[dict]) -> list[MCPImage]:
    """
    處理圖片資料，轉換為 MCP 圖片對象

    Args:
        images_data: 圖片資料列表

    Returns:
        List[MCPImage]: MCP 圖片對象列表
    """
    mcp_images = []

    for i, img in enumerate(images_data, 1):
        try:
            if not img.get("data"):
                debug_log(f"圖片 {i} 沒有資料，跳過")
                continue

            # 檢查數據類型並相應處理
            if isinstance(img["data"], bytes):
                # 如果是原始 bytes 數據，直接使用
                image_bytes = img["data"]
                debug_log(
                    f"圖片 {i} 使用原始 bytes 數據，大小: {len(image_bytes)} bytes"
                )
            elif isinstance(img["data"], str):
                # 如果是 base64 字符串，進行解碼
                image_bytes = base64.b64decode(img["data"])
                debug_log(f"圖片 {i} 從 base64 解碼，大小: {len(image_bytes)} bytes")
            else:
                debug_log(f"圖片 {i} 數據類型不支援: {type(img['data'])}")
                continue

            if len(image_bytes) == 0:
                debug_log(f"圖片 {i} 數據為空，跳過")
                continue

            # 根據文件名推斷格式
            file_name = img.get("name", "image.png")
            if file_name.lower().endswith((".jpg", ".jpeg")):
                image_format = "jpeg"
            elif file_name.lower().endswith(".gif"):
                image_format = "gif"
            else:
                image_format = "png"  # 默認使用 PNG

            # 創建 MCPImage 對象
            mcp_image = MCPImage(data=image_bytes, format=image_format)
            mcp_images.append(mcp_image)

            debug_log(f"圖片 {i} ({file_name}) 處理成功，格式: {image_format}")

        except Exception as e:
            # 使用統一錯誤處理（不影響 JSON RPC）
            error_id = ErrorHandler.log_error_with_context(
                e,
                context={"operation": "圖片處理", "image_index": i},
                error_type=ErrorType.FILE_IO,
            )
            debug_log(f"圖片 {i} 處理失敗 [錯誤ID: {error_id}]: {e}")

    debug_log(f"共處理 {len(mcp_images)} 張圖片")
    return mcp_images


# ===== MCP 工具定義 =====
@mcp.tool()
async def lab_meeting(
    project_directory: Annotated[str, Field(description="Project directory path")] = ".",
    summary: Annotated[
        str, Field(description="Summary of AI work completed for user review")
    ] = "Completed the requested task.",
    timeout: Annotated[int, Field(description="Pause duration in seconds (default: 15)")] = 15,
) -> str:
    """Pause briefly between lab turns to show the monitoring UI and collect optional feedback.

    This opens/refreshes the game-style lab monitoring UI in the browser.
    It pauses for a short period (default 15s) during which the user can
    optionally type feedback through the lab UI's feedback panel.

    Call this after EVERY autolab_record call during normal lab turns.

    LOOP PROTOCOL (normal turns):
      autolab_next -> (act as role) -> autolab_record -> lab_meeting -> autolab_next -> ...

    DO NOT use this for editorial submission waits. Use autolab_editorial instead.

    After this tool returns, ALWAYS call autolab_next to continue.
    NEVER stop the loop.

    Args:
        project_directory: Project directory path
        summary: Summary of AI work completed (displayed in the status bar)
        timeout: How long to pause for user feedback (seconds, default 15)

    Returns:
        str: Any user feedback collected, or a message to continue the loop
    """
    import asyncio

    try:
        if not os.path.exists(project_directory):
            project_directory = os.getcwd()
        project_directory = os.path.abspath(project_directory)

        # Ensure the game-style lab UI is open
        _ensure_lab_ui(project_directory)

        # Clear any stale feedback from the previous turn
        try:
            from .lab.state import load_state, save_state
            state = load_state(project_directory)
            old_feedback = state.get("user_feedback", "")
            if old_feedback:
                state["user_feedback"] = ""
                save_state(project_directory, state)
        except Exception:
            old_feedback = ""

        debug_log(f"lab_meeting: pausing {timeout}s for optional feedback")

        # Poll for user feedback during the pause window
        poll_interval = 3  # seconds
        elapsed = 0
        user_feedback = ""
        while elapsed < timeout:
            await asyncio.sleep(min(poll_interval, timeout - elapsed))
            elapsed += poll_interval
            try:
                state = load_state(project_directory)
                fb = state.get("user_feedback", "").strip()
                if fb:
                    user_feedback = fb
                    # Clear the feedback so the next turn starts fresh
                    state["user_feedback"] = ""
                    save_state(project_directory, state)
                    break  # Got feedback, no need to wait longer
            except Exception:
                continue

        if user_feedback:
            debug_log(f"lab_meeting: got user feedback: {user_feedback[:80]}")
            return (
                f"User feedback received:\n{user_feedback}\n\n"
                f"Incorporate this feedback into the next turn.\n"
                f"Continue the lab loop: call autolab_next now."
            )

        debug_log("lab_meeting: no feedback, continuing")
        return "No user feedback. Continue the lab loop: call autolab_next now."

    except Exception as e:
        debug_log(f"lab_meeting error: {e}")
        return f"Lab meeting error (non-fatal): {e}\n\nContinue the lab loop: call autolab_next now."


@mcp.tool()
def get_system_info() -> str:
    """
    獲取系統環境資訊

    Returns:
        str: JSON 格式的系統資訊
    """
    is_remote = is_remote_environment()
    is_wsl = is_wsl_environment()

    system_info = {
        "平台": sys.platform,
        "Python 版本": sys.version.split()[0],
        "WSL 環境": is_wsl,
        "遠端環境": is_remote,
        "介面類型": "Web UI",
        "環境變數": {
            "SSH_CONNECTION": os.getenv("SSH_CONNECTION"),
            "SSH_CLIENT": os.getenv("SSH_CLIENT"),
            "DISPLAY": os.getenv("DISPLAY"),
            "VSCODE_INJECTION": os.getenv("VSCODE_INJECTION"),
            "SESSIONNAME": os.getenv("SESSIONNAME"),
            "WSL_DISTRO_NAME": os.getenv("WSL_DISTRO_NAME"),
            "WSL_INTEROP": os.getenv("WSL_INTEROP"),
            "WSLENV": os.getenv("WSLENV"),
        },
    }

    return json.dumps(system_info, ensure_ascii=False, indent=2)


# ===== Autonomous Lab MCP Tools =====

# Mandatory loop instruction for NORMAL turns (PI/Trainee/Reviewer)
_LOOP_INSTRUCTION = (
    "\n\n---\n"
    "MANDATORY NEXT STEPS (you MUST do both, in order):\n"
    "1. Call lab_meeting NOW with a brief summary of what just happened.\n"
    "2. After lab_meeting returns, call autolab_next to continue the lab session.\n"
    "NEVER stop or end your turn. The lab loop runs until the user explicitly says 'stop'."
)

# Editorial instruction for submission / waiting-for-editor scenarios
_EDITORIAL_INSTRUCTION = (
    "\n\n---\n"
    "MANDATORY NEXT STEP:\n"
    "Call autolab_editorial NOW. It will wait (with no timeout) for the Editor's action.\n"
    "When it returns, follow its instructions to continue.\n"
    "NEVER stop or end your turn."
)


_lab_browser_opened: set[str] = set()  # tracks which project dirs already have a browser tab


def _ensure_lab_ui(project_directory: str, force_open: bool = False) -> None:
    """Ensure the monitoring web UI is running and points to this project.

    Only opens a new browser tab ONCE per project per MCP session.
    Subsequent calls just update lab_project_dir without spawning tabs.
    Set force_open=True to open the browser even if already opened (e.g., autolab_resume).
    """
    try:
        import urllib.parse

        from .web.main import get_web_ui_manager

        manager = get_web_ui_manager()
        manager.lab_project_dir = project_directory

        # Start the server if not running
        if manager.server_thread is None or not manager.server_thread.is_alive():
            import time
            manager.start_server()
            time.sleep(1)

        # Only open browser once per project (or on force_open)
        if project_directory not in _lab_browser_opened or force_open:
            encoded = urllib.parse.quote(project_directory, safe="")
            lab_url = f"{manager.get_server_url()}/lab?project={encoded}"
            manager.open_browser(lab_url)
            _lab_browser_opened.add(project_directory)
            debug_log(f"Lab UI browser opened for {project_directory}")
        else:
            debug_log("Lab UI already open, skipping browser launch")
    except Exception as e:
        debug_log(f"Lab UI ensure failed (non-fatal): {e}")


@mcp.tool()
async def autolab_init(
    project_directory: Annotated[str, Field(description="Project directory path")] = ".",
    idea: Annotated[str, Field(description="Research project idea and context")] = "",
    title: Annotated[str, Field(description="Paper title")] = "TITLE",
    force_new: Annotated[bool, Field(description="Force fresh start even if project exists")] = False,
) -> str:
    """Initialize an Autonomous Lab project.

    Creates .autolab/ directory with state, profiles, and meeting log.
    Creates paper/ directory with LaTeX template.
    Creates working directories: data/, scripts/, figures/, results/.
    Opens the game-style monitoring UI in the browser.

    If an existing project is found, returns a warning and suggests
    calling autolab_resume instead. Set force_new=True to overwrite.

    Args:
        project_directory: Path to the project directory
        idea: The research idea, context, and any preliminary data description
        title: Working title for the paper
        force_new: If True, overwrite existing project. Default False.

    Returns:
        str: Confirmation message with instructions to call autolab_next
    """
    from .lab.latex_template import create_paper_structure
    from .lab.state import init_project

    project_directory = os.path.abspath(project_directory)

    if not idea.strip():
        return "ERROR: 'idea' parameter is required. Provide the research idea and context."

    # Check if an existing project is present — suggest resume instead of overwriting
    autolab_dir = os.path.join(project_directory, ".autolab")
    state_file = os.path.join(autolab_dir, "state.json")
    if os.path.exists(state_file) and not force_new:
        try:
            with open(state_file, "r") as f:
                existing = json.load(f)
            iteration = existing.get("iteration", 0)
            status = existing.get("status", "active")
            role = existing.get("next_role", "pi")
            progress = existing.get("progress", 0)
            return (
                f"WARNING: An existing Autonomous Lab project was found in {project_directory}\n\n"
                f"  Iteration: {iteration}, Status: {status}, Next role: {role}, Progress: {progress}%\n\n"
                f"To RESUME the existing project, call autolab_resume(project_directory='{project_directory}').\n"
                f"To START FRESH (this will overwrite the existing project), call autolab_init again "
                f"with force_new=True.\n\n"
                f"Do NOT reinitialize unless the user explicitly asks for a fresh start."
            )
        except Exception:
            pass  # Corrupted state — safe to reinitialize

    try:
        state = init_project(project_directory, idea)
        create_paper_structure(project_directory, title)

        # Start the monitoring web UI (non-blocking)
        lab_url = ""
        try:
            _ensure_lab_ui(project_directory)
            from .web.main import get_web_ui_manager
            lab_url = f"{get_web_ui_manager().get_server_url()}/lab"
        except Exception as e:
            debug_log(f"Failed to start lab UI: {e}")
            lab_url = "(failed to start)"

        return (
            f"Autonomous Lab initialized in {project_directory}\n\n"
            f"Created:\n"
            f"  .autolab/ -- state, profiles, meeting log\n"
            f"  paper/ -- LaTeX template ({title})\n"
            f"  scripts/, figures/, results/ -- working directories\n\n"
            f"State: iteration={state['iteration']}, next_role={state['next_role']}\n"
            f"Monitoring UI: {lab_url}"
            + _LOOP_INSTRUCTION
        )
    except Exception as e:
        return f"ERROR initializing project: {e}"


@mcp.tool()
async def autolab_resume(
    project_directory: Annotated[str, Field(description="Project directory path")] = ".",
) -> str:
    """Resume an interrupted Autonomous Lab session.

    Call this when the AI agent was disconnected, the process crashed,
    or a new conversation is started for an existing project.

    This tool:
    1. Loads the full saved state (iteration, role, status, progress, editorial)
    2. Reads the last several meeting log entries to reconstruct context
    3. Lists all existing files (scripts, figures, paper sections)
    4. Opens the monitoring UI
    5. Returns a comprehensive briefing so the AI can continue seamlessly

    After this returns, follow the MANDATORY NEXT STEPS to re-enter the loop.

    Args:
        project_directory: Path to the project directory with .autolab/

    Returns:
        str: Full recovery briefing + instructions to continue
    """
    from .lab.state import (
        get_editorial,
        get_meeting_summaries,
        get_paper_progress,
        get_recent_meetings,
        load_idea,
        load_state,
        scan_project_files,
    )

    project_directory = os.path.abspath(project_directory)

    # Verify the project exists
    autolab_dir = os.path.join(project_directory, ".autolab")
    if not os.path.exists(autolab_dir):
        return (
            f"ERROR: No Autonomous Lab project found in {project_directory}\n"
            f"(.autolab/ directory does not exist)\n\n"
            f"To start a new project, call autolab_init instead."
        )

    try:
        state = load_state(project_directory)
    except Exception as e:
        return f"ERROR loading state: {e}"

    # Open the monitoring UI (force open since this is an explicit resume)
    _ensure_lab_ui(project_directory, force_open=True)

    # Gather context
    idea = load_idea(project_directory)
    iteration = state.get("iteration", 0)
    role = state.get("next_role", "pi")
    status = state.get("status", "active")
    progress = state.get("progress", 0)
    experts = state.get("experts", [])
    editorial = state.get("editorial", {})
    created_at = state.get("created_at", "unknown")

    # Recent meeting history (more than normal — 5 entries for recovery context)
    recent_meetings = get_recent_meetings(project_directory, n=5)
    summaries = get_meeting_summaries(project_directory)
    file_listings = scan_project_files(project_directory)
    paper_progress = get_paper_progress(project_directory)

    # Build file listing summary
    file_summary_parts = []
    for category, files in file_listings.items():
        if files:
            file_summary_parts.append(f"  {category}/: {len(files)} file(s) — {', '.join(files[:5])}")
    file_summary = "\n".join(file_summary_parts) if file_summary_parts else "  (no files yet)"

    # Build paper progress summary
    paper_parts = []
    for section, info in paper_progress.items():
        word_count = info.get("words", 0)
        paper_parts.append(f"  {section}: {word_count} words" if word_count > 0 else f"  {section}: not written")
    paper_summary = "\n".join(paper_parts) if paper_parts else "  (no paper sections yet)"

    # Build editorial status
    editorial_phase = editorial.get("phase", "none")
    editorial_summary = ""
    if editorial_phase != "none":
        editorial_summary = (
            f"\n\nEditorial Status:\n"
            f"  Phase: {editorial_phase}\n"
            f"  Round: {editorial.get('round', 1)}\n"
            f"  Decision: {editorial.get('decision', 'pending')}\n"
        )
        if editorial.get("reviewers"):
            rev_names = [r.get("name", "?") for r in editorial["reviewers"]]
            editorial_summary += f"  Reviewers: {', '.join(rev_names)}\n"
        if editorial.get("reviews"):
            done_reviews = list(editorial["reviews"].keys())
            editorial_summary += f"  Reviews completed: {', '.join(done_reviews)}\n"

    # Expert consultants
    expert_summary = ""
    if experts:
        expert_summary = "\n\nActive Consultants:\n"
        for ex in experts:
            expert_summary += f"  - {ex.get('name', '?')} ({ex.get('role', '?')})\n"

    # Determine the right instruction based on current state
    if status in ("submitted_to_editor", "reviews_complete"):
        next_instruction = _EDITORIAL_INSTRUCTION
    else:
        next_instruction = _LOOP_INSTRUCTION

    briefing = (
        f"[AUTOLAB] SESSION RECOVERED — Resuming from saved state\n"
        f"{'=' * 60}\n\n"
        f"Project: {project_directory}\n"
        f"Started: {created_at}\n"
        f"Iteration: {iteration}\n"
        f"Status: {status}\n"
        f"Progress: {progress}%\n"
        f"Next role: {role}\n"
        f"{editorial_summary}"
        f"{expert_summary}\n\n"
        f"Project Idea:\n{idea[:500]}{'...' if len(idea) > 500 else ''}\n\n"
        f"Files:\n{file_summary}\n\n"
        f"Paper Progress:\n{paper_summary}\n\n"
        f"Recent Meeting History (last 5 turns):\n"
        f"{'─' * 40}\n"
        f"{recent_meetings if recent_meetings.strip() else '(no meeting history yet)'}\n"
        f"{'─' * 40}\n\n"
        f"Meeting Summaries:\n"
        f"{summaries[:1000] if summaries.strip() else '(none yet)'}\n\n"
        f"{'=' * 60}\n"
        f"Recovery complete. The monitoring UI should be open.\n"
        f"Continue the lab from where it left off."
        f"{next_instruction}"
    )

    return briefing


@mcp.tool()
async def autolab_next(
    project_directory: Annotated[str, Field(description="Project directory path")] = ".",
) -> str:
    """Get the next role prompt for the Autonomous Lab session.

    Reads the current state to determine whose turn it is (PI or Trainee),
    assembles context (idea, profiles, meeting history, file listings),
    and returns a detailed role prompt.

    The AI should then follow the returned prompt exactly, acting as
    the specified role (PI or Trainee).

    MANDATORY LOOP PROTOCOL — the lab NEVER stops:

    Normal turns:
      autolab_next -> (act as role) -> autolab_record -> lab_meeting -> autolab_next -> ...

    Editorial submission:
      autolab_record(status='ready_for_review') -> autolab_editorial (blocks) -> follow instructions

    After EVERY autolab_record call, follow the MANDATORY NEXT STEPS in the return message.
    NEVER stop. NEVER end your turn.

    Args:
        project_directory: Path to the project directory

    Returns:
        str: The role prompt to follow (MUST act on it and continue the loop)
    """
    from .lab.prompts import (
        build_pi_prompt,
        build_reviewer_prompt,
        build_revision_prompt,
        build_submission_prompt,
        build_trainee_prompt,
    )
    from .lab.state import (
        get_editorial,
        get_meeting_summaries,
        get_paper_progress,
        get_recent_meetings,
        load_idea,
        load_profile,
        load_state,
        scan_project_files,
    )

    project_directory = os.path.abspath(project_directory)

    # Ensure the monitoring UI is visible on every new turn
    _ensure_lab_ui(project_directory)

    try:
        state = load_state(project_directory)
    except FileNotFoundError as e:
        return f"ERROR: {e}"

    idea = load_idea(project_directory)
    role = state["next_role"]
    iteration = state["iteration"]
    user_feedback = state.get("user_feedback", "")
    meeting_history = get_recent_meetings(project_directory, n=3)
    summaries = get_meeting_summaries(project_directory)
    file_listings = scan_project_files(project_directory)
    editorial = state.get("editorial", {})

    # --- Editorial workflow routing ---

    # 1. PI is preparing submission (ready_for_review)
    if state.get("status") == "ready_for_review" and editorial.get("phase", "none") == "none":
        paper_progress = get_paper_progress(project_directory)
        prompt = build_submission_prompt(paper_progress, file_listings, meeting_history)
        return (
            f"[AUTOLAB] SUBMISSION MODE -- Iteration {iteration}\n\n"
            f"{prompt}\n\n"
            f"After completing your submission, call autolab_record, then follow the MANDATORY NEXT STEPS."
        )

    # 2. Waiting for editor (submitted / reviews_complete)
    if state.get("status") in ("submitted_to_editor", "reviews_complete"):
        phase = editorial.get("phase", "none")
        if phase == "submitted":
            return (
                f"[AUTOLAB] AWAITING EDITOR -- Iteration {iteration}\n\n"
                f"The manuscript has been submitted to the Editor.\n"
                f"The Editor (user) will decide: Desk Reject or Invite Reviewers."
                + _EDITORIAL_INSTRUCTION
            )
        if phase == "reviews_complete":
            return (
                f"[AUTOLAB] AWAITING EDITORIAL DECISION -- Iteration {iteration}\n\n"
                f"All reviewer reports are in.\n"
                f"The Editor (user) will decide: Accept / Minor Revision / Major Revision / Reject."
                + _EDITORIAL_INSTRUCTION
            )
        # Catch-all for transitional editorial states
        return (
            f"[AUTOLAB] EDITORIAL IN PROGRESS -- Iteration {iteration}\n\n"
            f"Editorial phase: {phase}, Status: {state.get('status')}\n"
            f"The editorial workflow is in progress."
            + _EDITORIAL_INSTRUCTION
        )

    # 3. Reviewer turn
    if role.startswith("reviewer_"):
        paper_progress = get_paper_progress(project_directory)
        reviewers = editorial.get("reviewers", [])
        reviewer = next((r for r in reviewers if r["id"] == role), None)
        if not reviewer:
            return f"ERROR: Reviewer {role} not found in editorial state."

        cover_letter = editorial.get("cover_letter", "")
        round_num = editorial.get("round", 1)
        prompt = build_reviewer_prompt(
            reviewer=reviewer,
            paper_progress=paper_progress,
            file_listings=file_listings,
            cover_letter=cover_letter,
            round_number=round_num,
        )
        return (
            f"[AUTOLAB] REVIEWER TURN -- {reviewer['name']} ({reviewer['role']})\n"
            f"Iteration: {iteration}, Round: {round_num}\n\n"
            f"{prompt}\n\n"
            f"After completing your review, call autolab_record, then follow the MANDATORY NEXT STEPS."
        )

    # 4. Paper accepted — PI wraps up
    if state.get("status") == "accepted":
        return (
            f"[AUTOLAB] PAPER ACCEPTED -- Iteration {iteration}\n\n"
            f"Congratulations! The manuscript has been accepted by the Editor.\n"
            f"As PI, write a brief acknowledgment, finalize the paper, and celebrate.\n"
            f"If there is nothing left to do, inform the user the project is complete."
            + _LOOP_INSTRUCTION
        )

    # 5. PI handling revision
    if state.get("status") in ("revision_requested", "rejected"):
        profile = load_profile(project_directory, "pi")
        reviews = editorial.get("reviews", {})
        decision = editorial.get("decision", "")
        decision_fb = editorial.get("decision_feedback", "")
        round_num = editorial.get("round", 1)
        prompt = build_revision_prompt(
            idea=idea,
            profile=profile,
            reviews=reviews,
            editorial_decision=decision,
            editorial_feedback=decision_fb,
            meeting_history=meeting_history,
            file_listings=file_listings,
            round_number=round_num,
        )
        return (
            f"[AUTOLAB] REVISION MODE -- Round {round_num}\n"
            f"Decision: {decision.upper().replace('_', ' ')}\n"
            f"Iteration: {iteration}\n\n"
            f"{prompt}\n\n"
            f"After completing revision actions, call autolab_record, then follow the MANDATORY NEXT STEPS."
        )

    # --- Normal lab workflow ---
    profile = load_profile(project_directory, role)

    if role == "pi":
        prompt = build_pi_prompt(
            idea=idea,
            profile=profile,
            meeting_history=meeting_history,
            summaries=summaries,
            file_listings=file_listings,
            user_feedback=user_feedback,
            iteration=iteration,
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
        )

    return (
        f"[AUTOLAB] {role.upper()} TURN -- Iteration {iteration}\n\n"
        f"{prompt}\n\n"
        f"After completing your actions as {role.upper()}, call autolab_record, then follow the MANDATORY NEXT STEPS."
    )


@mcp.tool()
async def autolab_record(
    project_directory: Annotated[str, Field(description="Project directory path")] = ".",
    role: Annotated[str, Field(description="Role that just acted: 'pi' or 'trainee'")] = "",
    summary: Annotated[str, Field(description="Brief summary of what was done this turn")] = "",
    content: Annotated[str, Field(description="Full content of the turn (review, agenda, results, etc.)")] = "",
    status: Annotated[str, Field(description="Status: 'continue' or 'ready_for_review'")] = "continue",
    progress: Annotated[int, Field(description="PI overall project progress 0-100 (set only by PI)")] = -1,
) -> str:
    """Record a completed turn. Non-blocking — updates state and returns immediately.

    The monitoring UI (opened by autolab_init) auto-refreshes to show progress.
    User feedback is optional and picked up in the next autolab_next call.

    This tool:
    1. Saves the turn to meeting_log.md
    2. Updates state.json (advance iteration, flip role, store status)
    3. Compresses old meetings every 6 turns

    Call this AFTER completing your role actions (as PI or Trainee).

    CRITICAL LOOP PROTOCOL — you MUST follow the MANDATORY NEXT STEPS in the return message:
      - Normal turns: call lab_meeting, then autolab_next.
      - Submission/editorial: call autolab_editorial (blocks until editor acts).
      - NEVER stop. NEVER end your turn. The loop runs indefinitely.

    Args:
        project_directory: Path to the project directory
        role: Which role just completed: 'pi' or 'trainee'
        summary: Brief summary of actions taken
        content: Full structured content of the turn
        status: 'continue' to keep iterating, 'ready_for_review' when paper is done
        progress: PI's assessment of overall project progress 0-100 (only PI should set this)

    Returns:
        str: Instructions to continue the loop (MUST follow them)
    """
    from .lab.state import (
        append_meeting_log,
        compress_old_meetings,
        load_state,
        record_review,
        save_state,
        submit_manuscript,
        COMPRESS_EVERY,
    )

    project_directory = os.path.abspath(project_directory)

    # Validate role — allow pi, trainee, and reviewer_N
    valid_roles = role in ("pi", "trainee") or role.startswith("reviewer_")
    if not valid_roles:
        return "ERROR: role must be 'pi', 'trainee', or 'reviewer_N'"

    try:
        state = load_state(project_directory)
    except FileNotFoundError as e:
        return f"ERROR: {e}"

    iteration = state["iteration"]

    # 1. Append to meeting log
    append_meeting_log(project_directory, role, iteration, summary, content)

    # --- Reviewer turn handling ---
    if role.startswith("reviewer_"):
        # Extract recommendation from content (flexible matching)
        import re
        rec_match = re.search(
            r"(?:RECOMMENDATION|OVERALL\s+RECOMMENDATION|DECISION)[:\s]*([A-Za-z][A-Za-z _]*)",
            content, re.IGNORECASE,
        )
        conf_match = re.search(r"CONFIDENCE[^:]*:\s*(\d)", content, re.IGNORECASE)
        recommendation = rec_match.group(1).strip().lower().replace(" ", "_") if rec_match else "major_revision"
        confidence = int(conf_match.group(1)) if conf_match else 3

        review_data = {
            "recommendation": recommendation,
            "confidence": confidence,
            "report": content,
        }
        editorial = record_review(project_directory, role, review_data)
        remaining = [r["id"] for r in editorial.get("reviewers", [])
                     if r["id"] not in editorial.get("reviews", {})]

        if remaining:
            next_reviewer = remaining[0]
            # Find reviewer name for clarity
            rev_name = next((r.get("name", next_reviewer) for r in editorial.get("reviewers", [])
                           if r["id"] == next_reviewer), next_reviewer)
            return (
                f"Review by {role} recorded successfully.\n"
                f"Next reviewer: {rev_name} ({next_reviewer}).\n\n"
                f"---\n"
                f"MANDATORY NEXT STEPS:\n"
                f"1. Call autolab_next NOW — it will give you {rev_name}'s reviewer prompt.\n"
                f"2. Act as that reviewer and write the review.\n"
                f"3. Call autolab_record with role='{next_reviewer}'.\n"
                f"4. Then repeat: autolab_next for the next reviewer (or autolab_editorial if done).\n"
                f"NEVER stop. NEVER call lab_meeting between reviewer turns — go straight to autolab_next."
            )
        else:
            return (
                f"Review by {role} recorded. ALL REVIEWS ARE COMPLETE.\n\n"
                f"The Editor (user) will now make a final decision."
                + _EDITORIAL_INSTRUCTION
            )

    # --- PI submission handling ---
    if role == "pi" and status == "ready_for_review":
        submit_manuscript(project_directory, content)
        if progress >= 0:
            st = load_state(project_directory)
            st["progress"] = max(0, min(100, progress))
            save_state(project_directory, st)
        return (
            f"Manuscript submitted to Editor.\n"
            f"Cover letter and manuscript summary recorded.\n"
            f"The Editor (user) will now review via the monitoring UI."
            + _EDITORIAL_INSTRUCTION
        )

    # --- Normal turn handling ---
    next_role = "trainee" if role == "pi" else "pi"
    new_iteration = iteration + (1 if role == "trainee" else 0)
    state = load_state(project_directory)  # reload in case reviewer updated it
    state["iteration"] = new_iteration
    state["next_role"] = next_role
    state["status"] = status if status != "ready_for_review" else "active"
    # PI can set overall progress (0-100)
    if progress >= 0 and role == "pi":
        state["progress"] = max(0, min(100, progress))
    save_state(project_directory, state)

    # Compress old meetings periodically
    if new_iteration > 0 and new_iteration % COMPRESS_EVERY == 0:
        compress_old_meetings(project_directory)

    # Check if user left feedback via the web UI
    user_feedback = state.get("user_feedback", "").strip()

    if user_feedback:
        return (
            f"Turn recorded. Iteration {new_iteration}, next: {next_role.upper()}.\n\n"
            f"User feedback pending:\n{user_feedback}"
            + _LOOP_INSTRUCTION
        )

    return (
        f"Turn recorded. Iteration {new_iteration}, next: {next_role.upper()}."
        + _LOOP_INSTRUCTION
    )


@mcp.tool()
async def autolab_consult(
    project_directory: Annotated[str, Field(description="Project directory path")] = ".",
    expert_name: Annotated[str, Field(description="Expert's name (e.g., 'Dr. Sarah Chen')")] = "",
    expert_role: Annotated[str, Field(description="Expert's specialty (e.g., 'Statistician', 'Immunologist')")] = "",
    expert_avatar: Annotated[str, Field(
        description="Avatar key for the UI sprite: reviewer, bioethicist, science_writer, grant_reviewer, "
                    "immunologist, oncologist, neuroscientist, geneticist, cell_biologist, microbiologist, "
                    "pathologist, pharmacologist, structural_bio, systems_biologist, epidemiologist, "
                    "statistician, bioinformatician, data_scientist, ml_engineer, comp_biologist, "
                    "clinician, radiologist, surgeon, chemist, physicist, engineer, psychologist, ecologist, generic"
    )] = "generic",
    question: Annotated[str, Field(description="The specific question or topic to consult about")] = "",
) -> str:
    """Invite a domain expert for a one-time consultation during a PI turn.

    The PI can call this to get advice from a specialist on a specific topic.
    The AI should then role-play as that expert and provide their insight,
    then continue the PI turn normally.

    This does NOT interrupt the loop — it's a synchronous consultation within
    the current PI turn. The expert appears in the monitoring UI sidebar.

    Args:
        project_directory: Path to the project directory
        expert_name: Name of the expert
        expert_role: Their specialty/domain
        expert_avatar: Avatar sprite key for the UI
        question: The specific question being asked

    Returns:
        str: Expert consultation prompt — the AI should role-play as the expert
    """
    from .lab.state import load_state, save_state

    project_directory = os.path.abspath(project_directory)

    if not expert_name or not expert_role or not question:
        return "ERROR: expert_name, expert_role, and question are all required."

    try:
        state = load_state(project_directory)
    except FileNotFoundError as e:
        return f"ERROR: {e}"

    # Add expert to the list (dedup by name)
    experts = state.get("experts", [])
    existing = next((e for e in experts if e["name"] == expert_name), None)
    if not existing:
        experts.append({
            "name": expert_name,
            "role": expert_role,
            "avatar": expert_avatar,
        })
        state["experts"] = experts
        save_state(project_directory, state)

    # Build the expert consultation prompt
    idea = ""
    try:
        from .lab.state import load_idea
        idea = load_idea(project_directory)
    except Exception:
        pass

    return (
        f"[CONSULTATION] Expert: {expert_name} ({expert_role})\n"
        f"{'=' * 50}\n\n"
        f"You are now briefly acting as **{expert_name}**, a specialist in **{expert_role}**.\n"
        f"The PI has asked you the following question:\n\n"
        f"> {question}\n\n"
        f"Project context (brief):\n{idea[:300]}{'...' if len(idea) > 300 else ''}\n\n"
        f"Provide a concise, expert response (2-4 paragraphs). Be direct and specific.\n"
        f"Include:\n"
        f"- Your expert opinion on the question\n"
        f"- Key considerations the PI should be aware of\n"
        f"- Any methodological recommendations\n"
        f"- Potential pitfalls or caveats in your domain\n\n"
        f"After providing the consultation, IMMEDIATELY return to acting as the PI.\n"
        f"Do NOT call autolab_record for this consultation — it's part of the PI's turn.\n"
        f"Continue with whatever the PI was doing before the consultation."
    )


@mcp.tool()
async def autolab_editorial(
    project_directory: Annotated[str, Field(description="Project directory path")] = ".",
) -> str:
    """Wait for the Editor (user) to act on the submitted manuscript.

    Call this in two scenarios:
    1. After PI submits a manuscript (autolab_record with status='ready_for_review').
       Waits for the editor to either Desk Reject or Invite Reviewers.
    2. After all reviewer reports are complete.
       Waits for the editor to make a final decision (Accept/Minor/Major/Reject).

    This tool BLOCKS with NO TIMEOUT — it can wait minutes or hours.
    The editor acts via the monitoring UI (Editor's Desk overlay).

    Returns:
      - If reviewers invited: instructions to run the reviewer loop
        (call autolab_next for each reviewer, then autolab_editorial again)
      - If decision made (accept/minor/major/reject): the decision + instructions
        to call lab_meeting then autolab_next to continue PI/Trainee loop.
    """
    import asyncio

    from .lab.state import load_state

    project_directory = os.path.abspath(project_directory)

    try:
        state = load_state(project_directory)
    except FileNotFoundError as e:
        return f"ERROR: {e}"

    editorial = state.get("editorial", {})
    starting_phase = editorial.get("phase", "none")
    iteration = state.get("iteration", 0)
    current_status = state.get("status", "")
    next_role = state.get("next_role", "")

    debug_log(f"autolab_editorial: starting phase={starting_phase}, status={current_status}, next_role={next_role}")

    # ── Immediate returns: if the state already advanced, don't poll ──

    # Already in reviewer phase → return immediately with reviewer instructions
    if starting_phase in ("reviewers_invited", "under_review") and next_role.startswith("reviewer_"):
        reviewers = editorial.get("reviewers", [])
        reviewer_names = ", ".join(r.get("name", r["id"]) for r in reviewers)
        return (
            f"[EDITORIAL] Reviewers invited: {reviewer_names}\n\n"
            f"The editor has invited {len(reviewers)} reviewer(s) and reviews are in progress.\n"
            f"You must now play each reviewer role in sequence.\n\n"
            f"MANDATORY NEXT STEPS:\n"
            f"1. Call autolab_next — it will return the current reviewer's prompt.\n"
            f"2. Act as that reviewer, then call autolab_record with role='reviewer_N'.\n"
            f"3. Call autolab_next immediately for the next reviewer (do NOT call lab_meeting between reviewers).\n"
            f"4. Repeat until all reviews are done.\n"
            f"5. After the LAST review, call autolab_editorial to wait for the editor's final decision.\n"
            f"NEVER stop. Continue the loop."
        )

    # Already have a decision → return immediately
    if starting_phase == "decision_made":
        decision = editorial.get("decision", "")
        feedback = editorial.get("decision_feedback", "")
        dec_display = decision.upper().replace("_", " ")
        return (
            f"[EDITORIAL] Decision: {dec_display}\n\n"
            f"Editor feedback: {feedback}\n\n"
            f"The editorial process for this round is complete.\n"
            f"Decision: {dec_display}"
            + _LOOP_INSTRUCTION
        )

    # Reviews already complete → wait for decision only
    if starting_phase == "reviews_complete":
        debug_log("autolab_editorial: all reviews done, waiting for editor decision")

    # ── Poll loop: wait for the editorial phase to advance ──

    poll_interval = 5  # seconds
    while True:
        await asyncio.sleep(poll_interval)

        try:
            state = load_state(project_directory)
        except Exception:
            continue  # Retry on transient read errors

        editorial = state.get("editorial", {})
        current_phase = editorial.get("phase", "none")
        current_status = state.get("status", "")
        next_role = state.get("next_role", "")

        # --- Case 1: We were waiting for editor after submission ---
        if starting_phase == "submitted":
            if current_phase in ("reviewers_invited", "under_review"):
                # Editor invited reviewers → AI must now play reviewer roles
                reviewers = editorial.get("reviewers", [])
                reviewer_names = ", ".join(r.get("name", r["id"]) for r in reviewers)
                return (
                    f"[EDITORIAL] Reviewers invited: {reviewer_names}\n\n"
                    f"The editor has invited {len(reviewers)} reviewer(s).\n"
                    f"You must now play each reviewer role in sequence.\n\n"
                    f"MANDATORY NEXT STEPS:\n"
                    f"1. Call autolab_next — it will return the first reviewer's prompt.\n"
                    f"2. Act as that reviewer, then call autolab_record with role='reviewer_N'.\n"
                    f"3. Call autolab_next immediately for the next reviewer (do NOT call lab_meeting between reviewers).\n"
                    f"4. Repeat until all reviews are done.\n"
                    f"5. After the LAST review, call autolab_editorial to wait for the editor's final decision.\n"
                    f"NEVER stop. Continue the loop."
                )
            if current_phase == "decision_made":
                # Desk reject (editor decided without reviewers)
                decision = editorial.get("decision", "reject")
                feedback = editorial.get("decision_feedback", "")
                return (
                    f"[EDITORIAL] Decision: {decision.upper().replace('_', ' ')}\n\n"
                    f"The editor has desk-rejected the manuscript.\n"
                    f"Feedback: {feedback}\n"
                    + _LOOP_INSTRUCTION
                )

        # --- Case 2: We were waiting for final decision after reviews ---
        if starting_phase in ("reviews_complete", "under_review", "reviewers_invited"):
            if current_phase == "decision_made":
                decision = editorial.get("decision", "")
                feedback = editorial.get("decision_feedback", "")
                dec_display = decision.upper().replace("_", " ")
                return (
                    f"[EDITORIAL] Decision: {dec_display}\n\n"
                    f"Editor feedback: {feedback}\n\n"
                    f"The editorial process for this round is complete.\n"
                    f"Decision: {dec_display}"
                    + _LOOP_INSTRUCTION
                )
            # Reviewer role assigned while we're polling → break out to let AI play reviewer
            if next_role.startswith("reviewer_") and current_phase in ("reviewers_invited", "under_review"):
                reviewers = editorial.get("reviewers", [])
                reviewer_names = ", ".join(r.get("name", r["id"]) for r in reviewers)
                return (
                    f"[EDITORIAL] Reviewer turn ready: {next_role}\n\n"
                    f"Reviewers: {reviewer_names}\n\n"
                    f"MANDATORY NEXT STEPS:\n"
                    f"1. Call autolab_next to get the reviewer prompt.\n"
                    f"2. Act as the reviewer, then call autolab_record with role='{next_role}'.\n"
                    f"3. Call autolab_next immediately for the next reviewer (skip lab_meeting between reviewers).\n"
                    f"4. After the LAST review, call autolab_editorial again.\n"
                    f"NEVER stop."
                )

        # --- Generic: any phase change from what we started with ---
        if current_phase != starting_phase:
            return (
                f"[EDITORIAL] Phase changed: {starting_phase} → {current_phase}\n"
                f"Status: {current_status}\n"
                + _LOOP_INSTRUCTION
            )

        # Still waiting — continue polling
        # Gradually increase poll interval up to 15s to reduce disk I/O
        if poll_interval < 15:
            poll_interval = min(poll_interval + 1, 15)


@mcp.tool()
async def autolab_status(
    project_directory: Annotated[str, Field(description="Project directory path")] = ".",
) -> str:
    """Get the current status of an Autonomous Lab project.

    Returns the current state, file listings, paper progress,
    and a summary of recent meetings.

    Args:
        project_directory: Path to the project directory

    Returns:
        str: Formatted status report
    """
    from .lab.state import (
        get_paper_progress,
        get_recent_meetings,
        load_state,
        scan_project_files,
    )

    project_directory = os.path.abspath(project_directory)

    try:
        state = load_state(project_directory)
    except FileNotFoundError as e:
        return f"ERROR: {e}"

    file_listings = scan_project_files(project_directory)
    paper_progress = get_paper_progress(project_directory)
    recent = get_recent_meetings(project_directory, n=2)

    # Format file counts
    file_counts = {k: len(v) for k, v in file_listings.items()}

    # Format paper progress
    paper_lines = []
    for section, info in paper_progress.items():
        if info["exists"] and info["words"] > 0:
            paper_lines.append(f"  {section}: {info['words']} words")
        elif info["exists"]:
            paper_lines.append(f"  {section}: placeholder only")
        else:
            paper_lines.append(f"  {section}: not created")

    report = (
        f"# Autonomous Lab Status\n\n"
        f"**Iteration:** {state['iteration']}\n"
        f"**Next role:** {state['next_role']}\n"
        f"**Status:** {state['status']}\n"
        f"**Last updated:** {state.get('last_updated', 'unknown')}\n\n"
        f"## File Counts\n"
        f"  data: {file_counts.get('data', 0)} files\n"
        f"  scripts: {file_counts.get('scripts', 0)} files\n"
        f"  figures: {file_counts.get('figures', 0)} files\n"
        f"  results: {file_counts.get('results', 0)} files\n"
        f"  paper: {file_counts.get('paper', 0)} files\n\n"
        f"## Paper Progress\n"
        + "\n".join(paper_lines)
        + "\n\n"
        f"## Recent Meetings\n\n{recent if recent.strip() else 'No meetings yet.'}\n"
    )

    if state.get("user_feedback"):
        report += f"\n## Pending User Feedback\n\n{state['user_feedback']}\n"

    return report


# ===== Biomni Integration MCP Tools =====

@mcp.tool()
async def autolab_biomni_status(
    project_directory: Annotated[str, Field(description="Project directory path")] = ".",
) -> str:
    """Check Biomni integration status.

    Returns whether Biomni is installed, its version, whether it is enabled
    for this project, and datalake availability.

    Biomni (https://github.com/snap-stanford/Biomni) is an optional
    integration that provides 100+ biomedical tools and workflows.
    It is NOT bundled with Autonomous Lab — users must opt in.

    Args:
        project_directory: Path to the project directory

    Returns:
        str: JSON status report
    """
    from .integrations.biomni import get_status
    project_directory = os.path.abspath(project_directory)
    status = get_status(project_directory)
    return json.dumps(status, indent=2)


@mcp.tool()
async def autolab_biomni_install(
    from_source: Annotated[bool, Field(description="Install from GitHub (True) or PyPI (False)")] = True,
) -> str:
    """Install Biomni from GitHub or PyPI.

    This downloads and installs the Biomni package at runtime.
    Biomni is Apache 2.0 licensed, but some of its integrated tools
    may have more restrictive licenses — review before commercial use.

    NOTE: The full datalake (~11 GB) is downloaded on first agent use.
    Set skip_datalake=True in config to skip it.

    Args:
        from_source: If True, install latest from GitHub main. If False, use PyPI.

    Returns:
        str: Installation result
    """
    from .integrations.biomni import install_biomni
    result = install_biomni(from_source=from_source)
    if result["success"]:
        return f"SUCCESS: {result['message']}"
    else:
        return f"FAILED: {result['message']}"


@mcp.tool()
async def autolab_biomni_configure(
    project_directory: Annotated[str, Field(description="Project directory path")] = ".",
    enabled: Annotated[bool, Field(description="Enable Biomni integration")] = True,
    data_path: Annotated[str, Field(description="Path to Biomni datalake")] = "./data",
    skip_datalake: Annotated[bool, Field(description="Skip downloading the 11GB datalake")] = False,
) -> str:
    """Configure Biomni integration for this project.

    Saves settings to .autolab/biomni_config.json.
    Biomni must be installed separately (use autolab_biomni_install).

    NOTE: Biomni is used as a tool/database library only — we do NOT
    instantiate the Biomni A1 agent or use their LLM pipeline.

    Args:
        project_directory: Path to the project directory
        enabled: Whether to enable Biomni for this project
        data_path: Path for Biomni's data/datalake
        skip_datalake: If True, skip the ~11GB datalake download

    Returns:
        str: Saved configuration
    """
    from .integrations.biomni import save_biomni_config
    project_directory = os.path.abspath(project_directory)
    cfg = save_biomni_config(
        project_directory,
        enabled=enabled,
        data_path=data_path,
        skip_datalake=skip_datalake,
    )
    return f"Biomni config saved:\n{json.dumps(cfg, indent=2)}"


@mcp.tool()
async def autolab_biomni_tools(
    project_directory: Annotated[str, Field(description="Project directory path")] = ".",
) -> str:
    """List available Biomni tools and databases.

    Biomni provides curated biomedical tools (ADMET prediction,
    scRNA-seq analysis, CRISPR screen planning, etc.) and database
    connectors.  This lists what's available for import in scripts.

    NOTE: We use Biomni as a tool/database library, NOT as an agent.
    The Trainee imports and calls individual functions directly in
    their Python scripts.

    Returns:
        str: Available tools and databases, or installation instructions
    """
    from .integrations.biomni import (
        is_biomni_available,
        list_available_databases,
        list_available_tools,
    )
    if not is_biomni_available():
        return (
            "Biomni is not installed.\n"
            "Install with: autolab_biomni_install()\n"
            "See: https://github.com/snap-stanford/Biomni"
        )
    tools = list_available_tools()
    dbs = list_available_databases()

    lines = ["## Available Biomni Tools\n"]
    if tools:
        for t in tools:
            desc = f" — {t['description']}" if t.get('description') else ""
            lines.append(f"  - **{t['name']}** (`{t['module']}`){desc}")
    else:
        lines.append("  (no tools detected)")

    lines.append("\n## Available Biomni Databases\n")
    if dbs:
        for d in dbs:
            desc = f" — {d['description']}" if d.get('description') else ""
            lines.append(f"  - **{d['name']}** (`{d['module']}`){desc}")
    else:
        lines.append("  (no databases detected)")

    lines.append(
        "\n## Usage\n"
        "Import tools directly in scripts:\n"
        "```python\n"
        "from biomni.tools.<tool_name> import *\n"
        "```\n"
        "Do NOT instantiate the Biomni A1 agent."
    )
    return "\n".join(lines)


@mcp.tool()
async def autolab_biomni_env(
) -> str:
    """Check Biomni environment details.

    Shows which Biomni sub-packages are available, the active conda
    environment, datalake status, and Python path.

    Returns:
        str: Environment report
    """
    from .integrations.biomni import check_environment
    env = check_environment()
    return json.dumps(env, indent=2)


# ===== Main entry point =====
def main():
    """主要入口點，用於套件執行
    收集用戶的互動回饋，支援文字和圖片
    此工具使用 Web UI 介面收集用戶回饋，支援智能環境檢測。

    用戶可以：
    1. 執行命令來驗證結果
    2. 提供文字回饋
    3. 上傳圖片作為回饋
    4. 查看 AI 的工作摘要

    調試模式：
    - 設置環境變數 MCP_DEBUG=true 可啟用詳細調試輸出
    - 生產環境建議關閉調試模式以避免輸出干擾


    """
    # 檢查是否啟用調試模式
    debug_enabled = os.getenv("MCP_DEBUG", "").lower() in ("true", "1", "yes", "on")

    if debug_enabled:
        debug_log("🚀 啟動互動式回饋收集 MCP 服務器")
        debug_log(f"   服務器名稱: {SERVER_NAME}")
        debug_log(f"   版本: {__version__}")
        debug_log(f"   平台: {sys.platform}")
        debug_log(f"   編碼初始化: {'成功' if _encoding_initialized else '失敗'}")
        debug_log(f"   遠端環境: {is_remote_environment()}")
        debug_log(f"   WSL 環境: {is_wsl_environment()}")
        debug_log("   Interface: Web UI")
        debug_log("   等待來自 AI 助手的調用...")
        debug_log("準備啟動 MCP 伺服器...")
        debug_log("調用 mcp.run()...")

    try:
        # 使用正確的 FastMCP API
        mcp.run()
    except KeyboardInterrupt:
        if debug_enabled:
            debug_log("收到中斷信號，正常退出")
        sys.exit(0)
    except Exception as e:
        if debug_enabled:
            debug_log(f"MCP 服務器啟動失敗: {e}")
            import traceback

            debug_log(f"詳細錯誤: {traceback.format_exc()}")
        sys.exit(1)


if __name__ == "__main__":
    main()
