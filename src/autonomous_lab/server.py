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

import asyncio
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
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
    summary: Annotated[
        str, Field(description="Summary of AI work completed for user review")
    ] = "Completed the requested task.",
    timeout: Annotated[
        int, Field(description="Pause duration in seconds (default: 15)")
    ] = 15,
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

        from .lab.state import drain_feedback_queue, has_pending_feedback

        debug_log(f"lab_meeting: pausing {timeout}s for optional feedback")

        # Wait for the pause window, checking periodically for new feedback.
        # Messages accumulate in the queue — we don't clear anything yet.
        poll_interval = 3  # seconds
        elapsed = 0
        while elapsed < timeout:
            await asyncio.sleep(min(poll_interval, timeout - elapsed))
            elapsed += poll_interval
            try:
                if has_pending_feedback(project_directory):
                    # Give user a moment to finish typing follow-up messages
                    await asyncio.sleep(2)
                    break
            except Exception:
                continue

        # Drain ALL accumulated messages at once
        user_feedback = drain_feedback_queue(project_directory)

        if user_feedback:
            debug_log(f"lab_meeting: got user feedback ({len(user_feedback)} chars)")
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


_lab_browser_opened: set[str] = (
    set()
)  # tracks which project dirs already have a browser tab


def _ensure_lab_ui(project_directory: str, force_open: bool = False) -> None:
    """Ensure the monitoring web UI is running and points to this project.

    Only opens a new browser tab ONCE per project per MCP session.
    Subsequent calls just update lab_project_dir without spawning tabs.
    Set force_open=True to open the browser even if already opened (e.g., autolab_resume).

    If the web server thread died (crash, port conflict, etc.) this will
    restart it automatically.
    """
    try:
        import time
        import urllib.parse

        from .web.main import get_web_ui_manager

        manager = get_web_ui_manager()
        manager.lab_project_dir = project_directory

        # Start or restart the server if not running
        server_was_dead = (
            manager.server_thread is None or not manager.server_thread.is_alive()
        )
        if server_was_dead:
            debug_log("Lab UI server not running — (re)starting...")
            manager.start_server()
            # Wait with retries for the server to become healthy
            for attempt in range(5):
                time.sleep(1)
                try:
                    import urllib.request

                    health_url = f"{manager.get_server_url()}/api/autolab/state"
                    urllib.request.urlopen(health_url, timeout=2)
                    debug_log(f"Lab UI server healthy after {attempt + 1}s")
                    break
                except Exception:
                    if attempt == 4:
                        debug_log("Lab UI server failed health check after 5s")

        # Only open browser once per project (or on force_open, or after restart)
        if (
            project_directory not in _lab_browser_opened
            or force_open
            or server_was_dead
        ):
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
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
    idea: Annotated[str, Field(description="Research project idea and context")] = "",
    title: Annotated[str, Field(description="Paper title")] = "TITLE",
    force_new: Annotated[
        bool, Field(description="Force fresh start even if project exists")
    ] = False,
    domain: Annotated[
        str,
        Field(
            description="Project domain: 'research' (default), 'software', 'consulting', 'legal', 'medical', or 'creative'. "
            "Controls role labels (e.g., PI/Trainee/Editor vs Tech Lead/Developer/Code Reviewer) "
            "and artifact names (Paper vs Pull Request vs Report)."
        ),
    ] = "research",
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
        domain: Project domain controlling role labels and artifact names.

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

    # Validate domain
    valid_domains = (
        "research",
        "software",
        "consulting",
        "legal",
        "medical",
        "creative",
        "task",
    )
    if domain not in valid_domains:
        return f"ERROR: Invalid domain '{domain}'. Must be one of: {', '.join(valid_domains)}"

    try:
        state = init_project(
            project_directory, idea, domain=domain,
            config={"title": title},
        )
        if domain != "task":
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

        from .lab.state import get_domain_config

        dcfg = get_domain_config(project_directory)

        if domain == "task":
            dirs_line = "  scripts/, results/, data/ -- working directories"
        else:
            dirs_line = (
                f"  paper/ -- LaTeX template ({title})\n"
                f"  scripts/, figures/, results/ -- working directories"
            )

        return (
            f"Autonomous Lab initialized in {project_directory}\n\n"
            f"Domain: {domain} — roles: {dcfg['senior_label']}, {dcfg['junior_label']}, {dcfg['overseer_label']}\n"
            f"Artifact: {dcfg['artifact']}\n\n"
            f"Created:\n"
            f"  .autolab/ -- state, profiles, meeting log\n"
            f"{dirs_line}\n\n"
            f"State: iteration={state['iteration']}, next_role={state['next_role']}\n"
            f"Monitoring UI: {lab_url}" + _LOOP_INSTRUCTION
        )
    except Exception as e:
        return f"ERROR initializing project: {e}"


@mcp.tool()
async def autolab_resume(
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
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
        get_domain_config,
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

    dcfg = get_domain_config(project_directory)

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
            file_summary_parts.append(
                f"  {category}/: {len(files)} file(s) — {', '.join(files[:5])}"
            )
    file_summary = (
        "\n".join(file_summary_parts) if file_summary_parts else "  (no files yet)"
    )

    # Build paper progress summary
    paper_parts = []
    for section, info in paper_progress.items():
        word_count = info.get("words", 0)
        paper_parts.append(
            f"  {section}: {word_count} words"
            if word_count > 0
            else f"  {section}: not written"
        )
    paper_summary = (
        "\n".join(paper_parts) if paper_parts else "  (no paper sections yet)"
    )

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
    consultant_label = dcfg.get("consultant_label", "Consultant")
    expert_summary = ""
    if experts:
        expert_summary = f"\n\nActive {consultant_label}s:\n"
        for ex in experts:
            expert_summary += f"  - {ex.get('name', '?')} ({ex.get('role', '?')})\n"

    # Map internal role id to display label
    if role == "pi":
        role_display = dcfg["senior_label"]
    elif role == "trainee":
        role_display = dcfg["junior_label"]
    else:
        role_display = role  # reviewer_N etc.

    # Determine the right instruction based on current state
    if status == "completed":
        next_instruction = "\n\n---\nTask COMPLETED. No further action needed."
    elif status in ("submitted_to_editor", "reviews_complete"):
        next_instruction = _EDITORIAL_INSTRUCTION
    else:
        next_instruction = _LOOP_INSTRUCTION

    briefing = (
        f"[AUTOLAB] SESSION RECOVERED — Resuming from saved state\n"
        f"{'=' * 60}\n\n"
        f"Project: {project_directory}\n"
        f"Domain: {dcfg['senior_label']}/{dcfg['junior_label']}/{dcfg['overseer_label']}\n"
        f"Started: {created_at}\n"
        f"Iteration: {iteration}\n"
        f"Status: {status}\n"
        f"Progress: {progress}%\n"
        f"Next role: {role_display}\n"
        f"{editorial_summary}"
        f"{expert_summary}\n\n"
        f"Project Idea:\n{idea[:500]}{'...' if len(idea) > 500 else ''}\n\n"
        f"Files:\n{file_summary}\n\n"
        f"{dcfg['artifact']} Progress:\n{paper_summary}\n\n"
        f"Recent {dcfg['meeting_log_name']} (last 5 turns):\n"
        f"{'─' * 40}\n"
        f"{recent_meetings if recent_meetings.strip() else '(no meeting history yet)'}\n"
        f"{'─' * 40}\n\n"
        f"{dcfg['meeting_log_name']} Summaries:\n"
        f"{summaries[:1000] if summaries.strip() else '(none yet)'}\n\n"
        f"{'=' * 60}\n"
        f"Recovery complete. The monitoring UI should be open.\n"
        f"Continue the lab from where it left off."
        f"{next_instruction}"
    )

    return briefing


@mcp.tool()
async def autolab_next(
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
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
        get_domain_config,
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

    # --- Phase-aware routing for multi-character recruitment ---
    # Phase: recruiting — PI turn 0
    if state.get("phase") == "recruiting" and state.get("iteration", 0) == 0:
        idea = load_idea(project_directory)
        dc = get_domain_config(project_directory)
        from .lab.prompts import build_recruiting_prompt
        return build_recruiting_prompt(idea=idea, domain_config=dc)

    # Phase: recruiting but iteration > 0 — check readiness
    if state.get("phase") == "recruiting":
        from .lab.state import check_recruitment_ready
        if check_recruitment_ready(project_directory):
            state = load_state(project_directory)  # reload after phase change
        else:
            return (
                "## Recruiting Phase — Waiting\n\n"
                "Some required skills are still being learned. "
                "Call `autolab_skill_status` to check progress.\n\n"
                "Once all required skills are certified, the loop will begin automatically."
            )

    # --- Task mode: simplified routing ---
    from .lab.state import is_task_mode, load_config as _load_config_tm

    if is_task_mode(project_directory):
        from .lab.prompts import build_executor_prompt, build_verifier_prompt

        if state.get("status") == "completed":
            return (
                "[AUTOLAB] TASK COMPLETED\n\n"
                "All deliverables have been verified and the task is done.\n"
                "Inform the user the task is complete."
            )

        _tm_config = _load_config_tm(project_directory)
        _tm_max_iter = _tm_config.get("max_iterations", 3)
        role = state["next_role"]
        iteration = state["iteration"]
        dcfg = get_domain_config(project_directory)
        from .lab.state import drain_feedback_queue as _drain_fb_tm
        user_feedback = _drain_fb_tm(project_directory)
        meeting_history = get_recent_meetings(project_directory, n=3)
        summaries = get_meeting_summaries(project_directory)
        file_listings = scan_project_files(project_directory)
        profile = load_profile(project_directory, role)

        if role == "pi":  # Verifier turn
            prompt = build_verifier_prompt(
                idea=load_idea(project_directory),
                profile=profile,
                meeting_history=meeting_history,
                summaries=summaries,
                file_listings=file_listings,
                user_feedback=user_feedback,
                iteration=iteration,
                max_iterations=_tm_max_iter,
            )
            return (
                f"[AUTOLAB] VERIFIER TURN -- Iteration {iteration}\n\n"
                f"{prompt}\n\n"
                f"After completing your verification, call autolab_record, then follow the MANDATORY NEXT STEPS."
            )
        else:  # Executor turn
            # Build skill context if character has certified skills
            skill_context = ""
            from pathlib import Path as _Path
            from .lab.skills import build_skill_context

            char_slug = state.get("next_role", "trainee")
            char_dir = _Path(project_directory) / ".autolab" / "characters" / char_slug
            skills_dir = char_dir if (char_dir / "skills").exists() else _Path(project_directory)

            for char in state.get("recruitment", {}).get("characters", []):
                if char.get("slug") == char_slug or char.get("role") == "trainee":
                    certified = [
                        name for name, info in char.get("skills", {}).items()
                        if info.get("status") == "certified"
                    ]
                    if certified:
                        skill_context = build_skill_context(
                            skills_dir=skills_dir,
                            certified_skills=certified,
                            pi_agenda=meeting_history,
                        )
                    break

            prompt = build_executor_prompt(
                idea=load_idea(project_directory),
                profile=profile,
                meeting_history=meeting_history,
                summaries=summaries,
                file_listings=file_listings,
                user_feedback=user_feedback,
                iteration=iteration,
                max_iterations=_tm_max_iter,
                skill_context=skill_context,
            )
            return (
                f"[AUTOLAB] EXECUTOR TURN -- Iteration {iteration}\n\n"
                f"{prompt}\n\n"
                f"After completing your actions as Executor, call autolab_record, then follow the MANDATORY NEXT STEPS."
            )

    # --- Multi-agent orchestration (opt-in) ---
    # If orchestration: multi is set, delegate to dedicated agent subprocesses.
    # Auto-detects Claude Code native agent teams for richer collaboration.
    # On any failure or timeout, falls through to the single-agent path below.
    _orchestration_fallback_notice = ""
    from .lab.state import load_config as _load_config

    _config = _load_config(project_directory)

    async def _try_orchestrated_turn(
        _proj_dir: str, _state: dict, _cfg: dict
    ) -> str | None:
        """Attempt multi-agent orchestration. Returns result string or None to fall through."""
        _ma_role = _state["next_role"]
        if _ma_role not in ("pi", "trainee"):
            return None

        from .lab.orchestrator import Orchestrator

        _orch = Orchestrator(_proj_dir)

        if _ma_role == "pi":
            if _orch.is_available("pi"):
                _summary = await _orch.run_pi_turn()
                _dcfg = get_domain_config(_proj_dir)
                return (
                    f"[AUTOLAB MULTI-AGENT] {_dcfg['senior_label']} turn completed.\n\n"
                    f"**Summary:** {_summary}"
                )

        elif _ma_role == "trainee":
            from .lab.team_builder import TeamBuilder, is_agent_team_available

            if is_agent_team_available(_cfg):
                _tb = TeamBuilder(_proj_dir)
                _meetings = get_recent_meetings(_proj_dir, n=1)
                _team_prompt = _tb.build_team_prompt(
                    pi_agenda=_meetings,
                    idea=load_idea(_proj_dir),
                    file_listings=scan_project_files(_proj_dir),
                    iteration=_state["iteration"],
                )
                from .lab.event_log import log_event as _log_ev
                _log_ev(
                    _proj_dir,
                    "agent_team_mode",
                    iteration=_state["iteration"],
                )
                return _team_prompt

            elif _orch.is_available("trainee"):
                _summary = await _orch.run_trainee_turn()
                _dcfg = get_domain_config(_proj_dir)
                return (
                    f"[AUTOLAB MULTI-AGENT] {_dcfg['junior_label']} turn completed.\n\n"
                    f"**Summary:** {_summary}"
                )

        return None  # Fall through to single-agent

    if _config.get("orchestration") == "multi":
        _orch_timeout = _config.get("orchestrator_timeout", 900)
        try:
            _orch_result = await asyncio.wait_for(
                _try_orchestrated_turn(project_directory, state, _config),
                timeout=_orch_timeout,
            )
            if _orch_result is not None:
                return _orch_result + _LOOP_INSTRUCTION
            # None means fall through to single-agent
            from .lab.event_log import log_event as _log_ev
            _log_ev(
                project_directory,
                "orchestrator_fallback",
                role=state["next_role"],
                reason="agent_not_available",
            )
            _orchestration_fallback_notice = (
                "[NOTICE] Multi-agent orchestration is configured but the "
                f"{state['next_role']} agent is not available. "
                "Falling back to single-agent mode for this turn.\n\n"
            )
        except asyncio.TimeoutError:
            from .lab.event_log import log_event as _log_ev
            _log_ev(
                project_directory,
                "orchestrator_timeout",
                role=state["next_role"],
                timeout=_orch_timeout,
            )
            _orchestration_fallback_notice = (
                f"[NOTICE] Multi-agent orchestration timed out after {_orch_timeout}s. "
                "Falling back to single-agent mode for this turn.\n\n"
            )
        except Exception as _ma_err:
            from .lab.event_log import log_event as _log_ev
            _log_ev(
                project_directory,
                "orchestrator_error",
                role=state["next_role"],
                error=str(_ma_err),
            )
            _orchestration_fallback_notice = (
                f"[NOTICE] Multi-agent orchestration failed ({_ma_err}). "
                "Falling back to single-agent mode for this turn.\n\n"
            )

    idea = load_idea(project_directory)
    role = state["next_role"]
    iteration = state["iteration"]
    dcfg = get_domain_config(project_directory)
    # Drain accumulated user feedback queue (combines all pending messages)
    from .lab.state import drain_feedback_queue

    user_feedback = drain_feedback_queue(project_directory)
    meeting_history = get_recent_meetings(project_directory, n=3)
    summaries = get_meeting_summaries(project_directory)
    file_listings = scan_project_files(project_directory)
    editorial = state.get("editorial", {})

    # --- Editorial workflow routing ---

    # 1. PI is preparing submission (ready_for_review)
    if (
        state.get("status") == "ready_for_review"
        and editorial.get("phase", "none") == "none"
    ):
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
                f"[AUTOLAB] AWAITING {dcfg['overseer_label'].upper()} -- Iteration {iteration}\n\n"
                f"The {dcfg['artifact'].lower()} has been submitted to the {dcfg['overseer_label']}.\n"
                f"The {dcfg['overseer_label']} (user) will decide: Desk Reject or Invite {dcfg['reviewer_label']}s."
                + _EDITORIAL_INSTRUCTION
            )
        if phase == "reviews_complete":
            return (
                f"[AUTOLAB] AWAITING {dcfg['overseer_label'].upper()} DECISION -- Iteration {iteration}\n\n"
                f"All {dcfg['reviewer_label'].lower()} reports are in.\n"
                f"The {dcfg['overseer_label']} (user) will decide: {' / '.join(dcfg['decisions'])}."
                + _EDITORIAL_INSTRUCTION
            )
        # Catch-all for transitional editorial states
        return (
            f"[AUTOLAB] {dcfg['review_process'].upper()} IN PROGRESS -- Iteration {iteration}\n\n"
            f"Editorial phase: {phase}, Status: {state.get('status')}\n"
            f"The {dcfg['review_process'].lower()} is in progress."
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
            idea=load_idea(project_directory),
        )
        return (
            f"[AUTOLAB] REVIEWER TURN -- {reviewer['name']} ({reviewer['role']})\n"
            f"Iteration: {iteration}, Round: {round_num}\n\n"
            f"{prompt}\n\n"
            f"After completing your review, call autolab_record, then follow the MANDATORY NEXT STEPS."
        )

    # 4. Artifact accepted — senior wraps up
    if state.get("status") == "accepted":
        return (
            f"[AUTOLAB] {dcfg['artifact'].upper()} ACCEPTED -- Iteration {iteration}\n\n"
            f"Congratulations! The {dcfg['artifact'].lower()} has been accepted by the {dcfg['overseer_label']}.\n"
            f"As {dcfg['senior_label']}, write a brief acknowledgment, finalize the {dcfg['artifact'].lower()}, and celebrate.\n"
            f"If there is nothing left to do, inform the user the project is complete."
            + _LOOP_INSTRUCTION
        )

    # 5. Senior handling revision
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
            domain_config=dcfg,
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

    # Map internal role ids to display labels
    role_label = dcfg["senior_label"] if role == "pi" else dcfg["junior_label"]

    if role == "pi":
        from .lab.state import load_config as _load_cfg_sa
        _sa_config = _load_cfg_sa(project_directory)
        _sa_orch_mode = _sa_config.get("orchestration", "single")
        # Build trainee info for multi-agent mode prompt section
        _sa_trainees = None
        if _sa_orch_mode == "multi":
            _sa_agents = _sa_config.get("agents", {})
            if "trainees" in _sa_agents:
                _sa_trainees = [
                    {"name": t.get("name", "Trainee"), "focus": t.get("focus", "")}
                    for t in _sa_agents["trainees"]
                ]

        prompt = build_pi_prompt(
            idea=idea,
            profile=profile,
            meeting_history=meeting_history,
            summaries=summaries,
            file_listings=file_listings,
            user_feedback=user_feedback,
            iteration=iteration,
            domain_config=dcfg,
            orchestration_mode=_sa_orch_mode,
            current_trainees=_sa_trainees,
        )
    else:
        # Build skill context if character has certified skills
        skill_context = ""
        from pathlib import Path as _Path
        from .lab.skills import build_skill_context

        char_slug = state.get("next_role", "trainee")
        char_dir = _Path(project_directory) / ".autolab" / "characters" / char_slug
        skills_dir = char_dir if (char_dir / "skills").exists() else _Path(project_directory)

        for char in state.get("recruitment", {}).get("characters", []):
            if char.get("slug") == char_slug or char.get("role") == "trainee":
                certified = [
                    name for name, info in char.get("skills", {}).items()
                    if info.get("status") == "certified"
                ]
                if certified:
                    skill_context = build_skill_context(
                        skills_dir=skills_dir,
                        certified_skills=certified,
                        pi_agenda=meeting_history,
                    )
                break

        prompt = build_trainee_prompt(
            idea=idea,
            profile=profile,
            meeting_history=meeting_history,
            summaries=summaries,
            file_listings=file_listings,
            user_feedback=user_feedback,
            iteration=iteration,
            domain_config=dcfg,
            skill_context=skill_context,
        )

    return (
        f"{_orchestration_fallback_notice}"
        f"[AUTOLAB] {role_label.upper()} TURN -- Iteration {iteration}\n\n"
        f"{prompt}\n\n"
        f"After completing your actions as {role_label}, call autolab_record, then follow the MANDATORY NEXT STEPS."
    )


@mcp.tool()
async def autolab_record(
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
    role: Annotated[
        str, Field(description="Role that just acted: 'pi' or 'trainee'")
    ] = "",
    summary: Annotated[
        str, Field(description="Brief summary of what was done this turn")
    ] = "",
    content: Annotated[
        str,
        Field(description="Full content of the turn (review, agenda, results, etc.)"),
    ] = "",
    status: Annotated[
        str, Field(description="Status: 'continue' or 'ready_for_review'")
    ] = "continue",
    progress: Annotated[
        int, Field(description="PI overall project progress 0-100 (set only by PI)")
    ] = -1,
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
        get_domain_config,
        load_state,
        record_review,
        save_state,
        submit_manuscript,
        COMPRESS_EVERY,
    )

    project_directory = os.path.abspath(project_directory)

    # Normalize task-mode role aliases
    role_aliases = {"executor": "trainee", "verifier": "pi"}
    role = role_aliases.get(role, role)

    # Validate role — allow pi, trainee, and reviewer_N
    valid_roles = role in ("pi", "trainee") or role.startswith("reviewer_")
    if not valid_roles:
        return "ERROR: role must be 'pi', 'trainee', 'executor', 'verifier', or 'reviewer_N'"

    try:
        state = load_state(project_directory)
    except FileNotFoundError as e:
        return f"ERROR: {e}"

    iteration = state["iteration"]
    dcfg = get_domain_config(project_directory)

    # 1. Append to meeting log
    append_meeting_log(project_directory, role, iteration, summary, content)

    # --- Reviewer turn handling ---
    if role.startswith("reviewer_"):
        # Extract recommendation from content (flexible matching)
        import re

        rec_match = re.search(
            r"(?:RECOMMENDATION|OVERALL\s+RECOMMENDATION|DECISION)[:\s]*([A-Za-z][A-Za-z _]*)",
            content,
            re.IGNORECASE,
        )
        conf_match = re.search(r"CONFIDENCE[^:]*:\s*(\d)", content, re.IGNORECASE)
        recommendation = (
            rec_match.group(1).strip().lower().replace(" ", "_")
            if rec_match
            else "major_revision"
        )
        confidence = int(conf_match.group(1)) if conf_match else 3

        review_data = {
            "recommendation": recommendation,
            "confidence": confidence,
            "report": content,
        }
        editorial = record_review(project_directory, role, review_data)
        remaining = [
            r["id"]
            for r in editorial.get("reviewers", [])
            if r["id"] not in editorial.get("reviews", {})
        ]

        if remaining:
            next_reviewer = remaining[0]
            rev_name = next(
                (
                    r.get("name", next_reviewer)
                    for r in editorial.get("reviewers", [])
                    if r["id"] == next_reviewer
                ),
                next_reviewer,
            )
            return (
                f"{dcfg['reviewer_label']} review by {role} recorded successfully.\n"
                f"Next {dcfg['reviewer_label'].lower()}: {rev_name} ({next_reviewer}).\n\n"
                f"---\n"
                f"MANDATORY NEXT STEPS:\n"
                f"1. Call autolab_next NOW — it will give you {rev_name}'s {dcfg['reviewer_label'].lower()} prompt.\n"
                f"2. Act as that {dcfg['reviewer_label'].lower()} and write the review.\n"
                f"3. Call autolab_record with role='{next_reviewer}'.\n"
                f"4. Then repeat: autolab_next for the next {dcfg['reviewer_label'].lower()} (or autolab_editorial if done).\n"
                f"NEVER stop. NEVER call lab_meeting between {dcfg['reviewer_label'].lower()} turns — go straight to autolab_next."
            )
        else:
            return (
                f"{dcfg['reviewer_label']} review by {role} recorded. ALL REVIEWS ARE COMPLETE.\n\n"
                f"The {dcfg['overseer_label']} (user) will now make a final decision."
                + _EDITORIAL_INSTRUCTION
            )

    # --- Task mode special handling ---
    from .lab.state import is_task_mode as _is_task_mode_rec

    if _is_task_mode_rec(project_directory):
        if status == "ready_for_review":
            return "ERROR: Task mode does not support editorial. Use status='completed' or 'continue'."

        if status == "completed":
            # Verifier says all deliverables PASS — task is done
            state = load_state(project_directory)
            state["status"] = "completed"
            if progress >= 0:
                state["progress"] = max(0, min(100, progress))
            save_state(project_directory, state)
            append_meeting_log(project_directory, role, iteration, summary, content)
            return (
                "Task COMPLETED. All deliverables verified.\n\n"
                "Inform the user the task is done."
            )

        # Normal continue — flip roles, check iteration cap
        next_role = "trainee" if role == "pi" else "pi"
        new_iteration = iteration + (1 if role == "trainee" else 0)
        state = load_state(project_directory)
        from .lab.state import load_config as _load_cfg_task
        _task_cfg = _load_cfg_task(project_directory)
        _max_iter = _task_cfg.get("max_iterations", 3)

        # --- Strict cap enforcement ---
        # After executor turn: if new_iteration >= max, allow ONE final verifier
        if role == "trainee" and new_iteration >= _max_iter:
            state["iteration"] = new_iteration
            state["next_role"] = "pi"
            state["status"] = "active"
            if progress >= 0:
                state["progress"] = max(0, min(100, progress))
            save_state(project_directory, state)
            if new_iteration > 0 and new_iteration % COMPRESS_EVERY == 0:
                compress_old_meetings(project_directory)
            return (
                f"Turn recorded. Iteration {new_iteration} (MAX REACHED).\n"
                f"Final Verifier pass will run next."
                + _LOOP_INSTRUCTION
            )

        # After verifier turn: if at or past cap, force-complete — no more executor turns
        if role == "pi" and iteration >= _max_iter:
            state["status"] = "completed"
            if progress >= 0:
                state["progress"] = max(0, min(100, progress))
            save_state(project_directory, state)
            return (
                f"Turn recorded. Iteration cap ({_max_iter}) reached — task force-completed.\n\n"
                "Inform the user the task is done (iteration limit reached)."
            )

        state["iteration"] = new_iteration
        state["next_role"] = next_role
        state["status"] = "active"
        if progress >= 0 and role == "pi":
            state["progress"] = max(0, min(100, progress))
        save_state(project_directory, state)
        if new_iteration > 0 and new_iteration % COMPRESS_EVERY == 0:
            compress_old_meetings(project_directory)

        next_role_label = dcfg["junior_label"] if role == "pi" else dcfg["senior_label"]
        return (
            f"Turn recorded. Iteration {new_iteration}, next: {next_role_label.upper()}."
            + _LOOP_INSTRUCTION
        )

    # --- Senior submission handling ---
    if role == "pi" and status == "ready_for_review":
        submit_manuscript(project_directory, content)
        if progress >= 0:
            st = load_state(project_directory)
            st["progress"] = max(0, min(100, progress))
            save_state(project_directory, st)
        return (
            f"{dcfg['artifact']} submitted to {dcfg['overseer_label']}.\n"
            f"Cover letter and {dcfg['artifact'].lower()} summary recorded.\n"
            f"The {dcfg['overseer_label']} (user) will now review via the monitoring UI."
            + _EDITORIAL_INSTRUCTION
        )

    # --- Normal turn handling ---
    next_role = "trainee" if role == "pi" else "pi"
    next_role_label = dcfg["junior_label"] if role == "pi" else dcfg["senior_label"]
    new_iteration = iteration + (1 if role == "trainee" else 0)
    state = load_state(project_directory)  # reload in case reviewer updated it
    state["iteration"] = new_iteration
    state["next_role"] = next_role
    state["status"] = status if status != "ready_for_review" else "active"
    # Senior can set overall progress (0-100)
    if progress >= 0 and role == "pi":
        state["progress"] = max(0, min(100, progress))
    save_state(project_directory, state)

    # Compress old meetings periodically
    if new_iteration > 0 and new_iteration % COMPRESS_EVERY == 0:
        compress_old_meetings(project_directory)

    # Check if user left feedback via the web UI (peek, don't drain —
    # the feedback will be properly drained in autolab_next or lab_meeting)
    from .lab.state import has_pending_feedback

    pending_fb = has_pending_feedback(project_directory)

    if pending_fb:
        return (
            f"Turn recorded. Iteration {new_iteration}, next: {next_role_label.upper()}.\n\n"
            f"User feedback is queued — it will be included in the next autolab_next call."
            + _LOOP_INSTRUCTION
        )

    return (
        f"Turn recorded. Iteration {new_iteration}, next: {next_role_label.upper()}."
        + _LOOP_INSTRUCTION
    )


@mcp.tool()
async def autolab_consult(
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
    expert_name: Annotated[
        str, Field(description="Expert's name (e.g., 'Dr. Sarah Chen')")
    ] = "",
    expert_role: Annotated[
        str,
        Field(description="Expert's specialty (e.g., 'Statistician', 'Immunologist')"),
    ] = "",
    expert_avatar: Annotated[
        str,
        Field(
            description="Avatar key for the UI sprite: reviewer, bioethicist, science_writer, grant_reviewer, "
            "immunologist, oncologist, neuroscientist, geneticist, cell_biologist, microbiologist, "
            "pathologist, pharmacologist, structural_bio, systems_biologist, epidemiologist, "
            "statistician, bioinformatician, data_scientist, ml_engineer, comp_biologist, "
            "clinician, radiologist, surgeon, chemist, physicist, engineer, psychologist, ecologist, generic"
        ),
    ] = "generic",
    question: Annotated[
        str, Field(description="The specific question or topic to consult about")
    ] = "",
    character_repo: Annotated[
        str,
        Field(
            description="Optional GitHub repo of a marketplace character to use as the expert "
            "(e.g., 'albert-ying/autolab-char-biostatistician'). The character's skills "
            "and personality will define the expert's domain knowledge. "
            "If not provided, uses built-in domain knowledge."
        ),
    ] = "",
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
    from .lab.state import get_domain_config, load_state, save_state

    project_directory = os.path.abspath(project_directory)

    if not expert_name or not expert_role or not question:
        return "ERROR: expert_name, expert_role, and question are all required."

    try:
        state = load_state(project_directory)
    except FileNotFoundError as e:
        return f"ERROR: {e}"

    dcfg = get_domain_config(project_directory)
    senior = dcfg["senior_label"]
    consultant = dcfg.get("consultant_label", "Consultant")

    # Add expert to the list (dedup by name)
    experts = state.get("experts", [])
    existing = next((e for e in experts if e["name"] == expert_name), None)
    if not existing:
        experts.append(
            {
                "name": expert_name,
                "role": expert_role,
                "avatar": expert_avatar,
            }
        )
        state["experts"] = experts
        save_state(project_directory, state)

    # Build the expert consultation prompt with domain-specific knowledge
    idea = ""
    try:
        from .lab.state import load_idea

        idea = load_idea(project_directory)
    except Exception:
        pass

    # Look up domain-specific resources
    from .lab.prompts import CONSULTANT_DOMAINS, CONSULTANT_GENERIC

    # Load domain knowledge: marketplace character > built-in domains > generic
    domain_knowledge = ""
    character_skills_content = ""

    if character_repo:
        # Try to fetch character.yaml + skills from marketplace repo
        import urllib.request
        import json as _json

        clean_repo = character_repo.strip().strip("/")
        for branch in ("master", "main"):
            try:
                char_url = (
                    f"https://raw.githubusercontent.com/{clean_repo}"
                    f"/{branch}/character.yaml"
                )
                req = urllib.request.Request(char_url, headers={
                    "User-Agent": "AutonomousLab/0.5",
                })
                with urllib.request.urlopen(req, timeout=10) as resp:
                    import yaml as _yaml
                    char_data = _yaml.safe_load(resp.read().decode("utf-8"))

                # Use character's expertise and personality as domain knowledge
                char_expertise = char_data.get("expertise", "")
                char_personality = char_data.get("personality", [])
                char_skills = char_data.get("skills", [])

                # Override expert name/role from character if not provided
                if not expert_name or expert_name == "":
                    expert_name = char_data.get("name", expert_name)
                if not expert_role or expert_role == "":
                    expert_role = char_data.get("title", expert_role)

                domain_knowledge = (
                    f"Expertise: {char_expertise}\n"
                    f"Skills: {', '.join(char_skills)}\n"
                    f"Personality: {'; '.join(char_personality)}\n"
                )

                # Try to fetch SKILL.md files for richer context
                for skill_name in char_skills[:3]:  # Cap at 3 to limit context
                    try:
                        skill_url = (
                            f"https://raw.githubusercontent.com/{clean_repo}"
                            f"/{branch}/skills/{skill_name}/SKILL.md"
                        )
                        req2 = urllib.request.Request(skill_url, headers={
                            "User-Agent": "AutonomousLab/0.5",
                        })
                        with urllib.request.urlopen(req2, timeout=5) as resp2:
                            skill_content = resp2.read().decode("utf-8")
                            character_skills_content += (
                                f"\n### Skill: {skill_name}\n{skill_content[:2000]}\n"
                            )
                    except Exception:
                        pass
                break
            except Exception:
                continue

    if not domain_knowledge:
        # Fall back to built-in domain lookup
        avatar_key = expert_avatar.lower().replace(" ", "_")
        role_key = expert_role.lower().replace(" ", "_")
        domain_knowledge = (
            CONSULTANT_DOMAINS.get(avatar_key)
            or CONSULTANT_DOMAINS.get(role_key)
            or CONSULTANT_GENERIC
        )

    return (
        f"[{consultant.upper()}] Expert: {expert_name} ({expert_role})\n"
        f"{'=' * 50}\n\n"
        f"You are now briefly acting as **{expert_name}**, a specialist in **{expert_role}**.\n\n"
        f"**Your domain knowledge and standards:**\n{domain_knowledge}\n\n"
        + (f"**Your specialized skills:**\n{character_skills_content}\n\n"
           if character_skills_content else "")
        + f"The {senior} has asked you the following question:\n\n"
        f"> {question}\n\n"
        f"Project context (brief):\n{idea[:500]}{'...' if len(idea) > 500 else ''}\n\n"
        f"Provide a concise, expert response (2-4 paragraphs). Be direct and specific.\n"
        f"Ground your advice in the specific frameworks, guidelines, and standards listed above.\n"
        f"Include:\n"
        f"- Your expert opinion on the question, citing relevant standards or guidelines\n"
        f"- Key considerations the {senior} should be aware of\n"
        f"- Specific methodological recommendations with justification\n"
        f"- Potential pitfalls or caveats in your domain\n\n"
        f"After providing the consultation, IMMEDIATELY return to acting as the {senior}.\n"
        f"The {senior} should evaluate the {consultant.lower()}'s advice and decide what to adopt.\n"
        f"Do NOT call autolab_record for this consultation — it's part of the {senior}'s turn.\n"
        f"Continue with whatever the {senior} was doing before the consultation."
    )


@mcp.tool()
async def autolab_recruit(
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
    characters: Annotated[
        str,
        Field(
            description="JSON array of character specs. Each: "
            '{"slug": "bio-t", "role": "trainee", "name": "Dr. X", '
            '"avatar": "bioinformatician", '
            '"skills": [{"name": "scanpy", "required": true}]}'
        ),
    ] = "[]",
) -> str:
    """Recruit characters for the project team.

    Called by PI during the recruiting phase (iteration 0).
    Searches marketplace for matching characters, clones them,
    and initiates skill learning for missing skills.
    """
    import json as _json
    from autonomous_lab.lab.recruitment import execute_recruitment

    project_directory = os.path.abspath(project_directory)

    try:
        specs = _json.loads(characters)
    except _json.JSONDecodeError:
        return "ERROR: Invalid JSON in characters parameter."

    if not isinstance(specs, list) or not specs:
        return "ERROR: characters must be a non-empty JSON array."

    result = await execute_recruitment(project_directory, specs)

    return (
        result + "\n\n"
        "## MANDATORY NEXT STEPS\n\n"
        "1. Call `autolab_record` with role='pi', status='continue'\n"
        "2. Then call `lab_meeting`\n"
        "3. Then call `autolab_next` to continue\n"
    )


@mcp.tool()
async def autolab_skill_status(
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
) -> str:
    """Return current skill certification status for all characters."""
    from .lab.state import load_state

    project_directory = os.path.abspath(project_directory)

    try:
        state = load_state(project_directory)
    except FileNotFoundError as e:
        return f"ERROR: {e}"

    chars = state.get("recruitment", {}).get("characters", [])

    if not chars:
        return "No characters recruited yet."

    lines = ["# Skill Status\n"]
    for char in chars:
        slug = char.get("slug", "?")
        lines.append(f"## {slug} ({'READY' if char.get('ready') else 'NOT READY'})")
        for skill_name, info in char.get("skills", {}).items():
            status = info.get("status", "unknown")
            req = " (required)" if info.get("required") else ""
            lines.append(f"  - {skill_name}: **{status}**{req}")
        lines.append("")

    return "\n".join(lines)


@mcp.tool()
async def autolab_cite(
    action: Annotated[
        str,
        Field(
            description="Action: 'search' (find papers by topic), 'doi' (get BibTeX from DOI), "
            "'validate' (check references.bib for errors)"
        ),
    ] = "search",
    query: Annotated[
        str,
        Field(
            description="For 'search': topic description. For 'doi': the DOI string. "
            "For 'validate': path to .bib file (default: paper/references.bib)"
        ),
    ] = "",
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
    count: Annotated[
        int, Field(description="Number of results for search (max 10)")
    ] = 5,
    from_year: Annotated[
        int, Field(description="Only papers from this year onward (0 = no filter)")
    ] = 0,
) -> str:
    """Look up, search, and validate citations using CrossRef API.

    This tool ensures that ALL citations in the paper are real, verified papers
    with correct metadata. Use it for:

    1. **search**: Find real papers on a topic → returns verified BibTeX entries
       ready to paste into references.bib.
    2. **doi**: Convert a known DOI to a properly formatted BibTeX entry.
    3. **validate**: Check an existing .bib file — verify DOIs resolve, find
       entries missing DOIs, and suggest corrections.

    CRITICAL RULE FOR CITATIONS:
    - NEVER invent or fabricate citation entries. ALWAYS use this tool to get
      verified BibTeX from CrossRef.
    - EVERY \cite{key} in the paper MUST have a corresponding entry in
      references.bib that was retrieved via this tool.
    - When writing a paper section, first search for relevant papers, add
      their BibTeX to references.bib, THEN cite them in the text.

    Args:
        action: 'search', 'doi', or 'validate'
        query: Topic description, DOI string, or .bib file path
        project_directory: Project directory path
        count: Number of search results (max 10, default 5)
        from_year: Filter papers from this year onward (0 = no filter)

    Returns:
        BibTeX entries (search/doi) or validation report (validate)
    """
    from .lab.citations import doi_to_bibtex, search_papers, validate_bibtex_file

    project_directory = os.path.abspath(project_directory)

    try:
        if action == "doi":
            if not query:
                return "ERROR: 'query' must contain a DOI string (e.g., '10.1038/s41586-021-03819-2')"
            bibtex = doi_to_bibtex(query)
            return (
                f"Verified BibTeX entry from CrossRef:\n\n{bibtex}\n\n"
                f"Add this to paper/references.bib, then use the citation key in your LaTeX."
            )

        elif action == "search":
            if not query:
                return "ERROR: 'query' must describe what you need to cite (e.g., 'CRISPR efficiency in human cells')"
            year_filter = from_year if from_year > 0 else None
            results = search_papers(
                query, rows=min(count, 10), filter_from_year=year_filter
            )
            if not results or (len(results) == 1 and "error" in results[0]):
                return f"No papers found for: {query}\nTry broader search terms or remove the year filter."

            output_parts = [f"Found {len(results)} verified paper(s) for: {query}\n"]
            for i, r in enumerate(results, 1):
                if "error" in r:
                    continue
                output_parts.append(
                    f"--- Paper {i} ---\n"
                    f"Title: {r['title']}\n"
                    f"Authors: {r['authors']}\n"
                    f"Year: {r['year']} | Journal: {r['journal']}\n"
                    f"DOI: {r['doi']}\n\n"
                    f"{r['bibtex']}\n"
                )
            output_parts.append(
                "\nAdd the relevant entries to paper/references.bib.\n"
                "Use the citation keys (e.g., \\cite{Smith2024keyword}) in your LaTeX text."
            )
            return "\n".join(output_parts)

        elif action == "validate":
            bib_path = (
                query
                if query
                else os.path.join(project_directory, "paper", "references.bib")
            )
            if not os.path.exists(bib_path):
                return f"ERROR: BibTeX file not found: {bib_path}"
            report = validate_bibtex_file(bib_path)
            lines = [
                f"Citation Validation Report for {bib_path}",
                f"{'=' * 50}",
                f"Total entries: {report['total']}",
                f"Valid (DOI verified or match found): {report['valid']}",
                f"Errors: {len(report['errors'])}",
                f"Warnings: {len(report['warnings'])}",
            ]
            if report["errors"]:
                lines.append("\n--- ERRORS (must fix) ---")
                for err in report["errors"]:
                    lines.append(f"  [{err['key']}] {err['error']}")
            if report["warnings"]:
                lines.append("\n--- WARNINGS (review) ---")
                for w in report["warnings"]:
                    lines.append(f"  [{w['key']}] {w['warning']}")
                    if "suggested_doi" in w:
                        lines.append(f"    Suggested DOI: {w['suggested_doi']}")
            if report["corrected_entries"]:
                lines.append("\n--- CORRECTED ENTRIES (verified via CrossRef) ---")
                for c in report["corrected_entries"][:5]:  # Show max 5
                    lines.append(f"\n  Key: {c['key']}")
                    lines.append(f"  {c['corrected_bibtex']}")
            return "\n".join(lines)

        else:
            return (
                f"ERROR: Unknown action '{action}'. Use 'search', 'doi', or 'validate'."
            )

    except Exception as e:
        return f"Citation tool error: {e}\n\nYou can still manually look up papers at https://doi.org/ or https://scholar.google.com/"


@mcp.tool()
async def autolab_editorial(
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
) -> str:
    """Wait for the Editor (user) to act on the submitted manuscript.

    Call this in two scenarios:
    1. After PI submits a manuscript (autolab_record with status='ready_for_review').
       Waits for the editor to either Desk Reject or Invite Reviewers.
    2. After all reviewer reports are complete.
       Waits for the editor to make a final decision (Accept/Minor/Major/Reject).

    This tool BLOCKS until the editor acts OR the configured timeout expires
    (default 30 min, configurable via editor_timeout_minutes in config.yaml,
    set to 0 to wait forever). When timeout expires, returns an AI Editor
    prompt — you must then act as editor and call autolab_editor_act.

    Returns:
      - If reviewers invited: instructions to run the reviewer loop
        (call autolab_next for each reviewer, then autolab_editorial again)
      - If decision made (accept/minor/major/reject): the decision + instructions
        to call lab_meeting then autolab_next to continue PI/Trainee loop.
      - If timeout: AI Editor prompt + instructions to call autolab_editor_act.
    """
    import asyncio

    from .lab.state import load_state, get_domain_config, get_editor_timeout_seconds, is_task_mode

    project_directory = os.path.abspath(project_directory)

    if is_task_mode(project_directory):
        return "ERROR: Task mode does not use editorial review."

    # Ensure the web UI is running — the user needs it to make editorial decisions
    _ensure_lab_ui(project_directory)

    try:
        state = load_state(project_directory)
    except FileNotFoundError as e:
        return f"ERROR: {e}"

    editorial = state.get("editorial", {})
    starting_phase = editorial.get("phase", "none")
    iteration = state.get("iteration", 0)
    current_status = state.get("status", "")
    next_role = state.get("next_role", "")
    dcfg = get_domain_config(project_directory)

    debug_log(
        f"autolab_editorial: starting phase={starting_phase}, status={current_status}, next_role={next_role}"
    )

    # ── Immediate returns: if the state already advanced, don't poll ──

    # Already in reviewer phase → return immediately with reviewer instructions
    if starting_phase in ("reviewers_invited", "under_review") and next_role.startswith(
        "reviewer_"
    ):
        reviewers = editorial.get("reviewers", [])
        reviewer_names = ", ".join(r.get("name", r["id"]) for r in reviewers)
        rl = dcfg["reviewer_label"].lower()
        return (
            f"[EDITORIAL] {dcfg['reviewer_label']}s invited: {reviewer_names}\n\n"
            f"The {dcfg['overseer_label'].lower()} has invited {len(reviewers)} {rl}(s) and reviews are in progress.\n"
            f"You must now play each {rl} role in sequence.\n\n"
            f"MANDATORY NEXT STEPS:\n"
            f"1. Call autolab_next — it will return the current {rl}'s prompt.\n"
            f"2. Act as that {rl}, then call autolab_record with role='reviewer_N'.\n"
            f"3. Call autolab_next immediately for the next {rl} (do NOT call lab_meeting between {rl}s).\n"
            f"4. Repeat until all reviews are done.\n"
            f"5. After the LAST review, call autolab_editorial to wait for the {dcfg['overseer_label'].lower()}'s final decision.\n"
            f"NEVER stop. Continue the loop."
        )

    # Already have a decision → return immediately
    if starting_phase == "decision_made":
        decision = editorial.get("decision", "")
        feedback = editorial.get("decision_feedback", "")
        dec_display = decision.upper().replace("_", " ")
        return (
            f"[EDITORIAL] Decision: {dec_display}\n\n"
            f"{dcfg['overseer_label']} feedback: {feedback}\n\n"
            f"The {dcfg['review_process'].lower()} for this round is complete.\n"
            f"Decision: {dec_display}" + _LOOP_INSTRUCTION
        )

    # Reviews already complete → wait for decision only
    if starting_phase == "reviews_complete":
        debug_log("autolab_editorial: all reviews done, waiting for editor decision")

    # ── Timeout configuration ──
    timeout_seconds = get_editor_timeout_seconds(project_directory)
    elapsed_seconds = 0
    debug_log(f"autolab_editorial: timeout={timeout_seconds}s (0=disabled)")

    # ── Poll loop: wait for the editorial phase to advance ──

    poll_interval = 5  # seconds
    _health_check_counter = 0
    while True:
        await asyncio.sleep(poll_interval)
        elapsed_seconds += poll_interval
        _health_check_counter += 1

        # Periodic health check — restart the web server if it died (every ~60s)
        if _health_check_counter % 12 == 0:
            try:
                _ensure_lab_ui(project_directory)
            except Exception:
                pass  # non-fatal

        try:
            state = load_state(project_directory)
        except Exception:
            continue  # Retry on transient read errors

        editorial = state.get("editorial", {})
        current_phase = editorial.get("phase", "none")
        current_status = state.get("status", "")
        next_role = state.get("next_role", "")

        # --- Case 1: We were waiting for overseer after submission ---
        if starting_phase == "submitted":
            if current_phase in ("reviewers_invited", "under_review"):
                # Overseer invited reviewers → AI must now play reviewer roles
                reviewers = editorial.get("reviewers", [])
                reviewer_names = ", ".join(r.get("name", r["id"]) for r in reviewers)
                rl = dcfg["reviewer_label"].lower()
                return (
                    f"[EDITORIAL] {dcfg['reviewer_label']}s invited: {reviewer_names}\n\n"
                    f"The {dcfg['overseer_label'].lower()} has invited {len(reviewers)} {rl}(s).\n"
                    f"You must now play each {rl} role in sequence.\n\n"
                    f"MANDATORY NEXT STEPS:\n"
                    f"1. Call autolab_next — it will return the first {rl}'s prompt.\n"
                    f"2. Act as that {rl}, then call autolab_record with role='reviewer_N'.\n"
                    f"3. Call autolab_next immediately for the next {rl} (do NOT call lab_meeting between {rl}s).\n"
                    f"4. Repeat until all reviews are done.\n"
                    f"5. After the LAST review, call autolab_editorial to wait for the {dcfg['overseer_label'].lower()}'s final decision.\n"
                    f"NEVER stop. Continue the loop."
                )
            if current_phase == "decision_made":
                # Desk reject (overseer decided without reviewers)
                decision = editorial.get("decision", "reject")
                feedback = editorial.get("decision_feedback", "")
                return (
                    f"[EDITORIAL] Decision: {decision.upper().replace('_', ' ')}\n\n"
                    f"The {dcfg['overseer_label'].lower()} has desk-rejected the {dcfg['artifact'].lower()}.\n"
                    f"Feedback: {feedback}\n" + _LOOP_INSTRUCTION
                )

        # --- Case 2: We were waiting for final decision after reviews ---
        if starting_phase in ("reviews_complete", "under_review", "reviewers_invited"):
            if current_phase == "decision_made":
                decision = editorial.get("decision", "")
                feedback = editorial.get("decision_feedback", "")
                dec_display = decision.upper().replace("_", " ")
                return (
                    f"[EDITORIAL] Decision: {dec_display}\n\n"
                    f"{dcfg['overseer_label']} feedback: {feedback}\n\n"
                    f"The {dcfg['review_process'].lower()} for this round is complete.\n"
                    f"Decision: {dec_display}" + _LOOP_INSTRUCTION
                )
            # Reviewer role assigned while we're polling → break out to let AI play reviewer
            if next_role.startswith("reviewer_") and current_phase in (
                "reviewers_invited",
                "under_review",
            ):
                reviewers = editorial.get("reviewers", [])
                reviewer_names = ", ".join(r.get("name", r["id"]) for r in reviewers)
                rl = dcfg["reviewer_label"].lower()
                return (
                    f"[EDITORIAL] {dcfg['reviewer_label']} turn ready: {next_role}\n\n"
                    f"{dcfg['reviewer_label']}s: {reviewer_names}\n\n"
                    f"MANDATORY NEXT STEPS:\n"
                    f"1. Call autolab_next to get the {rl} prompt.\n"
                    f"2. Act as the {rl}, then call autolab_record with role='{next_role}'.\n"
                    f"3. Call autolab_next immediately for the next {rl} (skip lab_meeting between {rl}s).\n"
                    f"4. After the LAST review, call autolab_editorial again.\n"
                    f"NEVER stop."
                )

        # --- Generic: any phase change from what we started with ---
        if current_phase != starting_phase:
            return (
                f"[EDITORIAL] Phase changed: {starting_phase} → {current_phase}\n"
                f"Status: {current_status}\n" + _LOOP_INSTRUCTION
            )

        # ── Timeout check: if enabled and expired, switch to AI editor ──
        if timeout_seconds > 0 and elapsed_seconds >= timeout_seconds:
            debug_log(
                f"autolab_editorial: editor timeout after {elapsed_seconds}s — "
                f"switching to AI editor for phase={current_phase}"
            )
            return _build_ai_editor_return(
                project_directory, editorial, current_phase, elapsed_seconds
            )

        # Still waiting — continue polling
        # Gradually increase poll interval up to 15s to reduce disk I/O
        if poll_interval < 15:
            poll_interval = min(poll_interval + 1, 15)


def _build_ai_editor_return(
    project_directory: str,
    editorial: dict,
    phase: str,
    elapsed_seconds: int,
) -> str:
    """Construct the return value when editor timeout expires."""
    from .lab.prompts import build_ai_editor_prompt
    from .lab.state import get_domain_config, scan_project_files

    dcfg = get_domain_config(project_directory)
    file_listings = scan_project_files(project_directory)
    cover_letter = editorial.get("cover_letter", "")
    reviews = editorial.get("reviews", {})

    editor_prompt = build_ai_editor_prompt(
        editorial=editorial,
        phase=phase,
        cover_letter=cover_letter,
        reviews=reviews,
        file_listings=file_listings,
    )

    mins = elapsed_seconds // 60
    rl = dcfg["reviewer_label"].lower()
    return (
        f"[EDITORIAL — AI {dcfg['overseer_label'].upper()} ACTIVATED]\n\n"
        f"The human {dcfg['overseer_label'].lower()} did not respond within {mins} minutes.\n"
        f"You must now act as the {dcfg['overseer_label']} and make a decision.\n\n"
        f"{editor_prompt}\n\n"
        f"MANDATORY NEXT STEPS after making your decision:\n"
        f"1. Call `autolab_editor_act` with your action and feedback.\n"
        f"   - For desk reject: action='desk_reject', feedback='...'\n"
        f"   - For invite {rl}s: action='invite_reviewers', feedback='...', "
        f'reviewers=[{{"name":"Dr. X","role":"Specialty","avatar":"sprite_key"}}]\n'
        f"   - For final decision: action='accept'|'minor_revision'|'major_revision'|'reject', feedback='...'\n"
        f"2. Then follow the MANDATORY NEXT STEPS returned by autolab_editor_act.\n"
        f"NEVER stop. Continue the loop."
    )


@mcp.tool()
async def autolab_editor_act(
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
    action: Annotated[
        str,
        Field(
            description="Editorial action: 'desk_reject', 'invite_reviewers', "
            "'accept', 'minor_revision', 'major_revision', 'reject'"
        ),
    ] = "",
    feedback: Annotated[
        str, Field(description="Editor's feedback to the authors")
    ] = "",
    reviewers: Annotated[
        str,
        Field(
            description="JSON array of reviewer dicts for invite_reviewers, "
            'e.g. [{"name":"Dr. X","role":"Statistician","avatar":"statistician"}]. '
            "Ignored for other actions."
        ),
    ] = "[]",
) -> str:
    """Execute an editorial decision made by the AI editor (after timeout).

    Called by the AI after autolab_editorial returns an AI-editor prompt
    due to the human editor not responding within the configured timeout.

    Actions:
    - desk_reject: reject without review, feedback required
    - invite_reviewers: send to peer review, reviewers JSON required
    - accept / minor_revision / major_revision / reject: final decision after reviews

    Returns instructions to continue the loop.
    """
    from .lab.state import (
        get_domain_config,
        load_state,
        record_editorial_decision,
        invite_reviewers as _invite_reviewers,
    )

    project_directory = os.path.abspath(project_directory)

    try:
        state = load_state(project_directory)
    except FileNotFoundError as e:
        return f"ERROR: {e}"

    dcfg = get_domain_config(project_directory)
    action = action.strip().lower()
    feedback = (
        feedback.strip()
        or f"(No feedback provided by AI {dcfg['overseer_label'].lower()})"
    )

    debug_log(f"autolab_editor_act: action={action}")

    if action == "desk_reject":
        record_editorial_decision(project_directory, "reject", feedback)
        return (
            f"[AI {dcfg['overseer_label'].upper()}] {dcfg['artifact']} desk-rejected.\n\n"
            f"Feedback sent to {dcfg['senior_label']}: {feedback[:200]}...\n"
            + _LOOP_INSTRUCTION
        )

    elif action == "invite_reviewers":
        # Parse reviewer list
        import json as _json

        try:
            reviewer_list = _json.loads(reviewers)
        except (ValueError, TypeError):
            reviewer_list = []

        if not reviewer_list or not isinstance(reviewer_list, list):
            # Fallback: generate 3 default reviewers
            reviewer_list = [
                {
                    "id": "reviewer_1",
                    "name": "Dr. Methods",
                    "role": "Statistician",
                    "avatar": "statistician",
                },
                {
                    "id": "reviewer_2",
                    "name": "Dr. Domain",
                    "role": "Domain Expert",
                    "avatar": "reviewer",
                },
                {
                    "id": "reviewer_3",
                    "name": "Dr. Scope",
                    "role": "Generalist",
                    "avatar": "generic",
                },
            ]

        # Ensure each reviewer has an id
        for i, r in enumerate(reviewer_list):
            if "id" not in r:
                r["id"] = f"reviewer_{i + 1}"

        _invite_reviewers(project_directory, reviewer_list)
        reviewer_names = ", ".join(r.get("name", r["id"]) for r in reviewer_list)
        rl = dcfg["reviewer_label"].lower()
        return (
            f"[AI {dcfg['overseer_label'].upper()}] {dcfg['reviewer_label']}s invited: {reviewer_names}\n\n"
            f"{dcfg['overseer_label']} note: {feedback[:200]}\n\n"
            f"You must now play each {rl} role in sequence.\n\n"
            f"MANDATORY NEXT STEPS:\n"
            f"1. Call autolab_next — it will return the first {rl}'s prompt.\n"
            f"2. Act as that {rl}, then call autolab_record with role='reviewer_N'.\n"
            f"3. Call autolab_next immediately for the next {rl} (skip lab_meeting between {rl}s).\n"
            f"4. After the LAST review, call autolab_editorial again.\n"
            f"NEVER stop."
        )

    elif action in ("accept", "minor_revision", "major_revision", "reject"):
        record_editorial_decision(project_directory, action, feedback)
        dec_display = action.upper().replace("_", " ")
        return (
            f"[AI {dcfg['overseer_label'].upper()}] Decision: {dec_display}\n\n"
            f"Feedback: {feedback[:200]}...\n" + _LOOP_INSTRUCTION
        )

    else:
        return (
            f"ERROR: Unknown action '{action}'. Must be one of: "
            f"desk_reject, invite_reviewers, accept, minor_revision, major_revision, reject.\n"
            f"Call autolab_editor_act again with a valid action."
        )


@mcp.tool()
async def autolab_status(
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
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
        get_domain_config,
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

    dcfg = get_domain_config(project_directory)
    file_listings = scan_project_files(project_directory)
    paper_progress = get_paper_progress(project_directory)
    recent = get_recent_meetings(project_directory, n=2)

    # Map internal role id to display label
    next_role_id = state["next_role"]
    if next_role_id == "pi":
        next_role_display = dcfg["senior_label"]
    elif next_role_id == "trainee":
        next_role_display = dcfg["junior_label"]
    else:
        next_role_display = next_role_id  # reviewer_N etc.

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
        f"**Domain:** {dcfg.get('senior_label', 'PI')}/{dcfg.get('junior_label', 'Trainee')}/{dcfg.get('overseer_label', 'Editor')}\n"
        f"**Iteration:** {state['iteration']}\n"
        f"**Next role:** {next_role_display}\n"
        f"**Status:** {state['status']}\n"
        f"**Last updated:** {state.get('last_updated', 'unknown')}\n\n"
        f"## File Counts\n"
        f"  data: {file_counts.get('data', 0)} files\n"
        f"  scripts: {file_counts.get('scripts', 0)} files\n"
        f"  figures: {file_counts.get('figures', 0)} files\n"
        f"  results: {file_counts.get('results', 0)} files\n"
        f"  paper: {file_counts.get('paper', 0)} files\n\n"
        f"## {dcfg['artifact']} Progress\n" + "\n".join(paper_lines) + "\n\n"
        f"## Recent {dcfg['meeting_log_name']}\n\n{recent if recent.strip() else 'No meetings yet.'}\n"
    )

    # Show queued feedback count (don't drain it — just peek)
    queue = state.get("feedback_queue", [])
    if not isinstance(queue, list):
        queue = [queue] if queue else []
    if queue:
        report += f"\n## Pending User Feedback ({len(queue)} message(s) queued)\n\n"
        for i, msg in enumerate(queue, 1):
            report += f"  {i}. {msg[:100]}{'...' if len(msg) > 100 else ''}\n"

    return report


# ===== ToolUniverse Integration =====
# Replaces the former Biomni biotools section.
# ToolUniverse (aiscientist.tools) provides 1000+ scientific tools.


@mcp.tool()
async def autolab_tools_status(
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
) -> str:
    """Check whether the ToolUniverse SDK is installed and available.

    ToolUniverse provides 1000+ scientific tools (ML models, datasets,
    APIs, analysis packages) from https://aiscientist.tools.

    Returns JSON with: installed (bool), version, enabled, api_url.
    """
    from .integrations.tooluniverse import get_status

    project_directory = os.path.abspath(project_directory)
    status = get_status(project_directory)
    return json.dumps(status, indent=2)


@mcp.tool()
async def autolab_tools_install(
    from_source: Annotated[
        bool,
        Field(description="Install from GitHub (True) or PyPI (False)"),
    ] = False,
) -> str:
    """Install the ToolUniverse SDK (one-time setup).

    Installs the tooluniverse Python package for local tool execution.
    Without it, tools still work via the HTTP API at aiscientist.tools.
    """
    from .integrations.tooluniverse import install_sdk

    result = install_sdk(from_source=from_source)
    if result["success"]:
        return f"SUCCESS: {result['message']}"
    else:
        return f"FAILED: {result['message']}"


@mcp.tool()
async def autolab_tools_configure(
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
    enabled: Annotated[
        bool, Field(description="Enable ToolUniverse for this project")
    ] = True,
    prefer_local: Annotated[
        bool, Field(description="Prefer local SDK over HTTP API")
    ] = True,
    api_url: Annotated[
        str, Field(description="ToolUniverse API URL")
    ] = "https://aiscientist.tools",
    api_key: Annotated[
        str, Field(description="Optional API key")
    ] = "",
) -> str:
    """Configure ToolUniverse for this project.

    Saves settings to .autolab/tools_config.json.
    """
    from .integrations.tooluniverse import save_config

    project_directory = os.path.abspath(project_directory)
    cfg = save_config(
        project_directory,
        enabled=enabled,
        prefer_local=prefer_local,
        api_url=api_url,
        api_key=api_key or None,
    )
    return f"ToolUniverse config saved:\n{json.dumps(cfg, indent=2)}"


@mcp.tool()
async def autolab_tools_search(
    query: Annotated[
        str,
        Field(description="Natural language search query (e.g., 'protein docking', "
              "'gene expression analysis', 'drug-target interaction')"),
    ] = "",
    limit: Annotated[
        int, Field(description="Maximum results to return")
    ] = 10,
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
) -> str:
    """Search ToolUniverse for scientific tools matching a query.

    Searches the catalog of 1000+ scientific tools including ML models,
    datasets, APIs, and analysis packages. Works with local SDK or
    via the HTTP API at aiscientist.tools.

    Use this to discover tools during research planning or when
    the PI/Trainee needs specialized capabilities.
    """
    from .integrations.tooluniverse import search_tools

    if not query:
        return "ERROR: query is required (e.g., 'protein structure prediction')"

    results = search_tools(query.strip(), limit=limit)
    if not results:
        return f"No tools found matching '{query}'. Try a broader query."

    lines = [f"## ToolUniverse Search: '{query}'\n"]
    for i, tool in enumerate(results, 1):
        name = tool.get("name", "unknown")
        desc = tool.get("description", "")
        lines.append(f"{i}. **{name}** — {desc}")

    lines.append(
        f"\n*{len(results)} tool(s) found. "
        "Tools are automatically integrated during character recruitment.*"
    )
    return "\n".join(lines)


# ===== Character Creation =====

# Available avatar sprite keys (must match EXPERT_DEFS in sprites.js)
_AVATAR_KEYS = [
    "reviewer",
    "bioethicist",
    "science_writer",
    "grant_reviewer",
    "immunologist",
    "oncologist",
    "neuroscientist",
    "geneticist",
    "cell_biologist",
    "microbiologist",
    "pathologist",
    "pharmacologist",
    "structural_bio",
    "systems_biologist",
    "epidemiologist",
    "statistician",
    "bioinformatician",
    "data_scientist",
    "ml_engineer",
    "comp_biologist",
    "clinician",
    "radiologist",
    "surgeon",
    "chemist",
    "physicist",
    "engineer",
    "psychologist",
    "ecologist",
    "generic",
]

# Built-in skill catalog for reference (users can also create their own)
_SKILL_CATALOG = {
    "Bioinformatics & Omics": [
        "scanpy",
        "scvi-tools",
        "pydeseq2",
        "pysam",
        "deeptools",
        "anndata",
        "cellxgene-census",
        "geniml",
        "gtars",
        "arboreto",
    ],
    "Machine Learning & AI": [
        "scikit-learn",
        "pytorch-lightning",
        "transformers",
        "accelerate",
        "shap",
        "torch_geometric",
        "deepchem",
        "torchdrug",
        "umap-learn",
        "stable-baselines3",
    ],
    "Statistics & Modeling": [
        "statistical-analysis",
        "statsmodels",
        "pymc",
        "scikit-survival",
        "sympy",
        "pymoo",
    ],
    "Chemistry & Drug Discovery": [
        "rdkit",
        "datamol",
        "medchem",
        "molfeat",
        "pytdc",
        "diffdock",
        "rowan",
        "matchms",
    ],
    "Visualization": [
        "scientific-visualization",
        "matplotlib",
        "seaborn",
        "plotly",
    ],
    "Writing & Communication": [
        "scientific-writing",
        "literature-review",
        "peer-review",
        "scientific-analysis-review",
        "scientific-brainstorming",
        "hypothesis-generation",
        "venue-templates",
        "latex-posters",
        "scientific-slides",
    ],
    "Clinical & Databases": [
        "clinical-reports",
        "clinicaltrials-database",
        "clinvar-database",
        "pubmed-database",
        "geo-database",
        "opentargets-database",
        "string-database",
        "kegg-database",
        "reactome-database",
        "uniprot-database",
        "pdb-database",
        "pubchem-database",
        "chembl-database",
        "ensembl-database",
    ],
    "Infrastructure & MLOps": [
        "weights-and-biases",
        "mlflow",
        "tensorboard",
        "modal",
        "vllm",
        "sglang",
    ],
    "Data Processing": [
        "polars",
        "dask",
        "vaex",
        "zarr-python",
        "geopandas",
        "exploratory-data-analysis",
    ],
}

# Flatten catalog skills for reference (not for validation — users can add custom skills)
_KNOWN_SKILLS = sorted({s for skills in _SKILL_CATALOG.values() for s in skills})


@mcp.tool()
async def autolab_create_character(
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
    name: Annotated[
        str, Field(description="Character's display name (e.g., 'Dr. Maria Chen')")
    ] = "",
    role: Annotated[
        str, Field(description="Role type: 'pi', 'trainee', or 'collaborator'")
    ] = "",
    title: Annotated[
        str, Field(description="Professional title (e.g., 'Computational Biology PI')")
    ] = "",
    expertise: Annotated[
        str,
        Field(
            description="Areas of expertise (e.g., 'Single-cell genomics, machine learning')"
        ),
    ] = "",
    goal: Annotated[str, Field(description="Character's research goal")] = "",
    skills: Annotated[
        str,
        Field(
            description="Comma-separated skill names "
            "(e.g., 'scanpy,scvi-tools,scientific-writing,my-custom-tool'). "
            "Can be built-in Cursor skills or your own custom SKILL.md names. "
            "Use list_skills=True to browse the built-in catalog."
        ),
    ] = "",
    personality: Annotated[
        str,
        Field(
            description="Pipe-separated personality traits "
            "(e.g., 'Rigorous: demands statistical reproducibility|"
            "Visionary: spots novel research directions')"
        ),
    ] = "",
    avatar: Annotated[
        str,
        Field(
            description="Avatar sprite key for the game UI: reviewer, bioethicist, "
            "science_writer, grant_reviewer, immunologist, oncologist, "
            "neuroscientist, geneticist, cell_biologist, microbiologist, "
            "pathologist, pharmacologist, structural_bio, systems_biologist, "
            "epidemiologist, statistician, bioinformatician, data_scientist, "
            "ml_engineer, comp_biologist, clinician, radiologist, surgeon, "
            "chemist, physicist, engineer, psychologist, ecologist, generic"
        ),
    ] = "generic",
    deploy_as: Annotated[
        str,
        Field(
            description="Deploy locally as 'pi' or 'trainee' profile "
            "(saves to .autolab/profiles/). Leave empty to only generate the YAML."
        ),
    ] = "",
    github_repo: Annotated[
        str,
        Field(
            description="Optional GitHub repo path for marketplace listing "
            "(e.g., 'username/autolab-char-name')"
        ),
    ] = "",
    list_skills: Annotated[
        bool,
        Field(
            description="If True, ignore other params and return the built-in skill catalog "
            "(you can also use any custom skill name not listed here)"
        ),
    ] = False,
    list_avatars: Annotated[
        bool,
        Field(
            description="If True, ignore other params and return all available avatar keys"
        ),
    ] = False,
) -> str:
    """Create a character profile for the Autonomous Lab research team.

    A character is a YAML persona backed by Cursor agent skills that controls
    how the AI behaves during research sessions. Characters can be:
    - Deployed locally as a PI or Trainee profile
    - Published to GitHub for the marketplace
    - Used as consultants or reviewers during editorial review

    Each character is a **self-contained folder** with this structure:

        autolab-char-<name>/
        ├── character.yaml          # Profile (name, role, skills, personality)
        ├── README.md               # Description for marketplace
        └── skills/                 # Skill definitions (SKILL.md files)
            ├── <skill-name>/
            │   └── SKILL.md        # Detailed instructions for this skill
            └── ...

    Use list_skills=True or list_avatars=True to discover options first.

    Workflow:
    1. Call with list_skills=True to see available skills by domain.
    2. Call with all fields filled to generate the character folder.
    3. Optionally set github_repo to initialize a GitHub repo and push.

    Args:
        project_directory: Project directory path
        name: Character display name
        role: pi, trainee, or collaborator
        title: Professional title
        expertise: Areas of expertise
        goal: Research goal
        skills: Comma-separated Cursor skill names
        personality: Pipe-separated personality traits (Trait: description)
        avatar: Sprite key for the game UI
        deploy_as: Deploy as 'pi' or 'trainee' profile locally
        github_repo: GitHub repo for marketplace listing
        list_skills: Return skill catalog instead of creating
        list_avatars: Return avatar list instead of creating

    Returns:
        str: Generated YAML, deployment confirmation, or catalog listing
    """
    import yaml as _yaml
    from pathlib import Path

    # ── Discovery modes ──

    if list_skills:
        lines = ["# Built-in Cursor Skills (Reference Catalog)\n"]
        lines.append(
            f"Total: {len(_KNOWN_SKILLS)} built-in skills across "
            f"{len(_SKILL_CATALOG)} domains\n"
        )
        for domain, skill_list in _SKILL_CATALOG.items():
            lines.append(f"\n## {domain}")
            for s in skill_list:
                lines.append(f"  - {s}")
        lines.append(
            "\n---\n"
            "## Custom Skills\n"
            "You are NOT limited to this catalog. You can use ANY custom skill name.\n"
            "Custom skills just need a SKILL.md file in ~/.cursor/skills/<name>/SKILL.md.\n"
            "Example: skills='scanpy,my-protein-docking-tool,lab-specific-pipeline'\n"
            "\n## Tips\n"
            "- Choose 3-6 skills per character for focus\n"
            "- Group related tools into a single skill (e.g., 'ngs-pipeline' "
            "covering pysam + deeptools + samtools)\n"
            "- PI characters should include a writing skill\n"
            "- Trainee characters should include hands-on analysis tools\n"
        )
        return "\n".join(lines)

    if list_avatars:
        lines = ["# Available Avatar Sprite Keys\n"]
        lines.append("These control the pixel-art character in the game UI:\n")
        for a in _AVATAR_KEYS:
            lines.append(f"  - {a}")
        lines.append(
            "\n---\n" "Pass your chosen avatar key:\n" "  avatar='neuroscientist'"
        )
        return "\n".join(lines)

    # ── Validation ──

    errors = []
    if not name:
        errors.append("'name' is required (e.g., 'Dr. Maria Chen')")
    if role not in ("pi", "trainee", "collaborator"):
        errors.append("'role' must be 'pi', 'trainee', or 'collaborator'")
    if not title:
        errors.append("'title' is required (e.g., 'Computational Biology PI')")
    if not expertise:
        errors.append("'expertise' is required")
    if not goal:
        errors.append("'goal' is required")
    if not skills:
        errors.append("'skills' is required (comma-separated skill names)")
    if not personality:
        errors.append(
            "'personality' is required (pipe-separated 'Trait: description' entries)"
        )
    if avatar not in _AVATAR_KEYS:
        errors.append(
            f"'avatar' must be one of: {', '.join(_AVATAR_KEYS[:5])}... "
            "(use list_avatars=True)"
        )
    if deploy_as and deploy_as not in ("pi", "trainee"):
        errors.append("'deploy_as' must be 'pi' or 'trainee' (or empty)")

    if errors:
        return (
            "ERROR: Missing or invalid fields:\n"
            + "\n".join(f"  - {e}" for e in errors)
            + "\n\nUse list_skills=True or list_avatars=True to discover options."
        )

    # ── Parse inputs ──

    skill_list = [s.strip() for s in skills.split(",") if s.strip()]
    personality_list = [p.strip() for p in personality.split("|") if p.strip()]

    # ── Build character folder ──

    project_dir = Path(os.path.abspath(project_directory))
    char_slug = name.lower().replace(" ", "-").replace(".", "").replace("'", "")
    folder_name = f"autolab-char-{char_slug}"
    char_dir = project_dir / folder_name
    char_dir.mkdir(parents=True, exist_ok=True)

    # -- character.yaml --
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
    char_yaml_path = char_dir / "character.yaml"
    char_yaml_path.write_text(yaml_body, encoding="utf-8")

    # -- skills/ directory with SKILL.md skeletons --
    skills_dir = char_dir / "skills"
    skills_dir.mkdir(exist_ok=True)

    for skill_name in skill_list:
        skill_folder = skills_dir / skill_name
        skill_folder.mkdir(exist_ok=True)
        skill_md_path = skill_folder / "SKILL.md"
        if not skill_md_path.exists():
            skill_md_path.write_text(
                f"---\n"
                f"name: {skill_name}\n"
                f"description: TODO — describe what this skill does\n"
                f"metadata:\n"
                f"  skill-author: {name}\n"
                f"---\n\n"
                f"# {skill_name.replace('-', ' ').title()}\n\n"
                f"## When to use\n\n"
                f"- TODO: describe when the AI should apply this skill\n\n"
                f"## Standard workflow\n\n"
                f"```python\n"
                f"# TODO: add code examples and instructions\n"
                f"```\n\n"
                f"## Key decisions\n\n"
                f"- TODO: document important choices and defaults\n",
                encoding="utf-8",
            )
        # Write meta.yaml for self-learning skill tracking
        meta_path = skill_folder / "meta.yaml"
        if not meta_path.exists():
            meta_path.write_text(
                "status: uncertified\n"
                "validation_type: content-only\n"
                "source: created\n",
                encoding="utf-8",
            )

    # -- README.md --
    skills_tree = "\n".join(
        f"    ├── {s}/\n    │   ├── SKILL.md\n    │   └── meta.yaml" for s in skill_list[:-1]
    )
    if skill_list:
        skills_tree += f"\n    └── {skill_list[-1]}/\n        ├── SKILL.md\n        └── meta.yaml"

    readme_content = (
        f"# {name} — {title}\n\n"
        f"Character for [Autonomous Lab](https://autolab.kejunying.com). "
        f"{expertise.rstrip('.')}.\n\n"
        f"Skills are auto-learned and validated. Each skill has a `meta.yaml` "
        f"tracking certification status.\n\n"
        f"## Structure\n\n"
        f"```\n"
        f"{folder_name}/\n"
        f"├── character.yaml\n"
        f"├── README.md\n"
        f"└── skills/\n"
        f"{skills_tree}\n"
        f"```\n\n"
        f"## Install\n\n"
        f"```shell\n"
        f"git clone https://github.com/USERNAME/{folder_name} "
        f".autolab/characters/{char_slug}\n"
        f"```\n\n"
        f"## License\n\n"
        f"Apache 2.0\n"
    )
    readme_path = char_dir / "README.md"
    readme_path.write_text(readme_content, encoding="utf-8")

    # ── Output summary ──

    output_lines = [
        f"Character folder created: {char_dir}/\n",
        f"  {folder_name}/",
        f"  ├── character.yaml",
        f"  ├── README.md",
        f"  └── skills/",
    ]
    for s in skill_list:
        output_lines.append(f"      └── {s}/")
        output_lines.append(f"          ├── SKILL.md")
        output_lines.append(f"          └── meta.yaml")

    output_lines.append(f"\n--- character.yaml ---\n{yaml_body}--- End YAML ---\n")

    # ── IMPORTANT: remind about self-learning skills ──
    output_lines.append(
        "NEXT STEPS:\n"
        "  1. Run `autolab_skill_status` to see certification status of each skill.\n"
        "  2. Skills are auto-learned during the PI-Trainee loop. Each skill has a\n"
        "     meta.yaml tracking validation status (uncertified → certified).\n"
        "  3. To manually validate, run `validation.py` in the skill folder.\n"
        "  4. Fill in each skills/<name>/SKILL.md with detailed instructions:\n"
        "     - When to use the skill\n"
        "     - Standard workflow with code examples\n"
        "     - Key decisions and defaults\n"
        "The AI reads these instructions when acting as this character.\n"
    )

    # ── Deploy locally ──

    if deploy_as:
        autolab_dir = project_dir / ".autolab" / "profiles"
        autolab_dir.mkdir(parents=True, exist_ok=True)
        target = autolab_dir / f"{deploy_as}.yaml"
        target.write_text(yaml_body, encoding="utf-8")
        output_lines.append(
            f"Deployed as {deploy_as} profile: {target}\n"
            f"This character is now active. Run autolab_next to use it."
        )

    if len(skill_list) > 8:
        output_lines.append(
            "\nWARNING: More than 8 skills can dilute focus. "
            "Consider grouping related tools into fewer composite skills."
        )

    # ── GitHub publishing ──

    # Build topic commands for discoverability
    topic_cmds = "  gh repo edit --add-topic autolab-character\n"
    for s in skill_list:
        topic_cmds += f"  gh repo edit --add-topic autolab-skill-{s}\n"

    if github_repo:
        clean_repo = github_repo.strip().strip("/")
        output_lines.append(
            f"\n--- GitHub Publishing ---\n"
            f"To publish this character to the marketplace, run these commands "
            f"inside the character folder ({char_dir}):\n\n"
            f"  cd {char_dir}\n"
            f"  git init\n"
            f"  git add .\n"
            f"  git commit -m 'Character: {name}'\n"
            f"  gh repo create {clean_repo} --public --source=. --push\n\n"
            f"Add GitHub topics for auto-discovery by autolab_recruit:\n"
            f"{topic_cmds}\n"
            f"Then submit to the marketplace:\n"
            f"  Open an issue at https://github.com/albert-ying/autonomous-lab/issues/new\n"
            f"  Title: New Character: {name}\n"
            f"  Body: https://github.com/{clean_repo}\n"
        )
    else:
        output_lines.append(
            "\nTo publish to the marketplace, call again with "
            "github_repo='username/autolab-char-name', or run:\n"
            f"  cd {char_dir} && git init && git add . && "
            f"git commit -m 'Character: {name}'\n"
            f"  gh repo create username/{folder_name} --public --source=. --push\n\n"
            f"Add GitHub topics for auto-discovery by autolab_recruit:\n"
            f"{topic_cmds}"
        )

    return "\n".join(output_lines)


# ── Marketplace registry URL ──
_REGISTRY_URL = (
    "https://albert-ying.github.io/autonomous-lab/characters/index.json"
)


@mcp.tool()
async def autolab_acquire_skill(
    skill_name: Annotated[
        str,
        Field(
            description="The skill to acquire (e.g., 'scanpy', 'survival-analysis', "
            "'pytorch-lightning'). The tool first searches the character marketplace "
            "for an existing character that has this skill. If found, it downloads "
            "the SKILL.md. If not, it returns instructions to train the skill."
        ),
    ] = "",
    project_directory: Annotated[
        str, Field(description="Project directory path")
    ] = ".",
) -> str:
    """Search the character marketplace for a skill and acquire it.

    When the PI or Trainee encounters a skill gap during the loop, call this
    tool BEFORE attempting to build it from scratch. It checks the published
    character registry for any character that already has the needed skill.

    If a matching skill is found:
      - Downloads the SKILL.md from the character's GitHub repo
      - Saves it to .autolab/acquired_skills/<skill>/SKILL.md
      - Returns the skill content so you can use it immediately

    If no match is found:
      - Returns instructions to create a new character with the skill
      - Or to manually write the SKILL.md

    This enables in-loop skill acquisition: the research never stops because
    of a missing capability.

    Args:
        skill_name: Name of the skill to find
        project_directory: Project directory path

    Returns:
        str: The SKILL.md content if found, or instructions to create it
    """
    import json
    import urllib.parse
    import urllib.request
    from pathlib import Path

    if not skill_name:
        return "ERROR: skill_name is required (e.g., 'scanpy', 'survival-analysis')"

    skill_name = skill_name.strip().lower()
    project_dir = Path(os.path.abspath(project_directory))
    acquired_dir = project_dir / ".autolab" / "acquired_skills" / skill_name
    skill_md_path = acquired_dir / "SKILL.md"

    # Check if already acquired
    if skill_md_path.exists():
        content = skill_md_path.read_text(encoding="utf-8")
        return (
            f"Skill '{skill_name}' already acquired locally.\n\n"
            f"--- SKILL.md ---\n{content}\n--- End ---\n\n"
            f"Read and follow the instructions above."
        )

    output = [f"Searching marketplace for skill: {skill_name}...\n"]

    # Step 1: Fetch the character registry
    matches = []
    registry = {}
    try:
        req = urllib.request.Request(_REGISTRY_URL, headers={
            "User-Agent": "AutonomousLab/0.5",
            "Accept": "application/json",
        })
        with urllib.request.urlopen(req, timeout=10) as resp:
            registry = json.loads(resp.read().decode("utf-8"))

        for char in registry.get("characters", []):
            raw_skills = char.get("skills", [])
            if isinstance(raw_skills, dict):
                char_skills = [s.lower() for s in raw_skills.keys()]
            elif raw_skills and isinstance(raw_skills[0], dict):
                char_skills = [s.get("name", "").lower() for s in raw_skills]
            else:
                char_skills = [s.lower() for s in raw_skills]
            # Exact match or substring match
            if skill_name in char_skills or any(
                skill_name in s or s in skill_name for s in char_skills
            ):
                matches.append(char)
    except Exception as e:
        output.append(f"Registry fetch failed ({e}), trying GitHub search...\n")

    # Step 1b: No exact match — suggest closest available skills
    if not matches and registry:
        all_skills: dict[str, list[str]] = {}  # skill -> [char names]
        for char in registry.get("characters", []):
            raw_sk = char.get("skills", [])
            if isinstance(raw_sk, dict):
                sk_names = [s.lower() for s in raw_sk.keys()]
            elif raw_sk and isinstance(raw_sk[0], dict):
                sk_names = [s.get("name", "").lower() for s in raw_sk]
            else:
                sk_names = [s.lower() for s in raw_sk]
            for s in sk_names:
                all_skills.setdefault(s, []).append(char.get("name", char.get("id", "?")))

        # Score by keyword overlap between requested skill and available names
        def _score(available: str, query: str) -> int:
            q_words = set(query.replace("-", " ").replace("_", " ").split())
            a_words = set(available.replace("-", " ").replace("_", " ").split())
            return len(q_words & a_words) * 3 + (
                2 if query in available or available in query else 0
            )

        scored = [(s, _score(s, skill_name)) for s in all_skills]
        scored.sort(key=lambda x: -x[1])
        suggestions = [(s, sc) for s, sc in scored if sc > 0][:5]

        if suggestions:
            lines = [f"No exact match for '{skill_name}'. Similar available skills:"]
            for s, _ in suggestions:
                owners = ", ".join(all_skills[s][:2])
                lines.append(f"  - {s}  (characters: {owners})")
            lines.append(
                f"\nCall autolab_acquire_skill with one of these exact names "
                f"to download the SKILL.md."
            )
            output.append("\n".join(lines) + "\n")
            return "\n".join(output)
        else:
            # No keyword overlap — list all available skills
            output.append(
                f"No match for '{skill_name}'. Available skills in registry:\n"
                f"  {', '.join(sorted(all_skills.keys()))}\n"
                f"\nCall autolab_acquire_skill with one of these exact names.\n"
            )
            return "\n".join(output)

    # Step 2: Try GitHub search as fallback
    if not matches:
        try:
            query = urllib.parse.quote(
                f"autolab-char {skill_name} in:readme,name"
            )
            gh_url = (
                f"https://api.github.com/search/repositories"
                f"?q={query}+topic:autolab-character&per_page=5"
            )
            req = urllib.request.Request(gh_url, headers={
                "User-Agent": "AutonomousLab/0.5",
                "Accept": "application/vnd.github+json",
            })
            with urllib.request.urlopen(req, timeout=10) as resp:
                gh_data = json.loads(resp.read().decode("utf-8"))
            for repo in gh_data.get("items", []):
                matches.append({
                    "name": repo.get("name", ""),
                    "github": repo.get("full_name", ""),
                    "title": repo.get("description", ""),
                    "stars": repo.get("stargazers_count", 0),
                    "skills": [],  # Unknown until we fetch character.yaml
                })
        except Exception:
            pass

    # Step 3: If found, try to download the SKILL.md
    if matches:
        best = max(matches, key=lambda c: c.get("stars", 0))
        github_repo = best.get("github", "")
        output.append(
            f"Found matching character: {best.get('name', github_repo)}\n"
            f"  Repo: https://github.com/{github_repo}\n"
            f"  Stars: {best.get('stars', '?')}\n"
        )

        # Try to download the specific SKILL.md
        skill_content = None
        for try_name in [skill_name, skill_name.replace("-", "_")]:
            try:
                raw_url = (
                    f"https://raw.githubusercontent.com/{github_repo}"
                    f"/master/skills/{try_name}/SKILL.md"
                )
                req = urllib.request.Request(raw_url, headers={
                    "User-Agent": "AutonomousLab/0.5",
                })
                with urllib.request.urlopen(req, timeout=10) as resp:
                    skill_content = resp.read().decode("utf-8")
                break
            except Exception:
                # Try main branch
                try:
                    raw_url = (
                        f"https://raw.githubusercontent.com/{github_repo}"
                        f"/main/skills/{try_name}/SKILL.md"
                    )
                    req = urllib.request.Request(raw_url, headers={
                        "User-Agent": "AutonomousLab/0.5",
                    })
                    with urllib.request.urlopen(req, timeout=10) as resp:
                        skill_content = resp.read().decode("utf-8")
                    break
                except Exception:
                    continue

        if skill_content:
            # Save locally
            acquired_dir.mkdir(parents=True, exist_ok=True)
            skill_md_path.write_text(skill_content, encoding="utf-8")

            output.append(
                f"Downloaded SKILL.md → {skill_md_path}\n\n"
                f"--- SKILL.md ---\n{skill_content}\n--- End ---\n\n"
                f"Skill acquired. Read and follow the instructions above."
            )
            return "\n".join(output)
        else:
            # Character found but skill folder not matching
            output.append(
                f"Character exists but skill '{skill_name}' not found as a "
                f"separate SKILL.md in that repo.\n"
                f"The character's skills: {best.get('skills', [])}\n"
                f"You can still clone the full character:\n"
                f"  git clone https://github.com/{github_repo} "
                f".autolab/characters/{best.get('id', 'char')}\n"
            )

    # Step 2b: Try ToolUniverse catalog
    if not matches and not skill_md_path.exists():
        try:
            from autonomous_lab.integrations.tooluniverse import (
                search_tools, get_tool_spec, generate_skill_md,
            )
            tu_results = search_tools(skill_name, limit=3)
            if tu_results:
                best_tool = tu_results[0]
                spec = get_tool_spec(best_tool["name"]) or best_tool
                skill_content = generate_skill_md(spec)

                # Save locally
                acquired_dir.mkdir(parents=True, exist_ok=True)
                skill_md_path.write_text(skill_content, encoding="utf-8")
                # Write meta.yaml
                (acquired_dir / "meta.yaml").write_text(
                    "validation_type: content-only\nsource: tooluniverse\n",
                    encoding="utf-8",
                )

                output.append(
                    f"Found in ToolUniverse catalog: {best_tool.get('name', '')}\n"
                    f"  Description: {best_tool.get('description', '')}\n\n"
                    f"Generated SKILL.md → {skill_md_path}\n\n"
                    f"--- SKILL.md ---\n{skill_content}\n--- End ---\n\n"
                    f"Skill acquired from ToolUniverse. Read and follow the instructions above."
                )
                return "\n".join(output)
        except Exception as e:
            output.append(f"ToolUniverse search failed: {e}\n")

    # Step 4: Not found — provide instructions to create
    if not matches:
        output.append(f"No marketplace character found with skill '{skill_name}'.\n")

    output.append(
        f"\n--- Create This Skill ---\n"
        f"Option 1: Use the Character Builder UI at /character-builder\n"
        f"Option 2: Call autolab_create_character with skills='{skill_name}'\n"
        f"Option 3: Manually create the skill:\n"
        f"  mkdir -p .autolab/acquired_skills/{skill_name}\n"
        f"  # Write a SKILL.md with: When to use, Standard workflow, Key decisions\n"
        f"\nAfter creating, the skill will be available for all future sessions."
    )
    return "\n".join(output)


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
