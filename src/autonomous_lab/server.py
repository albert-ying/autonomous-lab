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
async def interactive_feedback(
    project_directory: Annotated[str, Field(description="專案目錄路徑")] = ".",
    summary: Annotated[
        str, Field(description="AI 工作完成的摘要說明")
    ] = "我已完成了您請求的任務。",
    timeout: Annotated[int, Field(description="等待用戶回饋的超時時間（秒）")] = 600,
) -> list:
    """Interactive feedback collection tool for LLM agents.

    USAGE RULES:
    1. During any process, task, or conversation, whether asking, replying, or completing phased tasks, you must call this tool to ask for feedback.
    2. Unless receiving termination instructions, all steps must repeatedly call this tool.
    3. Whenever user feedback is received, if the feedback content is not empty, you must call this tool again and adjust behavior based on the feedback content.
    4. Only when the user explicitly indicates "end" or "no more interaction needed" can you stop calling this tool, and the process is considered complete.
    5. You should summarize what have done, and provide project directory through args to let user know what you have done to provide feedback for next step.

    Args:
        project_directory: Project directory path for context
        summary: Summary of AI work completed for user review
        timeout: Timeout in seconds for waiting user feedback (default: 600 seconds)

    Returns:
        list: List containing TextContent and MCPImage objects representing user feedback
    """
    # 環境偵測
    is_remote = is_remote_environment()
    is_wsl = is_wsl_environment()

    debug_log(f"環境偵測結果 - 遠端: {is_remote}, WSL: {is_wsl}")
    debug_log("使用介面: Web UI")

    try:
        # 確保專案目錄存在
        if not os.path.exists(project_directory):
            project_directory = os.getcwd()
        project_directory = os.path.abspath(project_directory)

        # 使用 Web 模式
        debug_log("回饋模式: web")

        result = await launch_web_feedback_ui(project_directory, summary, timeout)

        # 處理取消情況
        if not result:
            return [TextContent(type="text", text="用戶取消了回饋。")]

        # 儲存詳細結果
        save_feedback_to_file(result)

        # 建立回饋項目列表
        feedback_items = []

        # 添加文字回饋
        if (
            result.get("interactive_feedback")
            or result.get("command_logs")
            or result.get("images")
        ):
            feedback_text = create_feedback_text(result)
            feedback_items.append(TextContent(type="text", text=feedback_text))
            debug_log("文字回饋已添加")

        # 添加圖片回饋
        if result.get("images"):
            mcp_images = process_images(result["images"])
            # 修復 arg-type 錯誤 - 直接擴展列表
            feedback_items.extend(mcp_images)
            debug_log(f"已添加 {len(mcp_images)} 張圖片")

        # 確保至少有一個回饋項目
        if not feedback_items:
            feedback_items.append(
                TextContent(type="text", text="用戶未提供任何回饋內容。")
            )

        debug_log(f"回饋收集完成，共 {len(feedback_items)} 個項目")
        return feedback_items

    except Exception as e:
        # 使用統一錯誤處理，但不影響 JSON RPC 響應
        error_id = ErrorHandler.log_error_with_context(
            e,
            context={"operation": "回饋收集", "project_dir": project_directory},
            error_type=ErrorType.SYSTEM,
        )

        # 生成用戶友好的錯誤信息
        user_error_msg = ErrorHandler.format_user_error(e, include_technical=False)
        debug_log(f"回饋收集錯誤 [錯誤ID: {error_id}]: {e!s}")

        return [TextContent(type="text", text=user_error_msg)]


async def launch_web_feedback_ui(project_dir: str, summary: str, timeout: int) -> dict:
    """
    啟動 Web UI 收集回饋，支援自訂超時時間

    Args:
        project_dir: 專案目錄路徑
        summary: AI 工作摘要
        timeout: 超時時間（秒）

    Returns:
        dict: 收集到的回饋資料
    """
    debug_log(f"啟動 Web UI 介面，超時時間: {timeout} 秒")

    try:
        # 使用新的 web 模組
        from .web import launch_web_feedback_ui as web_launch

        # 傳遞 timeout 參數給 Web UI
        return await web_launch(project_dir, summary, timeout)
    except ImportError as e:
        # 使用統一錯誤處理
        error_id = ErrorHandler.log_error_with_context(
            e,
            context={"operation": "Web UI 模組導入", "module": "web"},
            error_type=ErrorType.DEPENDENCY,
        )
        user_error_msg = ErrorHandler.format_user_error(
            e, ErrorType.DEPENDENCY, include_technical=False
        )
        debug_log(f"Web UI 模組導入失敗 [錯誤ID: {error_id}]: {e}")

        return {
            "command_logs": "",
            "interactive_feedback": user_error_msg,
            "images": [],
        }


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

@mcp.tool()
async def autolab_init(
    project_directory: Annotated[str, Field(description="Project directory path")] = ".",
    idea: Annotated[str, Field(description="Research project idea and context")] = "",
    title: Annotated[str, Field(description="Paper title")] = "TITLE",
) -> str:
    """Initialize an Autonomous Lab project.

    Creates .autolab/ directory with state, profiles, and meeting log.
    Creates paper/ directory with LaTeX template.
    Creates working directories: data/, scripts/, figures/, results/.

    Call this once at the start of a new research project.

    Args:
        project_directory: Path to the project directory
        idea: The research idea, context, and any preliminary data description
        title: Working title for the paper

    Returns:
        str: Confirmation message with instructions to call autolab_next
    """
    from .lab.latex_template import create_paper_structure
    from .lab.state import init_project

    project_directory = os.path.abspath(project_directory)

    if not idea.strip():
        return "ERROR: 'idea' parameter is required. Provide the research idea and context."

    try:
        state = init_project(project_directory, idea)
        create_paper_structure(project_directory, title)

        return (
            f"Autonomous Lab initialized in {project_directory}\n\n"
            f"Created:\n"
            f"  .autolab/ -- state, profiles, meeting log\n"
            f"  paper/ -- LaTeX template ({title})\n"
            f"  scripts/, figures/, results/ -- working directories\n\n"
            f"State: iteration={state['iteration']}, next_role={state['next_role']}\n\n"
            f"Next step: Call autolab_next to get the PI's first prompt."
        )
    except Exception as e:
        return f"ERROR initializing project: {e}"


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

    Call flow: autolab_next -> (AI acts as role) -> autolab_record -> autolab_next -> ...

    Args:
        project_directory: Path to the project directory

    Returns:
        str: The role prompt to follow
    """
    from .lab.prompts import build_pi_prompt, build_reviewer_prompt, build_trainee_prompt
    from .lab.state import (
        get_meeting_summaries,
        get_paper_progress,
        get_recent_meetings,
        load_idea,
        load_profile,
        load_state,
        scan_project_files,
    )

    project_directory = os.path.abspath(project_directory)

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

    # Check if paper review is needed
    if state.get("status") == "ready_for_review":
        paper_progress = get_paper_progress(project_directory)
        prompt = build_reviewer_prompt(paper_progress, file_listings, meeting_history)
        return (
            f"[AUTOLAB] REVIEWER MODE -- Paper Review\n"
            f"Iteration: {iteration}\n\n"
            f"{prompt}"
        )

    # Build role-specific prompt
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
        f"{prompt}"
    )


@mcp.tool()
async def autolab_record(
    project_directory: Annotated[str, Field(description="Project directory path")] = ".",
    role: Annotated[str, Field(description="Role that just acted: 'pi' or 'trainee'")] = "",
    summary: Annotated[str, Field(description="Brief summary of what was done this turn")] = "",
    content: Annotated[str, Field(description="Full content of the turn (review, agenda, results, etc.)")] = "",
    status: Annotated[str, Field(description="Status: 'continue' or 'ready_for_review'")] = "continue",
    timeout: Annotated[int, Field(description="Timeout in seconds for user feedback")] = 600,
) -> str:
    """Record a completed turn and show the monitoring UI for user feedback.

    This tool:
    1. Saves the turn to meeting_log.md
    2. Updates state.json (advance iteration, flip role, store status)
    3. Compresses old meetings every 6 turns
    4. Opens the web feedback UI with a rich dashboard summary
    5. Waits for user feedback (or timeout = auto-continue)
    6. Stores user feedback in state.json for the next turn

    Call this AFTER completing your role actions (as PI or Trainee).

    Args:
        project_directory: Path to the project directory
        role: Which role just completed: 'pi' or 'trainee'
        summary: Brief summary of actions taken
        content: Full structured content of the turn
        status: 'continue' to keep iterating, 'ready_for_review' when paper is done
        timeout: Seconds to wait for user feedback (default 600)

    Returns:
        str: User feedback text, or empty string if auto-continued
    """
    from .lab.state import (
        append_meeting_log,
        compress_old_meetings,
        get_paper_progress,
        load_state,
        save_state,
        scan_project_files,
        COMPRESS_EVERY,
    )

    project_directory = os.path.abspath(project_directory)

    if role not in ("pi", "trainee"):
        return "ERROR: role must be 'pi' or 'trainee'"

    try:
        state = load_state(project_directory)
    except FileNotFoundError as e:
        return f"ERROR: {e}"

    iteration = state["iteration"]

    # 1. Append to meeting log
    append_meeting_log(project_directory, role, iteration, summary, content)

    # 2. Update state
    next_role = "trainee" if role == "pi" else "pi"
    new_iteration = iteration + (1 if role == "trainee" else 0)
    state["iteration"] = new_iteration
    state["next_role"] = next_role
    state["status"] = status
    state["user_feedback"] = ""  # Clear previous feedback
    save_state(project_directory, state)

    # 3. Compress old meetings periodically
    if new_iteration > 0 and new_iteration % COMPRESS_EVERY == 0:
        compress_old_meetings(project_directory)

    # 4. Build dashboard summary for the web UI
    file_listings = scan_project_files(project_directory)
    paper_progress = get_paper_progress(project_directory)

    # Build a rich markdown summary for the feedback UI
    figures_list = file_listings.get("figures", [])
    figures_str = "\n".join(f"- {f}" for f in figures_list[:20]) if figures_list else "None yet"

    paper_str_parts = []
    for section, info in paper_progress.items():
        if info["exists"] and info["words"] > 0:
            paper_str_parts.append(f"  {section}: {info['words']} words")
        else:
            paper_str_parts.append(f"  {section}: --")
    paper_str = "\n".join(paper_str_parts)

    dashboard_summary = (
        f"[AUTOLAB]\n\n"
        f"# Autonomous Lab — Iteration {new_iteration}\n\n"
        f"**Completed:** {role.upper()} turn\n"
        f"**Next:** {next_role.upper()} turn\n"
        f"**Status:** {status}\n\n"
        f"## Turn Summary\n\n{summary}\n\n"
        f"## Figures\n\n{figures_str}\n\n"
        f"## Paper Progress\n\n{paper_str}\n\n"
        f"---\n\n"
        f"*Provide feedback below, or leave empty and submit to continue.*"
    )

    # 5. Launch the feedback UI (reusing mcp-feedback-enhanced's web UI)
    try:
        result = await launch_web_feedback_ui(project_directory, dashboard_summary, timeout)

        user_feedback_text = ""
        if result and result.get("interactive_feedback"):
            user_feedback_text = result["interactive_feedback"].strip()

        # 6. Store feedback in state
        state = load_state(project_directory)
        state["user_feedback"] = user_feedback_text
        save_state(project_directory, state)

        if user_feedback_text:
            return f"User feedback received:\n\n{user_feedback_text}"
        else:
            return "No user feedback. Continuing to next turn."

    except TimeoutError:
        return "Feedback timeout. Continuing to next turn."
    except Exception as e:
        debug_log(f"Feedback UI error: {e}")
        return f"Feedback UI error: {e}. Continuing to next turn."


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
