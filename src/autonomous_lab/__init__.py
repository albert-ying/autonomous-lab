"""
Autonomous Lab
==============

MCP server for iterative PI/Trainee research sessions.
Fork of mcp-feedback-enhanced with virtual lab extensions.
"""

__version__ = "0.5.2"
__author__ = "Albert Ying"

from .server import main as run_server
from .web import WebUIManager, get_web_ui_manager, launch_web_feedback_ui, stop_web_ui

__all__ = [
    "WebUIManager",
    "__author__",
    "__version__",
    "get_web_ui_manager",
    "launch_web_feedback_ui",
    "run_server",
    "stop_web_ui",
]


def main():
    """Main entry point"""
    from .__main__ import main as cli_main

    return cli_main()
