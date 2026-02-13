"""
Autonomous Lab - Core lab module

State management, prompt engineering, and LaTeX templates
for iterative PI/Trainee research sessions.
"""

from .state import (
    append_meeting_log,
    compress_old_meetings,
    drain_feedback_queue,
    get_paper_progress,
    has_pending_feedback,
    init_project,
    load_state,
    parse_meeting_log,
    save_state,
    scan_project_files,
    store_user_feedback,
)

__all__ = [
    "append_meeting_log",
    "compress_old_meetings",
    "drain_feedback_queue",
    "get_paper_progress",
    "has_pending_feedback",
    "init_project",
    "load_state",
    "parse_meeting_log",
    "save_state",
    "scan_project_files",
    "store_user_feedback",
]
