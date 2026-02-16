"""
Autonomous Lab - Structured Event Log

Append-only JSONL event log for comprehensive audit trail.
Every significant action produces a timestamped JSON line in
.autolab/event_log.jsonl.
"""

import json
import os
from datetime import datetime, timezone
from pathlib import Path

EVENT_LOG_FILE = "event_log.jsonl"


def _event_log_path(project_dir: str) -> Path:
    return Path(project_dir) / ".autolab" / EVENT_LOG_FILE


def log_event(project_dir: str, event: str, **kwargs: object) -> None:
    """Append a timestamped event to the JSONL log.

    Args:
        project_dir: Project directory path.
        event: Event type string (e.g. 'turn_start', 'consultation').
        **kwargs: Arbitrary event payload fields.
    """
    path = _event_log_path(project_dir)
    path.parent.mkdir(parents=True, exist_ok=True)

    entry = {
        "ts": datetime.now(timezone.utc).isoformat(),
        "event": event,
        **kwargs,
    }
    # Ensure JSON-serialisable: convert Path objects etc.
    safe = _make_serialisable(entry)
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps(safe, ensure_ascii=False) + "\n")


def query_events(
    project_dir: str,
    event_type: str | None = None,
    since: str | None = None,
    until: str | None = None,
    keyword: str | None = None,
    limit: int = 0,
) -> list[dict]:
    """Read and filter events from the log.

    Args:
        project_dir: Project directory.
        event_type: Filter by event type (exact match).
        since: ISO datetime lower bound (inclusive).
        until: ISO datetime upper bound (inclusive).
        keyword: Case-insensitive substring search across all values.
        limit: Max events to return (0 = unlimited). Returns most recent.

    Returns:
        List of matching event dicts, ordered chronologically.
    """
    path = _event_log_path(project_dir)
    if not path.exists():
        return []

    events: list[dict] = []
    kw_lower = keyword.lower() if keyword else None

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                entry = json.loads(line)
            except json.JSONDecodeError:
                continue

            if event_type and entry.get("event") != event_type:
                continue
            ts = entry.get("ts", "")
            if since and ts < since:
                continue
            if until and ts > until:
                continue
            if kw_lower:
                flat = json.dumps(entry, ensure_ascii=False).lower()
                if kw_lower not in flat:
                    continue
            events.append(entry)

    if limit > 0:
        events = events[-limit:]
    return events


def export_session_report(project_dir: str) -> str:
    """Generate a human-readable session summary from the event log.

    Returns a Markdown-formatted report suitable for display.
    """
    events = query_events(project_dir)
    if not events:
        return "No events recorded yet."

    lines = ["# Session Report\n"]

    # Summary stats
    event_counts: dict[str, int] = {}
    for e in events:
        et = e.get("event", "unknown")
        event_counts[et] = event_counts.get(et, 0) + 1

    lines.append(f"**Total events:** {len(events)}\n")
    lines.append("**Event breakdown:**")
    for et, count in sorted(event_counts.items()):
        lines.append(f"  - {et}: {count}")
    lines.append("")

    # Timeline of key events
    lines.append("## Timeline\n")
    for e in events:
        ts = e.get("ts", "?")[:19].replace("T", " ")
        et = e.get("event", "?")
        detail = _event_one_liner(e)
        lines.append(f"- `{ts}` **{et}** {detail}")

    return "\n".join(lines)


def diff_file_listings(
    before: dict[str, list[str]],
    after: dict[str, list[str]],
) -> dict[str, list[str]]:
    """Compute file changes between two scan_project_files() results.

    Returns dict with keys 'added', 'modified' (present in both but
    we can't detect content changes without hashing, so this is
    just set difference), and 'deleted'.
    """
    before_set: set[str] = set()
    after_set: set[str] = set()
    for files in before.values():
        before_set.update(files)
    for files in after.values():
        after_set.update(files)

    return {
        "added": sorted(after_set - before_set),
        "deleted": sorted(before_set - after_set),
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _event_one_liner(event: dict) -> str:
    """Produce a short human-readable description of an event."""
    et = event.get("event", "")
    if et == "turn_start":
        return f"Iteration {event.get('iteration', '?')}, role={event.get('role', '?')}"
    if et == "turn_end":
        summary = event.get("summary", "")[:80]
        return f"Iteration {event.get('iteration', '?')}, role={event.get('role', '?')}: {summary}"
    if et == "consultation":
        return f"{event.get('expert', '?')} ({event.get('expert_role', '?')}): {event.get('question', '')[:60]}"
    if et == "user_feedback":
        content = event.get("content", "")[:60]
        return f"from {event.get('source', 'ui')}: {content}"
    if et == "editorial_decision":
        return f"decision={event.get('decision', '?')}, round={event.get('round', '?')}"
    if et == "skill_acquired":
        return f"skill={event.get('skill', '?')}, source={event.get('source', '?')}"
    if et == "file_change":
        added = event.get("added", [])
        deleted = event.get("deleted", [])
        parts = []
        if added:
            parts.append(f"+{len(added)} files")
        if deleted:
            parts.append(f"-{len(deleted)} files")
        return ", ".join(parts) if parts else "no changes"
    if et == "error":
        return event.get("message", "")[:80]
    # Generic fallback
    keys = [k for k in event if k not in ("ts", "event")]
    if keys:
        return ", ".join(f"{k}={event[k]}" for k in keys[:3])
    return ""


def _make_serialisable(obj: object) -> object:
    """Recursively convert non-serialisable types to strings."""
    if isinstance(obj, dict):
        return {k: _make_serialisable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_make_serialisable(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return str(obj)
