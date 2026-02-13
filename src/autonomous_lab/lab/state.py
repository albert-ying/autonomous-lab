"""
Autonomous Lab - State Management

Manages the .autolab/ directory, state.json, meeting_log.md,
file scanning, and meeting summary compression.
"""

import json
import os
from datetime import datetime, timezone
from pathlib import Path

import yaml


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
AUTOLAB_DIR = ".autolab"
STATE_FILE = "state.json"
MEETING_LOG = "meeting_log.md"
MEETING_SUMMARIES = "meeting_summaries.md"
IDEA_FILE = "idea.md"
CONFIG_FILE = "config.yaml"
PROFILES_DIR = "profiles"
PI_PROFILE = "pi.yaml"
TRAINEE_PROFILE = "trainee.yaml"

# Compress meetings every N turns
COMPRESS_EVERY = 6

# Default profiles
DEFAULT_PI_PROFILE = {
    "title": "Principal Investigator",
    "expertise": "running a world-class science research lab",
    "goal": "perform research that maximizes scientific impact and produces "
    "publication-ready papers for top-tier journals (Nature, Science, Cell)",
    "personality": [
        "Visionary: identifies the 2-3 analyses that tell a compelling story",
        "Critical: demands statistical rigor and reproducibility",
        "Publication-focused: every figure must meet top-journal standards",
        "Efficient: prioritizes high-impact work over exhaustive analysis",
    ],
}

DEFAULT_TRAINEE_PROFILE = {
    "title": "Postdoctoral Researcher",
    "expertise": "technical implementation, data analysis, and scientific writing",
    "goal": "execute the PI's research vision with technical excellence, "
    "produce clean code, publication-quality figures, and clear LaTeX sections",
    "personality": [
        "Dedicated: completes tasks thoroughly and on time",
        "Technical: writes self-contained, reproducible code",
        "Proactive: identifies additional analyses aligned with PI vision",
        "Clear communicator: reports results with proper statistics",
    ],
}


# ---------------------------------------------------------------------------
# Project initialization
# ---------------------------------------------------------------------------
def init_project(
    project_dir: str,
    idea: str,
    pi_profile: dict | None = None,
    trainee_profile: dict | None = None,
    config: dict | None = None,
) -> dict:
    """
    Initialize a new Autonomous Lab project.

    Creates .autolab/ directory tree, paper/ structure, and working dirs.
    Returns the initial state dict.
    """
    base = Path(project_dir)
    autolab = base / AUTOLAB_DIR

    # Create directory structure
    (autolab / PROFILES_DIR).mkdir(parents=True, exist_ok=True)
    for d in ("data", "scripts", "figures", "results"):
        (base / d).mkdir(exist_ok=True)

    # Write idea
    (autolab / IDEA_FILE).write_text(idea.strip() + "\n", encoding="utf-8")

    # Write profiles
    pi = pi_profile or DEFAULT_PI_PROFILE
    trainee = trainee_profile or DEFAULT_TRAINEE_PROFILE
    with open(autolab / PROFILES_DIR / PI_PROFILE, "w", encoding="utf-8") as f:
        yaml.dump(pi, f, default_flow_style=False, allow_unicode=True)
    with open(autolab / PROFILES_DIR / TRAINEE_PROFILE, "w", encoding="utf-8") as f:
        yaml.dump(trainee, f, default_flow_style=False, allow_unicode=True)

    # Write config
    cfg = config or {"max_iterations": 50, "target_journal": "Nature"}
    with open(autolab / CONFIG_FILE, "w", encoding="utf-8") as f:
        yaml.dump(cfg, f, default_flow_style=False)

    # Initialize state
    now = datetime.now(timezone.utc).isoformat()
    state = {
        "iteration": 0,
        "next_role": "pi",
        "status": "active",
        "user_feedback": "",          # legacy compat (read by older code paths)
        "feedback_queue": [],         # accumulates messages until drained
        "progress": 0,
        "experts": [],
        "created_at": now,
        "last_updated": now,
    }
    save_state(project_dir, state)

    # Initialize empty meeting log and summaries
    (autolab / MEETING_LOG).write_text(
        "# Autonomous Lab Meeting Log\n\n", encoding="utf-8"
    )
    (autolab / MEETING_SUMMARIES).write_text(
        "# Meeting Summaries\n\n", encoding="utf-8"
    )

    return state


# ---------------------------------------------------------------------------
# State I/O
# ---------------------------------------------------------------------------
def _autolab_path(project_dir: str) -> Path:
    return Path(project_dir) / AUTOLAB_DIR


def load_state(project_dir: str) -> dict:
    """Load state.json from .autolab/"""
    state_path = _autolab_path(project_dir) / STATE_FILE
    if not state_path.exists():
        raise FileNotFoundError(
            f"No .autolab/state.json found in {project_dir}. "
            "Run autolab_init first."
        )
    with open(state_path, encoding="utf-8") as f:
        return json.load(f)


def save_state(project_dir: str, state: dict) -> None:
    """Write state.json to .autolab/"""
    state["last_updated"] = datetime.now(timezone.utc).isoformat()
    state_path = _autolab_path(project_dir) / STATE_FILE
    with open(state_path, "w", encoding="utf-8") as f:
        json.dump(state, f, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------------------
# Meeting log
# ---------------------------------------------------------------------------
def append_meeting_log(
    project_dir: str,
    role: str,
    iteration: int,
    summary: str,
    content: str,
) -> None:
    """Append a structured entry to meeting_log.md"""
    log_path = _autolab_path(project_dir) / MEETING_LOG
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    entry = (
        f"\n---\n\n"
        f"## Iteration {iteration} — {role.upper()} Turn\n"
        f"**Date:** {now}\n\n"
        f"### Summary\n{summary}\n\n"
        f"### Details\n{content}\n"
    )
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(entry)


def get_recent_meetings(project_dir: str, n: int = 3) -> str:
    """Return the last n meeting entries from meeting_log.md"""
    log_path = _autolab_path(project_dir) / MEETING_LOG
    if not log_path.exists():
        return ""
    text = log_path.read_text(encoding="utf-8")
    # Split on the separator
    entries = text.split("\n---\n")
    recent = entries[-n:] if len(entries) > n else entries
    return "\n---\n".join(recent)


def get_meeting_summaries(project_dir: str) -> str:
    """Return compressed meeting summaries"""
    path = _autolab_path(project_dir) / MEETING_SUMMARIES
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8")


def compress_old_meetings(project_dir: str, keep_recent: int = 3) -> None:
    """
    Compress older meeting entries into meeting_summaries.md.

    Keeps the most recent `keep_recent` entries intact in meeting_log.md.
    Older entries are condensed to one-line summaries and appended to
    meeting_summaries.md. The full log is then truncated to recent entries.
    """
    log_path = _autolab_path(project_dir) / MEETING_LOG
    summaries_path = _autolab_path(project_dir) / MEETING_SUMMARIES

    if not log_path.exists():
        return

    text = log_path.read_text(encoding="utf-8")
    parts = text.split("\n---\n")

    # First part is the header; rest are entries
    header = parts[0] if parts else "# Autonomous Lab Meeting Log\n\n"
    entries = parts[1:] if len(parts) > 1 else []

    if len(entries) <= keep_recent:
        return  # Nothing to compress

    old_entries = entries[: -keep_recent]
    recent_entries = entries[-keep_recent:]

    # Build compressed summaries from old entries
    compressed_lines = []
    for entry in old_entries:
        lines = entry.strip().split("\n")
        # Extract iteration header and summary
        iter_line = ""
        summary_line = ""
        for line in lines:
            if line.startswith("## Iteration"):
                iter_line = line.replace("## ", "").strip()
            if line.startswith("### Summary"):
                # Next non-empty line is the summary
                idx = lines.index(line)
                for subsequent in lines[idx + 1 :]:
                    if subsequent.strip():
                        summary_line = subsequent.strip()
                        break
        if iter_line:
            compressed_lines.append(f"- **{iter_line}**: {summary_line}")

    # Append to summaries file
    if compressed_lines:
        with open(summaries_path, "a", encoding="utf-8") as f:
            f.write("\n".join(compressed_lines) + "\n\n")

    # Rewrite meeting log with only recent entries
    with open(log_path, "w", encoding="utf-8") as f:
        f.write(header)
        for entry in recent_entries:
            f.write(f"\n---\n{entry}")


# ---------------------------------------------------------------------------
# File scanning
# ---------------------------------------------------------------------------
def scan_project_files(project_dir: str) -> dict:
    """
    Scan project directories and return file listings.

    Returns dict with keys: data, scripts, figures, results, paper
    Each value is a list of relative file paths.
    """
    base = Path(project_dir)
    result = {}
    for dirname in ("data", "scripts", "figures", "results", "paper"):
        dirpath = base / dirname
        if dirpath.exists():
            files = []
            for p in sorted(dirpath.rglob("*")):
                if p.is_file() and not p.name.startswith("."):
                    files.append(str(p.relative_to(base)))
            result[dirname] = files
        else:
            result[dirname] = []
    return result


def get_paper_progress(project_dir: str) -> dict:
    """
    Check which paper sections have content and their word counts.

    Returns dict like:
    {
        "abstract": {"exists": True, "words": 150},
        "introduction": {"exists": True, "words": 0},
        ...
    }
    """
    sections_dir = Path(project_dir) / "paper" / "sections"
    progress = {}
    for section_name in (
        "abstract",
        "introduction",
        "methods",
        "results",
        "discussion",
    ):
        tex_file = sections_dir / f"{section_name}.tex"
        if tex_file.exists():
            text = tex_file.read_text(encoding="utf-8")
            # Strip LaTeX comments and count words
            lines = [
                line
                for line in text.split("\n")
                if line.strip() and not line.strip().startswith("%")
            ]
            word_count = sum(len(line.split()) for line in lines)
            progress[section_name] = {"exists": True, "words": word_count}
        else:
            progress[section_name] = {"exists": False, "words": 0}
    return progress


# ---------------------------------------------------------------------------
# Profile loading
# ---------------------------------------------------------------------------
def load_idea(project_dir: str) -> str:
    """Load the project idea text"""
    path = _autolab_path(project_dir) / IDEA_FILE
    if path.exists():
        return path.read_text(encoding="utf-8")
    return ""


def load_profile(project_dir: str, role: str) -> dict:
    """Load a role profile (pi or trainee)"""
    filename = PI_PROFILE if role == "pi" else TRAINEE_PROFILE
    path = _autolab_path(project_dir) / PROFILES_DIR / filename
    if path.exists():
        with open(path, encoding="utf-8") as f:
            return yaml.safe_load(f) or {}
    return DEFAULT_PI_PROFILE if role == "pi" else DEFAULT_TRAINEE_PROFILE


def load_config(project_dir: str) -> dict:
    """Load project config"""
    path = _autolab_path(project_dir) / CONFIG_FILE
    if path.exists():
        with open(path, encoding="utf-8") as f:
            return yaml.safe_load(f) or {}
    return {}


# ---------------------------------------------------------------------------
# Meeting log parsing (for frontend)
# ---------------------------------------------------------------------------
def parse_meeting_log(project_dir: str) -> list[dict]:
    """
    Parse meeting_log.md into structured turns for the game UI.

    Returns a list of dicts:
        [{"iteration": 0, "role": "pi", "date": "...", "summary": "...", "content": "..."}, ...]
    """
    import re

    log_path = _autolab_path(project_dir) / MEETING_LOG
    if not log_path.exists():
        return []

    text = log_path.read_text(encoding="utf-8")
    sections = text.split("\n---\n")

    turns: list[dict] = []
    header_re = re.compile(
        r"##\s+Iteration\s+(\d+)\s+.+?\s+(PI|TRAINEE|REVIEWER\S*)\s+Turn",
        re.IGNORECASE,
    )
    date_re = re.compile(r"\*\*Date:\*\*\s*(.+)")
    summary_re = re.compile(r"###\s+Summary\s*\n([\s\S]*?)(?=###|\Z)")
    details_re = re.compile(r"###\s+Details\s*\n([\s\S]*)")

    for section in sections:
        hdr = header_re.search(section)
        if not hdr:
            continue
        iteration = int(hdr.group(1))
        role = hdr.group(2).lower()

        date_match = date_re.search(section)
        date_str = date_match.group(1).strip() if date_match else ""

        summary_match = summary_re.search(section)
        summary = summary_match.group(1).strip() if summary_match else ""

        details_match = details_re.search(section)
        content = details_match.group(1).strip() if details_match else ""

        turns.append(
            {
                "iteration": iteration,
                "role": role,
                "date": date_str,
                "summary": summary,
                "content": content,
            }
        )

    return turns


def store_user_feedback(project_dir: str, feedback: str, target: str = "") -> None:
    """Append user feedback to the queue in state.json.

    Messages accumulate in ``feedback_queue`` (a list) and are drained
    by ``drain_feedback_queue`` when the AI is ready to process them.
    This prevents messages from being lost when the user sends multiple
    messages before the current cycle completes.
    """
    state = load_state(project_dir)
    queue = state.get("feedback_queue", [])
    if not isinstance(queue, list):
        # Migrate from old string format
        queue = [queue] if queue else []

    entry = f"[To {target.upper()}] {feedback}" if target else feedback
    queue.append(entry)
    state["feedback_queue"] = queue
    # Also keep legacy field updated (last message) for compat
    state["user_feedback"] = entry
    save_state(project_dir, state)


def drain_feedback_queue(project_dir: str) -> str:
    """Drain and return all queued user feedback as a single string.

    Returns the combined feedback (empty string if none) and clears
    both ``feedback_queue`` and the legacy ``user_feedback`` field.
    """
    state = load_state(project_dir)
    queue = state.get("feedback_queue", [])
    if not isinstance(queue, list):
        queue = [queue] if queue else []

    # Also grab legacy field in case something wrote there directly
    legacy = state.get("user_feedback", "").strip()
    if legacy and legacy not in queue:
        queue.append(legacy)

    combined = "\n\n".join(msg.strip() for msg in queue if msg.strip())

    # Clear everything
    state["feedback_queue"] = []
    state["user_feedback"] = ""
    save_state(project_dir, state)
    return combined


def has_pending_feedback(project_dir: str) -> bool:
    """Check if there's any pending user feedback without draining it."""
    state = load_state(project_dir)
    queue = state.get("feedback_queue", [])
    if not isinstance(queue, list):
        return bool(queue)
    if queue:
        return True
    return bool(state.get("user_feedback", "").strip())


# ---------------------------------------------------------------------------
# Editorial workflow state
# ---------------------------------------------------------------------------
EDITORIAL_PHASES = (
    "none",               # Normal lab mode
    "submitted",          # PI submitted manuscript + cover letter
    "reviewers_invited",  # Editor invited reviewers
    "under_review",       # Reviewers generating reports
    "reviews_complete",   # All reviews in, awaiting editorial decision
    "decision_made",      # Editor made decision, returning to PI
)


def get_editorial(project_dir: str) -> dict:
    """Return the editorial state (initialising if absent)."""
    state = load_state(project_dir)
    return state.get("editorial", _default_editorial())


def _default_editorial() -> dict:
    return {
        "phase": "none",
        "cover_letter": "",
        "reviewers": [],
        "reviews": {},
        "decision": "",
        "decision_feedback": "",
        "round": 0,          # revision round counter
    }


def update_editorial(project_dir: str, **fields: object) -> dict:
    """Merge *fields* into editorial state and save."""
    state = load_state(project_dir)
    editorial = state.get("editorial", _default_editorial())
    editorial.update(fields)
    state["editorial"] = editorial
    save_state(project_dir, state)
    return editorial


def submit_manuscript(project_dir: str, cover_letter: str) -> dict:
    """PI submits manuscript.  Moves phase to 'submitted'."""
    state = load_state(project_dir)
    editorial = state.get("editorial", _default_editorial())
    editorial["phase"] = "submitted"
    editorial["cover_letter"] = cover_letter
    editorial["reviewers"] = []
    editorial["reviews"] = {}
    editorial["decision"] = ""
    editorial["decision_feedback"] = ""
    editorial["round"] = editorial.get("round", 0) + 1
    state["editorial"] = editorial
    state["status"] = "submitted_to_editor"
    save_state(project_dir, state)
    return editorial


def invite_reviewers(project_dir: str, reviewers: list[dict]) -> dict:
    """
    Editor invites reviewers.

    Each reviewer dict: {"id": "reviewer_1", "name": "Dr. X",
                         "role": "Immunologist", "avatar": "immunologist"}
    """
    state = load_state(project_dir)
    editorial = state.get("editorial", _default_editorial())
    editorial["phase"] = "reviewers_invited"
    editorial["reviewers"] = reviewers
    # Set next_role to first reviewer
    if reviewers:
        state["next_role"] = reviewers[0]["id"]
    state["editorial"] = editorial
    state["status"] = "under_review"
    save_state(project_dir, state)
    return editorial


def record_review(project_dir: str, reviewer_id: str, report: dict) -> dict:
    """
    Record a reviewer's report and advance to next reviewer or completion.

    report: {"score": 7, "recommendation": "minor_revision", "report": "..."}
    """
    state = load_state(project_dir)
    editorial = state.get("editorial", _default_editorial())
    editorial["reviews"][reviewer_id] = report

    # Determine next reviewer or complete
    reviewer_ids = [r["id"] for r in editorial.get("reviewers", [])]
    reviewed = set(editorial["reviews"].keys())
    remaining = [rid for rid in reviewer_ids if rid not in reviewed]

    if remaining:
        state["next_role"] = remaining[0]
        editorial["phase"] = "under_review"
    else:
        editorial["phase"] = "reviews_complete"
        state["status"] = "reviews_complete"
        state["next_role"] = "editor"  # signal UI to show decision panel

    state["editorial"] = editorial
    save_state(project_dir, state)
    return editorial


def record_editorial_decision(
    project_dir: str,
    decision: str,
    feedback: str,
) -> dict:
    """
    Editor makes a decision.

    decision: "accept" | "minor_revision" | "major_revision" | "reject"
    feedback: free-form editorial feedback text
    """
    state = load_state(project_dir)
    editorial = state.get("editorial", _default_editorial())
    editorial["phase"] = "decision_made"
    editorial["decision"] = decision
    editorial["decision_feedback"] = feedback
    state["editorial"] = editorial

    if decision == "reject":
        state["status"] = "rejected"
        state["next_role"] = "pi"
    elif decision == "accept":
        state["status"] = "accepted"
        state["next_role"] = "pi"
    else:
        # Revision requested → back to PI
        state["status"] = "revision_requested"
        state["next_role"] = "pi"

    save_state(project_dir, state)
    return editorial
