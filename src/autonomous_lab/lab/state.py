"""
Autonomous Lab - State Management

Manages the .autolab/ directory, state.json, meeting_log.md,
file scanning, and meeting summary compression.
"""

import importlib.resources
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
MEETING_PREAMBLE = "meeting_preamble.md"
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
    "skills": ["scientific-analysis-review"],
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
    domain: str = "research",
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

    # Write config (include domain)
    domain_cfg = DOMAIN_CONFIGS.get(domain, DOMAIN_CONFIGS["research"])
    default_cfg = {
        "domain": domain,
        "max_iterations": 50,
        "target_venue": domain_cfg.get("target_venue", "Nature"),
        "editor_timeout_minutes": 30,  # 0 = no timeout (wait forever)
    }
    if config:
        default_cfg.update(config)
    with open(autolab / CONFIG_FILE, "w", encoding="utf-8") as f:
        yaml.dump(default_cfg, f, default_flow_style=False)

    # Store title in config for Editorial Center display
    title_str = default_cfg.get("title", "")
    if not title_str or title_str == "TITLE":
        # Derive title from first line of idea (up to 80 chars)
        first_line = idea.strip().split("\n")[0][:80].strip()
        title_str = first_line if first_line else Path(project_dir).name
        default_cfg["title"] = title_str
        # Rewrite config with derived title
        with open(autolab / CONFIG_FILE, "w", encoding="utf-8") as f:
            yaml.dump(default_cfg, f, default_flow_style=False)

    # Initialize state
    now = datetime.now(timezone.utc).isoformat()
    state = {
        "title": title_str,
        "iteration": 0,
        "next_role": "pi",
        "status": "active",
        "user_feedback": "",  # legacy compat (read by older code paths)
        "feedback_queue": [],  # accumulates messages until drained
        "progress": 0,
        "experts": [],
        "phase": "recruiting",
        "recruitment": {
            "characters": [],
            "started_at": now,
            "ready_at": None,
        },
        "characters": [],
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

    # Deploy bundled skills to acquired_skills/
    _deploy_bundled_skills(project_dir)

    return state


def _deploy_bundled_skills(project_dir: str) -> None:
    """Copy bundled SKILL.md files from the package to .autolab/acquired_skills/."""
    dest_root = Path(project_dir) / AUTOLAB_DIR / "acquired_skills"
    try:
        # Try importlib.resources first (works for installed packages)
        skills_pkg = importlib.resources.files("autonomous_lab") / "skills"
        if skills_pkg.is_dir():
            _copy_skills_from(skills_pkg, dest_root)
            return
    except Exception:
        pass
    try:
        # Fallback: resolve relative to this file (works for dev/editable installs)
        pkg_root = Path(__file__).resolve().parent.parent  # autonomous_lab/
        skills_dir = pkg_root / "skills"
        if skills_dir.is_dir():
            _copy_skills_from(skills_dir, dest_root)
    except Exception:
        pass  # Non-critical: skills can still be acquired at runtime


def _copy_skills_from(skills_dir: Path, dest_root: Path) -> None:
    """Copy SKILL.md files from a skills directory to dest_root."""
    for skill_dir in skills_dir.iterdir():
        if not skill_dir.is_dir():
            continue
        skill_md = skill_dir / "SKILL.md"
        if not skill_md.is_file():
            continue
        target = dest_root / skill_dir.name
        target.mkdir(parents=True, exist_ok=True)
        target_file = target / "SKILL.md"
        if not target_file.exists():
            target_file.write_text(
                skill_md.read_text(encoding="utf-8"), encoding="utf-8"
            )


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
            f"No .autolab/state.json found in {project_dir}. " "Run autolab_init first."
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


def get_meeting_summaries(project_dir: str, max_chars: int = 4000) -> str:
    """Return pinned project preamble + recent meeting summaries.

    The preamble (~500 chars) captures the project origin, key decisions,
    and current state so the model never loses the narrative thread.
    The remaining budget is filled with the most recent summary entries.
    """
    autolab = _autolab_path(project_dir)

    # Load preamble (pinned context)
    preamble_path = autolab / MEETING_PREAMBLE
    preamble = ""
    if preamble_path.exists():
        preamble = preamble_path.read_text(encoding="utf-8").strip()

    # Load summaries
    path = autolab / MEETING_SUMMARIES
    if not path.exists():
        if preamble:
            return f"## Project Narrative\n{preamble}"
        return ""
    text = path.read_text(encoding="utf-8")

    # Budget: reserve space for preamble + section headers + buffer
    preamble_overhead = len(preamble) + 60 if preamble else 0  # headers + newlines
    summary_budget = max_chars - preamble_overhead

    if len(text) <= summary_budget:
        capped_summaries = text
    else:
        # Keep the header + most recent summaries within budget
        lines = text.strip().split("\n")
        header = lines[0] if lines else ""
        entries = lines[1:]
        kept = []
        total_len = len(header) + 100  # buffer for truncation note
        for line in reversed(entries):
            if total_len + len(line) + 1 > summary_budget:
                break
            kept.insert(0, line)
            total_len += len(line) + 1

        truncated = len(entries) - len(kept)
        note = (
            f"\n*(Oldest {truncated} summary entries omitted for context budget)*\n\n"
            if truncated > 0
            else "\n"
        )
        capped_summaries = header + note + "\n".join(kept)

    if preamble:
        return f"## Project Narrative\n{preamble}\n\n## Recent History\n{capped_summaries}"
    return capped_summaries


def _generate_preamble(summaries_text: str, idea_text: str = "") -> str:
    """Build a short project narrative from summaries for context retention.

    Extracts project origin (first summary), key pivot/decision entries,
    and current state (last summary) into a ~300-500 char block.
    """
    lines = [
        ln.strip()
        for ln in summaries_text.strip().split("\n")
        if ln.strip().startswith("- **")
    ]
    if not lines:
        return ""

    # Project origin: first summary entry, truncated
    origin = lines[0][:150]

    # Key decisions: entries containing strategic keywords
    decision_keywords = ("pivot", "strategic decision", "major revision", "restructur")
    decisions = []
    for ln in lines:
        if any(kw in ln.lower() for kw in decision_keywords):
            decisions.append(ln[:100])
        if len(decisions) >= 3:
            break

    # Current state: last entry
    current = lines[-1][:150] if len(lines) > 1 else ""

    parts = [f"**Project origin**: {origin}"]
    if decisions:
        parts.append("**Key decisions**: " + " | ".join(decisions))
    if current:
        parts.append(f"**Current state**: {current}")

    return "\n".join(parts)


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

    old_entries = entries[:-keep_recent]
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

    # Generate/update the meeting preamble for context retention
    preamble_path = _autolab_path(project_dir) / MEETING_PREAMBLE
    idea_text = load_idea(project_dir)
    all_summaries = summaries_path.read_text(encoding="utf-8") if summaries_path.exists() else ""
    preamble = _generate_preamble(all_summaries, idea_text)
    if preamble:
        preamble_path.write_text(preamble + "\n", encoding="utf-8")


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
# Domain configuration — maps domain IDs to role/artifact labels
# ---------------------------------------------------------------------------
DOMAIN_CONFIGS = {
    "research": {
        "senior_label": "Professor",
        "senior_short": "PI",
        "junior_label": "Trainee",
        "junior_short": "Trainee",
        "overseer_label": "Editor",
        "reviewer_label": "Reviewer",
        "artifact": "Paper",
        "artifact_plural": "Papers",
        "workspace": "Lab",
        "meeting_log_name": "Meeting Log",
        "review_process": "Peer Review",
        "consultant_label": "Consultant",
        "decisions": ["Accept", "Minor Revision", "Major Revision", "Reject"],
        "target_venue": "Nature",
    },
    "software": {
        "senior_label": "Tech Lead",
        "senior_short": "Lead",
        "junior_label": "Developer",
        "junior_short": "Dev",
        "overseer_label": "Code Reviewer",
        "reviewer_label": "Reviewer",
        "artifact": "Pull Request",
        "artifact_plural": "Pull Requests",
        "workspace": "Sprint",
        "meeting_log_name": "Standup Log",
        "review_process": "Code Review",
        "consultant_label": "Specialist",
        "decisions": ["Approve", "Request Changes", "Needs Discussion", "Reject"],
        "target_venue": "Production",
    },
    "consulting": {
        "senior_label": "Partner",
        "senior_short": "Partner",
        "junior_label": "Associate",
        "junior_short": "Associate",
        "overseer_label": "Client",
        "reviewer_label": "Advisor",
        "artifact": "Strategy Report",
        "artifact_plural": "Reports",
        "workspace": "Engagement",
        "meeting_log_name": "Engagement Log",
        "review_process": "Client Review",
        "consultant_label": "Subject Expert",
        "decisions": ["Approve", "Revise", "Rework", "Reject"],
        "target_venue": "Board Presentation",
    },
    "legal": {
        "senior_label": "Senior Partner",
        "senior_short": "Partner",
        "junior_label": "Associate",
        "junior_short": "Associate",
        "overseer_label": "Reviewing Partner",
        "reviewer_label": "Counsel",
        "artifact": "Legal Brief",
        "artifact_plural": "Briefs",
        "workspace": "Case",
        "meeting_log_name": "Case Log",
        "review_process": "Partner Review",
        "consultant_label": "Specialist",
        "decisions": ["File", "Minor Edits", "Major Rewrite", "Withdraw"],
        "target_venue": "Court Filing",
    },
    "medical": {
        "senior_label": "Attending",
        "senior_short": "Attending",
        "junior_label": "Resident",
        "junior_short": "Resident",
        "overseer_label": "Chief",
        "reviewer_label": "Reviewer",
        "artifact": "Treatment Plan",
        "artifact_plural": "Plans",
        "workspace": "Rounds",
        "meeting_log_name": "Rounds Log",
        "review_process": "Peer Review",
        "consultant_label": "Specialist",
        "decisions": ["Approve", "Modify Plan", "Reassess", "Override"],
        "target_venue": "Clinical Practice",
    },
    "creative": {
        "senior_label": "Art Director",
        "senior_short": "AD",
        "junior_label": "Designer",
        "junior_short": "Designer",
        "overseer_label": "Creative Director",
        "reviewer_label": "Reviewer",
        "artifact": "Campaign",
        "artifact_plural": "Campaigns",
        "workspace": "Studio",
        "meeting_log_name": "Creative Brief",
        "review_process": "Creative Review",
        "consultant_label": "Specialist",
        "decisions": ["Approve", "Polish", "Rework", "Scrap"],
        "target_venue": "Client Presentation",
    },
}


def get_domain_config(project_dir: str) -> dict:
    """
    Get domain-specific labels for a project.

    Reads the 'domain' key from config.yaml and returns the corresponding
    label dictionary. Defaults to 'research' for backward compatibility.
    """
    cfg = load_config(project_dir)
    domain = cfg.get("domain", "research")
    return DOMAIN_CONFIGS.get(domain, DOMAIN_CONFIGS["research"])


def get_domain_id(project_dir: str) -> str:
    """Get the domain ID for a project (default: 'research')."""
    cfg = load_config(project_dir)
    return cfg.get("domain", "research")


def get_editor_timeout_seconds(project_dir: str) -> int:
    """Return the editor timeout in seconds (0 = wait forever)."""
    cfg = load_config(project_dir)
    minutes = cfg.get("editor_timeout_minutes", 30)
    return max(0, int(minutes)) * 60


# ---------------------------------------------------------------------------
# Character roster and recruitment
# ---------------------------------------------------------------------------
def add_character(project_dir: str, char_data: dict) -> None:
    """Add a character to the recruitment roster in state.json."""
    state = load_state(project_dir)
    roster = state.get("recruitment", {}).get("characters", [])
    roster = [c for c in roster if c.get("slug") != char_data.get("slug")]
    roster.append(char_data)
    state["recruitment"]["characters"] = roster
    save_state(project_dir, state)


def update_skill_status(
    project_dir: str, slug: str, skill: str, status: str, **meta
) -> None:
    """Update a skill's status for a character in the roster."""
    state = load_state(project_dir)
    for char in state.get("recruitment", {}).get("characters", []):
        if char.get("slug") == slug:
            if skill in char.get("skills", {}):
                char["skills"][skill]["status"] = status
                char["skills"][skill].update(meta)
            break
    save_state(project_dir, state)


def check_recruitment_ready(project_dir: str) -> bool:
    """Check if all required skills are certified. Transition phase if so."""
    state = load_state(project_dir)
    characters = state.get("recruitment", {}).get("characters", [])
    if not characters:
        return False
    all_ready = True
    for char in characters:
        required_certified = all(
            skill_info.get("status") == "certified"
            for skill_info in char.get("skills", {}).values()
            if skill_info.get("required", False)
        )
        char["ready"] = required_certified
        if not required_certified:
            all_ready = False
    if all_ready:
        state["phase"] = "active"
        state["recruitment"]["ready_at"] = datetime.now(timezone.utc).isoformat()
    save_state(project_dir, state)
    return all_ready


def get_active_characters(project_dir: str) -> list[dict]:
    """Return all characters that are ready to participate in the loop."""
    state = load_state(project_dir)
    return [
        c for c in state.get("recruitment", {}).get("characters", [])
        if c.get("ready", False)
    ]


def load_character_profile(project_dir: str, slug: str) -> dict:
    """Load a character's profile from .autolab/characters/{slug}/character.yaml."""
    path = Path(project_dir) / AUTOLAB_DIR / "characters" / slug / "character.yaml"
    if path.exists():
        with open(path, encoding="utf-8") as f:
            return yaml.safe_load(f) or {}
    state = load_state(project_dir)
    for char in state.get("recruitment", {}).get("characters", []):
        if char.get("slug") == slug:
            return char
    return DEFAULT_TRAINEE_PROFILE


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


def store_direction_change(project_dir: str, new_direction: str) -> None:
    """Store an editorial direction change with high priority.

    This is different from regular feedback: direction changes are
    prefixed with [EDITORIAL REDIRECT] and placed at the FRONT of the
    queue so the PI sees them first and treats them as mandatory.
    """
    state = load_state(project_dir)
    queue = state.get("feedback_queue", [])
    if not isinstance(queue, list):
        queue = [queue] if queue else []

    entry = (
        f"[EDITORIAL REDIRECT — MANDATORY]\n"
        f"The Editor has changed the direction of this project.\n"
        f"You MUST pivot to address this new direction:\n\n"
        f"{new_direction}\n\n"
        f"Acknowledge this redirect in your next response and adjust "
        f"your research plan accordingly."
    )
    # Insert at front so it's seen first
    queue.insert(0, entry)
    state["feedback_queue"] = queue
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
    "none",  # Normal lab mode
    "submitted",  # PI submitted manuscript + cover letter
    "reviewers_invited",  # Editor invited reviewers
    "under_review",  # Reviewers generating reports
    "reviews_complete",  # All reviews in, awaiting editorial decision
    "decision_made",  # Editor made decision, returning to PI
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
        "round": 0,  # revision round counter
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
    editorial["waiting_since"] = datetime.now(timezone.utc).isoformat()
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
        editorial["waiting_since"] = datetime.now(timezone.utc).isoformat()
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
