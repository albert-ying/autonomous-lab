"""
Autonomous Lab — Skill Engine

Handles skill loading, validation, tiered prompt injection,
and the self-learning lifecycle.
"""

import subprocess
import re
from dataclasses import dataclass, field
from pathlib import Path

import yaml


@dataclass
class ValidationResult:
    status: str  # "certified" | "failed"
    passed: int = 0
    failed: int = 0
    total: int = 0
    attempt: int = 1
    error: str = ""


def load_skill_content(skills_base: Path, skill_name: str) -> str:
    """Load SKILL.md content for a given skill. Returns '' if not found."""
    if not isinstance(skills_base, Path):
        skills_base = Path(skills_base)
    candidates = [
        skills_base / "skills" / skill_name / "SKILL.md",
        skills_base / skill_name / "SKILL.md",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.read_text(encoding="utf-8")
    return ""


def load_skill_meta(skills_base: Path, skill_name: str) -> dict:
    """Load meta.yaml for a skill. Returns empty dict if not found."""
    candidates = [
        skills_base / "skills" / skill_name / "meta.yaml",
        skills_base / skill_name / "meta.yaml",
    ]
    for candidate in candidates:
        if candidate.exists():
            with open(candidate, encoding="utf-8") as f:
                return yaml.safe_load(f) or {}
    return {}


MAX_ACTIVE_SKILLS = 3
MAX_SKILL_TOKENS = 1500
FALLBACK_ACTIVE = 2


def select_active_skills(
    certified_skills: list[str],
    pi_agenda: str,
    skill_summaries: dict[str, str] | None = None,
) -> list[str]:
    """Select which skills are relevant to the current task."""
    if not certified_skills:
        return []
    agenda_lower = pi_agenda.lower()
    active = []
    for skill_name in certified_skills:
        if skill_name.replace("-", " ") in agenda_lower or skill_name in agenda_lower:
            active.append(skill_name)
            continue
        if skill_summaries:
            summary = skill_summaries.get(skill_name, "")
            for word in summary.lower().split():
                if len(word) > 4 and word in agenda_lower:
                    active.append(skill_name)
                    break
    active = list(dict.fromkeys(active))
    if not active:
        active = certified_skills[:FALLBACK_ACTIVE]
    return active[:MAX_ACTIVE_SKILLS]


def _truncate_to_tokens(text: str, max_tokens: int) -> str:
    """Rough truncation by word count (1 token ~ 0.75 words)."""
    max_words = int(max_tokens * 0.75)
    words = text.split()
    if len(words) <= max_words:
        return text
    return " ".join(words[:max_words]) + "\n\n[... truncated — see full SKILL.md for details]"


def build_skill_context(
    skills_dir: Path,
    certified_skills: list[str],
    pi_agenda: str = "",
) -> str:
    """Build tiered skill context for prompt injection."""
    if not certified_skills:
        return ""

    summaries = {}
    for skill_name in certified_skills:
        content = load_skill_content(skills_dir, skill_name)
        if content:
            first_meaningful = ""
            for line in content.split("\n"):
                stripped = line.strip()
                if stripped and not stripped.startswith("#"):
                    first_meaningful = stripped
                    break
            summaries[skill_name] = first_meaningful

    active = select_active_skills(certified_skills, pi_agenda, summaries)
    parts = []

    for skill_name in active:
        content = load_skill_content(skills_dir, skill_name)
        if content:
            truncated = _truncate_to_tokens(content, MAX_SKILL_TOKENS)
            parts.append(f"### {skill_name} (active this turn)\n\n{truncated}")

    remaining = [s for s in certified_skills if s not in active]
    if remaining:
        summary_lines = []
        for skill_name in remaining:
            desc = summaries.get(skill_name, "Available")
            summary_lines.append(f"- **{skill_name}**: {desc}")
        parts.append("### Other available skills:\n" + "\n".join(summary_lines))

    if parts:
        parts.append("Use your active skills as your primary toolkit. Follow their instructions exactly.")

    return "\n\n".join(parts)


def run_validation(skill_path: Path) -> ValidationResult:
    """Run validation.py for a skill via pytest."""
    meta = {}
    meta_path = skill_path / "meta.yaml"
    if meta_path.exists():
        with open(meta_path, encoding="utf-8") as f:
            meta = yaml.safe_load(f) or {}

    if meta.get("validation_type") == "content-only":
        skill_md = skill_path / "SKILL.md"
        if skill_md.exists():
            text = skill_md.read_text(encoding="utf-8")
            if len(text.split()) >= 100:
                return ValidationResult(status="certified", passed=1, total=1)
        return ValidationResult(status="failed", error="SKILL.md missing or too short")

    validation_file = skill_path / "validation.py"
    if not validation_file.exists():
        return ValidationResult(status="failed", error="No validation.py found")

    try:
        result = subprocess.run(
            ["python", "-m", "pytest", str(validation_file), "-v", "--tb=short", "-q"],
            capture_output=True, timeout=120, text=True,
        )
        passed = 0
        failed = 0
        for line in result.stdout.split("\n"):
            match = re.search(r"(\d+) passed", line)
            if match:
                passed = int(match.group(1))
            match = re.search(r"(\d+) failed", line)
            if match:
                failed = int(match.group(1))
        total = passed + failed
        status = "certified" if result.returncode == 0 else "failed"
        error = result.stderr if result.returncode != 0 else ""
        return ValidationResult(status=status, passed=passed, failed=failed, total=total, error=error)
    except subprocess.TimeoutExpired:
        return ValidationResult(status="failed", error="Validation timed out (>120s)")
    except FileNotFoundError:
        return ValidationResult(status="failed", error="pytest not found in PATH")
