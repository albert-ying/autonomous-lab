"""
Autonomous Lab — Recruitment Engine

Orchestrates the recruiting pipeline: marketplace search,
character resolution, and background skill learning.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path

from autonomous_lab.integrations.tooluniverse import (
    search_tools as search_tu_tools,
    get_tool_spec,
    generate_skill_md,
)
from autonomous_lab.lab.state import (
    add_character,
    update_skill_status,
    check_recruitment_ready,
    AUTOLAB_DIR,
)
from autonomous_lab.integrations.github import (
    search_characters,
    clone_character,
    fetch_skill_files,
)
from autonomous_lab.lab.skills import run_validation, load_skill_meta

logger = logging.getLogger(__name__)


SCORE_FULL_MATCH = 0.8
SCORE_PARTIAL_MATCH = 0.5
MAX_LEARNING_RETRIES = 3


@dataclass
class ResolutionStrategy:
    kind: str  # "full_match" | "partial_match" | "fresh_create"
    repo: str = ""
    missing_skills: list[str] = field(default_factory=list)
    available_skills: list[str] = field(default_factory=list)


def resolve_character(
    spec: dict,
    marketplace_results: list[dict],
) -> ResolutionStrategy:
    """Decide resolution strategy for a character spec."""
    needed = [s["name"] for s in spec.get("skills", []) if isinstance(s, dict)]
    if not needed:
        needed = [s for s in spec.get("skills", []) if isinstance(s, str)]

    if not marketplace_results:
        return ResolutionStrategy(
            kind="fresh_create",
            missing_skills=needed,
        )

    best = marketplace_results[0]
    if best["score"] >= SCORE_FULL_MATCH:
        char_skills = set(best.get("skills", []))
        missing = [s for s in needed if s not in char_skills]
        return ResolutionStrategy(
            kind="full_match",
            repo=best["repo"],
            missing_skills=missing,
            available_skills=[s for s in needed if s in char_skills],
        )
    elif best["score"] >= SCORE_PARTIAL_MATCH:
        char_skills = set(best.get("skills", []))
        missing = [s for s in needed if s not in char_skills]
        return ResolutionStrategy(
            kind="partial_match",
            repo=best["repo"],
            missing_skills=missing,
            available_skills=[s for s in needed if s in char_skills],
        )
    else:
        return ResolutionStrategy(
            kind="fresh_create",
            missing_skills=needed,
        )


async def execute_recruitment(
    project_dir: str,
    character_specs: list[dict],
) -> str:
    """Main recruitment pipeline."""
    base = Path(project_dir) / AUTOLAB_DIR / "characters"
    base.mkdir(parents=True, exist_ok=True)
    summaries = []

    for spec in character_specs:
        slug = spec.get("slug", "unnamed")
        skill_specs = spec.get("skills", [])
        needed_names = []
        required_map = {}
        for s in skill_specs:
            if isinstance(s, dict):
                needed_names.append(s["name"])
                required_map[s["name"]] = s.get("required", False)
            elif isinstance(s, str):
                needed_names.append(s)
                required_map[s] = True

        marketplace = await search_characters(needed_names, spec.get("role", ""))
        strategy = resolve_character(spec, marketplace)

        char_dir = base / slug
        if strategy.repo and strategy.kind in ("full_match", "partial_match"):
            cloned = await clone_character(strategy.repo, char_dir)
            source = "marketplace" if cloned else "fresh"
        else:
            char_dir.mkdir(parents=True, exist_ok=True)
            source = "fresh"

        # ── ToolUniverse fallback for missing skills ──
        tu_resolved: dict[str, str] = {}  # skill_name -> generated SKILL.md content
        if strategy.missing_skills:
            for missing_skill in strategy.missing_skills:
                try:
                    tu_matches = search_tu_tools(missing_skill, limit=3)
                    if tu_matches:
                        best_name = tu_matches[0].get("name", missing_skill)
                        tool_spec = get_tool_spec(best_name)
                        if tool_spec is None:
                            tool_spec = tu_matches[0]
                        skill_md = generate_skill_md(tool_spec)
                        # Save SKILL.md to character's skills dir
                        skill_dir = char_dir / "skills" / missing_skill
                        skill_dir.mkdir(parents=True, exist_ok=True)
                        (skill_dir / "SKILL.md").write_text(skill_md, encoding="utf-8")
                        # Write meta.yaml for content-only validation
                        (skill_dir / "meta.yaml").write_text(
                            "validation_type: content-only\nsource: tooluniverse\n",
                            encoding="utf-8",
                        )
                        tu_resolved[missing_skill] = skill_md
                except Exception as exc:
                    logger.warning("ToolUniverse lookup failed for %s: %s", missing_skill, exc)

        skills_status = {}
        for skill_name in needed_names:
            skill_path = char_dir / "skills" / skill_name
            if skill_path.exists() and (skill_path / "SKILL.md").exists():
                result = run_validation(skill_path)
                skills_status[skill_name] = {
                    "status": result.status,
                    "required": required_map.get(skill_name, False),
                    "source": source,
                    "tests_passed": result.passed,
                    "tests_total": result.total,
                }
            else:
                if strategy.repo:
                    fetched = await fetch_skill_files(
                        strategy.repo, skill_name,
                        char_dir / "skills" / skill_name,
                    )
                    if fetched:
                        result = run_validation(char_dir / "skills" / skill_name)
                        skills_status[skill_name] = {
                            "status": result.status,
                            "required": required_map.get(skill_name, False),
                            "source": "marketplace",
                        }
                        continue

                # Check if ToolUniverse resolved this skill
                if skill_name in tu_resolved:
                    skills_status[skill_name] = {
                        "status": "certified",
                        "required": required_map.get(skill_name, False),
                        "source": "tooluniverse",
                    }
                    continue

                skills_status[skill_name] = {
                    "status": "queued",
                    "required": required_map.get(skill_name, False),
                    "source": "pending",
                }

        add_character(project_dir, {
            "slug": slug,
            "role": spec.get("role", "trainee"),
            "name": spec.get("name", slug),
            "avatar": spec.get("avatar", "generic"),
            "source": source,
            "repo": strategy.repo,
            "skills": skills_status,
            "ready": False,
        })

        certified = sum(1 for s in skills_status.values() if s["status"] == "certified")
        total = len(skills_status)
        summaries.append(f"- **{slug}**: {certified}/{total} skills certified (source: {source})")

    check_recruitment_ready(project_dir)

    return "## Recruitment Summary\n\n" + "\n".join(summaries)
