"""
Autonomous Lab — GitHub Integration

Search GitHub for character repos and skills.
Uses `gh` CLI when available, falls back to httpx.
"""

import json
import re
import shutil
import subprocess
from pathlib import Path

import yaml


def _gh_available() -> bool:
    """Check if gh CLI is available and authenticated."""
    try:
        result = subprocess.run(
            ["gh", "auth", "status"], capture_output=True, timeout=10
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def _gh_api(endpoint: str) -> dict | list | None:
    """Call GitHub API via gh CLI. Returns parsed JSON or None."""
    try:
        result = subprocess.run(
            ["gh", "api", endpoint],
            capture_output=True, text=True, timeout=30,
        )
        if result.returncode == 0:
            return json.loads(result.stdout)
    except (FileNotFoundError, subprocess.TimeoutExpired, json.JSONDecodeError):
        pass
    return None


def parse_readme_frontmatter(readme_text: str) -> dict:
    """Parse YAML front-matter from a character README.md."""
    match = re.match(r"^---\s*\n(.*?)\n---", readme_text, re.DOTALL)
    if not match:
        return {}
    try:
        return yaml.safe_load(match.group(1)) or {}
    except yaml.YAMLError:
        return {}


def score_character_match(
    character_skills: list[str],
    needed_skills: list[str],
) -> float:
    """Score how well a character's skills match the needed skills. 0.0 to 1.0."""
    if not needed_skills:
        return 0.0
    char_set = set(s.lower() for s in character_skills)
    needed_set = set(s.lower() for s in needed_skills)
    matched = char_set & needed_set
    return len(matched) / len(needed_set)


async def search_characters(
    skills: list[str],
    role: str = "",
) -> list[dict]:
    """Search GitHub for autolab-char repos matching required skills."""
    if not _gh_available():
        return []

    query_parts = ["topic:autolab-character"]
    for skill in skills[:3]:
        query_parts.append(f"topic:autolab-skill-{skill}")

    query = "+".join(query_parts)
    data = _gh_api(f"search/repositories?q={query}&sort=stars&per_page=10")

    if not data or "items" not in data:
        data = _gh_api(
            "search/repositories?q=autolab-char+in:name&sort=updated&per_page=20"
        )

    if not data or "items" not in data:
        return []

    import base64
    results = []
    for repo in data["items"]:
        full_name = repo["full_name"]
        readme_data = _gh_api(f"repos/{full_name}/readme")
        if not readme_data:
            continue

        readme_text = ""
        if "content" in readme_data:
            try:
                readme_text = base64.b64decode(readme_data["content"]).decode("utf-8")
            except Exception:
                continue

        meta = parse_readme_frontmatter(readme_text)
        if not meta:
            char_data = _gh_api(f"repos/{full_name}/contents/character.yaml")
            if char_data and "content" in char_data:
                try:
                    char_text = base64.b64decode(char_data["content"]).decode("utf-8")
                    meta = yaml.safe_load(char_text) or {}
                except Exception:
                    continue

        char_skills = []
        for s in meta.get("skills", []):
            if isinstance(s, dict):
                char_skills.append(s.get("name", ""))
            elif isinstance(s, str):
                char_skills.append(s)

        score = score_character_match(char_skills, skills)

        results.append({
            "repo": full_name,
            "score": score,
            "skills": char_skills,
            "meta": meta,
            "stars": repo.get("stargazers_count", 0),
            "updated": repo.get("pushed_at", ""),
        })

    results.sort(key=lambda x: (-x["score"], -x["stars"]))
    return results


async def fetch_skill_files(
    repo: str,
    skill_name: str,
    target_dir: Path,
) -> bool:
    """Fetch a single skill's files from a character repo."""
    import base64
    target_dir.mkdir(parents=True, exist_ok=True)

    for filename in ("SKILL.md", "validation.py", "meta.yaml"):
        path = f"skills/{skill_name}/{filename}"
        data = _gh_api(f"repos/{repo}/contents/{path}")
        if data and "content" in data:
            try:
                content = base64.b64decode(data["content"]).decode("utf-8")
                (target_dir / filename).write_text(content, encoding="utf-8")
            except Exception:
                pass

    return (target_dir / "SKILL.md").exists()


async def clone_character(repo: str, target_dir: Path) -> bool:
    """Clone a character repo into target_dir."""
    target_dir.parent.mkdir(parents=True, exist_ok=True)
    if target_dir.exists():
        shutil.rmtree(target_dir)
    try:
        result = subprocess.run(
            ["gh", "repo", "clone", repo, str(target_dir), "--", "--depth=1"],
            capture_output=True, timeout=60,
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False
