# tests/test_integration.py
import json
import pytest
from pathlib import Path
from autonomous_lab.lab.state import (
    init_project, load_state, add_character,
    update_skill_status, check_recruitment_ready,
    get_active_characters,
)
from autonomous_lab.lab.skills import build_skill_context, load_skill_content
from autonomous_lab.lab.prompts import build_trainee_prompt, build_recruiting_prompt


@pytest.fixture
def full_project(tmp_path):
    """Set up a project with a character that has one certified skill."""
    init_project(str(tmp_path), "Analyze scRNA-seq from aging mouse brain")

    # Create character directory with a skill
    char_dir = tmp_path / ".autolab" / "characters" / "bio-trainee" / "skills" / "scanpy"
    char_dir.mkdir(parents=True)
    (char_dir / "SKILL.md").write_text(
        "# Scanpy\nUse scanpy for single-cell RNA-seq analysis.\n"
        "Key functions: sc.pp.normalize_total, sc.tl.pca, sc.tl.leiden.\n"
        "Always call sc.settings.verbosity = 1 for clean output.\n"
    )
    (char_dir / "meta.yaml").write_text("status: certified\nvalidation_type: content-only\n")

    # Add character to recruitment roster
    add_character(str(tmp_path), {
        "slug": "bio-trainee",
        "role": "trainee",
        "name": "Dr. Wei Zhang",
        "avatar": "bioinformatician",
        "skills": {
            "scanpy": {"status": "certified", "required": True},
        },
        "ready": False,
    })

    return str(tmp_path)


class TestFullRecruitingFlow:
    def test_init_creates_recruiting_phase(self, full_project):
        state = load_state(full_project)
        assert state["phase"] == "recruiting"

    def test_recruitment_becomes_ready(self, full_project):
        ready = check_recruitment_ready(full_project)
        assert ready is True
        state = load_state(full_project)
        assert state["phase"] == "active"

    def test_active_characters_returned(self, full_project):
        check_recruitment_ready(full_project)
        active = get_active_characters(full_project)
        assert len(active) == 1
        assert active[0]["slug"] == "bio-trainee"

    def test_skill_injected_into_prompt(self, full_project):
        check_recruitment_ready(full_project)
        char_dir = Path(full_project) / ".autolab" / "characters" / "bio-trainee"
        context = build_skill_context(
            skills_dir=char_dir,
            certified_skills=["scanpy"],
            pi_agenda="Run scanpy preprocessing",
        )
        assert "sc.pp.normalize_total" in context

        prompt = build_trainee_prompt(
            idea="Analyze scRNA-seq",
            profile={"title": "Postdoc", "expertise": "bioinformatics",
                     "goal": "analyze", "personality": ["Dedicated"]},
            meeting_history="PI: run scanpy clustering",
            summaries="",
            file_listings={"scripts": [], "figures": [], "data": [], "results": [], "paper": []},
            user_feedback="",
            iteration=1,
            skill_context=context,
        )
        assert "sc.pp.normalize_total" in prompt
        assert "Certified Skills" in prompt

    def test_recruiting_prompt_generated(self, full_project):
        prompt = build_recruiting_prompt(
            idea="Analyze scRNA-seq from aging mouse brain",
            domain_config={"senior_label": "PI", "senior_short": "PI"},
        )
        assert "autolab_recruit" in prompt
        assert "aging" in prompt.lower()

    def test_multi_character_flow(self, tmp_path):
        """Test with multiple characters, mixed skill readiness."""
        init_project(str(tmp_path), "Multi-omics integration study")

        # Character 1: all skills certified
        add_character(str(tmp_path), {
            "slug": "bio-trainee",
            "role": "trainee",
            "skills": {
                "scanpy": {"status": "certified", "required": True},
                "scvi-tools": {"status": "certified", "required": False},
            },
            "ready": False,
        })

        # Character 2: required skill still learning
        add_character(str(tmp_path), {
            "slug": "stats-trainee",
            "role": "trainee",
            "skills": {
                "pymc": {"status": "learning", "required": True},
                "arviz": {"status": "queued", "required": False},
            },
            "ready": False,
        })

        # Not ready yet — stats-trainee has required skill still learning
        assert check_recruitment_ready(str(tmp_path)) is False
        state = load_state(str(tmp_path))
        assert state["phase"] == "recruiting"

        # Bio-trainee should be ready, stats-trainee should not
        chars = state["recruitment"]["characters"]
        bio = next(c for c in chars if c["slug"] == "bio-trainee")
        stats = next(c for c in chars if c["slug"] == "stats-trainee")
        assert bio["ready"] is True
        assert stats["ready"] is False

        # Now certify the required skill
        update_skill_status(str(tmp_path), "stats-trainee", "pymc", "certified")
        assert check_recruitment_ready(str(tmp_path)) is True

        state = load_state(str(tmp_path))
        assert state["phase"] == "active"
        active = get_active_characters(str(tmp_path))
        assert len(active) == 2
