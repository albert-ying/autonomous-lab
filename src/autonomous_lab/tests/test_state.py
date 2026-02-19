"""
Tests for state layer — phase fields, character roster, and skill tracking.
"""

import json
import pytest
from pathlib import Path
from autonomous_lab.lab.state import (
    init_project,
    load_state,
    save_state,
    add_character,
    update_skill_status,
    check_recruitment_ready,
    get_active_characters,
    load_character_profile,
    DEFAULT_TRAINEE_PROFILE,
    AUTOLAB_DIR,
)


@pytest.fixture
def project_dir(tmp_path):
    init_project(str(tmp_path), "Test idea about scRNA-seq")
    return str(tmp_path)


class TestPhaseField:
    def test_init_sets_recruiting_phase(self, project_dir):
        state = load_state(project_dir)
        assert state["phase"] == "recruiting"

    def test_init_has_recruitment_block(self, project_dir):
        state = load_state(project_dir)
        assert "recruitment" in state
        assert state["recruitment"]["characters"] == []
        assert state["recruitment"]["started_at"] is not None
        assert state["recruitment"]["ready_at"] is None

    def test_init_has_characters_list(self, project_dir):
        state = load_state(project_dir)
        assert "characters" in state
        assert state["characters"] == []


class TestCharacterRoster:
    def test_add_character(self, project_dir):
        add_character(project_dir, {
            "slug": "bioinfo-trainee",
            "role": "trainee",
            "name": "Dr. Wei Zhang",
            "avatar": "bioinformatician",
            "source": "marketplace",
            "skills": {
                "scanpy": {"status": "learning", "required": True},
                "scvi-tools": {"status": "learning", "required": False},
            },
            "ready": False,
        })
        state = load_state(project_dir)
        assert len(state["recruitment"]["characters"]) == 1
        assert state["recruitment"]["characters"][0]["slug"] == "bioinfo-trainee"

    def test_add_character_deduplicates_by_slug(self, project_dir):
        """Adding a character with the same slug should replace the old one."""
        add_character(project_dir, {
            "slug": "bio-t",
            "role": "trainee",
            "name": "Version 1",
            "skills": {},
            "ready": False,
        })
        add_character(project_dir, {
            "slug": "bio-t",
            "role": "trainee",
            "name": "Version 2",
            "skills": {},
            "ready": False,
        })
        state = load_state(project_dir)
        assert len(state["recruitment"]["characters"]) == 1
        assert state["recruitment"]["characters"][0]["name"] == "Version 2"

    def test_add_multiple_characters(self, project_dir):
        add_character(project_dir, {
            "slug": "bio-t", "role": "trainee", "skills": {}, "ready": False,
        })
        add_character(project_dir, {
            "slug": "stats-t", "role": "trainee", "skills": {}, "ready": False,
        })
        state = load_state(project_dir)
        assert len(state["recruitment"]["characters"]) == 2

    def test_update_skill_status(self, project_dir):
        add_character(project_dir, {
            "slug": "bio-t",
            "role": "trainee",
            "skills": {"scanpy": {"status": "learning", "required": True}},
            "ready": False,
        })
        update_skill_status(project_dir, "bio-t", "scanpy", "certified",
                            tests_passed=5, tests_total=5)
        state = load_state(project_dir)
        char = state["recruitment"]["characters"][0]
        assert char["skills"]["scanpy"]["status"] == "certified"
        assert char["skills"]["scanpy"]["tests_passed"] == 5

    def test_update_skill_status_nonexistent_slug(self, project_dir):
        """Updating a nonexistent slug should not raise an error."""
        add_character(project_dir, {
            "slug": "bio-t",
            "role": "trainee",
            "skills": {"scanpy": {"status": "learning", "required": True}},
            "ready": False,
        })
        # Should silently do nothing
        update_skill_status(project_dir, "nonexistent", "scanpy", "certified")
        state = load_state(project_dir)
        # Original character unchanged
        assert state["recruitment"]["characters"][0]["skills"]["scanpy"]["status"] == "learning"

    def test_update_skill_status_nonexistent_skill(self, project_dir):
        """Updating a nonexistent skill should not raise an error."""
        add_character(project_dir, {
            "slug": "bio-t",
            "role": "trainee",
            "skills": {"scanpy": {"status": "learning", "required": True}},
            "ready": False,
        })
        # Should silently do nothing
        update_skill_status(project_dir, "bio-t", "nonexistent-skill", "certified")
        state = load_state(project_dir)
        assert "nonexistent-skill" not in state["recruitment"]["characters"][0]["skills"]

    def test_check_recruitment_ready_false(self, project_dir):
        add_character(project_dir, {
            "slug": "bio-t", "role": "trainee",
            "skills": {"scanpy": {"status": "learning", "required": True}},
            "ready": False,
        })
        assert check_recruitment_ready(project_dir) is False

    def test_check_recruitment_ready_true(self, project_dir):
        add_character(project_dir, {
            "slug": "bio-t", "role": "trainee",
            "skills": {"scanpy": {"status": "certified", "required": True}},
            "ready": False,
        })
        ready = check_recruitment_ready(project_dir)
        assert ready is True
        state = load_state(project_dir)
        assert state["phase"] == "active"
        assert state["recruitment"]["ready_at"] is not None

    def test_check_recruitment_ready_empty_roster(self, project_dir):
        """An empty roster should never be ready."""
        assert check_recruitment_ready(project_dir) is False

    def test_check_recruitment_ready_optional_skills_ignored(self, project_dir):
        """Characters with only optional (non-required) skills should be ready."""
        add_character(project_dir, {
            "slug": "bio-t", "role": "trainee",
            "skills": {
                "scanpy": {"status": "certified", "required": True},
                "scvi-tools": {"status": "learning", "required": False},
            },
            "ready": False,
        })
        assert check_recruitment_ready(project_dir) is True

    def test_check_recruitment_ready_multiple_characters(self, project_dir):
        """All characters must have their required skills certified."""
        add_character(project_dir, {
            "slug": "bio-t", "role": "trainee",
            "skills": {"scanpy": {"status": "certified", "required": True}},
            "ready": False,
        })
        add_character(project_dir, {
            "slug": "stats-t", "role": "trainee",
            "skills": {"pymc": {"status": "learning", "required": True}},
            "ready": False,
        })
        assert check_recruitment_ready(project_dir) is False
        state = load_state(project_dir)
        # bio-t should be marked ready, stats-t should not
        chars = {c["slug"]: c for c in state["recruitment"]["characters"]}
        assert chars["bio-t"]["ready"] is True
        assert chars["stats-t"]["ready"] is False

    def test_get_active_characters(self, project_dir):
        add_character(project_dir, {
            "slug": "bio-t", "role": "trainee",
            "skills": {"scanpy": {"status": "certified", "required": True}},
            "ready": True,
        })
        add_character(project_dir, {
            "slug": "stats-t", "role": "trainee",
            "skills": {"pymc": {"status": "learning", "required": True}},
            "ready": False,
        })
        active = get_active_characters(project_dir)
        assert len(active) == 1
        assert active[0]["slug"] == "bio-t"

    def test_get_active_characters_empty(self, project_dir):
        active = get_active_characters(project_dir)
        assert active == []


class TestLoadCharacterProfile:
    def test_load_from_yaml_file(self, project_dir):
        """When a character.yaml exists on disk, load from it."""
        import yaml as _yaml
        char_dir = Path(project_dir) / AUTOLAB_DIR / "characters" / "bio-t"
        char_dir.mkdir(parents=True)
        profile = {"name": "Dr. Wei Zhang", "role": "trainee", "expertise": "scRNA-seq"}
        with open(char_dir / "character.yaml", "w") as f:
            _yaml.dump(profile, f)

        loaded = load_character_profile(project_dir, "bio-t")
        assert loaded["name"] == "Dr. Wei Zhang"
        assert loaded["expertise"] == "scRNA-seq"

    def test_load_from_state_fallback(self, project_dir):
        """When no YAML file exists, fall back to state roster."""
        add_character(project_dir, {
            "slug": "bio-t",
            "role": "trainee",
            "name": "Dr. Wei Zhang",
            "skills": {},
            "ready": False,
        })
        loaded = load_character_profile(project_dir, "bio-t")
        assert loaded["name"] == "Dr. Wei Zhang"

    def test_load_default_fallback(self, project_dir):
        """When nothing is found, return the default trainee profile."""
        loaded = load_character_profile(project_dir, "nonexistent-slug")
        assert loaded == DEFAULT_TRAINEE_PROFILE
