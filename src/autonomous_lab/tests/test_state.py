"""
Tests for state layer — phase fields, character roster, skill tracking,
and meeting preamble context retention.
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
    append_meeting_log,
    compress_old_meetings,
    get_meeting_summaries,
    _generate_preamble,
    DEFAULT_TRAINEE_PROFILE,
    AUTOLAB_DIR,
    MEETING_PREAMBLE,
    MEETING_SUMMARIES,
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


class TestGeneratePreamble:
    def test_empty_summaries(self):
        assert _generate_preamble("") == ""
        assert _generate_preamble("# Meeting Summaries\n\n") == ""

    def test_single_entry(self):
        text = "# Meeting Summaries\n\n- **Iteration 1 — PI Turn**: Set up project goals\n"
        result = _generate_preamble(text)
        assert "**Project origin**" in result
        # Only one entry, so no "Current state" section
        assert "**Current state**" not in result

    def test_multiple_entries(self):
        text = (
            "# Meeting Summaries\n\n"
            "- **Iteration 1 — PI Turn**: Defined research question on gene regulation\n"
            "- **Iteration 2 — TRAINEE Turn**: Implemented initial pipeline\n"
            "- **Iteration 3 — PI Turn**: Reviewed results and approved\n"
        )
        result = _generate_preamble(text)
        assert "**Project origin**" in result
        assert "**Current state**" in result
        assert "Iteration 1" in result
        assert "Iteration 3" in result

    def test_pivot_keywords_captured(self):
        text = (
            "# Meeting Summaries\n\n"
            "- **Iteration 1 — PI Turn**: Initial setup\n"
            "- **Iteration 3 — PI Turn**: Strategic decision to pivot to UMAP\n"
            "- **Iteration 5 — PI Turn**: Major revision of methods section\n"
            "- **Iteration 7 — TRAINEE Turn**: Final figures done\n"
        )
        result = _generate_preamble(text)
        assert "**Key decisions**" in result
        assert "pivot" in result.lower()
        assert "major revision" in result.lower()

    def test_max_three_decisions(self):
        text = (
            "# Meeting Summaries\n\n"
            "- **Iteration 1 — PI Turn**: Origin\n"
            "- **Iteration 2 — PI Turn**: Strategic decision A\n"
            "- **Iteration 3 — PI Turn**: Pivot to new approach\n"
            "- **Iteration 4 — PI Turn**: Major revision round 1\n"
            "- **Iteration 5 — PI Turn**: Restructured the paper\n"
            "- **Iteration 6 — TRAINEE Turn**: Done\n"
        )
        result = _generate_preamble(text)
        # Should capture at most 3 decision entries
        decisions_line = [l for l in result.split("\n") if "**Key decisions**" in l]
        assert len(decisions_line) == 1
        # Count pipe separators (max 2 for 3 items)
        assert decisions_line[0].count(" | ") <= 2

    def test_truncation(self):
        long_summary = "x" * 300
        text = f"# Meeting Summaries\n\n- **Iteration 1 — PI Turn**: {long_summary}\n"
        result = _generate_preamble(text)
        # Origin line should be truncated to 150 chars
        origin_line = [l for l in result.split("\n") if "**Project origin**" in l][0]
        assert len(origin_line) < 200  # 150 chars + prefix


class TestCompressOldMeetingsWithPreamble:
    def test_preamble_created_on_compression(self, project_dir):
        """Compressing meetings should create meeting_preamble.md."""
        # Add enough meetings to trigger compression (> keep_recent)
        for i in range(1, 7):
            role = "pi" if i % 2 else "trainee"
            append_meeting_log(project_dir, role, i, f"Summary for iter {i}", "Details")

        compress_old_meetings(project_dir, keep_recent=3)

        preamble_path = Path(project_dir) / AUTOLAB_DIR / MEETING_PREAMBLE
        assert preamble_path.exists()
        content = preamble_path.read_text(encoding="utf-8")
        assert "**Project origin**" in content

    def test_preamble_not_created_when_nothing_to_compress(self, project_dir):
        """If there aren't enough entries, no preamble should be created."""
        append_meeting_log(project_dir, "pi", 1, "Only one", "Details")
        compress_old_meetings(project_dir, keep_recent=3)

        preamble_path = Path(project_dir) / AUTOLAB_DIR / MEETING_PREAMBLE
        assert not preamble_path.exists()

    def test_preamble_updated_on_subsequent_compression(self, project_dir):
        """Second compression should update the preamble with new entries."""
        for i in range(1, 7):
            append_meeting_log(project_dir, "pi", i, f"Summary {i}", "Details")
        compress_old_meetings(project_dir, keep_recent=3)

        # Add more meetings and compress again
        for i in range(7, 13):
            append_meeting_log(project_dir, "trainee", i, f"Summary {i}", "Details")
        compress_old_meetings(project_dir, keep_recent=3)

        preamble_path = Path(project_dir) / AUTOLAB_DIR / MEETING_PREAMBLE
        content = preamble_path.read_text(encoding="utf-8")
        assert "**Current state**" in content


class TestGetMeetingSummariesWithPreamble:
    def test_no_preamble_no_summaries(self, project_dir):
        """With no data, returns empty string."""
        # Delete the auto-created summaries file
        summaries_path = Path(project_dir) / AUTOLAB_DIR / MEETING_SUMMARIES
        summaries_path.unlink()
        result = get_meeting_summaries(project_dir)
        assert result == ""

    def test_preamble_only(self, project_dir):
        """If only preamble exists (no summaries), return preamble."""
        summaries_path = Path(project_dir) / AUTOLAB_DIR / MEETING_SUMMARIES
        summaries_path.unlink()
        preamble_path = Path(project_dir) / AUTOLAB_DIR / MEETING_PREAMBLE
        preamble_path.write_text("**Project origin**: Test project\n", encoding="utf-8")

        result = get_meeting_summaries(project_dir)
        assert "## Project Narrative" in result
        assert "Test project" in result

    def test_preamble_plus_summaries(self, project_dir):
        """Both preamble and summaries should appear in output."""
        preamble_path = Path(project_dir) / AUTOLAB_DIR / MEETING_PREAMBLE
        preamble_path.write_text("**Project origin**: Studying gene regulation\n", encoding="utf-8")
        summaries_path = Path(project_dir) / AUTOLAB_DIR / MEETING_SUMMARIES
        summaries_path.write_text(
            "# Meeting Summaries\n\n- **Iteration 1**: Did stuff\n",
            encoding="utf-8",
        )

        result = get_meeting_summaries(project_dir)
        assert "## Project Narrative" in result
        assert "## Recent History" in result
        assert "gene regulation" in result
        assert "Did stuff" in result

    def test_budget_respected(self, project_dir):
        """Output should not exceed max_chars."""
        preamble_path = Path(project_dir) / AUTOLAB_DIR / MEETING_PREAMBLE
        preamble_path.write_text("**Project origin**: Short\n", encoding="utf-8")
        # Create a large summaries file
        summaries_path = Path(project_dir) / AUTOLAB_DIR / MEETING_SUMMARIES
        lines = ["# Meeting Summaries\n"]
        for i in range(100):
            lines.append(f"- **Iteration {i}**: Summary entry number {i} with padding text\n")
        summaries_path.write_text("\n".join(lines), encoding="utf-8")

        result = get_meeting_summaries(project_dir, max_chars=2000)
        assert len(result) <= 2200  # allow some slack for headers

    def test_summaries_only_no_preamble(self, project_dir):
        """Without preamble, just return summaries as before."""
        summaries_path = Path(project_dir) / AUTOLAB_DIR / MEETING_SUMMARIES
        summaries_path.write_text(
            "# Meeting Summaries\n\n- **Iteration 1**: Did stuff\n",
            encoding="utf-8",
        )

        result = get_meeting_summaries(project_dir)
        assert "## Project Narrative" not in result
        assert "Did stuff" in result

    def test_default_max_chars_is_4000(self, project_dir):
        """Verify the default budget was lowered from 8000 to 4000."""
        import inspect
        sig = inspect.signature(get_meeting_summaries)
        assert sig.parameters["max_chars"].default == 4000
