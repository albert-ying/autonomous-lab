"""
Tests for task mode — domain config, init, prompts, and state machine.
"""

import json
import pytest
from pathlib import Path
from autonomous_lab.lab.state import (
    init_project,
    load_state,
    save_state,
    load_config,
    get_domain_id,
    is_task_mode,
    DOMAIN_CONFIGS,
    AUTOLAB_DIR,
)
from autonomous_lab.lab.prompts import (
    build_executor_prompt,
    build_verifier_prompt,
    CODING_RULES,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture
def task_project(tmp_path):
    """Initialize a task-mode project."""
    init_project(str(tmp_path), "Build a CSV parser that outputs JSON", domain="task")
    return str(tmp_path)


@pytest.fixture
def research_project(tmp_path):
    """Initialize a default research-mode project for comparison."""
    init_project(str(tmp_path), "Analyze scRNA-seq data")
    return str(tmp_path)


@pytest.fixture
def default_profile():
    return {
        "title": "Executor",
        "expertise": "technical implementation",
        "goal": "execute tasks",
        "personality": ["Thorough: checks everything twice"],
    }


# ---------------------------------------------------------------------------
# TestTaskModeInit — 9 tests
# ---------------------------------------------------------------------------
class TestTaskModeInit:
    def test_domain_config_exists(self):
        assert "task" in DOMAIN_CONFIGS

    def test_domain_config_labels(self):
        cfg = DOMAIN_CONFIGS["task"]
        assert cfg["senior_label"] == "Verifier"
        assert cfg["junior_label"] == "Executor"
        assert cfg["artifact"] == "Deliverable"

    def test_domain_config_decisions(self):
        cfg = DOMAIN_CONFIGS["task"]
        assert cfg["decisions"] == ["PASS", "FAIL"]

    def test_phase_is_active(self, task_project):
        state = load_state(task_project)
        assert state["phase"] == "active"

    def test_max_iterations_is_3(self, task_project):
        cfg = load_config(task_project)
        assert cfg["max_iterations"] == 3

    def test_skip_editorial_is_true(self, task_project):
        cfg = load_config(task_project)
        assert cfg["skip_editorial"] is True

    def test_no_figures_dir(self, task_project):
        assert not (Path(task_project) / "figures").exists()

    def test_working_dirs_created(self, task_project):
        for d in ("scripts", "results", "data"):
            assert (Path(task_project) / d).exists()

    def test_is_task_mode_true(self, task_project):
        assert is_task_mode(task_project) is True

    def test_is_task_mode_false_for_research(self, research_project):
        assert is_task_mode(research_project) is False

    def test_recruitment_ready_at_set(self, task_project):
        state = load_state(task_project)
        assert state["recruitment"]["ready_at"] is not None


# ---------------------------------------------------------------------------
# TestExecutorPrompt — 4 tests
# ---------------------------------------------------------------------------
class TestExecutorPrompt:
    def test_no_paper_rules(self, default_profile):
        prompt = build_executor_prompt(
            idea="Build X",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=0,
        )
        assert "LaTeX" not in prompt
        assert "paper/sections" not in prompt
        assert "references.bib" not in prompt
        assert "INTERPRETATION" not in prompt

    def test_coding_rules_present(self, default_profile):
        prompt = build_executor_prompt(
            idea="Build X",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=0,
        )
        assert CODING_RULES[0] in prompt

    def test_iteration_display(self, default_profile):
        prompt = build_executor_prompt(
            idea="Build X",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
            max_iterations=3,
        )
        assert "1 of 3 max" in prompt

    def test_skill_context_rendered(self, default_profile):
        prompt = build_executor_prompt(
            idea="Build X",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=0,
            skill_context="## scanpy\nUse scanpy for analysis.",
        )
        assert "scanpy" in prompt
        assert "Your Certified Skills" in prompt


# ---------------------------------------------------------------------------
# TestVerifierPrompt — 4 tests
# ---------------------------------------------------------------------------
class TestVerifierPrompt:
    def test_pass_fail_present(self, default_profile):
        prompt = build_verifier_prompt(
            idea="Build X",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
        )
        assert "PASS" in prompt
        assert "FAIL" in prompt
        assert "STATUS: completed" in prompt
        assert "STATUS: continue" in prompt

    def test_no_paper_planning(self, default_profile):
        prompt = build_verifier_prompt(
            idea="Build X",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
        )
        assert "LaTeX" not in prompt
        assert "citation" not in prompt.lower() or "Citation" not in prompt
        assert "PAPER PLANNING" not in prompt

    def test_remaining_display(self, default_profile):
        prompt = build_verifier_prompt(
            idea="Build X",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
            max_iterations=3,
        )
        assert "2 remaining" in prompt

    def test_last_iteration_message(self, default_profile):
        prompt = build_verifier_prompt(
            idea="Build X",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=3,
            max_iterations=3,
        )
        assert "LAST ITERATION" in prompt
        assert "STATUS: completed" in prompt


# ---------------------------------------------------------------------------
# TestTaskModeStateMachine — 1 test
# ---------------------------------------------------------------------------
class TestTaskModeStateMachine:
    def test_completed_status_round_trips(self, task_project):
        state = load_state(task_project)
        state["status"] = "completed"
        state["progress"] = 100
        save_state(task_project, state)

        reloaded = load_state(task_project)
        assert reloaded["status"] == "completed"
        assert reloaded["progress"] == 100
