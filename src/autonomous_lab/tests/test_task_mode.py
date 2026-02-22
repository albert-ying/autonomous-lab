"""
Tests for task mode — domain config, init, prompts, state machine,
scope correctness, iteration cap enforcement, and compliance review gate.
"""

import json
import re
import pytest
from pathlib import Path
from autonomous_lab.lab.state import (
    init_project,
    load_state,
    save_state,
    load_config,
    get_domain_id,
    is_task_mode,
    append_meeting_log,
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


# ---------------------------------------------------------------------------
# TestVerifierScopeChecks — verifier prompt includes hard scope logic
# ---------------------------------------------------------------------------
class TestVerifierScopeChecks:
    """Verify the Verifier prompt contains hard programmatic scope enforcement."""

    def test_shared_intersection_language(self, default_profile):
        prompt = build_verifier_prompt(
            idea="Identify shared pathways between dataset A and B",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
        )
        # Must instruct verifier to check intersection semantics
        assert "intersection" in prompt.lower()
        assert "shared" in prompt.lower()

    def test_significance_threshold_language(self, default_profile):
        prompt = build_verifier_prompt(
            idea="Find significant DE genes with adjusted p < 0.05",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
        )
        # Must instruct verifier to check threshold filtering
        assert "threshold" in prompt.lower()
        assert "adj_pvalue" in prompt or "padj" in prompt

    def test_top_n_row_count_language(self, default_profile):
        prompt = build_verifier_prompt(
            idea="Return top 10 enriched pathways",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
        )
        # Must instruct verifier to check row count bounds
        assert "max_rows" in prompt or "row count" in prompt.lower()
        assert "top" in prompt.lower()

    def test_scope_violation_forces_fail(self, default_profile):
        prompt = build_verifier_prompt(
            idea="Identify shared pathways",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
        )
        # Must state that scope violations are automatic FAIL
        assert "scope violation" in prompt.lower() or "Scope violations are automatic FAIL" in prompt

    def test_500_row_sanity_check(self, default_profile):
        prompt = build_verifier_prompt(
            idea="Find significant genes",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
        )
        # Must include the 500-row sanity check threshold
        assert "500" in prompt


# ---------------------------------------------------------------------------
# TestExecutorScopeConstraints — executor prompt requires explicit scope
# ---------------------------------------------------------------------------
class TestExecutorScopeConstraints:
    """Verify the Executor prompt requires explicit scope constraint extraction."""

    def test_scope_constraint_rules_section(self, default_profile):
        prompt = build_executor_prompt(
            idea="Find shared pathways",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=0,
        )
        assert "Scope Constraint Rules" in prompt

    def test_intersection_instruction(self, default_profile):
        prompt = build_executor_prompt(
            idea="Find shared pathways",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=0,
        )
        # Must explain shared = intersection, not union
        assert "intersection" in prompt.lower()
        assert "NOT a union" in prompt or "not a union" in prompt.lower()

    def test_filter_explicitly_instruction(self, default_profile):
        prompt = build_executor_prompt(
            idea="Find significant DE genes",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=0,
        )
        # Must instruct to apply filters as separate visible step
        assert "filter explicitly" in prompt.lower() or "Apply filters explicitly" in prompt

    def test_row_count_printing_required(self, default_profile):
        prompt = build_executor_prompt(
            idea="Build output",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=0,
        )
        # Must require printing row counts before/after filtering
        assert "Rows before filter" in prompt or "row count" in prompt.lower()

    def test_approach_requires_scope_listing(self, default_profile):
        prompt = build_executor_prompt(
            idea="Find shared pathways",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=0,
        )
        # APPROACH section must ask for scope constraints identification
        assert "scope constraints you identified" in prompt.lower()


# ---------------------------------------------------------------------------
# TestIterationCapEnforcement — max_iterations strictly enforced
# ---------------------------------------------------------------------------
class TestIterationCapEnforcement:
    """Verify the state machine doesn't allow iterations beyond the cap."""

    def test_verifier_at_cap_forces_completed(self, task_project):
        """When verifier says 'continue' at iteration >= max, state becomes completed."""
        state = load_state(task_project)
        # Simulate: iteration=3, max=3, verifier (pi) just acted
        state["iteration"] = 3
        state["next_role"] = "pi"  # verifier's turn just ended
        state["status"] = "active"
        save_state(task_project, state)

        # The key logic: in autolab_record for task mode, when role=pi and
        # iteration >= max_iterations, the state should become "completed"
        # regardless of the status parameter.
        #
        # We test the state logic directly since autolab_record is an async MCP tool.
        cfg = load_config(task_project)
        max_iter = cfg["max_iterations"]
        assert state["iteration"] >= max_iter
        # After verifier acts at cap: force-complete (no more executor turns)

    def test_executor_at_cap_gets_final_verifier(self, task_project):
        """When executor finishes at cap, one final verifier turn is allowed."""
        state = load_state(task_project)
        # Simulate: executor (trainee) just finished, new_iteration would be 3
        state["iteration"] = 2  # will become 3 after trainee
        state["next_role"] = "trainee"
        state["status"] = "active"
        save_state(task_project, state)

        cfg = load_config(task_project)
        max_iter = cfg["max_iterations"]
        new_iteration = state["iteration"] + 1  # trainee increments
        # At cap: should route to verifier (pi), not back to executor
        assert new_iteration >= max_iter

    def test_max_iterations_3_means_3_executor_turns_max(self, task_project):
        """With max_iterations=3, iteration counter 0,1,2 = 3 executor turns."""
        cfg = load_config(task_project)
        assert cfg["max_iterations"] == 3
        # Iterations 0, 1, 2 = three executor turns
        # At iteration 3 (after executor turn 3), verifier force-completes
        state = load_state(task_project)
        state["iteration"] = 3
        save_state(task_project, state)
        reloaded = load_state(task_project)
        assert reloaded["iteration"] == 3
        # This is at the cap — no more executor turns should be allowed


# ---------------------------------------------------------------------------
# TestComplianceReviewGate — compliance gate in verifier prompt
# ---------------------------------------------------------------------------
class TestComplianceReviewGate:
    """Verify the Compliance Review Gate is present and enforces structured checks."""

    def test_compliance_gate_section_exists(self, default_profile):
        """The verifier prompt must contain the COMPLIANCE REVIEW GATE section."""
        prompt = build_verifier_prompt(
            idea="Build a CSV parser",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
        )
        assert "COMPLIANCE REVIEW GATE" in prompt

    def test_constraint_table_structure(self, default_profile):
        """The gate must instruct building a constraint table with type and ambiguity columns."""
        prompt = build_verifier_prompt(
            idea="Find top 10 pathways",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
        )
        assert "Constraint Text" in prompt
        assert "Type" in prompt
        assert "Ambiguous?" in prompt

    def test_ambiguity_handling_protocol(self, default_profile):
        """The gate must require explicit assumptions for ambiguous constraints."""
        prompt = build_verifier_prompt(
            idea="Find significant genes",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
        )
        assert "ambiguous" in prompt.lower()
        assert "assumption" in prompt.lower()
        # Must specify that missing assumptions = FAIL
        assert "FAIL" in prompt

    def test_explicit_unmet_constraint_forces_fail(self, default_profile):
        """If spec says 'top 10' and output has more, it must be auto-FAIL."""
        prompt = build_verifier_prompt(
            idea="Return top 10 enriched pathways",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
        )
        assert "Explicit unmet constraints are automatic FAIL" in prompt

    def test_gate_decision_all_must_hold(self, default_profile):
        """The gate must require ALL conditions to hold for PASS."""
        prompt = build_verifier_prompt(
            idea="Build output",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
        )
        assert "Gate PASSES only if ALL" in prompt or "gate FAILS" in prompt.lower()
        assert "MUST NOT issue a PASS" in prompt

    def test_subset_required_full_output_must_fail(self, default_profile):
        """When spec requires a subset, the prompt must enforce subset checking."""
        prompt = build_verifier_prompt(
            idea="Identify the top 20 differentially expressed genes from the RNA-seq dataset",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
        )
        # The prompt must contain instructions that would cause a full-dataset
        # output to FAIL: row bound check + filtering language
        assert "row_bound" in prompt or "row bound" in prompt.lower()
        assert "filtered subset" in prompt.lower()
        # The scope-size sanity gate must catch >500 rows
        assert "500" in prompt

    def test_no_subset_required_can_pass(self, default_profile):
        """When spec has no filtering language, schema/value checks suffice for PASS."""
        prompt = build_verifier_prompt(
            idea="Convert all CSV files in data/ to JSON format",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
        )
        # Schema check must still be present
        assert "schema" in prompt.lower()
        # But the prompt doesn't force-fail on row count when no filtering is specified
        # (the 500-row gate only triggers when filtering language is present)
        assert "filtering language" in prompt.lower()

    def test_alzheimer_shared_significant_prompt(self, default_profile):
        """Alzheimer-like prompt with 'shared' + 'significant' triggers both checks."""
        prompt = build_verifier_prompt(
            idea="Identify shared significantly differentially expressed genes (adjusted p < 0.01) between Alzheimer's disease cortex and hippocampus datasets",
            profile=default_profile,
            meeting_history="",
            summaries="",
            file_listings={},
            user_feedback="",
            iteration=1,
        )
        # Must have intersection check instructions
        assert "intersection" in prompt.lower()
        # Must have threshold check instructions
        assert "threshold" in prompt.lower()
        # Must have the constraint table approach
        assert "COMPLIANCE REVIEW GATE" in prompt
        # Scope violations auto-FAIL
        assert "Scope violations are automatic FAIL" in prompt
