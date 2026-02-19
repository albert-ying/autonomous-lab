import pytest
from pathlib import Path
from autonomous_lab.lab.skills import (
    load_skill_content,
    build_skill_context,
    run_validation,
    select_active_skills,
    ValidationResult,
)

@pytest.fixture
def skill_dir(tmp_path):
    scanpy_dir = tmp_path / "skills" / "scanpy"
    scanpy_dir.mkdir(parents=True)
    (scanpy_dir / "SKILL.md").write_text(
        "# Scanpy\nUse scanpy for single-cell RNA-seq analysis.\n"
        "Key functions: sc.pp.normalize_total, sc.tl.pca, sc.tl.leiden.\n"
    )
    (scanpy_dir / "validation.py").write_text(
        "def test_import():\n    assert True\n"
    )
    (scanpy_dir / "meta.yaml").write_text(
        "status: certified\nsource: self-learned\nvalidation_type: pytest\n"
    )

    writing_dir = tmp_path / "skills" / "scientific-writing"
    writing_dir.mkdir(parents=True)
    (writing_dir / "SKILL.md").write_text(
        "# Scientific Writing\nWrite clear, concise LaTeX sections. "
        "Use active voice. Cite with bibtex keys. " * 10
    )
    (writing_dir / "meta.yaml").write_text(
        "status: certified\nsource: marketplace\nvalidation_type: content-only\n"
    )
    return tmp_path

class TestLoadSkillContent:
    def test_load_existing_skill(self, skill_dir):
        content = load_skill_content(skill_dir, "scanpy")
        assert "Scanpy" in content
        assert "sc.pp.normalize_total" in content

    def test_load_missing_skill_returns_empty(self, skill_dir):
        content = load_skill_content(skill_dir, "nonexistent")
        assert content == ""

class TestSelectActiveSkills:
    def test_matches_skill_name_in_agenda(self):
        active = select_active_skills(
            certified_skills=["scanpy", "scvi-tools", "scientific-writing"],
            pi_agenda="Run scanpy preprocessing and clustering analysis",
            skill_summaries={"scanpy": "scRNA-seq analysis", "scvi-tools": "VAE models", "scientific-writing": "LaTeX"},
        )
        assert "scanpy" in active

    def test_fallback_to_first_two(self):
        active = select_active_skills(
            certified_skills=["scanpy", "scvi-tools", "pymc"],
            pi_agenda="Do something unrelated",
            skill_summaries={},
        )
        assert len(active) == 2

class TestBuildSkillContext:
    def test_builds_tiered_context(self, skill_dir):
        context = build_skill_context(
            skills_dir=skill_dir / "skills",
            certified_skills=["scanpy", "scientific-writing"],
            pi_agenda="Run scanpy clustering",
        )
        assert "sc.pp.normalize_total" in context
        assert "Scientific Writing" in context or "scientific-writing" in context

class TestRunValidation:
    def test_passing_validation(self, skill_dir):
        result = run_validation(skill_dir / "skills" / "scanpy")
        assert result.status == "certified"
        assert result.passed >= 1

    def test_content_only_skill(self, skill_dir):
        result = run_validation(skill_dir / "skills" / "scientific-writing")
        assert result.status == "certified"
