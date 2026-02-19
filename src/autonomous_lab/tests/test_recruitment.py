"""
Tests for recruitment engine — resolve_character strategy selection.
"""

import pytest
from autonomous_lab.lab.recruitment import (
    resolve_character,
    ResolutionStrategy,
)


class TestResolveCharacter:
    def test_full_match(self):
        spec = {"slug": "bio-t", "skills": [
            {"name": "scanpy", "required": True},
        ]}
        marketplace = [{"repo": "user/autolab-char-bio", "score": 1.0, "skills": ["scanpy"]}]
        strategy = resolve_character(spec, marketplace)
        assert strategy.kind == "full_match"
        assert strategy.repo == "user/autolab-char-bio"
        assert strategy.missing_skills == []

    def test_partial_match(self):
        spec = {"slug": "bio-t", "skills": [
            {"name": "scanpy", "required": True},
            {"name": "pymc", "required": False},
        ]}
        marketplace = [{"repo": "user/autolab-char-bio", "score": 0.5, "skills": ["scanpy"]}]
        strategy = resolve_character(spec, marketplace)
        assert strategy.kind == "partial_match"
        assert "pymc" in strategy.missing_skills

    def test_no_match(self):
        spec = {"slug": "bio-t", "skills": [
            {"name": "scanpy", "required": True},
        ]}
        marketplace = []
        strategy = resolve_character(spec, marketplace)
        assert strategy.kind == "fresh_create"
        assert "scanpy" in strategy.missing_skills

    def test_low_score_creates_fresh(self):
        spec = {"slug": "bio-t", "skills": [
            {"name": "scanpy", "required": True},
            {"name": "pymc", "required": True},
        ]}
        marketplace = [{"repo": "user/some-char", "score": 0.3, "skills": ["unrelated"]}]
        strategy = resolve_character(spec, marketplace)
        assert strategy.kind == "fresh_create"

    def test_string_skills(self):
        spec = {"slug": "bio-t", "skills": ["scanpy", "pymc"]}
        marketplace = []
        strategy = resolve_character(spec, marketplace)
        assert strategy.kind == "fresh_create"
        assert "scanpy" in strategy.missing_skills
        assert "pymc" in strategy.missing_skills

    def test_full_match_with_missing_skills(self):
        """Full match repo but missing some skills."""
        spec = {"slug": "bio-t", "skills": [
            {"name": "scanpy", "required": True},
            {"name": "pymc", "required": True},
        ]}
        marketplace = [{"repo": "user/autolab-char-bio", "score": 0.8, "skills": ["scanpy"]}]
        strategy = resolve_character(spec, marketplace)
        assert strategy.kind == "full_match"
        assert "pymc" in strategy.missing_skills
        assert "scanpy" in strategy.available_skills

    def test_partial_match_available_skills(self):
        """Partial match correctly reports available skills."""
        spec = {"slug": "bio-t", "skills": [
            {"name": "scanpy", "required": True},
            {"name": "pymc", "required": False},
            {"name": "seurat", "required": True},
        ]}
        marketplace = [{"repo": "user/char", "score": 0.6, "skills": ["scanpy", "seurat"]}]
        strategy = resolve_character(spec, marketplace)
        assert strategy.kind == "partial_match"
        assert "pymc" in strategy.missing_skills
        assert "scanpy" in strategy.available_skills
        assert "seurat" in strategy.available_skills

    def test_empty_skills_spec_no_marketplace(self):
        """Spec with no skills and no marketplace results -> fresh_create."""
        spec = {"slug": "bio-t", "skills": []}
        marketplace = []
        strategy = resolve_character(spec, marketplace)
        assert strategy.kind == "fresh_create"
        assert strategy.missing_skills == []

    def test_empty_skills_spec_with_marketplace(self):
        """Spec with no skills but high-score marketplace hit -> full_match."""
        spec = {"slug": "bio-t", "skills": []}
        marketplace = [{"repo": "user/char", "score": 1.0, "skills": ["scanpy"]}]
        strategy = resolve_character(spec, marketplace)
        assert strategy.kind == "full_match"
        assert strategy.missing_skills == []

    def test_resolution_strategy_defaults(self):
        """ResolutionStrategy dataclass has correct defaults."""
        rs = ResolutionStrategy(kind="fresh_create")
        assert rs.repo == ""
        assert rs.missing_skills == []
        assert rs.available_skills == []
