import pytest
from pathlib import Path
from autonomous_lab.integrations.github import (
    parse_readme_frontmatter,
    score_character_match,
)

class TestParseReadmeFrontmatter:
    def test_parses_yaml_block(self):
        readme = """---
autolab: character
name: Dr. Wei Zhang
role: trainee
skills:
  - name: scanpy
    status: certified
  - name: pymc
    status: certified
tags: [bioinformatics]
---

# Dr. Wei Zhang
A trainee specialized in...
"""
        meta = parse_readme_frontmatter(readme)
        assert meta["name"] == "Dr. Wei Zhang"
        assert meta["role"] == "trainee"
        assert len(meta["skills"]) == 2

    def test_no_frontmatter_returns_empty(self):
        readme = "# Just a README\nNo YAML here."
        meta = parse_readme_frontmatter(readme)
        assert meta == {}

    def test_invalid_yaml_returns_empty(self):
        readme = "---\n: invalid: yaml: [broken\n---\nContent"
        meta = parse_readme_frontmatter(readme)
        assert meta == {}

class TestScoreCharacterMatch:
    def test_perfect_match(self):
        score = score_character_match(["scanpy", "scvi-tools", "pymc"], ["scanpy", "scvi-tools", "pymc"])
        assert score == 1.0

    def test_partial_match(self):
        score = score_character_match(["scanpy"], ["scanpy", "scvi-tools"])
        assert score == 0.5

    def test_no_match(self):
        score = score_character_match(["pytorch"], ["scanpy"])
        assert score == 0.0

    def test_empty_needed(self):
        score = score_character_match(["scanpy"], [])
        assert score == 0.0

    def test_case_insensitive(self):
        score = score_character_match(["Scanpy"], ["scanpy"])
        assert score == 1.0
