import pytest
from autonomous_lab.lab.prompts import (
    build_trainee_prompt,
    build_recruiting_prompt,
)


class TestSkillInjection:
    def test_trainee_prompt_includes_skill_context(self):
        prompt = build_trainee_prompt(
            idea="Analyze scRNA-seq data",
            profile={"title": "Postdoc", "expertise": "bioinformatics",
                     "goal": "analyze data", "personality": ["Dedicated"]},
            meeting_history="PI asked for clustering",
            summaries="",
            file_listings={"scripts": [], "figures": [], "data": [], "results": [], "paper": []},
            user_feedback="",
            iteration=1,
            skill_context="### scanpy (active this turn)\nUse scanpy for scRNA-seq.\n",
        )
        assert "scanpy (active this turn)" in prompt
        assert "Certified Skills" in prompt

    def test_trainee_prompt_without_skills(self):
        prompt = build_trainee_prompt(
            idea="Analyze data",
            profile={"title": "Postdoc", "expertise": "stats",
                     "goal": "analyze", "personality": []},
            meeting_history="",
            summaries="",
            file_listings={"scripts": [], "figures": [], "data": [], "results": [], "paper": []},
            user_feedback="",
            iteration=1,
        )
        assert "Certified Skills" not in prompt


class TestRecruitingPrompt:
    def test_recruiting_prompt_includes_idea(self):
        prompt = build_recruiting_prompt(
            idea="Study aging in mouse brain scRNA-seq",
            domain_config={"senior_label": "PI", "senior_short": "PI"},
        )
        assert "aging" in prompt.lower()
        assert "autolab_recruit" in prompt
        assert "skill" in prompt.lower()

    def test_recruiting_prompt_default_config(self):
        prompt = build_recruiting_prompt(idea="Some research idea")
        assert "Principal Investigator" in prompt
