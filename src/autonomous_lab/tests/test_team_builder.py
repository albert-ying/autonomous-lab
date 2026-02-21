"""Tests for team_builder.is_agent_team_available opt-in behavior."""

import os

import pytest

from autonomous_lab.lab.team_builder import is_agent_team_available


class TestIsAgentTeamAvailable:
    """Verify that native agent teams require explicit opt-in."""

    def test_empty_config_returns_false(self):
        """Empty config should never enable agent teams."""
        assert is_agent_team_available({}) is False

    def test_subprocess_mode_returns_false(self):
        """Explicit subprocess mode should not enable agent teams."""
        assert is_agent_team_available({"team_mode": "subprocess"}) is False

    def test_native_mode_with_env_var(self, monkeypatch):
        """Native mode + env var should enable agent teams."""
        monkeypatch.setenv("CLAUDE_CODE_EXPERIMENTAL_AGENT_TEAMS", "1")
        assert is_agent_team_available({"team_mode": "native"}) is True

    def test_native_mode_without_env_var(self, monkeypatch):
        """Native mode without env var should not enable agent teams."""
        monkeypatch.delenv("CLAUDE_CODE_EXPERIMENTAL_AGENT_TEAMS", raising=False)
        assert is_agent_team_available({"team_mode": "native"}) is False

    def test_none_config_returns_false(self):
        """None config should never enable agent teams."""
        assert is_agent_team_available(None) is False

    def test_no_args_returns_false(self):
        """Default (no arguments) should never enable agent teams."""
        assert is_agent_team_available() is False
