"""Tests for ToolUniverse integration module."""

import json
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock


class TestDetection:
    def test_is_sdk_available_when_installed(self):
        from autonomous_lab.integrations.tooluniverse import is_sdk_available
        import autonomous_lab.integrations.tooluniverse as mod
        mod._sdk_available = None
        with patch("autonomous_lab.integrations.tooluniverse.importlib.import_module"):
            assert is_sdk_available() is True

    def test_is_sdk_available_when_not_installed(self):
        from autonomous_lab.integrations.tooluniverse import is_sdk_available
        import autonomous_lab.integrations.tooluniverse as mod
        mod._sdk_available = None
        with patch(
            "autonomous_lab.integrations.tooluniverse.importlib.import_module",
            side_effect=ImportError,
        ):
            assert is_sdk_available() is False

    def test_get_sdk_version_returns_version(self):
        from autonomous_lab.integrations.tooluniverse import get_sdk_version
        import autonomous_lab.integrations.tooluniverse as mod
        mod._sdk_available = None
        mock_mod = MagicMock(__version__="1.0.7")
        with patch(
            "autonomous_lab.integrations.tooluniverse.importlib.import_module",
            return_value=mock_mod,
        ):
            assert get_sdk_version() == "1.0.7"

    def test_get_sdk_version_not_installed(self):
        from autonomous_lab.integrations.tooluniverse import get_sdk_version
        import autonomous_lab.integrations.tooluniverse as mod
        mod._sdk_available = False
        assert get_sdk_version() is None


class TestConfig:
    def test_load_default_config(self, tmp_path):
        from autonomous_lab.integrations.tooluniverse import load_config
        cfg = load_config(str(tmp_path))
        assert cfg["enabled"] is True
        assert cfg["prefer_local"] is True
        assert "api_url" in cfg

    def test_save_and_load_config(self, tmp_path):
        from autonomous_lab.integrations.tooluniverse import save_config, load_config
        (tmp_path / ".autolab").mkdir()
        save_config(str(tmp_path), enabled=False, api_key="test-key")
        cfg = load_config(str(tmp_path))
        assert cfg["enabled"] is False
        assert cfg["api_key"] == "test-key"

    def test_save_creates_autolab_dir(self, tmp_path):
        from autonomous_lab.integrations.tooluniverse import save_config
        save_config(str(tmp_path), enabled=True)
        assert (tmp_path / ".autolab" / "tools_config.json").exists()


class TestGenerateSkillMd:
    def test_generate_from_tool_spec(self):
        from autonomous_lab.integrations.tooluniverse import generate_skill_md
        spec = {
            "name": "alphafold",
            "description": "Predict protein 3D structure from amino acid sequence",
            "parameters": [
                {"name": "sequence", "type": "string", "required": True,
                 "description": "Amino acid sequence"},
                {"name": "model_version", "type": "string", "required": False,
                 "description": "Model version to use"},
            ],
        }
        md = generate_skill_md(spec)
        assert "# alphafold" in md
        assert "ToolUniverse" in md
        assert "sequence" in md
        assert "tu.tools.alphafold" in md

    def test_generate_minimal_spec(self):
        from autonomous_lab.integrations.tooluniverse import generate_skill_md
        spec = {"name": "my_tool", "description": "Does something"}
        md = generate_skill_md(spec)
        assert "# my_tool" in md
        assert "Does something" in md

    def test_generate_has_enough_content(self):
        """Generated SKILL.md should be usable with content-only validation."""
        from autonomous_lab.integrations.tooluniverse import generate_skill_md
        spec = {"name": "test_tool", "description": "A test tool for unit testing"}
        md = generate_skill_md(spec)
        assert len(md.split()) >= 50


class TestSearchTools:
    def test_search_with_sdk(self):
        from autonomous_lab.integrations.tooluniverse import search_tools
        with patch(
            "autonomous_lab.integrations.tooluniverse.is_sdk_available",
            return_value=True,
        ):
            with patch(
                "autonomous_lab.integrations.tooluniverse._sdk_search_tools",
                return_value=[{"name": "alphafold", "description": "Protein structure prediction"}],
            ):
                results = search_tools("protein structure")
                assert isinstance(results, list)
                assert len(results) == 1

    def test_search_without_sdk_uses_http(self):
        from autonomous_lab.integrations.tooluniverse import search_tools
        with patch(
            "autonomous_lab.integrations.tooluniverse.is_sdk_available",
            return_value=False,
        ):
            with patch(
                "autonomous_lab.integrations.tooluniverse._http_search_tools",
                return_value=[{"name": "alphafold", "description": "Predict protein structure"}],
            ) as mock_http:
                results = search_tools("protein structure")
                mock_http.assert_called_once()
                assert len(results) == 1

    def test_search_graceful_failure(self):
        from autonomous_lab.integrations.tooluniverse import search_tools
        with patch(
            "autonomous_lab.integrations.tooluniverse.is_sdk_available",
            return_value=False,
        ):
            with patch(
                "autonomous_lab.integrations.tooluniverse._http_search_tools",
                side_effect=Exception("network error"),
            ):
                results = search_tools("anything")
                assert results == []


class TestGetToolSpec:
    def test_get_spec_with_sdk(self):
        from autonomous_lab.integrations.tooluniverse import get_tool_spec
        with patch(
            "autonomous_lab.integrations.tooluniverse.is_sdk_available",
            return_value=True,
        ):
            with patch(
                "autonomous_lab.integrations.tooluniverse._sdk_get_tool_spec",
                return_value={"name": "pubchem", "description": "Search PubChem", "parameters": []},
            ):
                spec = get_tool_spec("pubchem")
                assert spec["name"] == "pubchem"

    def test_get_spec_not_found(self):
        from autonomous_lab.integrations.tooluniverse import get_tool_spec
        with patch(
            "autonomous_lab.integrations.tooluniverse.is_sdk_available",
            return_value=False,
        ):
            with patch(
                "autonomous_lab.integrations.tooluniverse._http_get_tool_spec",
                return_value=None,
            ):
                spec = get_tool_spec("nonexistent")
                assert spec is None


class TestGetStatus:
    def test_status_with_sdk(self, tmp_path):
        from autonomous_lab.integrations.tooluniverse import get_status
        import autonomous_lab.integrations.tooluniverse as mod
        mod._sdk_available = None
        mock_mod = MagicMock(__version__="1.0.7")
        with patch(
            "autonomous_lab.integrations.tooluniverse.importlib.import_module",
            return_value=mock_mod,
        ):
            status = get_status(str(tmp_path))
            assert status["installed"] is True
            assert status["version"] == "1.0.7"

    def test_status_without_sdk(self, tmp_path):
        from autonomous_lab.integrations.tooluniverse import get_status
        import autonomous_lab.integrations.tooluniverse as mod
        mod._sdk_available = False
        status = get_status(str(tmp_path))
        assert status["installed"] is False
