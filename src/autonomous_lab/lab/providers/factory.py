"""
Autonomous Lab - Provider Factory

Maps provider configuration names to concrete LLMProvider instances.
"""

from .base import LLMProvider
from .cli_provider import CLIProvider

# CLI provider name → CLI binary key
CLI_PROVIDER_MAP = {
    "claude-cli": "claude",
    "codex-cli": "codex",
    "cursor-cli": "cursor-agent",
}


def create_provider(
    provider_name: str,
    model: str,
    project_dir: str,
    api_key: str = "",
    timeout: int = 600,
) -> LLMProvider:
    """Create an LLM provider instance from config values.

    Args:
        provider_name: One of 'claude-cli', 'codex-cli', 'cursor-cli',
                       'anthropic', or 'openai'.
        model: Model name or alias.
        project_dir: Project directory path.
        api_key: API key (only needed for 'anthropic' / 'openai').
        timeout: Subprocess timeout in seconds (CLI providers only).

    Returns:
        Configured LLMProvider instance.

    Raises:
        ValueError: If provider_name is unknown.
    """
    if provider_name in CLI_PROVIDER_MAP:
        return CLIProvider(
            cli=CLI_PROVIDER_MAP[provider_name],
            model=model,
            project_dir=project_dir,
            timeout=timeout,
        )

    if provider_name == "anthropic":
        from .anthropic_provider import AnthropicProvider

        return AnthropicProvider(model=model, api_key=api_key)

    if provider_name == "openai":
        from .openai_provider import OpenAIProvider

        return OpenAIProvider(model=model, api_key=api_key)

    raise ValueError(
        f"Unknown provider: {provider_name!r}. "
        f"Supported: {list(CLI_PROVIDER_MAP.keys()) + ['anthropic', 'openai']}"
    )
