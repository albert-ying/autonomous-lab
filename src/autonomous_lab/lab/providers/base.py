"""
Autonomous Lab - LLM Provider Interface

Abstract base class and shared dataclasses for all LLM providers
(CLI agents and API providers).
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field


@dataclass
class ToolCall:
    """A tool invocation requested by the model."""

    id: str
    name: str
    arguments: dict = field(default_factory=dict)


@dataclass
class ToolResult:
    """Result of executing a tool call."""

    tool_call_id: str
    content: str
    is_error: bool = False


@dataclass
class LLMResponse:
    """Unified response from any LLM provider."""

    content: str = ""
    tool_calls: list[ToolCall] = field(default_factory=list)
    stop_reason: str = ""
    input_tokens: int = 0
    output_tokens: int = 0
    model: str = ""


class LLMProvider(ABC):
    """Abstract LLM provider interface."""

    @abstractmethod
    async def complete(
        self,
        system_prompt: str,
        messages: list[dict],
        tools: list[dict] | None = None,
        max_tokens: int = 8192,
    ) -> LLMResponse:
        """Send a completion request and return the response.

        Args:
            system_prompt: System-level instructions.
            messages: Conversation messages in OpenAI-style format.
            tools: Tool schemas (for API providers). Ignored by CLI providers.
            max_tokens: Maximum tokens in response.

        Returns:
            Unified LLMResponse.
        """

    def format_tools(self, tool_specs: list[dict]) -> list[dict]:
        """Format tool specifications for this provider's API.

        Default implementation returns specs unchanged. Override for
        provider-specific formatting.
        """
        return tool_specs

    @property
    @abstractmethod
    def provider_name(self) -> str:
        """Human-readable provider identifier."""
