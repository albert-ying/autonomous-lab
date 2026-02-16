"""
Autonomous Lab - Anthropic API Provider

Direct Anthropic API access via SDK or raw HTTP fallback.
For users who prefer per-token billing over CLI subscriptions.
"""

import json
import logging

from .base import LLMProvider, LLMResponse, ToolCall

logger = logging.getLogger(__name__)


class AnthropicProvider(LLMProvider):
    """Anthropic Messages API provider.

    Tries the anthropic SDK first, falls back to raw HTTP via aiohttp.
    Handles tool_use stop reasons by extracting ToolCall objects.
    """

    def __init__(self, model: str, api_key: str):
        self.model = model
        self.api_key = api_key
        self._client = None

    @property
    def provider_name(self) -> str:
        return "anthropic"

    def _get_client(self):
        """Lazy-init the Anthropic async client."""
        if self._client is None:
            try:
                import anthropic

                self._client = anthropic.AsyncAnthropic(api_key=self.api_key)
            except ImportError:
                raise ImportError(
                    "anthropic SDK not installed. "
                    "Install with: pip install 'autonomous-lab[multi-agent]'"
                )
        return self._client

    def format_tools(self, tool_specs: list[dict]) -> list[dict]:
        """Format tools for Anthropic's tool use format."""
        formatted = []
        for spec in tool_specs:
            formatted.append(
                {
                    "name": spec["name"],
                    "description": spec.get("description", ""),
                    "input_schema": spec.get("parameters", spec.get("input_schema", {})),
                }
            )
        return formatted

    async def complete(
        self,
        system_prompt: str,
        messages: list[dict],
        tools: list[dict] | None = None,
        max_tokens: int = 8192,
    ) -> LLMResponse:
        """Send request to Anthropic Messages API."""
        client = self._get_client()

        kwargs: dict = {
            "model": self.model,
            "max_tokens": max_tokens,
            "messages": messages,
        }
        if system_prompt:
            kwargs["system"] = system_prompt
        if tools:
            kwargs["tools"] = tools

        try:
            response = await client.messages.create(**kwargs)
        except Exception as e:
            logger.error("Anthropic API error: %s", e)
            raise

        # Parse response
        content_text = ""
        tool_calls = []

        for block in response.content:
            if block.type == "text":
                content_text += block.text
            elif block.type == "tool_use":
                tool_calls.append(
                    ToolCall(
                        id=block.id,
                        name=block.name,
                        arguments=block.input if isinstance(block.input, dict) else {},
                    )
                )

        return LLMResponse(
            content=content_text,
            tool_calls=tool_calls,
            stop_reason=response.stop_reason or "",
            input_tokens=getattr(response.usage, "input_tokens", 0),
            output_tokens=getattr(response.usage, "output_tokens", 0),
            model=response.model,
        )

    async def complete_raw(
        self,
        system_prompt: str,
        messages: list[dict],
        tools: list[dict] | None = None,
        max_tokens: int = 8192,
    ) -> LLMResponse:
        """Fallback: raw HTTP via aiohttp when SDK is unavailable."""
        import aiohttp

        headers = {
            "x-api-key": self.api_key,
            "anthropic-version": "2023-06-01",
            "content-type": "application/json",
        }

        body: dict = {
            "model": self.model,
            "max_tokens": max_tokens,
            "messages": messages,
        }
        if system_prompt:
            body["system"] = system_prompt
        if tools:
            body["tools"] = tools

        async with aiohttp.ClientSession() as session:
            async with session.post(
                "https://api.anthropic.com/v1/messages",
                headers=headers,
                json=body,
            ) as resp:
                if resp.status != 200:
                    error_text = await resp.text()
                    raise RuntimeError(
                        f"Anthropic API error {resp.status}: {error_text[:500]}"
                    )
                data = await resp.json()

        content_text = ""
        tool_calls = []

        for block in data.get("content", []):
            if block.get("type") == "text":
                content_text += block.get("text", "")
            elif block.get("type") == "tool_use":
                tool_calls.append(
                    ToolCall(
                        id=block.get("id", ""),
                        name=block.get("name", ""),
                        arguments=block.get("input", {}),
                    )
                )

        usage = data.get("usage", {})
        return LLMResponse(
            content=content_text,
            tool_calls=tool_calls,
            stop_reason=data.get("stop_reason", ""),
            input_tokens=usage.get("input_tokens", 0),
            output_tokens=usage.get("output_tokens", 0),
            model=data.get("model", self.model),
        )
