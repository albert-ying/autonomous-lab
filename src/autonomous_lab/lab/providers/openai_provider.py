"""
Autonomous Lab - OpenAI API Provider

Direct OpenAI API access via SDK or raw HTTP fallback.
Supports o3, GPT-4o, and other OpenAI models.
"""

import json
import logging

from .base import LLMProvider, LLMResponse, ToolCall

logger = logging.getLogger(__name__)


class OpenAIProvider(LLMProvider):
    """OpenAI Chat Completions API provider.

    Tries the openai SDK first, falls back to raw HTTP via aiohttp.
    Handles function/tool calling by extracting ToolCall objects.
    """

    def __init__(self, model: str, api_key: str):
        self.model = model
        self.api_key = api_key
        self._client = None

    @property
    def provider_name(self) -> str:
        return "openai"

    def _get_client(self):
        """Lazy-init the OpenAI async client."""
        if self._client is None:
            try:
                import openai

                self._client = openai.AsyncOpenAI(api_key=self.api_key)
            except ImportError:
                raise ImportError(
                    "openai SDK not installed. "
                    "Install with: pip install 'autonomous-lab[multi-agent]'"
                )
        return self._client

    def format_tools(self, tool_specs: list[dict]) -> list[dict]:
        """Format tools for OpenAI's function calling format."""
        formatted = []
        for spec in tool_specs:
            formatted.append(
                {
                    "type": "function",
                    "function": {
                        "name": spec["name"],
                        "description": spec.get("description", ""),
                        "parameters": spec.get("parameters", spec.get("input_schema", {})),
                    },
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
        """Send request to OpenAI Chat Completions API."""
        client = self._get_client()

        # Prepend system message
        full_messages = []
        if system_prompt:
            full_messages.append({"role": "system", "content": system_prompt})
        full_messages.extend(messages)

        kwargs: dict = {
            "model": self.model,
            "max_tokens": max_tokens,
            "messages": full_messages,
        }
        if tools:
            kwargs["tools"] = tools

        try:
            response = await client.chat.completions.create(**kwargs)
        except Exception as e:
            logger.error("OpenAI API error: %s", e)
            raise

        choice = response.choices[0] if response.choices else None
        if not choice:
            return LLMResponse(model=self.model)

        message = choice.message
        content_text = message.content or ""
        tool_calls = []

        if message.tool_calls:
            for tc in message.tool_calls:
                arguments = {}
                if tc.function.arguments:
                    try:
                        arguments = json.loads(tc.function.arguments)
                    except json.JSONDecodeError:
                        arguments = {"raw": tc.function.arguments}

                tool_calls.append(
                    ToolCall(
                        id=tc.id,
                        name=tc.function.name,
                        arguments=arguments,
                    )
                )

        usage = response.usage
        return LLMResponse(
            content=content_text,
            tool_calls=tool_calls,
            stop_reason=choice.finish_reason or "",
            input_tokens=getattr(usage, "prompt_tokens", 0) if usage else 0,
            output_tokens=getattr(usage, "completion_tokens", 0) if usage else 0,
            model=response.model or self.model,
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
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json",
        }

        full_messages = []
        if system_prompt:
            full_messages.append({"role": "system", "content": system_prompt})
        full_messages.extend(messages)

        body: dict = {
            "model": self.model,
            "max_tokens": max_tokens,
            "messages": full_messages,
        }
        if tools:
            body["tools"] = tools

        async with aiohttp.ClientSession() as session:
            async with session.post(
                "https://api.openai.com/v1/chat/completions",
                headers=headers,
                json=body,
            ) as resp:
                if resp.status != 200:
                    error_text = await resp.text()
                    raise RuntimeError(
                        f"OpenAI API error {resp.status}: {error_text[:500]}"
                    )
                data = await resp.json()

        choices = data.get("choices", [])
        if not choices:
            return LLMResponse(model=self.model)

        message = choices[0].get("message", {})
        content_text = message.get("content", "") or ""
        tool_calls = []

        for tc in message.get("tool_calls", []):
            fn = tc.get("function", {})
            arguments = {}
            raw_args = fn.get("arguments", "")
            if raw_args:
                try:
                    arguments = json.loads(raw_args)
                except json.JSONDecodeError:
                    arguments = {"raw": raw_args}

            tool_calls.append(
                ToolCall(
                    id=tc.get("id", ""),
                    name=fn.get("name", ""),
                    arguments=arguments,
                )
            )

        usage = data.get("usage", {})
        return LLMResponse(
            content=content_text,
            tool_calls=tool_calls,
            stop_reason=choices[0].get("finish_reason", ""),
            input_tokens=usage.get("prompt_tokens", 0),
            output_tokens=usage.get("completion_tokens", 0),
            model=data.get("model", self.model),
        )
