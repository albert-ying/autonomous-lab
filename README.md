# Autonomous Lab

MCP server for iterative PI/Trainee research sessions with a monitoring web UI, producing submission-ready LaTeX papers.

## Usage

Add to your Cursor MCP settings:

```json
{
  "mcpServers": {
    "autonomous-lab": {
      "command": "uv",
      "args": ["--directory", "/path/to/autonomous-lab", "run", "autonomous-lab"],
      "timeout": 600,
      "env": {
        "MCP_WEB_PORT": "8766"
      }
    }
  }
}
```

## MCP Tools

- `autolab_init` -- Initialize a research project
- `autolab_next` -- Get the next role prompt (PI or Trainee)
- `autolab_record` -- Record a turn and show monitoring UI
- `autolab_status` -- Get project status
- `interactive_feedback` -- General feedback collection (inherited)
