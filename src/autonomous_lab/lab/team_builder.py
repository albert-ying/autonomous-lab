"""
Autonomous Lab - Agent Team Builder

Generates agent-team instructions for Claude Code's native teams feature.
Instead of spawning subprocesses, this builds a prompt that tells the
Claude Code session to create a native agent team where the PI is the
team lead and trainees are teammates with their own context windows.

Used automatically in `orchestration: multi` when Claude Code's agent
teams feature is detected (CLAUDE_CODE_EXPERIMENTAL_AGENT_TEAMS=1).
Falls back to CLI subprocess spawning when not available.
"""

import logging
import os

from .state import load_config, load_idea, scan_project_files

logger = logging.getLogger(__name__)


def is_agent_team_available() -> bool:
    """Check if Claude Code's native agent teams feature is available.

    Detects via environment variable set by Claude Code.
    """
    return os.environ.get("CLAUDE_CODE_EXPERIMENTAL_AGENT_TEAMS") == "1"


class TeamBuilder:
    """Builds agent-team instructions for Claude Code's native teams feature.

    Reads trainee configs from `agents.trainees` (same section as the
    Orchestrator) and team-specific settings from the `team` section.
    """

    def __init__(self, project_dir: str):
        self.project_dir = project_dir
        self.delegate_mode: bool = True
        self.require_plan_approval: bool = False
        self.teammate_model: str = "sonnet"
        self.trainees: list[dict] = []
        self._load_config()

    def _load_config(self) -> None:
        """Load config from config.yaml.

        Trainee definitions come from agents.trainees (shared with Orchestrator).
        Team-specific settings come from the team section.
        """
        config = load_config(self.project_dir)

        # Team-specific settings
        team_cfg = config.get("team", {})
        self.delegate_mode = team_cfg.get("delegate_mode", True)
        self.require_plan_approval = team_cfg.get(
            "require_plan_approval", False
        )
        self.teammate_model = team_cfg.get("teammate_model", "sonnet")

        # Trainee definitions — shared with Orchestrator
        agents_cfg = config.get("agents", {})
        if "trainees" in agents_cfg:
            self.trainees = [
                {
                    "name": t.get("name", "Trainee"),
                    "focus": t.get("focus", "general tasks"),
                    "model": t.get("model", self.teammate_model),
                }
                for t in agents_cfg["trainees"]
            ]
        elif "trainee" in agents_cfg:
            self.trainees = [
                {
                    "name": "Trainee",
                    "focus": "technical implementation, data analysis, and scientific writing",
                    "model": agents_cfg["trainee"].get(
                        "model", self.teammate_model
                    ),
                }
            ]
        else:
            # Default single trainee
            self.trainees = [
                {
                    "name": "Trainee",
                    "focus": "technical implementation, data analysis, and scientific writing",
                    "model": self.teammate_model,
                }
            ]

    def build_team_prompt(
        self,
        pi_agenda: str,
        idea: str,
        file_listings: str | dict,
        iteration: int,
    ) -> str:
        """Build instructions for the MCP client to create an agent team.

        Returns a prompt that the AI (in Claude Code) should follow
        to spawn teammates and coordinate work via the native agent
        teams feature.

        Args:
            pi_agenda: The PI's agenda/directives for this iteration.
            idea: The project idea description.
            file_listings: Current file listings (str or dict).
            iteration: Current iteration number.

        Returns:
            Formatted prompt string with team-building instructions.
        """
        # Format file listings if dict
        if isinstance(file_listings, dict):
            parts = []
            for category, files in file_listings.items():
                if files:
                    parts.append(f"**{category}/**: {', '.join(files)}")
            file_listings_str = "\n".join(parts) if parts else "(no files yet)"
        else:
            file_listings_str = str(file_listings)

        # Build trainee specs
        trainee_specs = []
        for t in self.trainees:
            model = t.get("model", self.teammate_model)
            spec = (
                f'- "{t["name"]}": Focus on {t.get("focus", "general tasks")}. '
                f"Model: {model}."
            )
            if self.require_plan_approval:
                spec += " Require plan approval before implementation."
            trainee_specs.append(spec)

        delegate_note = ""
        if self.delegate_mode:
            delegate_note = (
                "\nIMPORTANT: Use delegate mode (Shift+Tab). "
                "As PI, you ONLY coordinate -- do NOT implement yourself. "
                "Assign all implementation tasks to your teammates.\n"
            )

        approval_note = ""
        if self.require_plan_approval:
            approval_note = (
                "\nEach trainee must share their plan before implementing. "
                "Review and approve plans before they proceed.\n"
            )

        return f"""[AUTOLAB AGENT-TEAM MODE] -- Iteration {iteration}

You are the PI (team lead). Create an agent team to execute this iteration's work.

## PI Agenda
{pi_agenda}

## Project Idea
{idea}

## Current Files
{file_listings_str}

## Team Structure
Create an agent team with these teammates:
{chr(10).join(trainee_specs)}
{delegate_note}{approval_note}
## Instructions
1. Create the agent team with the teammates listed above
2. Break the PI agenda into tasks on the shared task list
3. Assign tasks to teammates based on their focus areas
4. Teammates will self-coordinate and claim unassigned tasks
5. Once all tasks complete, synthesize results
6. Call autolab_record with the combined summary and content
7. Then call lab_meeting, then autolab_next to continue

## Communication
- Teammates can message each other directly to coordinate
- If a trainee is stuck, they can ask another trainee for help
- Monitor progress and redirect approaches that aren't working

## Project Directory
{self.project_dir}

## MCP Tools Available
All teammates have access to the autolab MCP tools (autolab_cite, etc.)
and standard file/shell tools. They should write code to scripts/,
figures to figures/, results to results/, and paper sections to paper/sections/.
"""

    def build_team_prompt_from_state(self, pi_agenda: str) -> str:
        """Convenience: build team prompt loading idea and files from state.

        Args:
            pi_agenda: PI's agenda for this iteration.

        Returns:
            Formatted team prompt string.
        """
        from .state import load_state

        idea = load_idea(self.project_dir)
        file_listings = scan_project_files(self.project_dir)
        state = load_state(self.project_dir)
        iteration = state.get("iteration", 0)

        return self.build_team_prompt(
            pi_agenda=pi_agenda,
            idea=idea,
            file_listings=file_listings,
            iteration=iteration,
        )
