# Coding Notes for Agents

## Shared preferences and attribution

Codex and current Claude Code use this repository's `AGENTS.md` directly.
Personal defaults also live in `~/.codex/AGENTS.md`, imported by
`~/.claude/CLAUDE.md`. Keep repository instructions usable without those
personal files. More specific project rules take precedence.

- Never add AI attribution to commits, pull requests, or release notes: no
  `Co-Authored-By:` trailer for an AI, no "Generated with ..." line, and no
  other AI credit.
- Prefer focused changes, preserve useful TODO and debugging comments, and
  report validation performed and work that remains.
- Do not publish packages, approve PyPI deployments, or push release tags
  without an explicit request. Prepare artifacts and commands for the maintainer.
