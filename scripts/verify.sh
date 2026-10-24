#!/usr/bin/env bash
set -euo pipefail

task_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$task_dir"

lake build
scripts/check_axioms.sh
lake env lean scripts/CompletionContract.lean

if rg -n '\bsorry\b|\badmit\b|\bnative_decide\b|^\s*(private\s+|protected\s+)?axiom\s' --glob '*.lean' --glob '!.lake/**'; then
  echo "Lean placeholder, proof escape, or custom axiom declaration found" >&2
  exit 1
fi

uv run ruff check .
uv run ruff format --check .
uv run mypy kernelgecp
uv run python -m pytest
uv run pip-audit
