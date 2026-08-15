#!/usr/bin/env bash
set -euo pipefail

task_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$task_dir"

lake build
scripts/check_axioms.sh

if rg -n '\bsorry\b|\badmit\b' --glob '*.lean' --glob '!.lake/**'; then
  echo "Lean placeholder found" >&2
  exit 1
fi

uv run ruff check .
uv run ruff format --check .
uv run mypy kernelgecp
uv run pytest
uv run pip-audit
