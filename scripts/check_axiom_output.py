"""Fail-closed validation of Lean dependency reports for the public audit list."""

from __future__ import annotations

import re
import sys
from collections import Counter
from collections.abc import Sequence
from pathlib import Path

ALLOWED_AXIOMS = frozenset({"propext", "Classical.choice", "Quot.sound"})
RECORD = re.compile(
    r"'(?P<name>[^'\n]+)' (?:depends on axioms:\s*\[(?P<axioms>[^\]]*)\]"
    r"|does not depend on any axioms)"
)


def validate_axiom_output(output: str, expected: Sequence[str]) -> int:
    """Require every requested declaration and only the permitted axioms."""
    if not expected:
        raise ValueError("the audit declaration list is empty")
    matches = list(RECORD.finditer(output))
    observed = Counter(match.group("name") for match in matches)
    required = Counter(expected)
    if observed != required:
        missing = list((required - observed).elements())
        unexpected = list((observed - required).elements())
        raise ValueError(
            f"declaration coverage mismatch: missing={missing}, extra={unexpected}"
        )
    remainder = RECORD.sub("", output)
    if "depends on axioms" in remainder or "does not depend on any axioms" in remainder:
        raise ValueError("unparsed axiom dependency report")
    for match in matches:
        axioms = {
            value.strip()
            for value in (match.group("axioms") or "").split(",")
            if value.strip()
        }
        disallowed = axioms - ALLOWED_AXIOMS
        if disallowed:
            raise ValueError(
                f"{match.group('name')} uses disallowed axioms: {sorted(disallowed)}"
            )
    return len(matches)


def main() -> None:
    source = Path(sys.argv[1]).read_text(encoding="utf-8")
    expected = re.findall(r"^#print axioms (\S+)\s*$", source, re.MULTILINE)
    output = sys.stdin.read()
    try:
        count = validate_axiom_output(output, expected)
    except ValueError as error:
        print(output, end="")
        print(f"Axiom audit failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
    print(f"Axiom audit passed: {count} declarations; permitted axioms only.")


if __name__ == "__main__":
    main()
