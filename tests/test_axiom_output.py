import os
import subprocess
import sys
from pathlib import Path

import pytest

from scripts.check_axiom_output import validate_axiom_output


def test_allowed_multiline_and_axiom_free_reports() -> None:
    output = (
        "'A' depends on axioms: [propext,\n Classical.choice,\n Quot.sound]\n"
        "'B' does not depend on any axioms\n"
    )
    assert validate_axiom_output(output, ["A", "B"]) == 2


@pytest.mark.parametrize("axiom", ["unproved_claim", "Lean.ofReduceBool", "sorryAx"])
def test_unexpected_axioms_are_rejected(axiom: str) -> None:
    with pytest.raises(ValueError, match="disallowed axioms"):
        validate_axiom_output(f"'A' depends on axioms: [propext, {axiom}]", ["A"])


@pytest.mark.parametrize(
    "output",
    [
        "",
        "'A' depends on axioms: [propext",
        "'B' does not depend on any axioms",
        "'A' does not depend on any axioms\n'A' does not depend on any axioms",
    ],
)
def test_missing_malformed_or_extra_reports_fail_closed(output: str) -> None:
    with pytest.raises(ValueError, match="coverage mismatch"):
        validate_axiom_output(output, ["A"])


def test_unknown_report_format_and_empty_manifest_fail_closed() -> None:
    with pytest.raises(ValueError, match="unparsed"):
        validate_axiom_output(
            "'A' does not depend on any axioms\n'X' depends on axioms: INVALID", ["A"]
        )
    with pytest.raises(ValueError, match="empty"):
        validate_axiom_output("", [])


def test_cli_fails_on_disallowed_axiom(tmp_path: Path) -> None:
    manifest = tmp_path / "audit.lean"
    manifest.write_text("#print axioms A\n", encoding="utf-8")
    script = Path(__file__).parents[1] / "scripts/check_axiom_output.py"
    result = subprocess.run(
        [sys.executable, str(script), str(manifest)],
        input="'A' depends on axioms: [unproved_claim]",
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 1
    assert "disallowed axioms" in result.stderr


def test_shell_gate_preserves_lean_failure_even_with_valid_reports(
    tmp_path: Path,
) -> None:
    root = Path(__file__).parents[1]
    shim = tmp_path / "lake"
    shim.write_text(
        f"#!{sys.executable}\n"
        "import pathlib,re,sys\n"
        "source=pathlib.Path('scripts/AxiomAudit.lean').read_text()\n"
        "for name in re.findall(r'^#print axioms (\\S+)',source,re.MULTILINE):\n"
        "    print(f\"'{name}' does not depend on any axioms\")\n"
        "sys.exit(42)\n",
        encoding="utf-8",
    )
    shim.chmod(0o755)
    result = subprocess.run(
        ["bash", "scripts/check_axioms.sh"],
        cwd=root,
        env={**os.environ, "PATH": f"{tmp_path}{os.pathsep}{os.environ['PATH']}"},
        text=True,
        capture_output=True,
        check=False,
    )
    assert "Axiom audit passed" in result.stdout
    assert result.returncode == 42
