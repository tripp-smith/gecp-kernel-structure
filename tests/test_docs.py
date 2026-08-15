from pathlib import Path

ROOT = Path(__file__).parents[1]


def test_readme_progress_contract() -> None:
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    required = (
        "Current phase",
        "Phase state",
        "Last verification",
        "Verification command",
        "Claim level",
        "$phase-cadence",
    )
    assert all(field in text for field in required)
    for phase in ("0", "A", "B", "C", "D", "E", "F", "R"):
        assert f"| {phase} |" in text


def test_documentation_links_exist() -> None:
    for name in (
        "APPLICATION.md",
        "FINDINGS.md",
        "MATHLIB.md",
        "RESEARCH.md",
        "SPEC.md",
    ):
        assert (ROOT / name).is_file(), name


def test_research_blockers_are_consistent() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    findings = (ROOT / "FINDINGS.md").read_text(encoding="utf-8")
    research = (ROOT / "RESEARCH.md").read_text(encoding="utf-8")
    assert "| D | Separated approximation | complete |" in readme
    assert "| E | Structural GECP | verified |" in readme
    assert "endpoint-resolved 128-bit finite-grid" in findings
    assert "| C | Fermionic structure/census | complete |" in readme
    assert "expFamily_separatedApprox" in research
    assert "fermionicKernel_separatedApprox" in research
    assert "CrossRatioControl" in research
