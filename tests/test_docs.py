from pathlib import Path

ROOT = Path(__file__).parents[1]
PROJECT_DOCS = (
    "AGENTS.md",
    "APPLICATION.md",
    "FINDINGS.md",
    "MATHLIB.md",
    "README.md",
    "RESEARCH.md",
    "SPEC.md",
)


def test_readme_progress_contract() -> None:
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    required = (
        "Current phase",
        "Phase state",
        "Last verification",
        "Verification command",
        "Claim level",
        "$phase-cadence",
        "Python implementation",
        "Why this implementation demonstrates the project objectives",
        "Realistic synthetic applications",
    )
    assert all(field in text for field in required)
    for phase in ("0", "A", "B", "C", "D", "E", "F", "R", "S"):
        assert f"| {phase} |" in text


def test_documentation_links_exist() -> None:
    for name in (
        "APPLICATION.md",
        "FINAL_HANDOFF.md",
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
    assert "| E | Structural GECP | complete |" in readme
    assert "| R | Release readiness | complete |" in readme
    assert "| S | Post-v1 synthetic applications | complete |" in readme
    assert "endpoint-resolved 128-bit finite-grid" in findings
    assert "| C | Fermionic structure/census | complete |" in readme
    assert "expFamily_separatedApprox" in research
    assert "fermionicKernel_separatedApprox" in research
    assert "CrossRatioControl" in research
    assert "run_synthetic_application_suite" in research
    assert (ROOT / "experiments/data/synthetic_applications.json").is_file()
    assert (ROOT / "experiments/data/synthetic_applications.png").is_file()


def test_implementation_run_metadata_is_complete_and_linked() -> None:
    handoff = (ROOT / "FINAL_HANDOFF.md").read_text(encoding="utf-8")
    required = (
        "Implementation run metadata",
        "Declared model identity",
        "Exact serving model/version",
        "Collaboration mode",
        "Goal-accounting usage at completion",
        "1,801,132",
        "10,747 seconds (2:59:07)",
        "$18.01",
        "$9.01–$54.03",
        "API-equivalent estimates",
    )
    assert all(field in handoff for field in required)

    metadata_link = "FINAL_HANDOFF.md#implementation-run-metadata"
    for name in PROJECT_DOCS:
        text = (ROOT / name).read_text(encoding="utf-8")
        assert metadata_link in text, name
