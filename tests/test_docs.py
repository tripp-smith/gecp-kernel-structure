import re
import tomllib
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
    "COMPLETION_AUDIT.md",
    "FINAL_AUDIT.md",
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
        "COMPLETION_AUDIT.md",
        "FINAL_AUDIT.md",
    ):
        assert (ROOT / name).is_file(), name


def test_distributions_exclude_local_working_notes() -> None:
    config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    excluded = config["tool"]["hatch"]["build"]["exclude"]
    assert "/FINDINGS 2.md" in excluded
    assert "/MATHLIB 2.md" in excluded
    assert "/autonomous-implementation.md" in excluded


def test_completion_audit_accounts_for_every_requirement_group() -> None:
    text = (ROOT / "COMPLETION_AUDIT.md").read_text(encoding="utf-8")
    rows = re.findall(
        r"^\| (R\d{2}) \|.*? \| (Supported|Incomplete) \|", text, re.MULTILINE
    )
    assert {name for name, _ in rows} == {f"R{i:02d}" for i in range(1, 41)}
    summary = re.search(r"\*\*(\d+) supported and (\d+) incomplete\*\*", text)
    assert summary is not None
    assert sum(status == "Supported" for _, status in rows) == int(summary[1])
    assert sum(status == "Incomplete" for _, status in rows) == int(summary[2])


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


def test_independent_audit_covers_original_sections_not_only_inventory_labels() -> None:
    text = (ROOT / "FINAL_AUDIT.md").read_text(encoding="utf-8")
    sections = re.findall(r"^\| (\d+) \|", text, re.MULTILINE)
    groups = re.findall(r"^\| (R\d{2}) \|", text, re.MULTILINE)
    assert len(sections) == 30 and set(map(int, sections)) == set(range(1, 31))
    assert len(groups) == 40 and set(groups) == {f"R{i:02d}" for i in range(1, 41)}


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
