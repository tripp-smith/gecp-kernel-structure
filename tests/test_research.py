from decimal import Decimal
from pathlib import Path

import pytest

from kernelgecp.research import (
    analyze_census_file,
    analyze_census_record,
    certify_fermionic_block_contraction,
    dyadic_cutoff_scale,
    fermionic_first_step_ratio,
)


def test_first_step_formula_obstructs_uniform_contraction() -> None:
    assert dyadic_cutoff_scale(Decimal("1")) == 0
    assert dyadic_cutoff_scale(Decimal("10")) == 4
    assert fermionic_first_step_ratio(Decimal("1")) == Decimal(
        "0.86466471676338730810600050502751559659236845409042411853184112734592662589851231"
    )
    assert fermionic_first_step_ratio(Decimal("100")) > Decimal("0.999999999999999999")


def test_record_analysis_separates_complete_and_trailing_blocks() -> None:
    evidence = analyze_census_record(
        {
            "cutoff": "1",
            "tolerance": "1e-8",
            "pivots": ["1", "0.8", "0.4", "0.2", "0.1"],
            "residual": "0.05",
        }
    )
    assert evidence.block_size == 2
    assert evidence.complete_block_ratios == ("0.4", "0.25")
    assert evidence.trailing_steps == 1
    assert evidence.trailing_ratio == "0.5"
    assert evidence.continuous_certified is False


def test_canonical_census_supports_but_does_not_certify_block_conjecture() -> None:
    root = Path(__file__).parents[1]
    evidence = analyze_census_file(root / "experiments/data/gecp_census.jsonl")
    assert [item.cutoff for item in evidence] == [
        "1",
        "10",
        "100",
        "1000",
        "10000",
        "100000",
        "1000000",
    ]
    assert all(item.complete_block_ratios for item in evidence)
    assert all(
        Decimal(item.maximum_complete_block_ratio) < Decimal("0.5")
        for item in evidence
        if item.maximum_complete_block_ratio is not None
    )
    assert not any(item.continuous_certified for item in evidence)


def test_interval_certificate_closes_first_cutoff_block() -> None:
    evidence = certify_fermionic_block_contraction(
        1.0,
        2,
        precision_bits=128,
        abs_tol=1e-8,
        rel_tol=1e-6,
        max_cells=20_000,
    )
    assert evidence.trajectory_certified
    assert evidence.half_contraction_certified
    assert evidence.completed_steps == 2
    assert evidence.contraction_ratio_upper < 0.078
    assert len(evidence.pivot_certificates) == 3


@pytest.mark.parametrize("cutoff", [Decimal("0"), Decimal("Infinity")])
def test_research_diagnostics_reject_invalid_cutoffs(cutoff: Decimal) -> None:
    with pytest.raises(ValueError):
        fermionic_first_step_ratio(cutoff)
