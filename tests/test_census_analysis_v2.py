from decimal import Decimal

import pytest

from kernelgecp.research import analyze_census_record


def _record() -> dict[str, object]:
    return {
        "schema_version": 2,
        "cutoff": "1",
        "tolerance": "1e-12",
        "rank": 5,
        "pivots": ["0.999999", "0.799999", "0.399999", "0.199999", "0.099999"],
        "residual": "0.05",
        "residual_history": ["1", "0.8", "0.4", "0.2", "0.1", "0.05"],
    }


def test_schema_two_analysis_uses_sampled_maxima_not_nearmaximal_pivots() -> None:
    record = _record()
    evidence = analyze_census_record(record)
    assert Decimal(evidence.first_step_observed) == Decimal("0.8")
    assert evidence.complete_block_ratios == ("0.4", "0.25")
    assert evidence.trailing_steps == 1
    assert evidence.trailing_ratio == "0.5"
    assert evidence.maximum_complete_block_ratio == "0.4"
    assert evidence.continuous_certified is False
    legacy = analyze_census_record(record | {"schema_version": 1})
    assert legacy.first_step_observed != evidence.first_step_observed
    assert legacy.complete_block_ratios != evidence.complete_block_ratios


@pytest.mark.parametrize(
    "history",
    [
        None,
        ["1", "0.8", "0.4", "0.2", "0.05"],
        ["1", "0.8", "0.4", "0.2", "0.1", "0.04"],
        ["1", "0.8", "0", "0.2", "0.1", "0.05"],
        ["1", "0.8", "nan", "0.2", "0.1", "0.05"],
    ],
)
def test_schema_two_analysis_never_falls_back_to_pivot_proxy(history: object) -> None:
    with pytest.raises(ValueError, match="residual|preterminal"):
        analyze_census_record(_record() | {"residual_history": history})
