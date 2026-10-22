import json
import math
from pathlib import Path
from typing import Any

import pytest

from experiments import plot_census as plotting

plotting.plt.switch_backend("Agg")


def _record(
    cutoff: str = "4",
    tolerance: str = "1e-12",
    history: tuple[str, ...] = ("1", "0.1", "0.001"),
    precision: int = 128,
) -> dict[str, Any]:
    return {
        "schema_version": 2,
        "cutoff": cutoff,
        "tolerance": tolerance,
        "rank": len(history) - 1,
        "precision_bits": precision,
        "residual_history": list(history),
    }


def test_normalized_curves_include_every_prefix_and_true_log1p() -> None:
    record = _record(cutoff="1e-20", history=("1", "0.1", "1e-12", "0"))
    curve = plotting.normalized_curves([record])[0]
    assert curve.ranks == (0, 1, 2, 3)
    assert curve.normalized_ranks == pytest.approx(
        tuple(index / math.log1p(1e-20) for index in range(4))
    )
    assert curve.residuals == (1, 0.1, 1e-12, 0)
    assert curve.zero_indices == (3,)


def test_normalized_selects_one_deepest_compatible_trajectory_per_cutoff() -> None:
    deepest = _record()
    prefix = _record(tolerance="1e-6", history=("1", "0.1"))
    same_cutoff = _record(cutoff="4.00", tolerance="1e-8", history=("1", "0.1"))
    other = _record(cutoff="16", precision=192)
    records = [other, prefix, same_cutoff, deepest, deepest.copy()]
    curves = plotting.normalized_curves(records)
    assert [curve.cutoff for curve in curves] == [4, 16]
    assert [curve.tolerance for curve in curves] == [1e-12, 1e-12]
    assert [curve.precision_bits for curve in curves] == [128, 192]
    assert curves == plotting.normalized_curves(list(reversed(records)))


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"schema_version": 1}, "schema_version 2"),
        ({"residual_history": None}, r"rank \+ 1"),
        ({"residual_history": ["1"]}, r"rank \+ 1"),
        ({"residual_history": [1, 0.1, 0]}, "decimal strings"),
        ({"residual_history": ["1", "NaN", "0"]}, "finite"),
        ({"residual_history": ["1", "-0.1", "0"]}, "nonnegative"),
        ({"residual_history": ["1", "1e-400", "0"]}, "plotting range"),
        ({"cutoff": "0"}, "positive"),
        ({"cutoff": "-1"}, "positive"),
        ({"cutoff": "infinity"}, "finite"),
        ({"cutoff": "1e400"}, "plotting range"),
        ({"tolerance": "0"}, "positive"),
        ({"rank": True}, "integer"),
        ({"rank": 1.5}, "integer"),
        ({"rank": -1}, "nonnegative"),
        ({"precision_bits": 0}, "positive"),
    ],
)
def test_normalized_rejects_invalid_records(
    change: dict[str, Any], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        plotting.normalized_curves([_record() | change])


def test_normalized_rejects_missing_or_empty_input() -> None:
    with pytest.raises(ValueError, match="at least one"):
        plotting.normalized_curves([])
    record = _record()
    del record["residual_history"]
    with pytest.raises(ValueError, match="requires residual_history"):
        plotting.normalized_curves([record])


@pytest.mark.parametrize(
    "prefix",
    [
        _record(tolerance="1e-6", history=("1", "0.2")),
        _record(tolerance="1e-6", history=("1", "0.1"), precision=192),
        _record(tolerance="1e-6", history=("1", "0.1", "0.001", "0")),
    ],
)
def test_normalized_rejects_incompatible_prefixes(prefix: dict[str, Any]) -> None:
    with pytest.raises(ValueError, match="incompatible.*cutoff 4"):
        plotting.normalized_curves([_record(), prefix])


def test_normalized_plot_masks_and_discloses_zeros_and_closes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    observed: dict[str, Any] = {}
    original_save = plotting.plt.Figure.savefig

    def inspect_save(figure: Any, *args: Any, **kwargs: Any) -> None:
        axis = figure.axes[0]
        observed["title"] = axis.get_title()
        observed["label"] = axis.get_ylabel()
        observed["texts"] = [item.get_text() for item in axis.texts]
        observed["ys"] = axis.lines[0].get_ydata()
        observed["yscale"] = axis.get_yscale()
        original_save(figure, *args, **kwargs)

    monkeypatch.setattr(plotting.plt.Figure, "savefig", inspect_save)
    before = plotting.plt.get_fignums()
    output = tmp_path / "plots" / "normalized.png"
    summary = plotting.plot_normalized(
        [_record(history=("1", "0.1", "0"), precision=192)], output
    )
    assert output.read_bytes().startswith(b"\x89PNG")
    assert summary.curve_count == 1
    assert summary.omitted_zero_count == 1
    assert summary.precision_bits == (192,)
    assert "192 bits" in observed["title"]
    assert "not continuum certificates" in observed["title"]
    assert "sampled finite-grid" in observed["label"]
    assert "omitted" in observed["texts"][0]
    assert "no positive display floor" in observed["texts"][0]
    assert math.isnan(observed["ys"][-1])
    assert observed["yscale"] == "log"
    assert plotting.plt.get_fignums() == before


def test_normalized_all_zero_plot_reports_no_positive_samples(tmp_path: Path) -> None:
    output = tmp_path / "zeros.svg"
    summary = plotting.plot_normalized([_record(history=("0", "0"))], output)
    assert summary.omitted_zero_count == 2
    assert "no positive display floor" in output.read_text(encoding="utf-8")


def test_default_cli_accepts_legacy_rank_records(tmp_path: Path) -> None:
    source = tmp_path / "legacy.jsonl"
    record = _record() | {"schema_version": 1}
    del record["residual_history"]
    source.write_text(json.dumps(record) + "\n", encoding="utf-8")
    output = tmp_path / "ranks.svg"
    plotting.main([str(source), str(output)])
    text = output.read_text(encoding="utf-8")
    assert "finite-grid GECP rank" in text
    assert "observed precision: 128 bits" in text
    with pytest.raises(SystemExit, match="2"):
        plotting.main([str(source), str(output), "--mode", "normalized"])


def test_normalized_cli_and_jsonl_validation(tmp_path: Path) -> None:
    source = tmp_path / "records.jsonl"
    source.write_text(json.dumps(_record()) + "\n\n", encoding="utf-8")
    output = tmp_path / "normalized.svg"
    plotting.main([str(source), str(output), "--mode", "normalized"])
    assert "k / log(1 + cutoff)" in output.read_text(encoding="utf-8")
    for bad_text, message in [
        ("\n", "at least one"),
        ("[]\n", "must be an object"),
        ("{\n", "invalid census JSON at line 1"),
    ]:
        source.write_text(bad_text, encoding="utf-8")
        with pytest.raises(ValueError, match=message):
            plotting.read_census(source)


def test_failed_save_does_not_leak_a_figure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def fail_save(*args: Any, **kwargs: Any) -> None:
        raise OSError("test save failure")

    monkeypatch.setattr(plotting.plt.Figure, "savefig", fail_save)
    before = plotting.plt.get_fignums()
    with pytest.raises(OSError, match="test save failure"):
        plotting.plot_normalized([_record()], tmp_path / "normalized.png")
    with pytest.raises(OSError, match="test save failure"):
        plotting.plot_ranks([_record()], tmp_path / "ranks.png")
    assert plotting.plt.get_fignums() == before
