import copy
from fractions import Fraction

import mpmath as mp
import numpy as np
import pytest

from experiments.audit_epw_frequency import validate_pole_dictionary
from kernelgecp.integrations._frequency_cross import (
    CauchyFrequencyPlan,
    _cauchy_dictionary,
)
from kernelgecp.types import NumericalBreakdownError


def test_schur_cardinal_independent_and_analytic_pole():
    x = np.array([1.0, 2.0, 3.0, 5.0, 8.0])
    plan = CauchyFrequencyPlan(x, rank=3, dictionary_size=7, pole_cutoff_ev=2)
    poles = _cauchy_dictionary(7, 4)
    residual = 1 / (x[:, None] - poles)
    for i, j, score in zip(
        plan.pivot_rows,
        [np.flatnonzero(poles == p)[0] for p in plan.poles],
        plan.log_residual_history,
        strict=True,
    ):
        assert (i, j) == np.unravel_index(np.argmax(abs(residual)), residual.shape)
        assert np.exp(score) == pytest.approx(abs(residual[i, j]), rel=1e-10)
        residual -= np.outer(residual[:, j], residual[i, :]) / residual[i, j]
    q = np.array([1.25, 4.5, 10.0])
    kernel_nodes = 1 / (plan.nodes[:, None] - plan.poles)
    assert np.allclose(
        plan.design(q) @ kernel_nodes,
        1 / (q[:, None] - plan.poles),
        rtol=1e-12,
        atol=1e-12,
    )
    amplitudes = np.array([0.7, -0.4, 0.2])
    values = (1 / (x[:, None] - plan.poles)) @ amplitudes
    assert np.allclose(
        plan.design(q) @ plan.fit(values),
        (1 / (q[:, None] - plan.poles)) @ amplitudes,
        atol=1e-12,
    )


def test_exact_and_near_nodes_and_large_query():
    plan = CauchyFrequencyPlan([1, 2, 4, 8], rank=3, dictionary_size=8)
    assert np.array_equal(plan.design(plan.nodes), np.eye(3))
    for k, z in enumerate(plan.nodes):
        query = np.nextafter(z, np.inf)
        assert plan.design([query])[0, k] == pytest.approx(1, abs=2e-14)
    assert np.all(np.isfinite(plan.design([1e300])))
    assert (
        plan.evaluation_diagnostics([1.5])[1]["certification_status"] == "not_certified"
    )


def test_fraction_exact_kernel():
    plan = CauchyFrequencyPlan(
        [1.0, 2.0, 4.0], rank=2, dictionary_size=3, pole_cutoff_ev=1
    )
    q = Fraction(3, 2)
    actual = plan.design([float(q)]) @ (1 / (plan.nodes[:, None] - plan.poles))
    expected = [float(1 / (q - Fraction.from_float(float(p)))) for p in plan.poles]
    assert np.allclose(actual[0], expected, atol=1e-14)


def test_payload_ownership_and_corruption():
    x = np.array([1.0, 2.0, 4.0, 8.0])
    plan = CauchyFrequencyPlan(x, rank=3, dictionary_size=8)
    x[0] = 20
    metadata, arrays = plan.to_payload()
    loaded = CauchyFrequencyPlan.from_payload(metadata, arrays)
    arrays["weights"][0] *= 2
    assert not np.array_equal(loaded.weights, arrays["weights"])
    with pytest.raises(ValueError):
        CauchyFrequencyPlan.from_payload(metadata, arrays)
    metadata2, arrays2 = plan.to_payload()
    metadata2 = copy.deepcopy(metadata2)
    metadata2["pivot_rows"][0] = 99
    with pytest.raises(ValueError):
        CauchyFrequencyPlan.from_payload(metadata2, arrays2)
    with pytest.raises(ValueError):
        plan.nodes.setflags(write=True)


def test_complex_fits_cached_and_context_restored():
    old = mp.mp.prec
    mp.mp.prec = 77
    try:
        plan = CauchyFrequencyPlan(
            [1, 2, 3, 4], rank=3, dictionary_size=7, precision_bits=180
        )
        assert mp.mp.prec == 77
        assert plan.prepare() is plan.prepare()
        plan.prepare(rcond=0)
        with pytest.raises(ValueError):
            plan.prepare(rcond=False)
        assert plan.prepare(rcond=1e-8) is plan.prepare(rcond=1e-8)
        values = np.arange(24).reshape(2, 4, 3) * (1 + 2j)
        coef = plan.fit(values, axis=1)
        assert coef.shape == (2, 3, 3)
        assert np.allclose(plan.fit(1j * values, axis=1), 1j * coef)
    finally:
        mp.mp.prec = old


def test_rejections_and_condition_limit():
    for x in ([0, 1], [1, 1], [1, np.inf], [1, True]):
        with pytest.raises(ValueError):
            CauchyFrequencyPlan(x, rank=1)
    with pytest.raises(ValueError):
        CauchyFrequencyPlan([1, 2], rank=3)
    with pytest.raises(NumericalBreakdownError):
        CauchyFrequencyPlan([1e308], rank=1, pole_cutoff_ev=1e154)
    plan = CauchyFrequencyPlan(
        [1, 2, 3, 4, 6], rank=3, dictionary_size=7, max_evaluation_condition=1
    )
    with pytest.raises(NumericalBreakdownError):
        plan.design([2.5])


def test_portable_load_does_not_reselect(monkeypatch):
    plan = CauchyFrequencyPlan([1.0, 2.0, 4.0, 8.0], rank=3, dictionary_size=8)
    metadata, arrays = plan.to_payload()
    metadata["version"] = 1
    del metadata["dictionary_policy"]
    del metadata["selected_pivot_columns"]
    pole = next(i for i, p in enumerate(arrays["poles"]) if p < 0)
    arrays["poles"][pole] = np.nextafter(arrays["poles"][pole], 0.0)
    with mp.workprec(metadata["precision_bits"]):
        nodes = [mp.mpf(float(v)) for v in arrays["nodes"]]
        poles = [mp.mpf(float(v)) for v in arrays["poles"]]
        arrays["weights"] = np.array(
            [
                float(
                    mp.fprod(z - p for p in poles)
                    / mp.fprod(z - v for i, v in enumerate(nodes) if i != j)
                )
                for j, z in enumerate(nodes)
            ]
        )
    monkeypatch.setattr(np, "geomspace", lambda *a, **k: pytest.fail("reselected"))
    loaded = CauchyFrequencyPlan.from_payload(metadata, arrays)
    assert np.array_equal(loaded.poles, arrays["poles"])
    assert np.array_equal(loaded.weights, arrays["weights"])


def test_sealed_parameters_immutable_results_and_bounded_cache():
    plan = CauchyFrequencyPlan([1.0, 2.0, 4.0, 8.0], rank=3, dictionary_size=8)
    for name, value in (
        ("rank", 1),
        ("weights", np.ones(3)),
        ("nodes", np.ones(3)),
        ("max_evaluation_condition", 1),
    ):
        with pytest.raises(AttributeError):
            setattr(plan, name, value)
    with pytest.raises(AttributeError):
        del plan.rank
    for array in (plan.design([1.5]), plan.fit([1, 2, 3, 4])):
        with pytest.raises(ValueError):
            array.setflags(write=True)
    for i in range(12):
        plan.prepare(rcond=i * 1e-8)
    assert len(plan._plans) <= 8
    assert plan.prepare(rcond=11e-8) is plan.prepare(rcond=11e-8)


def test_reject_bad_saved_model():
    plan = CauchyFrequencyPlan([1.0, 2.0, 4.0, 8.0], rank=3, dictionary_size=8)
    for key, bad in (("nodes", -1.0), ("poles", 1.0), ("weights", 0.0)):
        meta, arrays = plan.to_payload()
        arrays[key][0] = bad
        with pytest.raises(ValueError):
            CauchyFrequencyPlan.from_payload(meta, arrays)
    meta, arrays = plan.to_payload()
    meta["near_tie_counts"][0] = 0
    with pytest.raises(ValueError):
        CauchyFrequencyPlan.from_payload(meta, arrays)


@pytest.mark.parametrize("size", [2, 3, 7, 101, 1201])
def test_v2_dictionary_independent_decimal_and_endpoints(size):
    dictionary = _cauchy_dictionary(size, 169.0)
    assert dictionary[0] == 0 and dictionary[1] == -1e-10
    if size > 2:
        assert dictionary[-1] == -169.0
    plan = {
        "version": 2,
        "rank": size,
        "precision_bits": 256,
        "max_evaluation_condition": 1e4,
        "pivot_rows": [],
        "log_residual_history": [],
        "tie_counts": [],
        "near_tie_counts": [],
        "certification_status": "not_certified",
        "dictionary_size": size,
        "pole_cutoff_ev": 13,
        "dictionary_policy": "binary64_endpoints_mpmath256_log_linear_v1",
        "selected_pivot_columns": list(range(size)),
    }
    validate_pole_dictionary(plan, dictionary)
    assert np.isfinite(dictionary).all()


def test_grid_duplicate_and_payload_v1_roundtrip():
    with pytest.raises(ValueError, match="duplicate"):
        _cauchy_dictionary(100, np.nextafter(1e-10, np.inf))
    plan = CauchyFrequencyPlan([1.0, 2.0, 3.0, 4.0], rank=3, dictionary_size=7)
    metadata, arrays = plan.to_payload()
    assert metadata["version"] == 2
    legacy = copy.deepcopy(metadata)
    legacy["version"] = 1
    del legacy["dictionary_policy"]
    del legacy["selected_pivot_columns"]
    loaded = CauchyFrequencyPlan.from_payload(legacy, arrays)
    actual_meta, actual_arrays = loaded.to_payload()
    assert actual_meta == legacy
    for key in arrays:
        assert np.array_equal(actual_arrays[key], arrays[key])
    assert loaded.selected_pivot_columns is None
    for key, change in (
        ("dictionary_policy", "wrong"),
        ("selected_pivot_columns", [0, 0, 0]),
        ("selected_pivot_columns", [True, 1, 2]),
        ("selected_pivot_columns", [0, 1, 2]),
    ):
        bad = copy.deepcopy(metadata)
        bad[key] = change
        if key == "selected_pivot_columns" and change == metadata[key]:
            continue
        with pytest.raises(ValueError):
            CauchyFrequencyPlan.from_payload(bad, arrays)
    corrupt = {k: v.copy() for k, v in arrays.items()}
    corrupt["poles"][0] = np.nextafter(corrupt["poles"][0], -np.inf)
    with pytest.raises(ValueError):
        CauchyFrequencyPlan.from_payload(metadata, corrupt)


def test_v2_dictionary_rejects_changed_zero_sign():
    plan = CauchyFrequencyPlan([1.0, 2.0, 3.0, 4.0], rank=3, dictionary_size=7)
    metadata, arrays = plan.to_payload()
    selected = metadata["selected_pivot_columns"]
    if 0 not in selected:
        pytest.skip("fixture did not select the zero pole")
    arrays["poles"][selected.index(0)] = 0.0
    with pytest.raises(ValueError, match="dictionary columns"):
        CauchyFrequencyPlan.from_payload(metadata, arrays)
    with pytest.raises(ValueError, match="dictionary ordinals"):
        validate_pole_dictionary(metadata, arrays["poles"])


def test_dictionary_restores_mpmath_context_on_success_and_failure():
    old = mp.mp.prec
    mp.mp.prec = 77
    try:
        _cauchy_dictionary(7, 169.0)
        assert mp.mp.prec == 77
        with pytest.raises(ValueError, match="duplicate"):
            _cauchy_dictionary(100, np.nextafter(1e-10, np.inf))
        assert mp.mp.prec == 77
    finally:
        mp.mp.prec = old
