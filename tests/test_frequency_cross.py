import copy
from fractions import Fraction

import mpmath as mp
import numpy as np
import pytest

from kernelgecp.integrations._frequency_cross import CauchyFrequencyPlan
from kernelgecp.types import NumericalBreakdownError


def test_schur_cardinal_independent_and_analytic_pole():
    x = np.array([1.0, 2.0, 3.0, 5.0, 8.0])
    plan = CauchyFrequencyPlan(x, rank=3, dictionary_size=7, pole_cutoff_ev=2)
    poles = -np.r_[0.0, np.geomspace(1e-10, 4, 6)]
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
