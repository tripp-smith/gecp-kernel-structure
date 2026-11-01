"""Explicit high-precision same-mesh projection contracts."""

import json

import mpmath as mp
import numpy as np
import pytest

from kernelgecp.integrations.maxent import SpectralDataset, prepare_spectral_projection


class Design:
    beta = 1.0

    def __init__(self, matrix):
        self.matrix = matrix
        self.calls = 0

    def design_tau(self, tau):
        self.calls += 1
        return self.matrix


def high_precision_svd_projector(a):
    with mp.workprec(240):
        matrix = mp.matrix([[mp.mpf(float(x)) for x in row] for row in a])
        u, _, _ = mp.svd(matrix, full_matrices=False)
        p = u * u.T
        return np.array(
            [[float(p[i, j]) for j in range(a.shape[0])] for i in range(a.shape[0])]
        )


def test_high_precision_projection_independent_oracle_and_propagation():
    rng = np.random.default_rng(912)
    q, _ = np.linalg.qr(rng.normal(size=(9, 4)))
    a = q @ np.diag([1, 1e-3, 1e-6, 1e-10])
    original = a.copy()
    design = Design(a)
    tau = np.linspace(0, 1, 9)
    plan = prepare_spectral_projection(
        design,
        tau,
        tau.copy(),
        rcond=1e-12,
        precision_bits=192,
        metadata={"caller": [1, "preserved"]},
    )
    assert design.calls == 1
    assert plan.diagnostics.effective_rank == 4
    assert plan.diagnostics.condition_number > 1e9
    assert plan.transformed_diagnostics.condition_number < 1.000001
    assert np.allclose(
        plan.matrix, high_precision_svd_projector(original), atol=2e-8, rtol=1e-8
    )
    assert np.allclose(plan.matrix, q @ q.T, atol=2e-8, rtol=1e-8)
    assert json.loads(plan.metadata_json) == {"caller": [1, "preserved"]}
    provenance = json.loads(plan.projection_provenance_json)
    assert provenance["precision_bits"] == 192
    assert provenance["method"] == "mpmath_qr_then_binary64_cached_svd"
    assert provenance["source_design_sha256"] != provenance["decomposed_design_sha256"]
    a[:] = 0
    assert np.array_equal(plan.source_design, original)
    for arr in (plan.source_design, plan.decomposed_design):
        with pytest.raises(ValueError):
            arr.setflags(write=True)
    covariance = np.diag(np.linspace(1, 2, 9))
    values = np.linspace(0.1, 0.9, 9)
    data = SpectralDataset(1, tau, values, covariance)
    projected = plan.apply(data)
    assert np.allclose(projected.values, plan.matrix @ values)
    assert np.allclose(projected.covariance, plan.matrix @ covariance @ plan.matrix.T)
    assert json.loads(projected.provenance_json)["projection_arithmetic"] == provenance
    kernel = rng.normal(size=(9, 5)) + 1j * rng.normal(size=(9, 5))
    assert np.allclose(plan.matrix @ (1j * kernel), 1j * (plan.matrix @ kernel))
    assert np.allclose(
        plan.matrix @ kernel,
        high_precision_svd_projector(original) @ kernel,
        atol=2e-8,
        rtol=1e-8,
    )


@pytest.mark.parametrize("precision", [True, 54, 127, 1025, 128.0, "192"])
def test_invalid_precision_rejected(precision):
    tau = np.linspace(0, 1, 3)
    with pytest.raises(ValueError, match="precision_bits"):
        prepare_spectral_projection(
            Design(np.eye(3)), tau, tau, precision_bits=precision
        )


def test_rank_truncation_ambiguity_and_cross_mesh_are_rejected():
    tau = np.linspace(0, 1, 3)
    for singular, message in [(1e-13, "full column rank"), (1.5e-12, "ambiguous rank")]:
        with pytest.raises(ValueError, match=message):
            prepare_spectral_projection(
                Design(np.diag([1.0, 0.5, singular])),
                tau,
                tau,
                rcond=1e-12,
                precision_bits=192,
            )
    with pytest.raises(ValueError, match="ambiguous rank"):
        prepare_spectral_projection(
            Design(np.diag([1.0, 0.5, 1e-16])),
            tau,
            tau,
            rcond=0.0,
            precision_bits=192,
        )
    with pytest.raises(ValueError, match="same mesh"):
        prepare_spectral_projection(
            Design(np.eye(3)), tau, np.array([0.0, 0.25, 1.0]), precision_bits=192
        )


def test_default_path_uses_single_owned_source_design():
    tau = np.linspace(0, 1, 3)
    design = Design(np.eye(3))
    plan = prepare_spectral_projection(design, tau, tau)
    assert design.calls == 1
    assert plan.transformed_diagnostics is None
    assert np.array_equal(plan.source_design, plan.decomposed_design)
    assert (
        json.loads(plan.projection_provenance_json)["method"] == "binary64_cached_svd"
    )
