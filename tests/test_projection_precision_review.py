"""Regression fragment for projection resource, precision and provenance contracts."""

import hashlib
import json

import mpmath as mp
import numpy as np
import pytest

import kernelgecp.integrations.maxent as maxent


class FixedDesign:
    beta = 1.0

    def __init__(self, matrix):
        self.matrix = matrix

    def design_tau(self, tau):
        return self.matrix


class PolynomialDesign:
    beta = 1.0

    def design_tau(self, tau):
        tau = np.asarray(tau)
        return np.column_stack((np.ones_like(tau), tau, tau**2))


@pytest.mark.parametrize("rows,columns", [(257, 2), (129, 129)])
def test_oversize_high_precision_rejected_before_svd(monkeypatch, rows, columns):
    def forbidden(*args, **kwargs):
        pytest.fail("SVD constructed before the high-precision resource guard")

    monkeypatch.setattr(maxent, "_SpectralLeastSquaresPlan", forbidden)
    tau = np.linspace(0, 1, rows)
    with pytest.raises(ValueError, match="resource limits"):
        maxent.prepare_spectral_projection(
            FixedDesign(np.ones((rows, columns))), tau, tau, precision_bits=192
        )


def test_nearly_collinear_columns_against_independent_mp_svd():
    rng = np.random.default_rng(314159)
    q, _ = np.linalg.qr(rng.normal(size=(12, 4)))
    # Columns share a dominant direction rather than merely differing in scale.
    a = np.column_stack(
        (
            q[:, 0],
            q[:, 0] + 1e-5 * q[:, 1],
            q[:, 0] + 1e-7 * q[:, 2],
            q[:, 0] + 1e-9 * q[:, 3],
        )
    )
    tau = np.linspace(0, 1, len(a))
    plan = maxent.prepare_spectral_projection(
        FixedDesign(a), tau, tau, precision_bits=256, rcond=1e-13
    )
    with mp.workprec(300):
        source = mp.matrix([[mp.mpf(float(x)) for x in row] for row in a])
        u, _, _ = mp.svd(source, full_matrices=False)
        oracle = u * u.T
        expected = np.array(
            [[float(oracle[i, j]) for j in range(len(a))] for i in range(len(a))]
        )
    assert np.max(np.abs(plan.matrix - expected)) <= 1e-13


def canonical_hash(array):
    return hashlib.sha256(
        np.ascontiguousarray(array, dtype="<f8").tobytes()
    ).hexdigest()


@pytest.mark.parametrize("same_mesh,precision", [(True, 192), (True, 53), (False, 53)])
def test_design_hashes_bind_exposed_arrays_and_shapes(same_mesh, precision):
    source = np.linspace(0, 1, 7)
    target = source.copy() if same_mesh else np.linspace(0, 1, 5)
    plan = maxent.prepare_spectral_projection(
        PolynomialDesign(), source, target, precision_bits=precision
    )
    provenance = json.loads(plan.projection_provenance_json)
    for label, array in (
        ("source_design", plan.source_design),
        ("decomposed_design", plan.decomposed_design),
    ):
        assert provenance[label + "_sha256"] == canonical_hash(array)
        assert provenance[label + "_shape"] == list(array.shape)
    if not same_mesh:
        np.testing.assert_array_equal(plan.decomposed_design, plan.source_design.T)


def test_projected_provenance_survives_archive_round_trip(tmp_path):
    tau = np.linspace(0, 1, 7)
    plan = maxent.prepare_spectral_projection(
        PolynomialDesign(), tau, tau, precision_bits=192, metadata={"basis": "poly"}
    )
    data = maxent.SpectralDataset(1, tau, np.linspace(1, 2, 7), np.eye(7))
    projected = plan.apply(data)
    path = tmp_path / "projected.npz"
    maxent.save_spectral_dataset(
        projected, path, units={"tau": "1"}, provider={"name": "test"}
    )
    restored = maxent.load_spectral_dataset(path)
    assert restored.provenance == projected.provenance
    assert restored.provenance["projection_arithmetic"] == json.loads(
        plan.projection_provenance_json
    )
    assert restored.fingerprint == projected.fingerprint
    np.testing.assert_array_equal(restored.projection, projected.projection)
