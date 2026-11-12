"""Focused tests; none establish installed-provider or campaign acceptance."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

PATH = Path(__file__).parents[1] / "scripts/providers/projection_precision_smoke.py"
spec = importlib.util.spec_from_file_location("projection_precision_smoke", PATH)
smoke = importlib.util.module_from_spec(spec)
spec.loader.exec_module(smoke)


def test_installed_path_rejects_checkout_and_symlink(tmp_path):
    root = tmp_path / "checkout"
    root.mkdir()
    inside = root / "maxent.py"
    inside.write_text("")
    outside = tmp_path / "installed.py"
    outside.write_text("")
    alias = tmp_path / "alias.py"
    alias.symlink_to(inside)
    with pytest.raises(RuntimeError, match="checkout"):
        smoke.validate_import_path(inside, root, True)
    with pytest.raises(RuntimeError, match="checkout"):
        smoke.validate_import_path(alias, root, True)
    assert smoke.validate_import_path(outside, root, True) == outside


def test_independent_oracle_ill_conditioned_known_span():
    q, _ = np.linalg.qr(np.random.default_rng(4).standard_normal((12, 3)))
    a = q @ np.diag([1.0, 1e-5, 1e-10])
    p, rank = smoke.independent_projector(a)
    assert rank == 3
    assert np.max(np.abs(p - q @ q.T)) < 1e-13
    assert np.allclose(p, p.T, atol=1e-15)


def test_canonical_hash_layout():
    a = np.arange(12, dtype=">f8").reshape(3, 4)
    assert smoke.canonical_sha(a) == smoke.canonical_sha(
        np.asfortranarray(a.astype("<f8"))
    )
    changed = a.copy()
    changed[0, 0] = 1
    assert smoke.canonical_sha(a) != smoke.canonical_sha(changed)


def test_matrix_checks_rejects_provider_mismatch():
    from kernelgecp.integrations.maxent import (
        SpectralDataset,
        prepare_spectral_projection,
    )

    class Basis:
        beta = 10

        def design_tau(self, tau):
            return np.column_stack((np.ones(len(tau)), tau, tau**2))

    tau = np.linspace(0, 10, 9)
    plan = prepare_spectral_projection(
        Basis(), tau, tau, rcond=1e-12, precision_bits=256
    )
    result = smoke.matrix_checks(plan, plan.source_design.copy(), SpectralDataset, tau)
    assert result["bitwise_equal"] and result["projector_max_abs"] < 1e-13
    changed = plan.source_design.copy()
    changed[0, 0] += 1e-10
    with pytest.raises(RuntimeError, match="repeated provider design"):
        smoke.matrix_checks(plan, changed, SpectralDataset, tau)
