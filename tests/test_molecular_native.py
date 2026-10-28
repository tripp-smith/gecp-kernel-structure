"""Mandatory native-workflow physics checks; ordinary CI retains explicit skips."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import numpy as np
import pytest

from kernelgecp.integrations.green import read_green_tau


@pytest.fixture(scope="module")
def native_root():
    declared = os.environ.get("GECP_MOLECULAR_NATIVE_ROOT")
    if declared is None:
        pytest.skip(
            "actual four-run native artifacts are required by the scientific workflow"
        )
    root = Path(declared)
    assert (root / "native-campaign.json").is_file(), (
        "declared native artifact root is incomplete"
    )
    return root


@pytest.mark.parametrize("bond", [0.74, 1.5])
@pytest.mark.parametrize("label", ["1e4", "1e5"])
def test_actual_native_basis_tables_weighted_fields_dyson_and_independent_quadrature(
    native_root, bond, label
):
    from experiments.molecular_native import NativeIRSolution

    folder = native_root / f"h2-{bond}-{label}"
    data = read_green_tau(
        folder / "input.h5", folder / "sim.h5", folder / "bound-manifest.json"
    )
    solution = NativeIRSolution(data, folder / "native-grid.h5", float(label))
    assert max(solution.table_errors.values()) < 2e-10
    assert max(solution.independent_field_reconstruction_errors.values()) < 2e-9
    # Includes the Lambda1e5 argument-roundtrip regression. Columnwise error
    # alone is not a field-error bound and cannot reject a correct weighted G.
    assert solution.physical_tau_roundtrip_matrix_error >= 0
    indices = np.array([-32, -1, 0, 31], dtype=np.int64)
    assert solution.dyson_residual(indices) < 1e-6
    x, weights = np.polynomial.legendre.leggauss(512)
    times = data.beta * (x + 1) / 2
    frequency = (2 * indices + 1) * np.pi / data.beta
    kernel = np.exp(1j * frequency[:, None] * times) * (weights * data.beta / 2)
    for field in ("green", "sigma"):
        # Direct exp(+iw*tau) quadrature does not call the library Fourier
        # basis or native frequency transform; it checks their sign/units.
        transformed = np.tensordot(
            kernel, solution.evaluate_tau(times, field=field), axes=(1, 0)
        )
        expected = solution.evaluate_matsubara(indices, field=field)
        assert np.max(abs(transformed - expected)) < 2e-8


def test_wrong_native_grid_metadata_or_nonfinite_table_is_rejected(
    native_root, tmp_path
):
    import h5py

    from experiments.molecular_native import NativeIRSolution

    folder = native_root / "h2-0.74-1e4"
    data = read_green_tau(
        folder / "input.h5", folder / "sim.h5", folder / "bound-manifest.json"
    )
    grid = tmp_path / "grid.h5"
    shutil.copyfile(folder / "native-grid.h5", grid)
    with h5py.File(grid, "r+") as source:
        source.attrs["__grids_version__"] = "wrong"
    with pytest.raises(ValueError, match="versions differ"):
        NativeIRSolution(data, grid, 1e4)
    shutil.copyfile(folder / "native-grid.h5", grid)
    with h5py.File(grid, "r+") as source:
        source["fermi/uwl"][0, 0] = np.nan
    with pytest.raises(ValueError, match="must be finite"):
        NativeIRSolution(data, grid, 1e4)
