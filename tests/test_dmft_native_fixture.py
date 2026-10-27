"""Small attributed native TRIQS data exercised in ordinary package CI."""

import hashlib
import json
from pathlib import Path

import numpy as np

from kernelgecp import FermionicBasis


def test_actual_native_interacting_green_and_dynamic_sigma_compress():
    source = Path(__file__).with_name("data") / "dmft_native_ipt_v1.npz"
    manifest = json.loads(source.with_suffix(".json").read_text())
    assert hashlib.sha256(source.read_bytes()).hexdigest() == manifest["archive_sha256"]
    assert manifest["provider_version"] == "4.0.2"
    assert (
        manifest["chemical_potential_convention"]
        == "mu_U_over_2_cancels_static_Hartree"
    )
    with np.load(source, allow_pickle=False) as native:
        for field, cutoff in (("g", 4), ("sigma", 12)):
            basis = FermionicBasis(
                10, cutoff, 1e-8, precision="reference", reference_backend="python"
            )
            fit = basis.fit_tau(native[field + "_tau"], native["tau"], rcond=1e-13)
            time_error = np.max(
                np.abs(
                    fit.evaluate_tau(native["validation_tau"])
                    - native["validation_" + field + "_tau"]
                )
            )
            frequency_error = np.max(
                np.abs(
                    fit.evaluate_matsubara(native["validation_indices"])
                    - native["validation_" + field + "_iw"]
                )
            )
            assert time_error <= 1e-6
            assert frequency_error <= 1e-6
