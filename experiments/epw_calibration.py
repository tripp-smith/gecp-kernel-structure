"""Reproduce a 17 K EPW calibration miss; final rows remain reserved.

This scouting run is not an acceptance or a native-restart test. It permits a
large condition cap to document failure modes; production defaults are stricter.
"""

from __future__ import annotations

import argparse
import gzip
import json
import tempfile
from pathlib import Path

import numpy as np

from kernelgecp import FermionicBasis
from kernelgecp.integrations.epw import (
    fit_epw_nambu,
    nambu_green,
    read_epw_anisotropic,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tol", type=float, default=1e-3)
    parser.add_argument("--wmax", type=float, default=13.0)
    parser.add_argument("--manifest", default="manifest.json")
    args = parser.parse_args()
    fixture = Path(__file__).resolve().parents[1] / "tests/data/epw"
    records = json.loads((fixture / args.manifest).read_text())
    entry = next(
        e for e in records["files"] if e["context"]["temperature_kelvin"] == 17
    )
    with tempfile.TemporaryDirectory() as directory:
        source = Path(directory) / entry["context"]["native_filename"]
        source.write_bytes(
            gzip.decompress((fixture / entry["compressed_file"]).read_bytes())
        )
        data = read_epw_anisotropic(source, entry["context"])
    calibration = np.array([36, 40, 44])
    final = {45, 49, 53}
    rows = np.array(
        [i for i in range(len(data.indices)) if i not in set(calibration) | final]
    )
    basis = FermionicBasis(
        beta=data.beta, wmax=args.wmax, tol=args.tol, precision="auto"
    )
    truth = nambu_green(data)
    print(
        json.dumps(
            {
                "basis_rank": basis.rank,
                "calibration_rows": calibration.tolist(),
                "reserved_final_rows": sorted(final),
                "fit_rows": rows.tolist(),
                "basis_tolerance": args.tol,
                "spectral_cutoff_ev": args.wmax,
            }
        )
    )
    for rcond in (1e-10, 1e-12, 1e-14):
        fit = fit_epw_nambu(
            data,
            basis,
            fit_rows=rows,
            rcond=rcond,
            complete_adjoint=True,
            max_retained_condition=1e15,
        )
        values = fit.evaluate_matsubara(data.indices)
        recovered = fit.reconstruct_native(real_tolerance=0.1)
        ls = fit.fit.diagnostics.least_squares
        metrics = {
            "rcond": rcond,
            "effective_rank": ls.effective_rank,
            "retained_condition": ls.retained_condition_number,
            "scaled_nambu_error": float(
                np.max(np.abs(values[calibration] - truth[calibration]))
                / np.max(np.abs(truth[calibration]))
            ),
            "z_error": float(
                np.max(np.abs(recovered.data.z[calibration] - data.z[calibration]))
            ),
            "gap_error_ev": float(
                np.max(
                    np.abs(
                        recovered.data.delta_ev[calibration]
                        - data.delta_ev[calibration]
                    )
                )
            ),
            "chi_error_ev": float(
                np.max(
                    np.abs(
                        recovered.data.chi_ev[calibration] - data.chi_ev[calibration]
                    )
                )
            ),
        }
        print(json.dumps(metrics, sort_keys=True))


if __name__ == "__main__":
    main()
