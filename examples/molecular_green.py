"""Read actual bound Green GF2 output and archive its matrix compression."""

from __future__ import annotations

import argparse
from pathlib import Path

from kernelgecp import FermionicBasis
from kernelgecp.integrations.green import read_green_tau


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("results", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--wmax", type=float, default=4.0)
    args = parser.parse_args()
    data = read_green_tau(args.input, args.results, args.manifest)
    basis = FermionicBasis(data.beta, args.wmax, tol=1e-9, grid_order=8)
    fit = basis.fit_tau(data.g_tau, data.tau, axis=0, rcond=1e-13)
    fit.save(args.archive)
    print(f"Native iteration {data.iteration}, axes {data.axes}, rank {basis.rank}")
    print(f"Total GF2 energy: {data.energies.total_hartree:.12g} Hartree")
    print(f"Observed fit residual: {fit.diagnostics.training_residual_max:.3g}")
    print("Static Sigma, overlap and mu remain in the bound native record.")
    print("The training residual does not certify unsampled representation error.")


if __name__ == "__main__":
    main()
