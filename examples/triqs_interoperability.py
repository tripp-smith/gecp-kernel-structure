"""Actual TRIQS 4 complex-matrix fitting and ordered block export.

Run from this checkout inside the isolated environments/triqs.yml runtime::

    python -m examples.triqs_interoperability

This validates an analytic noninteracting fixture, not an interacting DMFT run.
"""

from __future__ import annotations

import json

import numpy as np

from kernelgecp import FermionicBasis
from kernelgecp.integrations.triqs import (
    fit_triqs_tau,
    to_triqs_iw,
    to_triqs_tau,
)


def run_example() -> dict[str, object]:
    from triqs.gfs import BlockGf, Gf
    from triqs.mesh import MeshImFreq, MeshImTime

    beta = 4.0
    basis = FermionicBasis(beta, 1.0, 1e-10, grid_order=8, max_rank=64)
    mesh = MeshImTime(beta=beta, statistic="Fermion", n_tau=81)
    tau = np.asarray([point.value for point in mesh])
    energies = np.array([-0.4, 0.25])
    residues = np.array([[[0.6, 0.1j], [-0.1j, 0.3]], [[0.4, -0.1j], [0.1j, 0.7]]])
    matrix = Gf(mesh=mesh, target_shape=[2, 2])
    matrix.data[:] = -np.einsum(
        "te,eij->tij",
        np.exp(-tau[:, None] * energies) / (1 + np.exp(-beta * energies)),
        residues,
    )
    single = Gf(mesh=mesh, target_shape=[1, 1])
    single.data[:, 0, 0] = -np.exp(-tau * 0.2) / (1 + np.exp(-beta * 0.2))
    source = BlockGf(
        name_list=["zeta", "alpha"], block_list=[matrix, single], make_copies=True
    )
    names = [name for name, _ in source]
    fits = [fit_triqs_tau(basis, block) for _, block in source]
    iw_mesh = MeshImFreq(beta=beta, statistic="Fermion", n_iw=32)
    converted = BlockGf(
        name_list=names,
        block_list=[to_triqs_iw(fit, iw_mesh) for fit in fits],
        make_copies=True,
    )
    output_mesh = MeshImTime(beta=beta, statistic="Fermion", n_tau=127)
    time_output = to_triqs_tau(fits[0], output_mesh)
    query_tau = np.asarray([point.value for point in output_mesh])
    indices = np.asarray([point.index for point in iw_mesh])
    truth_tau = -np.einsum(
        "te,eij->tij",
        np.exp(-query_tau[:, None] * energies) / (1 + np.exp(-beta * energies)),
        residues,
    )
    truth_iw = np.einsum(
        "ne,eij->nij",
        1 / (1j * (2 * indices[:, None] + 1) * np.pi / beta - energies),
        residues,
    )
    time_error = float(np.max(np.abs(time_output.data - truth_tau)))
    frequency_error = float(np.max(np.abs(converted["zeta"].data - truth_iw)))
    if (
        [name for name, _ in converted] != names
        or time_error > 1e-8
        or frequency_error > 1e-8
    ):
        raise RuntimeError("actual TRIQS block export or physical accuracy failed")
    return {
        "provider": "TRIQS 4.0.2",
        "rank": basis.rank,
        "ordered_blocks": names,
        "time_max_error": time_error,
        "frequency_max_error": frequency_error,
        "accuracy_passed": True,
        "scope": "analytic_complex_matrix_fixture_not_interacting_DMFT",
        "certification_status": "not_certified_sampled_checks",
    }


if __name__ == "__main__":
    print(json.dumps(run_example(), indent=2, sort_keys=True))
