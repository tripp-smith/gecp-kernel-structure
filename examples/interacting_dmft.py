"""Bounded half-filled Bethe IPT; no optional native solver is needed."""

from kernelgecp import FermionicBasis
from kernelgecp.dmft import IPTConfig, IPTPlan


def main() -> None:
    config = IPTConfig(interaction=2.0)
    green = FermionicBasis(
        10, 4, 1e-8, precision="reference", reference_backend="python"
    )
    sigma = FermionicBasis(
        10, 12, 1e-8, precision="reference", reference_backend="python"
    )
    result = IPTPlan(green, sigma, config=config).solve()
    if not result.converged or result.final_update is None:
        raise RuntimeError(
            f"IPT did not converge: {result.stop_reason}: {result.failure_detail}"
        )
    print(
        f"Iterations: {len(result.history)}; "
        f"unmixed defect: {result.unmixed_defect:.3e}"
    )
    density = -(green.design_tau([config.beta]) @ result.g_coefficients)[0]
    print(f"Per-spin density: {density}; dynamic-Sigma rank: {sigma.rank}")
    print("Convergence is discrete IPT evidence; IPT is an approximate impurity model.")


if __name__ == "__main__":
    main()
