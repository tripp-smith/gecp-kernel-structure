"""Continue a covariance-aware scalar dataset with actual optional TauMaxEnt."""

from __future__ import annotations

import argparse
import dataclasses
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path, help="bound spectral-dataset archive")
    parser.add_argument("--mu0", type=float, default=1e-8)
    parser.add_argument(
        "--default-model", choices=["flat", "broad_gaussian"], default="flat"
    )
    args = parser.parse_args()
    from kernelgecp.integrations.maxent import (
        MaxEntConfig,
        continue_tau_maxent,
        load_spectral_dataset,
    )

    data = load_spectral_dataset(args.dataset)
    config = MaxEntConfig(
        n_alpha=33, minimizer_mu0=args.mu0, default_model=args.default_model
    )
    result = continue_tau_maxent(data, config)
    print(
        json.dumps(
            {
                "dataset_fingerprint": data.fingerprint,
                "effective_rank": result.effective_rank,
                "config": dataclasses.asdict(config),
                "mass": result.mass,
                "selected_alpha": float(result.alphas[result.selected_index]),
                "normalized_training_whitened_rms": (
                    result.normalized_training_whitened_rms
                ),
                "scope": (
                    "new analysis with declared settings; "
                    "no held-out or unique-peak claim"
                ),
            },
            indent=2,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
