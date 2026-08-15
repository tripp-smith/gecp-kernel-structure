"""Run and plot the realistic synthetic application suite."""

from __future__ import annotations

import argparse
import subprocess
import tomllib
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from kernelgecp import (
    SyntheticApplicationConfig,
    SyntheticApplicationSuite,
    quadrature_weights,
    run_synthetic_application_suite,
    synthetic_sensor_points,
    synthetic_spectral_density,
)


def _load_config(path: Path) -> SyntheticApplicationConfig:
    with path.open("rb") as stream:
        raw: dict[str, Any] = tomllib.load(stream)
    if raw.pop("schema_version", None) != 1:
        raise ValueError("unsupported synthetic application schema_version")
    return SyntheticApplicationConfig(**raw)


def _git_revision(root: Path) -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _plot_suite(
    suite: SyntheticApplicationSuite,
    config: SyntheticApplicationConfig,
    output: Path,
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)

    omega = np.linspace(-config.cutoff, config.cutoff, 2_001)
    integration_weights = quadrature_weights(omega)
    for fixture, label in (
        ("hubbard_three_peak", "Hubbard-like three peak"),
        ("gapped_two_band", "Gapped two band"),
    ):
        density = synthetic_spectral_density(fixture, omega)  # type: ignore[arg-type]
        density /= float(integration_weights @ density)
        axes[0, 0].plot(omega, density, label=label)
    axes[0, 0].set_title("Synthetic fermionic spectra")
    axes[0, 0].set_xlabel("scaled frequency ω")
    axes[0, 0].set_ylabel("normalized density ρ(ω)")
    axes[0, 0].legend(frameon=False)

    names = [record.fixture.replace("_", "\n") for record in suite.green_compression]
    positions = np.arange(len(names), dtype=np.float64)
    width = 0.36
    axes[0, 1].bar(
        positions - width / 2,
        [record.green_validation_error for record in suite.green_compression],
        width,
        label="Green error",
    )
    axes[0, 1].bar(
        positions + width / 2,
        [record.transfer_bound for record in suite.green_compression],
        width,
        label="L¹ transfer bound",
    )
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_xticks(positions, names)
    axes[0, 1].set_title("Universal GECP compression")
    axes[0, 1].set_ylabel("held-out maximum error")
    axes[0, 1].legend(frameon=False)

    for row, record in enumerate(suite.sparse_recovery):
        axes[1, 0].scatter(
            record.true_frequencies,
            np.full(record.true_atom_count, row + 0.12),
            s=300 * np.asarray(record.true_weights),
            marker="o",
            facecolors="none",
            edgecolors="black",
            label="true" if row == 0 else None,
        )
        recovered_sizes = 300 * np.minimum(1.0, np.abs(record.recovered_weights))
        axes[1, 0].scatter(
            record.recovered_frequencies,
            np.full(record.recovered_atom_count, row - 0.12),
            s=recovered_sizes,
            marker="x",
            label="recovered" if row == 0 else None,
        )
        axes[1, 0].text(
            config.cutoff,
            row,
            f"held-out {record.held_out_error:.1e}",
            ha="right",
            va="center",
            fontsize=8,
        )
    axes[1, 0].set_xlim(-config.cutoff, config.cutoff)
    axes[1, 0].set_yticks(
        range(len(suite.sparse_recovery)),
        [record.fixture.replace("_", "\n") for record in suite.sparse_recovery],
    )
    axes[1, 0].set_xlabel("scaled frequency ω")
    axes[1, 0].set_title("Sparse spectrum recovery")
    axes[1, 0].legend(frameon=False, loc="lower left")

    points = synthetic_sensor_points(config.psd_point_count, config.seed + 1)
    selected = np.asarray(suite.psd_landmarks.selected_indices, dtype=int)
    axes[1, 1].scatter(points[:, 0], points[:, 1], s=18, alpha=0.35, label="sites")
    axes[1, 1].scatter(
        points[selected, 0],
        points[selected, 1],
        c=np.arange(selected.size),
        cmap="viridis",
        s=45,
        label="GECP/Cholesky landmarks",
    )
    axes[1, 1].set_title(
        f"PSD covariance landmarks: {selected.size}/{points.shape[0]} sites"
    )
    axes[1, 1].set_xlabel("x")
    axes[1, 1].set_ylabel("y")
    axes[1, 1].set_aspect("equal")
    axes[1, 1].legend(frameon=False, fontsize=8)

    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("experiments/synthetic_applications.toml"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("experiments/data/synthetic_applications.json"),
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=Path("experiments/data/synthetic_applications.png"),
    )
    arguments = parser.parse_args()
    root = arguments.config.resolve().parents[1]
    config = _load_config(arguments.config)
    suite = run_synthetic_application_suite(config, git_commit=_git_revision(root))
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = arguments.output.with_suffix(arguments.output.suffix + ".tmp")
    temporary.write_text(suite.canonical_json() + "\n", encoding="utf-8")
    temporary.replace(arguments.output)
    _plot_suite(suite, config, arguments.plot)


if __name__ == "__main__":
    main()
