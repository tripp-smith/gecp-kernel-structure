"""Independently audit one authenticated molecular payload; not release eligibility.

The caller must authenticate and safely materialize the artifact before calling.
The supplied Python must be the pinned molecular science environment.
"""

import argparse
import json
import os
import re
import stat
import subprocess
import sys
from pathlib import Path

from scripts.providers.verify_receipt import read_receipt
from scripts.release_provider_audits import digest, require, safe, source_snapshot

ROOT = Path(__file__).resolve().parents[1]
CASES = {f"h2-{bond}-{grid}" for bond in ("0.74", "1.5") for grid in ("1e4", "1e5")}
GATES = {
    "representation",
    "density",
    "matrix_causality",
    "coefficient_prediction_consistency",
    "native_Dyson_and_independent_Fourier",
    "native_grid_refinement",
    "doubled_cutoff_controls",
    "finite_model_and_DF_accounting",
}
THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)


def inventory(root):
    require(
        root.is_dir() and not any(p.is_symlink() for p in (root, *root.parents)),
        "unsafe payload root",
    )
    result = {}
    for path in root.rglob("*"):
        mode = path.lstat().st_mode
        require(not stat.S_ISLNK(mode), "symlink in payload")
        if stat.S_ISDIR(mode):
            continue
        require(
            stat.S_ISREG(mode) and path.stat().st_nlink == 1,
            "nonregular or hardlinked payload",
        )
        result[path.relative_to(root).as_posix()] = digest(path)
    return result


def checked_path(path):
    """Reject symlinks in the original traversal, before normalizing dot-dot."""
    path = Path(path).absolute()
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current = current / part
        require(not current.is_symlink(), "symlink in path")
    return path


def log_text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value or ""


def audit(
    payload,
    revision,
    python,
    output,
    *,
    repository=ROOT,
    runner=subprocess.run,
    snapshot=source_snapshot,
):
    """Write all results under a new exclusive attempt directory, including failure."""
    raw_root = checked_path(payload)
    raw_output = checked_path(output)
    root = raw_root.resolve(strict=True)
    require(root.is_dir(), "payload must be a directory")
    parent = raw_output.parent.resolve(strict=True)
    require(parent.is_dir(), "output parent must be a directory")
    output = raw_output.resolve(strict=False)
    require(
        output != root and root not in output.parents, "output cannot be inside payload"
    )
    require(not output.exists() and not output.is_symlink(), "output must be new")
    output.mkdir(parents=False, exist_ok=False)
    receipt = {
        "schema": "kernelgecp.molecular-payload-audit.v1",
        "status": "failed",
        "source_revision": revision,
        "scope": "molecular saved-payload audit only; not release eligibility",
    }
    try:
        require(
            isinstance(revision, str) and re.fullmatch(r"[0-9a-f]{40}", revision),
            "invalid revision",
        )
        before = inventory(root)
        report = read_receipt(safe(root, "comparison/report.json"))
        native = read_receipt(safe(root, "native/native-campaign.json"))
        require(isinstance(report, dict) and isinstance(native, dict), "invalid report")
        require(
            report.get("schema") == "kernelgecp.molecular-green.v2"
            and report.get("pilot") is False
            and report.get("source_revision") == revision
            and report.get("native_source_revision") == revision,
            "comparison source or schema mismatch",
        )
        require(
            native.get("schema") == "kernelgecp.green-molecular-native.v1"
            and native.get("status") == "passed"
            and native.get("source_revision") == revision
            and report.get("native_campaign_sha256")
            == before["native/native-campaign.json"],
            "native source or campaign mismatch",
        )
        rows = native["cases"]
        require(
            isinstance(rows, list)
            and len(rows) == 4
            and {r["directory"] for r in rows} == CASES
            and {(r["bond_length_angstrom"], r["native_grid_lambda"]) for r in rows}
            == {(b, g) for b in (0.74, 1.5) for g in (1e4, 1e5)},
            "incorrect native cohort",
        )
        require(
            {p.name for p in (root / "native").iterdir() if p.is_dir()} == CASES,
            "orphan native case",
        )
        for row in rows:
            case = row["directory"]
            label = f"{int(row['native_grid_lambda']):.0e}".replace("e+0", "e")
            require(
                case == f"h2-{row['bond_length_angstrom']}-{label}",
                "native case label mismatch",
            )
            base = f"native/{case}/"
            manifest = read_receipt(safe(root, base + "bound-manifest.json"))
            provider = read_receipt(safe(root, base + "provider-receipt.json"))
            require(
                provider.get("schema") == "kernelgecp.provider-smoke.v1"
                and provider.get("provider") == "green"
                and provider.get("status") == "passed"
                and provider.get("source_revision") == revision,
                "native provider identity mismatch",
            )
            require(
                before[base + "bound-manifest.json"] == row["manifest_sha256"]
                and manifest["source_revision"] == revision
                and before[base + "provider-receipt.json"]
                == row["provider_receipt_sha256"],
                "native binding mismatch",
            )
            for record in manifest["files"].values():
                require(
                    digest(safe(root, base + record["path"])) == record["sha256"],
                    "native file hash mismatch",
                )
        require(
            before["comparison/selection.json"] == report["selection_sha256"],
            "selection hash mismatch",
        )
        expected = {"report.json", "selection.json"}
        require(
            set(report["reference_artifacts"]) == CASES, "reference cohort mismatch"
        )
        for case, files in report["reference_artifacts"].items():
            require(
                set(files) == {"calibration.npz", "validation.npz"},
                "reference inventory mismatch",
            )
            for name, value in files.items():
                path = f"{case}/{name}"
                require(
                    digest(safe(root, "comparison/" + path)) == value,
                    "reference hash mismatch",
                )
                expected.add(path)
        for kind in ("calibration", "cutoff_refinement_attempts", "validation"):
            for row in report[kind]:
                directory = (
                    row["worker_directory"]
                    if kind != "validation"
                    else (
                        f"validation/repeat-{row['repeat']}-"
                        f"{row['method']}-{row['case']}"
                    )
                )
                require(
                    isinstance(row["artifacts"], dict) and row["artifacts"],
                    "missing worker artifacts",
                )
                for name, value in row["artifacts"].items():
                    path = f"{directory}/{name}"
                    require(
                        path not in expected
                        and digest(safe(root, "comparison/" + path)) == value,
                        "duplicate or altered worker artifact",
                    )
                    expected.add(path)
        actual = {
            p[len("comparison/") :] for p in before if p.startswith("comparison/")
        }
        require(actual == expected, "orphan or missing comparison artifact")
        env = dict(os.environ)
        env.update(dict.fromkeys(THREADS, "1"))
        env.pop("PYTHONPATH", None)
        with snapshot(repository, revision) as source:
            command = [
                str(Path(python).absolute()),
                "-I",
                "-c",
                (
                    "import sys,runpy; sys.path.insert(0,sys.argv.pop(1)); "
                    "runpy.run_module('experiments.audit_molecular_green',"
                    "run_name='__main__',alter_sys=True)"
                ),
                str(source),
                str(root / "comparison"),
                "--native-root",
                str(root / "native"),
                "--output",
                str(output / "independent.json"),
            ]
            result = runner(
                command,
                cwd=source,
                env=env,
                text=True,
                capture_output=True,
                check=False,
                timeout=300,
            )
        (output / "auditor.log").write_text(
            (result.stdout or "") + "\n" + (result.stderr or "")
        )
        require(result.returncode == 0, "independent molecular auditor failed")
        checked = read_receipt(output / "independent.json")
        require(
            isinstance(checked, dict)
            and checked.get("schema") == "kernelgecp.molecular-green-audit.v1"
            and checked.get("status") == "passed"
            and checked.get("generation_revision") == revision
            and checked.get("native_source_revision") == revision
            and type(checked.get("validated_identity_count")) is int
            and checked["validated_identity_count"] == 36
            and isinstance(checked.get("gates"), dict)
            and set(checked["gates"]) == GATES
            and all(v is True for v in checked["gates"].values()),
            "invalid independent audit result",
        )
        for name, count in (
            ("final_trials", 36),
            ("native", 4),
            ("native_grid_refinement", 2),
            ("cutoff_controls", 12),
            ("first_passing_calibration", 12),
            ("model_differences", 4),
        ):
            require(
                isinstance(checked.get(name), list) and len(checked[name]) == count,
                "independent audit cohort mismatch",
            )
        require(inventory(root) == before, "payload changed during audit")
        receipt.update(
            status="passed",
            payload_sha256=before,
            independent_sha256=digest(output / "independent.json"),
            execution_source="private detached Git checkout",
        )
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        AttributeError,
        subprocess.SubprocessError,
    ) as error:
        receipt["error"] = f"{type(error).__name__}: {error}"
        if not (output / "auditor.log").exists():
            if isinstance(error, subprocess.TimeoutExpired):
                details = log_text(error.stdout) + "\n" + log_text(error.stderr)
                (output / "auditor.log").write_text(
                    details + "\n" + receipt["error"] + "\n"
                )
            else:
                (output / "auditor.log").write_text(receipt["error"] + "\n")
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--molecular-python", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit(
            args.payload, args.source_revision, args.molecular_python, args.output
        )
    except (OSError, ValueError) as error:
        print(f"cannot create audit attempt: {error}", file=sys.stderr)
        return 1
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
