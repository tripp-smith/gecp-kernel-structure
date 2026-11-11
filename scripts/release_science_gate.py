"""Replay one authenticated science role. Caller authenticates metadata receipts."""

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from scripts import release_science_payloads as payloads
from scripts.providers.verify_receipt import read_receipt
from scripts.release_science_runs import REPOSITORY, V4_JOBS, V4_SCOPE

ROLES = ("covariance", "cthyb", "molecular", "epw", "providers", "spectroscopy")
SCHEMAS = {r: "kernelgecp." + r + "-payload-audit.v1" for r in ROLES[:-1]}
SCHEMAS["providers"] = "kernelgecp.provider-payload-audit.v1"
SCHEMAS["spectroscopy"] = "kernelgecp.spectroscopy-v4-saved-replay.v1"
FIELDS = {
    "covariance": (
        "controller_identity",
        "source_sha256",
        "dataset_sha256",
        "payload_sha256",
        "runtime",
        "independent_sha256",
        "provider_sha256",
    ),
    "cthyb": (
        "controller_identity",
        "source_sha256",
        "runtime",
        "payload_sha256",
        "independent_sha256",
    ),
    "molecular": ("payload_sha256", "independent_sha256", "execution_source"),
    "epw": (
        "auditor_identity",
        "native_auditor_sha256",
        "runtime",
        "archive_sha256",
        "native_sha256",
        "payload_sha256",
    ),
    "providers": ("lanes",),
    "spectroscopy": ("counts", "source_sha256", "receipt_sha256"),
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def inventory(root):
    result = {}
    for p in sorted(Path(root).rglob("*")):
        require(not p.is_symlink(), "symlink in input")
        if p.is_file():
            require(p.stat().st_nlink == 1, "hardlink in input")
            result[str(p.relative_to(root))] = (p.stat().st_size, sha(p))
        else:
            require(p.is_dir(), "nonregular input")
    return result


def origin_check(o):
    require(
        isinstance(o, dict)
        and o.get("scope") == V4_SCOPE
        and o.get("repository") == REPOSITORY
        and o.get("workflow") == "maxent-science.yml",
        "invalid spectroscopy origin",
    )
    require(
        isinstance(o.get("revision"), str)
        and re.fullmatch("[0-9a-f]{40}", o["revision"])
        and payloads.positive(o.get("id"))
        and type(o.get("attempt")) is int
        and o["attempt"] == 1,
        "invalid origin identity",
    )
    jobs = o.get("jobs")
    require(
        isinstance(jobs, list)
        and len(jobs) == len(V4_JOBS)
        and all(isinstance(j, dict) and payloads.positive(j.get("id")) for j in jobs)
        and {j.get("name") for j in jobs} == set(V4_JOBS)
        and len({j["id"] for j in jobs}) == len(jobs),
        "invalid origin jobs",
    )
    a = o.get("artifact")
    require(
        isinstance(a, dict)
        and a.get("name") == f"maxent-fixed-bath-final-{o['id']}-1"
        and payloads.positive(a.get("id"))
        and payloads.positive(a.get("size_in_bytes"))
        and a["size_in_bytes"] <= payloads.MAX_ZIP
        and isinstance(a.get("digest"), str)
        and re.fullmatch("sha256:[0-9a-f]{64}", a["digest"]),
        "invalid origin artifact",
    )
    return o


def validate_wrapper(role, value, revision, origin=None):
    require(
        isinstance(value, dict)
        and value.get("schema") == SCHEMAS[role]
        and value.get("status") == "passed"
        and value.get("source_revision") == revision,
        "wrapper receipt identity",
    )
    require(all(value.get(k) for k in FIELDS[role]), "wrapper evidence missing")
    if role in ("covariance", "cthyb", "epw"):
        require(
            value.get(
                "inputs_unchanged" if role == "covariance" else "payload_unchanged"
            )
            is True,
            "wrapper input changed",
        )
    if role == "providers":
        from scripts.release_provider_audits import LANES

        lanes = value["lanes"]
        require(
            isinstance(lanes, dict)
            and set(lanes) == set(LANES)
            and all(
                isinstance(x, dict)
                and x.get("receipt_valid") is True
                and x.get("receipt_sha256")
                and isinstance(x.get("bound_payload_sha256"), dict)
                for x in lanes.values()
            )
            and all(
                lanes[name].get("independent_dmft")
                for name, provider in LANES.items()
                if provider == "triqs"
            ),
            "provider lanes incomplete",
        )
    if role == "spectroscopy":
        require(
            value.get("counts")
            == {
                "synthetic_calibration": 24,
                "synthetic_final": 72,
                "measured_calibration": 8,
                "measured_final": 24,
                "raw_chains": 256,
            },
            "spectroscopy cohort incomplete",
        )
        hashes = value["receipt_sha256"]
        require(
            isinstance(hashes, dict)
            and set(hashes)
            == {
                "synthetic-cal",
                "synthetic-final",
                "measured-cal",
                "raw",
                "measured-final",
            }
            and all(
                isinstance(h, str) and re.fullmatch("[0-9a-f]{64}", h)
                for h in hashes.values()
            ),
            "spectroscopy independent receipt inventory",
        )
        require(
            value.get("archive_sha256") == origin["artifact"]["digest"][7:]
            and value.get("archive_size") == origin["artifact"]["size_in_bytes"]
            and value.get("run_id") == str(origin["id"]),
            "spectroscopy archive mismatch",
        )


def run_role(
    role,
    selected_file,
    origin_file,
    runtime_python,
    output,
    *,
    downloader=payloads.download,
    runner=subprocess.run,
    root=None,
):
    root = Path(root or Path(__file__).resolve().parents[1]).resolve()
    selected_file, origin_file = (
        Path(selected_file).resolve(),
        Path(origin_file).resolve(),
    )
    output = Path(output).absolute()
    require(role in ROLES, "unknown role")
    require(
        not output.exists() and not output.is_symlink() and output.parent.is_dir(),
        "output must be new",
    )
    require(root not in output.resolve().parents, "output must be outside source")
    output.mkdir()
    try:
        selected, origin = (
            read_receipt(selected_file),
            origin_check(read_receipt(origin_file)),
        )
        payloads.validate_selection(selected)
        require(
            selected.get("profile") == "production-components",
            "wrong selection profile",
        )
        p = selected["revision"]
        require(
            subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=root, text=True
            ).strip()
            == p,
            "source HEAD mismatch",
        )
        require(
            not subprocess.check_output(
                ["git", "status", "--porcelain", "--untracked-files=all"],
                cwd=root,
                text=True,
            ).strip(),
            "source checkout is dirty",
        )
        selection_hash, origin_hash = sha(selected_file), sha(origin_file)
        runtime = str(Path(runtime_python).absolute())
        require(Path(runtime).is_file(), "runtime python missing")
        source = origin["revision"] if role == "spectroscopy" else p
        module = {
            "covariance": "release_covariance_audits",
            "cthyb": "release_cthyb_audits",
            "molecular": "release_molecular_audits",
            "epw": "release_epw_audits",
            "providers": "release_provider_audits",
            "spectroscopy": "release_spectroscopy_v4",
        }[role]
        args = [sys.executable, "-m", "scripts." + module]
        records = []
        if role == "spectroscopy":
            a = origin["artifact"]
            archive = output / (a["name"] + ".zip")
            downloader(a["id"], archive, output / "download.log")
            require(
                archive.stat().st_size == a["size_in_bytes"]
                and sha(archive) == a["digest"][7:],
                "spectroscopy archive size or digest mismatch",
            )
            args += [
                "--archive",
                str(archive),
                "--archive-sha256",
                a["digest"][7:],
                "--archive-size",
                str(a["size_in_bytes"]),
                "--run-id",
                str(origin["id"]),
            ]
            before = (archive.stat().st_size, sha(archive))
            records = [a]
        else:
            mat = payloads.materialize(
                selected,
                output / "payloads",
                output / "materialized.json",
                downloader,
                role=role,
            )
            records = mat["artifacts"]
            byname = {x["name"]: Path(x["directory"]) for x in records}
            run = selected["runs"][role]
            if role == "providers":
                args += ["audit", "--artifacts-root", str(output / "payloads" / role)]
            else:
                name = {
                    "covariance": "real-data-covariance",
                    "cthyb": f"cthyb-evidence-replay-{run['id']}",
                    "molecular": "molecular-green-native-linux",
                    "epw": "epw-explicit-restart-archive-linux",
                }[role]
                args += ["--payload", str(byname[name])]
                if role == "covariance":
                    args += [
                        "--data",
                        str(
                            byname[
                                f"covariance-uci-dataset-{run['id']}-{run['attempt']}"
                            ]
                            / "uci-501.zip"
                        ),
                    ]
            before = inventory(output / "payloads")
        runtime_flag = {
            "covariance": "--gp-python",
            "molecular": "--molecular-python",
            "epw": "--epw-python",
        }.get(role, "--triqs-python")
        wrapper_out = output / ("provider.json" if role == "providers" else "replay")
        args += [
            "--source-revision",
            source,
            runtime_flag,
            runtime,
            "--output",
            str(wrapper_out),
        ]
        env = os.environ.copy()
        env.pop("GH_TOKEN", None)
        env.pop("GITHUB_TOKEN", None)
        with (output / "wrapper.log").open("xb") as log:
            proc = runner(
                args,
                cwd=root,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=315 * 60,
                check=False,
                env=env,
            )
        require(proc.returncode == 0, "wrapper command failed")
        actual = wrapper_out if role == "providers" else wrapper_out / "receipt.json"
        validate_wrapper(role, read_receipt(actual), source, origin)
        require(
            (archive.stat().st_size, sha(archive)) == before
            if role == "spectroscopy"
            else inventory(output / "payloads") == before,
            "inputs changed during replay",
        )
        require(
            sha(selected_file) == selection_hash and sha(origin_file) == origin_hash,
            "metadata changed during replay",
        )
        shutil.copyfile(actual, output / "wrapper-receipt.json")
        result = dict(
            schema="kernelgecp.science-role-replay.v1",
            status="passed",
            role=role,
            producer_revision=p,
            source_revision=source,
            selection_sha256=selection_hash,
            origin_sha256=origin_hash,
            wrapper_receipt_sha256=sha(output / "wrapper-receipt.json"),
            command_exit=proc.returncode,
            artifacts=records,
        )
        (output / "role-receipt.json").write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n"
        )
        return result
    except BaseException as exc:
        (output / "failure.log").write_text(f"{type(exc).__name__}: {exc}\n")
        raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="action", required=True)
    q = sub.add_parser("run-role")
    q.add_argument("--role", choices=ROLES, required=True)
    for arg in ("selected-runs", "spectroscopy-origin", "runtime-python", "output"):
        q.add_argument("--" + arg, required=True, type=Path)
    a = p.parse_args()
    try:
        run_role(
            a.role, a.selected_runs, a.spectroscopy_origin, a.runtime_python, a.output
        )
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        p.exit(1, f"role replay failed: {exc}\n")


if __name__ == "__main__":
    main()
