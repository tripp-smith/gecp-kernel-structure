"""DRAFT saved EPW comparator evidence replay; caller authenticates artifact.

No native executable is inspected or run. This is not scientific eligibility.
"""

import argparse
import contextlib
import gzip
import hashlib
import json
import math
import re
import signal
import stat
import subprocess
import threading
import time
import zlib
from decimal import Decimal, DecimalException
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
UPSTREAM_IN = "a67e51f747674319d3c5659e402529bcdf3d571b54f12d2802f2531d74ba8a82"
UPSTREAM_IO = "416ee21c65fd49b2a04894c2d5e18ffb1a5e8f7c223cab950a28f760cf56ec36"
SCOUT = "7a74f36be9a1604ec8b4184e1bb985f14d237835ec43b19c6ac7e9f4dc2937e1"
CHANNELS = (
    "Z_real",
    "Z_imag",
    "gap_real_ev",
    "gap_imag_ev",
    "shift_real_ev",
    "shift_imag_ev",
)
LIMIT = 300
MAX_FILE = 256 * 1024 * 1024
MAX_ROWS = 2_000_000
MAX_INVENTORY = 8 * 1024**3


def require(ok, text):
    if not ok:
        raise ValueError(text)


def safe(root, relative):
    path = root
    for piece in relative.split("/"):
        require(piece not in ("", ".", "..") and "\\" not in piece, "unsafe path")
        path = path / piece
        require(not path.is_symlink(), "symlink")
    require(path.is_file() and path.stat().st_nlink == 1, "missing or hardlinked input")
    require(path.stat().st_size <= MAX_FILE, "input byte limit")
    return path


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path):
    require(path.stat().st_size < 16 * 1024 * 1024, "JSON size")

    def pairs(seq):
        d = {}
        for k, v in seq:
            require(k not in d, "duplicate JSON key")
            d[k] = v
        return d

    def invalid(x):
        raise ValueError("nonfinite JSON: " + x)

    obj = json.loads(path.read_text(), object_pairs_hook=pairs, parse_constant=invalid)
    json.dumps(obj, allow_nan=False)
    require(isinstance(obj, dict), "invalid JSON")
    return obj


def blob_hash(revision, path, repository):
    data = subprocess.check_output(
        ["git", "-C", str(repository), "show", revision + ":" + path], timeout=20
    )
    return hashlib.sha256(data).hexdigest()


def checked_path(path):
    path = Path(path).absolute()
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current = current / part
        require(not current.is_symlink(), "symlink path")
    return path


def inventory(root):
    result = {}
    total = 0
    for path in root.rglob("*"):
        mode = path.lstat().st_mode
        require(not stat.S_ISLNK(mode), "symlink in payload")
        if stat.S_ISDIR(mode):
            continue
        require(stat.S_ISREG(mode) and path.stat().st_nlink == 1, "unsafe payload file")
        total += path.stat().st_size
        require(
            len(result) < 20000
            and path.stat().st_size <= MAX_FILE
            and total <= MAX_INVENTORY,
            "inventory budget",
        )
        result[path.relative_to(root).as_posix()] = sha(path)
    return result


def restart_axes(path, nfreq, nstates):
    require(
        type(nfreq) is int and type(nstates) is int and 0 < nfreq * nstates <= MAX_ROWS,
        "native row budget",
    )
    frequencies, energies = [], []
    with path.open("rb") as stream:
        require(b"shift(w)" in stream.readline(4096), "native header")
        for i in range(nfreq):
            for j in range(nstates):
                line = stream.readline(1024).rstrip(b"\r\n")
                require(len(line) == 100, "native row size")
                fields = [
                    float(line[k : k + 20].replace(b"D", b"E"))
                    for k in range(0, 100, 20)
                ]
                require(
                    all(math.isfinite(x) for x in fields)
                    and fields[0] > 0
                    and fields[2] > 0,
                    "invalid native fields",
                )
                if j == 0:
                    frequencies.append(fields[0])
                else:
                    require(fields[0] == frequencies[-1], "native frequency block")
                if i == 0:
                    energies.append(fields[1])
                else:
                    require(fields[1] == energies[j], "native state ordering")
        require(not stream.readline(1024), "extra native rows")
    require(
        all(a < b for a, b in zip(frequencies, frequencies[1:], strict=False)),
        "native frequency ordering",
    )
    return frequencies, energies


def compare_native_axes(first, second):
    require(first == second, "raw/restored native axes differ")


def expected_input(text):
    require(
        not re.search(r"(?im)^\s*imag_read\s*=", text), "unexpected source imag_read"
    )
    for pattern, replacement in (
        (r"(?m)^(\s*temps\s*=\s*)17\s+20\s+30(?=\s*!|\s*$)", r"\g<1>20"),
        (r"(?m)^(\s*iverbosity\s*=\s*)0(?=\s|$)", r"\g<1>2"),
        (
            r"(?im)^([ \t]*limag[ \t]*=[ \t]*\.true\.[ \t]*)$",
            r"\1\n  imag_read = .true.",
        ),
    ):
        text, count = re.subn(pattern, replacement, text)
        require(count == 1, "source setting")
    require(
        re.search(r"(?m)^\s*wscut\s*=\s*1\.0\s*!", text) is not None, "source wscut"
    )
    return text


def printed_fraction(value):
    number = Decimal(value.decode("ascii").replace("D", "E"))
    require(
        number.is_finite() and (number.is_zero() or -400 <= number.adjusted() <= 308),
        "printed decimal range",
    )
    return Fraction(number)


def exact_field_budget(raw, fitted):
    if raw == fitted:
        return 0, 0
    raw_number, fitted_number = printed_fraction(raw), printed_fraction(fitted)
    difference = abs(fitted_number - raw_number)
    return difference, difference / (
        Fraction(1, 100000) + Fraction(1, 1000) * abs(raw_number)
    )


def compare_recorded_metrics(actual, recorded):
    require(
        isinstance(recorded, dict) and set(actual) == set(recorded), "comparison schema"
    )
    for key in ("state_count", "frequency_count", "row_count"):
        require(
            type(recorded[key]) is int and recorded[key] == actual[key],
            "comparison count",
        )
    for key in ("atol", "rtol", "budget_reference", "passed"):
        require(
            recorded[key] == actual[key] and (key != "passed" or recorded[key] is True),
            "comparison policy",
        )
    require(
        isinstance(recorded["channels"], dict)
        and set(recorded["channels"]) == set(CHANNELS),
        "comparison channels",
    )
    for name in CHANNELS:
        found, expected = recorded["channels"][name], actual["channels"][name]
        require(
            isinstance(found, dict) and set(found) == set(expected),
            "comparison channel schema",
        )
        require(
            type(found["failures"]) is int
            and found["failures"] == expected["failures"]
            and found["unit"] == expected["unit"],
            "comparison channel count/unit",
        )
        for key in ("max_absolute_error", "max_budget_ratio"):
            value = found[key]
            require(
                type(value) in (int, float)
                and math.isfinite(value)
                and math.isclose(value, expected[key], rel_tol=1e-12, abs_tol=1e-15),
                "comparison recorded maximum",
            )


def table_streams(raw, fitted, states, nfreq, step, receipts, deadline):
    require(
        type(nfreq) is int
        and 2 <= nfreq <= MAX_ROWS
        and 0 < len(states) * nfreq <= MAX_ROWS,
        "row budget",
    )
    stats = {
        name: {
            "max_absolute_error": Fraction(0),
            "max_budget_ratio": Fraction(0),
            "failures": 0,
            "unit": "dimensionless" if name.startswith("Z_") else "eV",
        }
        for name in CHANNELS
    }
    h = [hashlib.sha256(), hashlib.sha256()]
    sizes = [0, 0]
    with gzip.open(raw, "rb") as a, gzip.open(fitted, "rb") as b:

        def read(i, f):
            require(time.monotonic() < deadline, "table replay timeout")
            line = f.readline(4096)
            sizes[i] += len(line)
            require(sizes[i] <= MAX_FILE, "gzip expansion limit")
            h[i].update(line)
            return line

        required = (
            b"Re[znorm(w)]",
            b"Im[znorm(w)]",
            b"Re[delta(w)]",
            b"Im[delta(w)]",
            b"Re[shift(w)]",
            b"Im[shift(w)]",
        )
        for i, f in enumerate((a, b)):
            header = read(i, f)
            require(
                header.lstrip().startswith(b"#") and all(x in header for x in required),
                "table header",
            )
        axis = []
        for state, energy in enumerate(states):
            for j in range(nfreq):
                lines = [read(0, a).rstrip(b"\r\n"), read(1, b).rstrip(b"\r\n")]
                require(all(len(x) == 160 for x in lines), "table row shape")
                fields = [
                    [line[k : k + 20] for k in range(0, 160, 20)] for line in lines
                ]
                require(fields[0][:2] == fields[1][:2], "table axis mismatch")
                rows = [[float(x.replace(b"D", b"E")) for x in row] for row in fields]
                require(
                    all(math.isfinite(x) for row in rows for x in row),
                    "nonfinite table",
                )
                require(rows[0][1] == energy, "state axis mismatch")
                if state == 0:
                    require(
                        abs(rows[0][0] - (j + 1) * step) <= 6e-10,
                        "frequency axis mismatch",
                    )
                    axis.append(fields[0][0])
                else:
                    require(fields[0][0] == axis[j], "repeated frequency mismatch")
                for k, name in enumerate(CHANNELS, 2):
                    diff, ratio = exact_field_budget(fields[0][k], fields[1][k])
                    q = stats[name]
                    q["max_absolute_error"] = max(q["max_absolute_error"], diff)
                    q["max_budget_ratio"] = max(q["max_budget_ratio"], ratio)
                    q["failures"] += int(ratio > 1)
        require(not read(0, a) and not read(1, b), "extra table rows")
    for i, r in enumerate(receipts):
        require(
            sizes[i] == r["bytes"] and h[i].hexdigest() == r["sha256"],
            "decompressed table hash/size",
        )
    for item in stats.values():
        item["max_absolute_error"] = float(item["max_absolute_error"])
        item["max_budget_ratio"] = float(item["max_budget_ratio"])
    return {
        "state_count": len(states),
        "frequency_count": nfreq,
        "row_count": len(states) * nfreq,
        "atol": 1e-5,
        "rtol": 1e-3,
        "budget_reference": "absolute raw value",
        "channels": stats,
        "passed": all(v["failures"] == 0 for v in stats.values()),
    }


def audit_native_comparator(
    payload_root, expected_generation_revision, output_dir, *, repository=ROOT
):
    require(
        threading.current_thread() is threading.main_thread(), "main thread required"
    )
    require(signal.getitimer(signal.ITIMER_REAL) == (0.0, 0.0), "active caller timer")
    old_handler = signal.getsignal(signal.SIGALRM)

    def timed_out(signum, frame):
        raise TimeoutError("whole audit deadline exceeded")

    signal.signal(signal.SIGALRM, timed_out)
    signal.setitimer(signal.ITIMER_REAL, LIMIT)
    deadline = time.monotonic() + LIMIT
    try:
        raw_root = checked_path(payload_root)
        raw_output = checked_path(output_dir)
        root = raw_root.resolve(strict=True)
        parent = raw_output.parent.resolve(strict=True)
        require(root.is_dir() and parent.is_dir(), "invalid directory")
        output = raw_output.resolve(strict=False)
        require(
            output != root and root not in output.parents and not output.exists(),
            "unsafe output",
        )
        # Check inputs before creating output; the caller guarantees no writers.
        before = inventory(root)
        output.mkdir()
    except BaseException:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old_handler)
        raise
    result = {
        "schema": "kernelgecp.epw-native-payload-audit.v1",
        "status": "failed",
        "scope": (
            "static caller-authenticated payload; recorded binary hash lineage "
            "only; binary bytes absent; case executable path/log/source are "
            "attestations; no process rerun or release eligibility"
        ),
    }
    audit_log = output / "audit.log"
    try:
        audit_log.write_text("started saved comparator audit\n")
        require(
            re.fullmatch(r"[0-9a-f]{40}", expected_generation_revision) is not None,
            "revision",
        )
        rec = read_json(safe(root, "restart20-native-comparison.json"))
        native = read_json(safe(root, "epw-anisotropic.json"))
        require(
            rec["schema"] == "kernelgecp.epw-restart-comparison.v1"
            and rec["status"] == "passed"
            and rec["errors"] == []
            and rec["source_revision"] == expected_generation_revision
            and rec["temperature_kelvin"] == 20
            and rec["provider"] == {"qe": "7.5", "epw": "6.0"}
            and rec["self_consistency_iteration"] is False,
            "comparator identity",
        )
        require(
            native["source_revision"] == expected_generation_revision
            and native["status"] == "passed"
            and native["schema"] == "kernelgecp.epw-native-anisotropic-fixture.v1"
            and native["provider"] == rec["provider"],
            "native identity",
        )
        source = rec["source"]
        exe = native["native_anisotropic_execution"]
        binaries = [
            x["executable_sha256"]
            for x in native["baseline"]["native_executions"]
            if x.get("input") == "epw1.in"
        ]
        require(
            len(binaries) == 1
            and binaries[0] == exe["executable_sha256"]
            and exe["returncode"] == 0,
            "native binary lineage",
        )
        require(
            source["binary_sha256"] == exe["executable_sha256"]
            and re.fullmatch("[0-9a-f]{64}", source["binary_sha256"])
            and source["io_source_sha256"]
            == native["upstream_io_sha256"]
            == UPSTREAM_IO
            and source["official_input_sha256"]
            == native["upstream_input_sha256"]
            == UPSTREAM_IN
            and source["script_sha256"]
            == blob_hash(
                expected_generation_revision,
                "scripts/providers/epw_restart.py",
                repository,
            )
            and native["script_sha256"]
            == blob_hash(
                expected_generation_revision,
                "scripts/providers/epw_anisotropic_fixture.py",
                repository,
            ),
            "source binding",
        )
        base = "epw-native/test-suite/epw_super/"
        source_in = safe(root, base + "epw6.in")
        source_log = safe(root, base + "epw6.out")
        a2f = safe(root, base + "MgB2.a2f")
        expected_restart_input = expected_input(source_in.read_text())
        require(
            sha(source_in)
            == exe["input_sha256"]
            == native["executed_input_sha256"]
            == SCOUT
            and sha(source_log) == exe["log_sha256"]
            and rec["input"]["original_sha256"] == SCOUT,
            "native source input/log",
        )
        lines = a2f.read_text().splitlines()
        require(len(lines) > 500 and "w[meV]" in lines[0], "a2f mesh")
        vals = [float(line.split()[0]) for line in lines[1:501]]
        require(
            all(math.isfinite(x) and x > 0 for x in vals)
            and all(y > x for x, y in zip(vals, vals[1:], strict=False)),
            "a2f values",
        )
        ratio = 1000 / vals[-1]
        require(abs(ratio - (math.floor(ratio) + 0.5)) > 1e-5, "ambiguous mesh")
        nfreq = 500 * math.floor(ratio + 0.5)
        step = vals[-1] / 500000
        grid = rec["expected_real_grid"]
        require(
            grid["a2f_sha256"] == sha(a2f)
            and grid["frequency_count"] == nfreq
            and grid["step_ev_derived_from_rounded_a2f"] == step,
            "a2f binding",
        )
        require(
            set(rec["cases"]) == {"raw", "reconstructed", "wrong-frequency"}, "case set"
        )
        states = rec["inputs"]["raw"]["state_energy_minus_ef0_ev"]
        require(
            isinstance(states, list)
            and len(states) == rec["inputs"]["raw"]["state_count"]
            and all(type(x) in (int, float) and math.isfinite(x) for x in states),
            "states",
        )
        raw = safe(root, "restart20/case-0/raw-replay.dat")
        restored = safe(root, "restart20/case-0/restored-native.dat")
        archive_case = read_json(safe(root, "restart20/case-0/case.json"))
        corpus = [x for x in native["files"] if x.get("temperature_kelvin") == 20]
        require(
            len(corpus) == 1
            and corpus[0]["filename"] == "MgB2.imag_aniso_020.00"
            and sha(safe(root, base + corpus[0]["filename"]))
            == corpus[0]["sha256"]
            == sha(raw)
            and archive_case["files_sha256"]["raw-replay.dat"] == sha(raw)
            and archive_case["files_sha256"]["restored-native.dat"] == sha(restored),
            "archive source binding",
        )
        parsed_axes = {}
        for label, path in (("raw", raw), ("reconstructed", restored)):
            meta = rec["inputs"][label]
            require(
                meta["sha256"] == sha(path)
                and meta["bytes"] == path.stat().st_size
                and meta["positive_only"] is True
                and meta["state_count"] == len(states)
                and meta["state_energy_minus_ef0_ev"] == states,
                "native input metadata",
            )
            frequency, energy = restart_axes(
                path, meta["frequency_count"], meta["state_count"]
            )
            require(
                frequency == meta["frequency_ev"]
                and energy == meta["state_energy_minus_ef0_ev"],
                "independent native axes",
            )
            parsed_axes[label] = (frequency, energy)
        compare_native_axes(parsed_axes["raw"], parsed_axes["reconstructed"])
        old = "/artifacts/restart20-native-comparison-restart-cases"
        require(rec["case_root"] == old, "case root binding")
        outmeta = []
        for name in ("raw", "reconstructed", "wrong-frequency"):
            c = rec["cases"][name]
            rel = "restart20-native-comparison-restart-cases/" + name + "/"
            oldcase = old + "/" + name
            log = safe(root, rel + "restart.out")
            inp = safe(root, rel + "epw6.in")
            restart = safe(root, rel + "MgB2.imag_aniso_020.00")
            ex = c["execution"]
            require(
                c["workdir"] == oldcase
                and ex["log"] == oldcase + "/restart.out"
                and ex["command"] == ["/opt/qe/EPW/bin/epw.x", "-in", "epw6.in"]
                and ex["log_sha256"] == sha(log)
                and ex["log_bytes"] == log.stat().st_size
                and ex["mpi_ranks"] == 1
                and ex["threads"] == 1
                and ex["timed_out"] is False
                and c["status"] == "passed",
                "execution binding",
            )
            text = log.read_text(errors="strict")
            require(
                inp.read_text() == expected_restart_input
                and c["input_sha256"] == sha(inp) == rec["input"]["modified_sha256"]
                and c["restart_file_sha256"] == sha(restart),
                "case hashes or modified input",
            )
            if name == "wrong-frequency":
                require(
                    type(ex["returncode"]) is int
                    and ex["returncode"] != 0
                    and "eliashberg_read_aniso_iaxis" in text.lower()
                    and "temperature not the same with the input" in text.lower(),
                    "negative execution",
                )
                original = raw.read_bytes().splitlines(keepends=True)
                damaged = restart.read_bytes().splitlines(keepends=True)
                ix = len(states)
                require(len(original) == len(damaged) > ix, "negative rows")
                f = float(original[ix][:20].replace(b"D", b"E"))
                replacement = f"{f * 1.1:20.10E}".encode()
                original[ix] = replacement + original[ix][20:]
                require(
                    original == damaged
                    and c["negative_control"]
                    == {
                        "line_1_based": ix + 1,
                        "state_ordinal_1_based": len(states),
                        "original_frequency_ev": f,
                        "replacement_frequency_ev": float(replacement),
                        "change": (
                            "last state of first frequency block, frequency times 1.1"
                        ),
                    },
                    "negative mutation",
                )
            else:
                require(
                    ex["returncode"] == 0
                    and c["negative_control"] is None
                    and all(
                        x in text
                        for x in (
                            "Pade approximant of full-bandwidth anisotropic",
                            "Convergence was reached for N =",
                            "raxis_pade",
                            "Total program execution",
                        )
                    )
                    and re.search(r"EPW\s+v\.6\.0\b", text)
                    and "error in routine" not in text.lower(),
                    "positive execution",
                )
                require(
                    sha(restart) == sha(raw if name == "raw" else restored),
                    "restart link",
                )
                m = c["output"]
                require(m["path"] == oldcase + "/MgB2.pade_aniso_020.00", "output path")
                outmeta.append(m)
        metrics = table_streams(
            safe(
                root,
                "restart20-native-comparison-restart-cases/raw/MgB2.pade_aniso_020.00.gz",
            ),
            safe(
                root,
                "restart20-native-comparison-restart-cases/reconstructed/MgB2.pade_aniso_020.00.gz",
            ),
            states,
            nfreq,
            step,
            outmeta,
            deadline,
        )
        require(metrics["passed"] is True, "exact printed-number budget failed")
        compare_recorded_metrics(metrics, rec["comparison"])
        result.update(
            status="passed",
            source_revision=expected_generation_revision,
            comparison=metrics,
            cases=["raw", "reconstructed", "wrong-frequency"],
            payload_sha256=before,
            comparison_scope=(
                "exact finite printed-number budgets, not continuous certification"
            ),
        )
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        AttributeError,
        subprocess.SubprocessError,
        EOFError,
        TimeoutError,
        zlib.error,
        DecimalException,
    ) as error:
        # A bounded one-second cleanup alarm follows an expired computation budget.
        if isinstance(error, TimeoutError):
            signal.setitimer(signal.ITIMER_REAL, 1.0)
        result["error"] = f"{type(error).__name__}: {error}"
    try:
        # Write all normal outputs before the final, independent input inventory.
        result["payload_unchanged"] = True
        with audit_log.open("a") as stream:
            stream.write(result["status"] + " " + result.get("error", "") + "\n")
        receipt_path = output / "receipt.json"
        receipt_path.write_text(
            json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        try:
            require(inventory(root) == before, "payload changed during audit")
        except (OSError, ValueError, TimeoutError) as error:
            result.update(
                status="failed",
                payload_unchanged=False,
                error=f"final inventory failed: {type(error).__name__}: {error}",
            )
            with audit_log.open("a") as stream:
                stream.write(result["error"] + "\n")
            receipt_path.write_text(
                json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
            )
            # The failure receipt is also followed by a final input inventory.
            with contextlib.suppress(OSError, ValueError, TimeoutError):
                inventory(root)
            # No concurrent writers are permitted; do not reclassify this as passed.
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old_handler)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload", required=True, type=Path)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--repository", type=Path, default=ROOT)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        result = audit_native_comparator(
            args.payload, args.source_revision, args.output, repository=args.repository
        )
    except (OSError, ValueError, TimeoutError) as error:
        parser.exit(1, f"native payload audit failed: {error}\n")
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
