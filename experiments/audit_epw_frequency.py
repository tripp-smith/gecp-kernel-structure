"""Independent audit of EPW frequency-cross calibration and frozen campaigns.

Calibration never reads spectral values from reserved final rows. This is finite
native-grid interoperability evidence, not an experimental or speedup claim.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CAL = [36, 40, 44]
FINAL = [45, 49, 53]
RANKS = [48, 60, 72]
BUDGET = {"green_scaled": 1e-6, "z": 1e-5, "gap_ev": 1e-6, "chi_ev": 1e-6}
CONTROLS = {
    "pole_cutoff_ev": 13,
    "dictionary_size": 1201,
    "precision_bits": 256,
    "rcond": 1e-12,
    "max_evaluation_condition": 1e4,
}
KB = 1.380649e-23 / 1.602176634e-19
THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)
SOURCES = {
    "experiments/epw_frequency.py",
    "experiments/audit_epw_frequency.py",
    "kernelgecp/integrations/epw.py",
    "kernelgecp/integrations/epw_frequency.py",
    "kernelgecp/integrations/_frequency_cross.py",
    "kernelgecp/_basis_io.py",
    "kernelgecp/basis.py",
    "kernelgecp/fitting.py",
    "kernelgecp/types.py",
    "kernelgecp/__init__.py",
    "kernelgecp/integrations/__init__.py",
    "kernelgecp/_basis_numerics.py",
    "kernelgecp/_reference_backend.py",
    "kernelgecp/gecp.py",
    "kernelgecp/kernels.py",
    "scripts/providers/epw_anisotropic_fixture.py",
    "environments/epw-science/pyproject.toml",
    "environments/epw-science/uv.lock",
}


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    obj = json.loads(Path(path).read_text())
    require(isinstance(obj, dict) and obj, f"empty metadata: {path}")
    return obj


def close(a, b, name, rtol=2e-10, atol=2e-12):
    require(
        np.shape(a) == np.shape(b) and np.allclose(a, b, rtol=rtol, atol=atol), name
    )


def blob(rev, path):
    require(re.fullmatch(r"[0-9a-f]{40}", rev) is not None, "invalid source revision")
    require(
        subprocess.run(
            ["git", "merge-base", "--is-ancestor", rev, "HEAD"],
            cwd=ROOT,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).returncode
        == 0,
        "generation source is not an ancestor",
    )
    return subprocess.check_output(["git", "show", f"{rev}:{path}"], cwd=ROOT)


def provenance(p):
    require(
        p.get("tracked_clean") is True
        and re.match(r"^3\.12\.", p.get("python", "")) is not None,
        "unclean or wrong Python provenance",
    )
    require(
        p.get("versions")
        == {"numpy": "2.2.6", "scipy": "1.18.1", "mpmath": "1.3.0", "pytest": "9.0.3"},
        "scientific version mismatch",
    )
    require(all(p.get("threads", {}).get(k) == "1" for k in THREADS), "thread metadata")
    require(
        p.get("platform")
        and p.get("python_executable")
        and p.get("python_prefix", "").endswith("environments/epw-science/.venv"),
        "runtime provenance missing",
    )
    hashes = p.get("source_hashes", {})
    require(set(hashes) == SOURCES, "source inventory mismatch")
    for name, digest in hashes.items():
        require(
            hashlib.sha256(blob(p["revision"], name)).hexdigest() == digest,
            f"generation source changed: {name}",
        )


def native(native_root, context, allowed):
    root = Path(native_root)
    receipt_path = root / "epw-anisotropic.json"
    rec = read(receipt_path)
    execution = rec.get("native_anisotropic_execution", {})
    require(
        rec.get("schema") == "kernelgecp.epw-native-anisotropic-fixture.v1"
        and rec.get("status") == "passed"
        and rec.get("provider") == {"qe": "7.5", "epw": "6.0"}
        and rec.get("wscut_ev") == 1.0,
        "native receipt",
    )
    require(
        hashlib.sha256(
            blob(rec["source_revision"], "scripts/providers/epw_anisotropic_fixture.py")
        ).hexdigest()
        == rec.get("script_sha256"),
        "native producer Git blob mismatch",
    )
    require(
        rec.get("upstream_input_sha256")
        == "a67e51f747674319d3c5659e402529bcdf3d571b54f12d2802f2531d74ba8a82"
        and rec.get("upstream_io_sha256")
        == "416ee21c65fd49b2a04894c2d5e18ffb1a5e8f7c223cab950a28f760cf56ec36",
        "QE source identity",
    )
    baseline = rec.get("baseline", {})
    checks = baseline.get("checks", {})
    require(
        bool(checks) and all(v is True for v in checks.values()),
        "native baseline checks",
    )
    binaries = [
        r.get("executable_sha256")
        for r in baseline.get("native_executions", [])
        if r.get("input") == "epw1.in"
    ]
    require(
        execution.get("returncode") == 0
        and len(binaries) == 1
        and binaries[0] == execution.get("executable_sha256")
        and re.fullmatch(r"[0-9a-f]{64}", binaries[0] or "") is not None,
        "native binary identity",
    )
    work = root / "epw-native/test-suite/epw_super"
    inp = work / "epw6.in"
    log = work / "epw6.out"
    require(
        sha(inp)
        == rec.get("executed_input_sha256")
        == execution.get("input_sha256")
        == "7a74f36be9a1604ec8b4184e1bb985f14d237835ec43b19c6ac7e9f4dc2937e1"
        and sha(log) == execution.get("log_sha256"),
        "native input/log hash",
    )
    logtext = log.read_text()
    number = r"([+-]?\d+\.\d+[Ee][+-]\d+)"
    ef = re.findall(r"Fermi level \(eV\)\s*=\s*" + number, logtext)
    mu = re.findall(r"Chemical potential \(itemp\s*=\s*1\)\s*=\s*" + number, logtext)
    require(
        len(ef) == len(mu) == 1
        and len(re.findall("Convergence was reached in nsiter", logtext)) >= 3
        and abs(float(mu[0]) - float(ef[0])) <= 1e-10,
        "native chemical potential/convergence",
    )
    require(
        context.get("reference_fermi_energy_ev") == float(ef[0])
        and context.get("chemical_potential_offset_ev") == 0
        and context.get("temperature_kelvin") == 17
        and context.get("state_order") == "frequency_major_ik_ibnd_ordinal"
        and context.get("provider") == rec["provider"],
        "native context",
    )
    files = rec.get("files", [])
    require(
        len(files) == 3
        and {x.get("temperature_kelvin") for x in files} == {17, 20, 30},
        "native corpus inventory",
    )
    for f in files:
        name = f["filename"]
        require(Path(name).name == name, "unsafe native filename")
        path = work / name
        require(
            sha(path) == f["sha256"] and path.stat().st_size == f["bytes"],
            "native corpus file hash",
        )
    f = next(x for x in files if x["temperature_kelvin"] == 17)
    path = work / f["filename"]
    require(
        context.get("native_sha256") == sha(path)
        and context.get("native_filename") == path.name
        and context.get("source_provenance", {}).get("receipt_sha256")
        == sha(receipt_path),
        "context receipt binding",
    )
    expected_provenance = {
        "receipt_sha256": sha(receipt_path),
        "input_sha256": sha(inp),
        "log_sha256": sha(log),
        "executable_sha256": execution["executable_sha256"],
        "source_revision": rec["source_revision"],
        "source_input_sha256": rec["upstream_input_sha256"],
        "source_io_sha256": rec["upstream_io_sha256"],
        "printed_fermi_ev": ef[0],
        "printed_mu_ev": mu[0],
        "corpus_sha256": {x["filename"]: x["sha256"] for x in files},
    }
    require(
        context.get("source_provenance") == expected_provenance,
        "native provenance binding",
    )
    require(
        context.get("schema") == "kernelgecp.epw-anisotropic-input.v1"
        and context.get("units") == "eV"
        and context.get("positive_only") is True
        and context.get("full_bandwidth") is True
        and context.get("chi_convention") == "native_fbw_column",
        "native conventions",
    )
    count = f["frequency_count"]
    states = f["state_count"]
    require(
        count == 109
        and states > 0
        and context.get("state_count") == states
        and context.get("matsubara_indices") == list(range(count)),
        "native mesh counts",
    )
    require(
        f["columns"]
        == ["frequency_ev", "energy_minus_ef0_ev", "Z", "gap_ev", "shift_ev"]
        and f["positive_only"] is True,
        "native field layout",
    )
    frequency, energies, values = parse_native_file(path, count, states, allowed)
    require(
        frequency == f["frequency_ev"] and energies == f["state_energy_minus_ef0_ev"],
        "native receipt axes",
    )
    return {
        "frequency": np.array(frequency),
        "energy": np.array(energies),
        "fields": values,
        "path": path,
        "context": context,
    }


def parse_native_file(path, count, states, allowed):
    frequency = []
    energies = []
    values = {}
    selected = set(allowed)
    with path.open() as stream:
        require("shift(w)" in stream.readline(), "native header")
        for i in range(count):
            rowvals = []
            for j in range(states):
                line = stream.readline().rstrip("\r\n")
                require(len(line) == 100, "truncated native row")
                w = float(line[:20].replace("D", "E"))
                e = float(line[20:40].replace("D", "E"))
                require(
                    math.isfinite(w) and math.isfinite(e) and w > 0,
                    "native axis nonfinite",
                )
                if j == 0:
                    frequency.append(w)
                else:
                    require(w == frequency[-1], "native frequency block order")
                if i == 0:
                    energies.append(e)
                else:
                    require(e == energies[j], "native state ordering")
                if i in selected:
                    fields = [
                        float(line[k : k + 20].replace("D", "E")) for k in (40, 60, 80)
                    ]
                    require(
                        all(math.isfinite(x) for x in fields) and fields[0] > 0,
                        "native field invalid",
                    )
                    rowvals.append(fields)
            if i in selected:
                values[i] = rowvals
        require(not stream.readline(), "extra native rows")
    return frequency, energies, values


def partitions():
    return {
        "train": [i for i in range(109) if i not in CAL + FINAL],
        "calibration": CAL,
        "final": FINAL,
    }


def truth(data, rows):
    w = data["frequency"][rows, None]
    xi = data["energy"][None, :] - data["context"]["chemical_potential_offset_ev"]
    fields = np.array([data["fields"][int(i)] for i in rows])
    z = fields[:, :, 0]
    gap = fields[:, :, 1]
    chi = fields[:, :, 2]
    xi = xi + chi
    phi = z * gap
    denominator = (w * z) ** 2 + xi**2 + phi**2
    require(
        np.all(denominator > 0) and np.all(np.isfinite(denominator)),
        "native Nambu singularity",
    )
    green = np.empty((*z.shape, 2, 2), complex)
    green[:, :, 0, 0] = -(1j * w * z + xi) / denominator
    green[:, :, 1, 1] = -(1j * w * z - xi) / denominator
    green[:, :, 0, 1] = green[:, :, 1, 0] = -phi / denominator
    return green, z, gap, chi


def direct_design(x, nodes, poles, weights):
    require(
        len(nodes) == len(poles) == len(weights) and len(nodes) > 0,
        "cardinal dimensions",
    )
    result = []
    conditions = []
    with mp.workprec(256):
        ns = list(map(lambda x: mp.mpf(float(x)), nodes))
        ps = list(map(lambda x: mp.mpf(float(x)), poles))
        ws = list(map(lambda x: mp.mpf(float(x)), weights))
        for value in x:
            v = mp.mpf(float(value))
            hit = [j for j, n in enumerate(ns) if n == v]
            if hit:
                require(len(hit) == 1, "duplicate nodes")
                row = [mp.mpf(int(j == hit[0])) for j in range(len(ns))]
                condition = mp.mpf(1)
            else:
                product = mp.fprod(v - n for n in ns) / mp.fprod(v - p for p in ps)
                row = [w * product / (v - n) for w, n in zip(ws, ns, strict=True)]
                terms = [w / (v - n) for w, n in zip(ws, ns, strict=True)]
                denominator = 1 + mp.fsum(terms)
                require(denominator != 0, "cardinal singularity")
                condition = (1 + mp.fsum(abs(t) for t in terms)) / abs(denominator)
            result.append([float(y) for y in row])
            conditions.append(float(condition))
    result = np.asarray(result)
    require(
        np.all(np.isfinite(result))
        and all(math.isfinite(c) and c <= 1e4 for c in conditions),
        "cardinal evaluation condition",
    )
    return result, max(conditions)


def recover(g, w, energy, offset):
    det = g[:, :, 0, 0] * g[:, :, 1, 1] - g[:, :, 0, 1] * g[:, :, 1, 0]
    require(np.all(np.isfinite(det)) and np.all(det != 0), "singular predicted Nambu")
    a = g[:, :, 1, 1] / det
    d = g[:, :, 0, 0] / det
    b = -g[:, :, 0, 1] / det
    c = -g[:, :, 1, 0] / det
    z = (a + d) / (2j * w[:, None])
    chi = (d - a) / 2 - energy[None, :] + offset
    gap = -(b + c) / (2 * z)
    require(
        np.all(np.isfinite(z))
        and np.all(np.isfinite(gap))
        and np.all(np.isfinite(chi)),
        "nonfinite recovered fields",
    )
    singular = np.linalg.svd(g, compute_uv=False)
    condition = singular[..., 0] / singular[..., -1]
    gain = 1 / (singular[..., -1] ** 2 * w[:, None])
    require(
        np.max(condition) <= 1e10 and np.max(gain) <= 1e10, "unsafe recovery condition"
    )
    return z, gap, chi, float(np.max(condition)), float(np.max(gain))


def archive(path, data, rank, parts):
    record = read(path)
    require(
        record.get("schema") == "kernelgecp.epw-parity-cross-application.v1"
        and record.get("representation") == "explicit_real_eliashberg_nambu_parity"
        and record.get("units") == "eV"
        and record.get("representation_units") == {
            "observation_x_nodes_poles": "eV^2",
            "cardinal_weights": "eV^2",
            "channels_and_coefficients": ["eV^-1", "eV^-2", "eV^-1"],
            "green_function": "eV^-1",
            "z": "dimensionless",
        },
        "fit manifest schema",
    )
    copy = {k: v for k, v in record.items() if k != "record_sha256"}
    encoded = json.dumps(
        copy, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    ).encode()
    require(
        hashlib.sha256(encoded).hexdigest() == record.get("record_sha256"),
        "fit manifest hash",
    )
    component = record.get("archive", {})
    filename = component.get("file")
    require(
        isinstance(filename, str)
        and filename == Path(path).with_suffix(".npz").name
        and sha(Path(path).parent / filename) == component.get("sha256"),
        "fit archive hash",
    )
    with np.load(Path(path).parent / filename, allow_pickle=False) as z:
        require(
            set(z.files)
            == {
                "__metadata__",
                "observation_x",
                "nodes",
                "poles",
                "weights",
                "coefficients",
                "fit_rows",
                "frequencies_ev",
                "state_energies_ev",
            },
            "fit payload inventory",
        )
        a = {k: z[k] for k in z.files}
        metadata = json.loads(a.pop("__metadata__").tobytes().decode())
    require(
        metadata.get("format") == "kernelgecp.epw_nambu_parity_cross"
        and metadata.get("version") == 2
        and metadata.get("completion_policy") == "explicit_real_eliashberg_nambu_parity"
        and metadata.get("context") == record.get("source_context") == data["context"],
        "fit metadata",
    )
    require(
        np.array_equal(a["fit_rows"], parts["train"])
        and np.array_equal(a["frequencies_ev"], data["frequency"])
        and np.array_equal(a["state_energies_ev"], data["energy"]),
        "fit axes/rows",
    )
    plan = metadata["plan"]
    require(
        plan["rank"] == rank
        and plan["dictionary_size"] == 1201
        and plan["precision_bits"] == 256
        and plan["pole_cutoff_ev"] == 13
        and plan["max_evaluation_condition"] == 1e4,
        "fit controls",
    )
    require(
        a["coefficients"].shape == (rank, len(data["energy"]), 3), "coefficient axes"
    )
    for k in a:
        require(np.all(np.isfinite(a[k])), f"fit array nonfinite: {k}")
    nodes = a["nodes"]
    poles = a["poles"]
    weights = a["weights"]
    require(
        nodes.shape == poles.shape == weights.shape == (rank,)
        and len(np.unique(nodes)) == len(np.unique(poles)) == rank
        and np.all(nodes > 0)
        and np.all(poles <= 0)
        and np.all(poles >= -169)
        and np.all(weights != 0),
        "fit cardinal structure",
    )
    dictionary = -np.r_[0.0, np.geomspace(1e-10, 169.0, 1200)]
    require(np.all(np.isin(poles, dictionary)), "poles not in declared dictionary")
    beta = 1 / (KB * 17)
    canonical = ((2 * np.array(parts["train"]) + 1) * (np.pi / beta)) ** 2
    close(a["observation_x"], canonical, "observation coordinates", rtol=1e-14, atol=0)
    pivots = plan["pivot_rows"]
    require(
        len(pivots) == rank
        and len(set(pivots)) == rank
        and all(type(i) is int and 0 <= i < len(canonical) for i in pivots)
        and np.array_equal(nodes, a["observation_x"][pivots]),
        "pivot rows",
    )
    with mp.workprec(256):
        ns = [mp.mpf(float(x)) for x in nodes]
        ps = [mp.mpf(float(x)) for x in poles]
        expected = np.array(
            [
                float(
                    mp.fprod(n - p for p in ps)
                    / mp.fprod(n - o for j, o in enumerate(ns) if j != i)
                )
                for i, n in enumerate(ns)
            ]
        )
    close(weights, expected, "stored cardinal weights", rtol=1e-14, atol=0)
    design, condition = direct_design(canonical, nodes, poles, weights)
    sv = np.linalg.svd(design, compute_uv=False)
    ls = metadata["diagnostics"]["least_squares"]
    require(
        ls["effective_rank"] == rank
        and ls["rcond"] == 1e-12
        and ls["converged"] is True
        and sv[-1] > 0
        and sv[0] / sv[-1] < 1e8,
        "retained fit rank/condition",
    )
    close(ls["singular_values"], sv, "fit singular values", rtol=1e-8, atol=1e-10)
    require(
        metadata["diagnostics"]["construction_precision_bits"] == 256
        and metadata["diagnostics"]["evaluation_precision_bits"] == 53
        and metadata["diagnostics"]["channel_order"]
        == ["ReG11", "ImG11_over_omega", "ReG12"],
        "fit numerical metadata",
    )
    for value in metadata["diagnostics"]["timing_seconds"].values():
        require(math.isfinite(value) and value >= 0, "fit timing")
    return a, metadata, condition


def metrics(data, a, rows, saved, prefix, report):
    truth_g, z, gap, chi = truth(data, rows)
    omega = (2 * np.asarray(rows) + 1) * (np.pi / (1 / (KB * 17)))
    design, condition = direct_design(omega**2, a["nodes"], a["poles"], a["weights"])
    channels = np.einsum("ij,jsc->isc", design, a["coefficients"])
    predicted = np.empty_like(truth_g)
    predicted[:, :, 0, 0] = channels[:, :, 0] + 1j * omega[:, None] * channels[:, :, 1]
    predicted[:, :, 1, 1] = -channels[:, :, 0] + 1j * omega[:, None] * channels[:, :, 1]
    predicted[:, :, 0, 1] = predicted[:, :, 1, 0] = channels[:, :, 2]
    scale = np.max(np.abs(truth_g), axis=(-2, -1))
    require(np.all(scale > 0), "zero Green scale")
    zz, dd, cc, cond, gain = recover(
        predicted,
        data["frequency"][rows],
        data["energy"],
        data["context"]["chemical_potential_offset_ev"],
    )
    dif = np.abs(predicted - truth_g)
    result = {
        "green_scaled": float(np.max(dif / scale[:, :, None, None])),
        "green_global_scaled": float(np.max(dif) / np.max(np.abs(truth_g))),
        "z": float(np.max(np.abs(zz.real - z))),
        "gap_ev": float(np.max(np.abs(dd.real - gap))),
        "chi_ev": float(np.max(np.abs(cc.real - chi))),
        "imaginary_defect": max(float(np.max(np.abs(q.imag))) for q in (zz, dd, cc)),
    }
    checks = {
        "truth_green": truth_g,
        "predicted_green": predicted,
        "truth_z": z,
        "truth_gap": gap,
        "truth_chi": chi,
        "predicted_z": zz.real,
        "predicted_gap": dd.real,
        "predicted_chi": cc.real,
    }
    for name, val in checks.items():
        require(prefix + name in saved, f"missing {prefix + name}")
        close(
            saved[prefix + name],
            val,
            f"independent {prefix + name}",
            rtol=3e-8,
            atol=2e-10,
        )
    # Recompute the generator's reported scalar metrics from its saved arrays,
    # then separately test direct-product values against the fixed budget.
    p = saved[prefix + "predicted_green"]
    dif0 = np.abs(p - truth_g)
    zz0, dd0, cc0, _, _ = recover(
        p,
        data["frequency"][rows],
        data["energy"],
        data["context"]["chemical_potential_offset_ev"],
    )
    values = {
        "green_scaled": float(np.max(dif0 / scale[:, :, None, None])),
        "green_global_scaled": float(np.max(dif0) / np.max(np.abs(truth_g))),
        "z": float(np.max(np.abs(zz0.real - z))),
        "gap_ev": float(np.max(np.abs(dd0.real - gap))),
        "chi_ev": float(np.max(np.abs(cc0.real - chi))),
        "imaginary_defect": max(float(np.max(np.abs(q.imag))) for q in (zz0, dd0, cc0)),
    }
    for k, v in values.items():
        close(report[k], v, f"reported {prefix}{k}", rtol=2e-5, atol=2e-12)
    require(result["imaginary_defect"] <= 1e-8, "imaginary recovered fields")
    return result, {
        "cardinal_condition_max": condition,
        "recovery_condition_max": cond,
        "z_derivative_gain_max": gain,
    }


def passed(m, factor):
    return (
        all(math.isfinite(m[k]) and m[k] <= factor * v for k, v in BUDGET.items())
        and m["imaginary_defect"] <= 1e-8
    )


def npz(path):
    with np.load(path, allow_pickle=False) as z:
        require(bool(z.files), "empty scientific archive")
        data = {k: z[k] for k in z.files}
    require(
        all(np.all(np.isfinite(v)) for v in data.values()), "nonfinite scientific array"
    )
    return data


def stable(p):
    return {
        k: v
        for k, v in p.items()
        if k
        not in ("revision", "tracked_clean", "development_lock_sha256_context_only")
    }


def audit_calibration(calibration_path, freeze_path, native_root):
    cal = read(calibration_path)
    freeze = read(freeze_path)
    require(
        cal.get("schema") == "kernelgecp.epw-frequency-calibration.v1"
        and cal.get("status") == "passed"
        and cal.get("pilot_noncanonical") is False,
        "noncanonical or failed calibration",
    )
    require(
        freeze.get("schema") == "kernelgecp.epw-frequency-freeze.v1"
        and freeze.get("status") == "passed"
        and freeze.get("pilot_noncanonical") is False,
        "noncanonical or failed freeze",
    )
    provenance(cal["provenance"])
    provenance(freeze["provenance"])
    require(
        stable(cal["provenance"]) == stable(freeze["provenance"])
        and cal["context"] == freeze["context"]
        and cal["partitions"] == freeze["partitions"] == partitions(),
        "freeze calibration binding",
    )
    require(
        sha(calibration_path) == freeze.get("calibration_sha256")
        and freeze.get("rank_ladder") == RANKS
        and freeze.get("budget") == BUDGET
        and freeze.get("controls") == CONTROLS
        and freeze.get("calibration_margin") == 0.1,
        "freeze controls/hash",
    )
    data = native(native_root, freeze["context"], partitions()["train"] + CAL)
    require(sha(data["path"]) == freeze["native_sha256"], "frozen native hash")
    selected = freeze["selected_rank"]
    require(selected in RANKS and cal.get("selected_rank") == selected, "selected rank")
    attempts = cal.get("attempts", [])
    require(
        [a.get("rank") for a in attempts] == RANKS[: RANKS.index(selected) + 1],
        "candidate inventory",
    )
    details = []
    for attempt in attempts:
        rank = attempt["rank"]
        require(
            "error" not in attempt and isinstance(attempt.get("passed"), bool),
            "unresolved candidate cannot establish selection",
        )
        folder = Path(calibration_path).parent
        payload = folder / f"rank-{rank}.npz"
        fit = folder / f"rank-{rank}-fit.json"
        require(
            sha(payload) == attempt.get("arrays_sha256")
            and sha(fit) == attempt.get("fit_json_sha256")
            and sha(fit.with_suffix(".npz")) == attempt.get("fit_npz_sha256"),
            "candidate archive hash",
        )
        a, meta, condition = archive(fit, data, rank, partitions())
        require(
            attempt.get("diagnostics") == meta["diagnostics"],
            "candidate diagnostics binding",
        )
        timing = attempt.get("timing_seconds", {})
        require(
            set(timing)
            == {
                "basis_construction",
                "nambu_observation_construction",
                "least_squares_preparation",
                "coefficient_fit",
                "training_replay",
                "fit_total",
            }
            and all(math.isfinite(v) and v >= 0 for v in timing.values())
            and math.isfinite(attempt.get("construction_and_fit_seconds", -1))
            and attempt["construction_and_fit_seconds"] >= 0,
            "candidate timing metadata",
        )
        samples = npz(payload)
        expected = {
            name + "_" + key
            for name in ("training", "calibration")
            for key in (
                "truth_green",
                "predicted_green",
                "truth_z",
                "truth_gap",
                "truth_chi",
                "predicted_z",
                "predicted_gap",
                "predicted_chi",
            )
        }
        require(set(samples) == expected, "candidate scientific inventory")
        from kernelgecp.integrations.epw_frequency import EPWFrequencyCrossFit

        loaded = EPWFrequencyCrossFit.load(fit)
        for name, rows in (("training", partitions()["train"]), ("calibration", CAL)):
            close(
                loaded.evaluate_matsubara(np.array(rows)),
                samples[name + "_predicted_green"],
                "public calibration archive replay",
                rtol=1e-12,
                atol=1e-13,
            )
        train, _ = metrics(
            data, a, partitions()["train"], samples, "training_", attempt["training"]
        )
        held, _ = metrics(data, a, CAL, samples, "calibration_", attempt["calibration"])
        success = passed(train, 0.1) and passed(held, 0.1)
        require(
            success == attempt["passed"] and success == (rank == selected),
            "first passing rank mismatch",
        )
        details.append(
            {
                "rank": rank,
                "passed": success,
                "training": train,
                "calibration": held,
                "design_condition_max": condition,
            }
        )
    return {"calibration": details, "selected_rank": selected}, data, freeze


def parse_fitted(path, data):
    n = len(data["frequency"])
    count = len(data["energy"])
    fields = []
    with Path(path).open() as stream:
        require("shift(w)" in stream.readline(), "fitted native header")
        for i in range(n):
            block = []
            for j in range(count):
                line = stream.readline().rstrip("\r\n")
                require(len(line) == 100, "fitted native truncation")
                row = [
                    float(line[k : k + 20].replace("D", "E"))
                    for k in (0, 20, 40, 60, 80)
                ]
                require(
                    all(math.isfinite(x) for x in row)
                    and row[0] == data["frequency"][i]
                    and row[1] == data["energy"][j]
                    and row[2] > 0,
                    "fitted native axes/fields",
                )
                block.append(row[2:])
            fields.append(block)
        require(not stream.readline(), "extra fitted native rows")
    return np.asarray(fields)


def audit_campaign(manifest_path, freeze_path, native_root, calibration_path):
    calibration, _, freeze = audit_calibration(
        calibration_path, freeze_path, native_root
    )
    manifest = read(manifest_path)
    require(
        manifest.get("schema") == "kernelgecp.epw-frequency-campaign.v1"
        and manifest.get("status") == "passed"
        and manifest.get("freeze_sha256") == sha(freeze_path),
        "campaign metadata",
    )
    provenance(manifest["provenance"])
    require(
        stable(manifest["provenance"]) == stable(freeze["provenance"]),
        "campaign source/runtime",
    )
    children = manifest.get("children", [])
    require(
        len(children) == 3 and [c.get("repeat") for c in children] == [0, 1, 2],
        "three fresh children required",
    )
    data = native(native_root, freeze["context"], list(range(109)))
    folder = Path(manifest_path).parent
    first = None
    archives = None
    details = []
    for i, child in enumerate(children):
        require(
            child.get("status") == "passed"
            and child.get("exit_code") == 0
            and math.isfinite(child.get("process_wall_seconds", -1))
            and child["process_wall_seconds"] >= 0,
            "failed child exit",
        )
        command = child.get("command", [])
        require(
            isinstance(command, list)
            and len(command) == 14
            and str(Path(command[0]).resolve())
            == manifest["provenance"]["python_executable"]
            and command[1:4] == ["-m", "experiments.epw_frequency", "case"]
            and command[4] == "--native-root"
            and Path(command[5]).resolve() == Path(native_root).resolve()
            and command[6] == "--freeze"
            and Path(command[7]).resolve() == Path(freeze_path).resolve()
            and command[8:10] == ["--chemical-potential-offset-ev", "0.0"]
            and command[10:12] == ["--repeat", str(i)]
            and command[12] == "--output"
            and Path(command[13]).resolve() == (folder / f"case-{i}").resolve(),
            "fresh child command",
        )
        log = folder / f"case-{i}.log"
        case = folder / f"case-{i}"
        recfile = case / "case.json"
        require(
            sha(log) == child.get("log_sha256")
            and sha(recfile) == child.get("receipt_sha256"),
            "child log/receipt hash",
        )
        rec = read(recfile)
        require(
            rec.get("schema") == "kernelgecp.epw-frequency-case.v1"
            and rec.get("status") == "passed"
            and rec.get("repeat") == i
            and rec.get("freeze_sha256") == sha(freeze_path)
            and rec.get("rank") == freeze["selected_rank"]
            and rec.get("partitions") == partitions(),
            "child metadata",
        )
        provenance(rec["provenance"])
        require(
            stable(rec["provenance"]) == stable(freeze["provenance"]),
            "child provenance",
        )
        files = rec.get("files_sha256", {})
        expected = {
            "fit.json",
            "fit.npz",
            "original-replay.dat",
            "fitted-native.dat",
            "scientific-arrays.npz",
        }
        require(
            set(files) == expected
            and all(sha(case / name) == digest for name, digest in files.items()),
            "child archive inventory/hash",
        )
        require(
            (case / "original-replay.dat").read_bytes() == data["path"].read_bytes(),
            "raw lossless replay",
        )
        a, meta, condition = archive(
            case / "fit.json", data, freeze["selected_rank"], partitions()
        )
        samples = npz(case / "scientific-arrays.npz")
        expect = {
            name + "_" + key
            for name in ("training", "calibration", "final", "full")
            for key in (
                "truth_green",
                "predicted_green",
                "truth_z",
                "truth_gap",
                "truth_chi",
                "predicted_z",
                "predicted_gap",
                "predicted_chi",
            )
        } | {"coefficients", "frequencies_ev", "state_energies_ev"}
        require(
            set(samples) == expect
            and np.array_equal(samples["coefficients"], a["coefficients"])
            and np.array_equal(samples["frequencies_ev"], data["frequency"])
            and np.array_equal(samples["state_energies_ev"], data["energy"]),
            "scientific inventory/axes",
        )
        if first is not None:
            require(
                all(np.array_equal(first[k], samples[k]) for k in first),
                "scientific repeats differ",
            )
            require(
                all(np.array_equal(archives[k], a[k]) for k in a),
                "fit numeric repeats differ",
            )
        else:
            first = samples
            archives = a
        # The public facade is checked for compatibility; it is not the oracle.
        from kernelgecp.integrations.epw_frequency import EPWFrequencyCrossFit

        model = EPWFrequencyCrossFit.load(case / "fit.json")
        results = {}
        for name, rows in (
            ("training", partitions()["train"]),
            ("calibration", CAL),
            ("final", FINAL),
            ("full", list(range(109))),
        ):
            m, diagnostics = metrics(
                data, a, rows, samples, name + "_", rec["metrics"][name]
            )
            require(passed(m, 1), f"{name} error exceeds fixed budget")
            close(
                model.evaluate_matsubara(np.array(rows)),
                samples[name + "_predicted_green"],
                f"public fit replay {name}",
                rtol=1e-12,
                atol=1e-13,
            )
            results[name] = {"metrics": m, **diagnostics}
        fitted_context = dict(data["context"])
        fitted_context["source_native_sha256"] = fitted_context["native_sha256"]
        fitted_context["native_sha256"] = sha(case / "fitted-native.dat")
        fitted_context["native_filename"] = "fitted-native.dat"
        require(rec.get("fitted_context") == fitted_context, "fitted export context")
        fitted = parse_fitted(case / "fitted-native.dat", data)
        expected_fields = np.stack(
            [
                samples["full_predicted_z"],
                samples["full_predicted_gap"],
                samples["full_predicted_chi"],
            ],
            axis=-1,
        )
        close(
            fitted,
            expected_fields,
            "native fitted field output",
            rtol=1e-10,
            atol=1e-10,
        )
        recovery = rec.get("recovery", {})
        _, _, _, cond, gain = recover(
            samples["full_predicted_green"],
            data["frequency"],
            data["energy"],
            data["context"]["chemical_potential_offset_ev"],
        )
        close(
            recovery.get("condition_max"),
            cond,
            "reported recovery condition",
            rtol=1e-8,
        )
        close(
            recovery.get("z_derivative_gain_max"),
            gain,
            "reported recovery gain",
            rtol=1e-8,
        )
        require(
            recovery.get("precision_bits") == 53
            and recovery.get("certification_status")
            == "finite_float64_diagnostics_not_a_certificate",
            "recovery precision metadata",
        )
        size = {
            "raw_ascii": data["path"].stat().st_size,
            "raw_binary_arrays": data["frequency"].nbytes
            + data["energy"].nbytes
            + 109 * len(data["energy"]) * 3 * 8,
            "gzip_ascii": len(gzip.compress(data["path"].read_bytes(), mtime=0)),
            "fit_npz_plus_json": (case / "fit.npz").stat().st_size
            + (case / "fit.json").stat().st_size,
        }
        require(rec.get("bytes") == size, "storage counts")
        times = rec.get("timing_seconds", {})
        keys = {
            "construction_and_fit",
            "basis_construction",
            "nambu_observation_construction",
            "least_squares_preparation",
            "coefficient_fit",
            "training_replay",
            "fit_total",
            "evaluation",
            "archive_and_export",
            "total",
        }
        require(
            set(times) == keys
            and all(math.isfinite(x) and x >= 0 for x in times.values())
            and rec.get("peak_rss_bytes", 0) > 0,
            "timing/resource metadata",
        )
        details.append(
            {
                "repeat": i,
                "results": results,
                "bytes": size,
                "training_design_condition_max": condition,
            }
        )
    return {
        "status": "passed",
        "kind": "campaign",
        "calibration": calibration,
        "cases": details,
        "scope": "finite native grid only; no speedup or experimental prediction",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--calibration", type=Path)
    parser.add_argument("--freeze", type=Path, required=True)
    parser.add_argument("--native-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        require(
            args.calibration is not None or args.manifest is not None,
            "--calibration or --manifest required",
        )
        freeze = read(args.freeze)
        cal = args.calibration or Path(freeze["calibration_path"])
        if args.manifest is None:
            details, _, _ = audit_calibration(cal, args.freeze, args.native_root)
            result = {
                "status": "passed",
                "kind": "calibration",
                "details": details,
                "scope": "training and calibration rows only; final rows unscored",
            }
        else:
            result = audit_campaign(args.manifest, args.freeze, args.native_root, cal)
        code = 0
    except Exception as error:
        result = {"status": "failed", "error": f"{type(error).__name__}: {error}"}
        code = 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    return code


if __name__ == "__main__":
    raise SystemExit(main())
