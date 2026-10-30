"""Independent finite-grid audit of the 20 K model-plus-XOR campaign.

The all-field corrections are truth-assisted, not fit predictions. Native solver
acceptance requires a separate actual 20 K restart comparator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import mpmath as mp
import numpy as np

from experiments.audit_epw_frequency import (
    KB,
    THREADS,
    blob,
    close,
    direct_design,
    parse_native_file,
    passed,
    recover,
    require,
    sha,
    truth,
)

HOLDOUT = [45, 49, 53]
BUDGET = {"green_scaled": 1e-6, "z": 1e-5, "gap_ev": 1e-6, "chi_ev": 1e-6}
CONTROLS = {
    "temperature_kelvin": 20,
    "rank": 60,
    "dictionary_size": 1201,
    "pole_cutoff_ev": 13,
    "precision_bits": 256,
    "rcond": 1e-12,
    "max_evaluation_condition": 1e4,
    "real_nambu_parity": True,
    "holdout_rows": HOLDOUT,
    "budget": BUDGET,
    "policy": "lossless_native_xor_v1",
}
KEYS = {
    "observation_x",
    "nodes",
    "poles",
    "weights",
    "coefficients",
    "fit_rows",
    "frequencies_ev",
    "state_energies_ev",
    "corrections",
}
FIELDS = {
    "truth_green",
    "predicted_green",
    "truth_z",
    "truth_gap",
    "truth_chi",
    "predicted_z",
    "predicted_gap",
    "predicted_chi",
}


def read(path):
    def unique(pairs):
        out = {}
        for k, v in pairs:
            require(k not in out, "duplicate JSON key")
            out[k] = v
        return out

    value = json.loads(
        Path(path).read_text(),
        object_pairs_hook=unique,
        parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)),
    )
    require(isinstance(value, dict), "JSON object required")
    return value


def digest_json(value):
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        ).encode()
    ).hexdigest()


def bits(a):
    return np.asarray(a, dtype="<f8", order="C").view("<u8")


def hash_bits(a):
    return hashlib.sha256(np.asarray(a, dtype="<u8", order="C").tobytes()).hexdigest()


def provenance(env):
    require(
        re.fullmatch(r"[0-9a-f]{40}", env["revision"]) is not None
        and env["versions"] == {"numpy": "2.2.6", "scipy": "1.18.1", "mpmath": "1.3.0"}
        and re.fullmatch(r"3\.12\.\d+", env["python"]) is not None
        and env["prefix"].endswith("environments/epw-science/.venv")
        and env["executable"]
        and env["platform"]
        and env["threads"] == {k: "1" for k in THREADS},
        "runtime provenance",
    )
    hashes = env["source_hashes"]
    require(
        {
            "experiments/audit_epw_restart_roundtrip.py",
            "tests/test_epw_restart_roundtrip_audit.py",
            "experiments/epw_restart_roundtrip.py",
            "kernelgecp/integrations/epw_restart_archive.py",
            "environments/epw-science/uv.lock",
            "environments/epw-science/pyproject.toml",
        }
        <= set(hashes),
        "source inventory",
    )
    for p, h in hashes.items():
        require(
            hashlib.sha256(blob(env["revision"], p)).hexdigest() == h,
            "source blob: " + p,
        )


def native(root, context):
    root = Path(root)
    receipt = root / "epw-anisotropic.json"
    rec = read(receipt)
    exe = rec["native_anisotropic_execution"]
    baseline = rec["baseline"]
    require(
        rec["schema"] == "kernelgecp.epw-native-anisotropic-fixture.v1"
        and rec["status"] == "passed"
        and rec["provider"] == {"qe": "7.5", "epw": "6.0"}
        and rec["wscut_ev"] == 1.0
        and exe["returncode"] == 0
        and baseline["checks"]
        and all(v is True for v in baseline["checks"].values()),
        "native receipt",
    )
    require(
        hashlib.sha256(
            blob(rec["source_revision"], "scripts/providers/epw_anisotropic_fixture.py")
        ).hexdigest()
        == rec["script_sha256"],
        "native producer blob",
    )
    binaries = [
        x["executable_sha256"]
        for x in baseline["native_executions"]
        if x["input"] == "epw1.in"
    ]
    require(
        len(binaries) == 1
        and binaries[0] == exe["executable_sha256"]
        and re.fullmatch(r"[0-9a-f]{64}", binaries[0]),
        "executable identity",
    )
    require(
        rec["upstream_input_sha256"]
        == "a67e51f747674319d3c5659e402529bcdf3d571b54f12d2802f2531d74ba8a82"
        and rec["upstream_io_sha256"]
        == "416ee21c65fd49b2a04894c2d5e18ffb1a5e8f7c223cab950a28f760cf56ec36"
        and rec["executed_input_sha256"]
        == "7a74f36be9a1604ec8b4184e1bb985f14d237835ec43b19c6ac7e9f4dc2937e1",
        "upstream identity",
    )
    work = root / "epw-native/test-suite/epw_super"
    inp, log = work / "epw6.in", work / "epw6.out"
    require(
        sha(inp) == exe["input_sha256"] == rec["executed_input_sha256"]
        and sha(log) == exe["log_sha256"],
        "input/log hashes",
    )
    text = log.read_text()
    num = r"([+-]?\d+\.\d+[Ee][+-]\d+)"
    ef = re.findall(r"Fermi level \(eV\)\s*=\s*" + num, text)
    mu = re.findall(r"Chemical potential \(itemp\s*=\s*2\)\s*=\s*" + num, text)
    require(
        len(ef) == len(mu) == 1
        and len(re.findall("Convergence was reached in nsiter", text)) >= 3,
        "convergence",
    )
    files = rec["files"]
    require(
        len(files) == 3 and {x["temperature_kelvin"] for x in files} == {17, 20, 30},
        "corpus inventory",
    )
    for f in files:
        name = f"MgB2.imag_aniso_{f['temperature_kelvin']:06.2f}"
        require(
            f["filename"] == name
            and sha(work / name) == f["sha256"]
            and (work / name).stat().st_size == f["bytes"],
            "corpus hash",
        )
    f = next(x for x in files if x["temperature_kelvin"] == 20)
    path = work / f["filename"]
    n, s = f["frequency_count"], f["state_count"]
    prov = {
        "receipt_sha256": sha(receipt),
        "input_sha256": sha(inp),
        "log_sha256": sha(log),
        "executable_sha256": exe["executable_sha256"],
        "source_revision": rec["source_revision"],
        "source_input_sha256": rec["upstream_input_sha256"],
        "source_io_sha256": rec["upstream_io_sha256"],
        "printed_fermi_ev": ef[0],
        "printed_mu_ev": mu[0],
        "corpus_sha256": {x["filename"]: x["sha256"] for x in files},
    }
    require(
        context
        == {
            "schema": "kernelgecp.epw-anisotropic-input.v1",
            "provider": rec["provider"],
            "temperature_kelvin": 20,
            "units": "eV",
            "state_order": "frequency_major_ik_ibnd_ordinal",
            "state_count": s,
            "matsubara_indices": list(range(n)),
            "positive_only": True,
            "full_bandwidth": True,
            "chi_convention": "native_fbw_column",
            "reference_fermi_energy_ev": float(ef[0]),
            "chemical_potential_offset_ev": context["chemical_potential_offset_ev"],
            "native_sha256": sha(path),
            "native_filename": path.name,
            "source_provenance": prov,
        }
        and abs(float(mu[0]) - float(ef[0]) - context["chemical_potential_offset_ev"])
        <= 1e-10
        and n > max(HOLDOUT)
        and s > 0
        and f["positive_only"] is True
        and f["state_ordering"] == "native_frequency_major_then_ik_then_ibnd_ordinal"
        and f["columns"]
        == ["frequency_ev", "energy_minus_ef0_ev", "Z", "gap_ev", "shift_ev"],
        "context",
    )
    frequency, energy, fields = parse_native_file(path, n, s, range(n))
    require(
        frequency == f["frequency_ev"] and energy == f["state_energy_minus_ef0_ev"],
        "native axes",
    )
    return {
        "frequency": np.array(frequency),
        "energy": np.array(energy),
        "fields": fields,
        "path": path,
        "context": context,
    }


def bundle(path, data):
    rec = read(path)
    require(
        set(rec)
        == {
            "schema",
            "policy",
            "units",
            "provider",
            "state_order",
            "extra_encoded_information",
            "archive",
            "record_sha256",
        }
        and rec["schema"] == "kernelgecp.epw-restart-xor-application.v1"
        and rec["policy"] == CONTROLS["policy"]
        and rec["units"] == "eV"
        and rec["provider"] == data["context"]["provider"]
        and rec["state_order"] == data["context"]["state_order"]
        and rec["extra_encoded_information"]
        == "uint64 XOR corrections for all native Z, gap, and chi fields"
        and digest_json({k: v for k, v in rec.items() if k != "record_sha256"})
        == rec["record_sha256"],
        "bundle manifest",
    )
    require(
        rec["archive"]
        == {"file": "bundle.npz", "sha256": sha(path.with_suffix(".npz"))},
        "bundle hash",
    )
    with np.load(path.with_suffix(".npz"), allow_pickle=False) as z:
        require(
            len(z.files) == len(set(z.files))
            and set(z.files) == KEYS | {"__metadata__"},
            "bundle inventory",
        )
        a = {k: z[k] for k in z.files}
    raw = a.pop("__metadata__")
    require(raw.dtype == np.uint8 and raw.ndim == 1, "metadata encoding")
    meta = json.loads(raw.tobytes().decode())
    require(
        set(meta)
        == {
            "format",
            "policy",
            "plan",
            "context",
            "diagnostics",
            "analysis_provenance",
            "baseline_sha256",
            "target_sha256",
        }
        and meta["format"] == "kernelgecp.epw-restart-model-xor.v1"
        and meta["policy"] == CONTROLS["policy"]
        and meta["context"] == data["context"],
        "bundle metadata",
    )
    for k, v in a.items():
        require(
            v.dtype
            == (
                np.dtype("<u8")
                if k == "corrections"
                else np.dtype("<i8")
                if k == "fit_rows"
                else np.dtype("<f8")
            )
            and np.all(np.isfinite(v)),
            "dtype/nonfinite " + k,
        )
    n, s = len(data["frequency"]), len(data["energy"])
    train = [i for i in range(n) if i not in HOLDOUT]
    require(
        np.array_equal(a["fit_rows"], train)
        and np.array_equal(bits(a["frequencies_ev"]), bits(data["frequency"]))
        and np.array_equal(bits(a["state_energies_ev"]), bits(data["energy"]))
        and a["coefficients"].shape == (60, s, 3)
        and a["corrections"].shape == (n, s, 3),
        "axes/shapes",
    )
    plan = meta["plan"]
    require(
        all(
            plan[k] == v
            for k, v in {
                "rank": 60,
                "dictionary_size": 1201,
                "precision_bits": 256,
                "pole_cutoff_ev": 13,
                "max_evaluation_condition": 1e4,
            }.items()
        ),
        "plan controls",
    )
    nodes, poles, weights = (a[k] for k in ("nodes", "poles", "weights"))
    require(
        nodes.shape == poles.shape == weights.shape == (60,)
        and len(np.unique(nodes)) == len(np.unique(poles)) == 60
        and np.all(nodes > 0)
        and np.all(poles <= 0)
        and np.all(weights != 0)
        and np.all(np.isin(poles, -np.r_[0.0, np.geomspace(1e-10, 169.0, 1200)])),
        "cardinals",
    )
    x = ((2 * np.array(train) + 1) * (np.pi * KB * 20)) ** 2
    close(a["observation_x"], x, "training coordinates", rtol=1e-14, atol=0)
    pivots = plan["pivot_rows"]
    require(
        len(pivots) == len(set(pivots)) == 60
        and all(type(i) is int and 0 <= i < len(train) for i in pivots)
        and np.array_equal(nodes, a["observation_x"][pivots]),
        "pivots",
    )
    with mp.workprec(256):
        ns, ps = [mp.mpf(float(v)) for v in nodes], [mp.mpf(float(v)) for v in poles]
        expected = [
            float(
                mp.fprod(n - p for p in ps)
                / mp.fprod(n - o for j, o in enumerate(ns) if j != i)
            )
            for i, n in enumerate(ns)
        ]
    close(weights, expected, "cardinal weights", rtol=1e-14, atol=0)
    design, condition = direct_design(x, nodes, poles, weights)
    green, _, _, _ = truth(data, train)
    omega = (2 * np.array(train) + 1) * (np.pi * KB * 20)
    channels = np.stack(
        (
            green[..., 0, 0].real,
            green[..., 0, 0].imag / omega[:, None],
            green[..., 0, 1].real,
        ),
        axis=-1,
    )
    coef, _, rank, sv = np.linalg.lstsq(
        design, channels.reshape(len(train), -1), rcond=1e-12
    )
    ls = meta["diagnostics"]["least_squares"]
    require(
        rank == ls["effective_rank"] == 60
        and ls["rcond"] == 1e-12
        and ls["converged"] is True
        and sv[-1] > 0
        and sv[0] / sv[-1] < 1e8,
        "fit rank/condition",
    )
    close(sv, ls["singular_values"], "singular values", rtol=1e-8, atol=1e-10)
    close(
        coef.reshape(60, s, 3),
        a["coefficients"],
        "training least squares",
        rtol=3e-6,
        atol=2e-8,
    )
    target = bits(np.array([data["fields"][i] for i in range(n)]))
    require(
        hash_bits(target) == meta["target_sha256"]
        and hash_bits(np.bitwise_xor(target, a["corrections"]))
        == meta["baseline_sha256"],
        "XOR hashes",
    )
    return a, meta, condition


def metrics(data, a, rows, saved, prefix, report):
    rows = np.asarray(rows, dtype=int)
    actual, z, gap, chi = truth(data, rows)
    omega = (2 * rows + 1) * (np.pi * KB * 20)
    design, condition = direct_design(omega**2, a["nodes"], a["poles"], a["weights"])
    c = np.einsum("ij,jsc->isc", design, a["coefficients"])
    pred = np.empty_like(actual)
    pred[..., 0, 0] = c[..., 0] + 1j * omega[:, None] * c[..., 1]
    pred[..., 1, 1] = -c[..., 0] + 1j * omega[:, None] * c[..., 1]
    pred[..., 0, 1] = pred[..., 1, 0] = c[..., 2]
    zz, dd, cc, rc, gain = recover(
        pred,
        data["frequency"][rows],
        data["energy"],
        data["context"]["chemical_potential_offset_ev"],
    )
    scale = np.max(np.abs(actual), axis=(-2, -1))
    require(np.all(scale > 0), "zero Green scale")

    def errors(p, z0, d0, c0):
        diff = np.abs(p - actual)
        return {
            "green_scaled": float(np.max(diff / scale[..., None, None])),
            "green_global_scaled": float(np.max(diff) / np.max(np.abs(actual))),
            "z": float(np.max(np.abs(z0.real - z))),
            "gap_ev": float(np.max(np.abs(d0.real - gap))),
            "chi_ev": float(np.max(np.abs(c0.real - chi))),
            "imaginary_defect": max(
                float(np.max(np.abs(v.imag))) for v in (z0, d0, c0)
            ),
        }

    vals = {
        "truth_green": actual,
        "predicted_green": pred,
        "truth_z": z,
        "truth_gap": gap,
        "truth_chi": chi,
        "predicted_z": zz.real,
        "predicted_gap": dd.real,
        "predicted_chi": cc.real,
    }
    for k, v in vals.items():
        close(saved[prefix + k], v, "independent " + prefix + k, rtol=3e-8, atol=2e-10)
    p = saved[prefix + "predicted_green"]
    zs, ds, cs, _, _ = recover(
        p,
        data["frequency"][rows],
        data["energy"],
        data["context"]["chemical_potential_offset_ev"],
    )
    for k, v in errors(p, zs, ds, cs).items():
        close(report[k], v, "reported " + prefix + k, rtol=2e-5, atol=2e-12)
    return errors(pred, zz, dd, cc), {
        "cardinal_condition_max": condition,
        "recovery_condition_max": rc,
        "z_derivative_gain_max": gain,
    }


def audit_campaign(manifest_path):
    folder = Path(manifest_path).resolve().parent
    campaign = read(manifest_path)
    freeze_path = folder / "freeze.json"
    freeze = read(freeze_path)
    require(
        campaign["schema"] == "kernelgecp.epw-restart-roundtrip-campaign.v1"
        and campaign["status"] == "passed"
        and campaign["native_acceptance"] == "not_tested"
        and campaign["repeat_arrays_identical"] is True
        and campaign["freeze_sha256"] == sha(freeze_path),
        "campaign",
    )
    require(
        freeze["schema"] == "kernelgecp.epw-restart-roundtrip-freeze.v1"
        and freeze["controls"] == CONTROLS
        and freeze["offset"] == freeze["context"]["chemical_potential_offset_ev"],
        "frozen protocol",
    )
    provenance(freeze["environment"])
    data = native(freeze["native_root"], freeze["context"])
    children = campaign["children"]
    require(
        len(children) == 3 and [x["repeat"] for x in children] == [0, 1, 2],
        "three children",
    )
    first = None
    first_bundle = None
    results = []
    n = len(data["frequency"])
    train = [i for i in range(n) if i not in HOLDOUT]
    for i, child in enumerate(children):
        case = folder / f"case-{i}"
        log = folder / f"case-{i}.log"
        receipt = case / "case.json"
        cmd = child["command"]
        expected = [
            freeze["environment"]["executable"],
            "-m",
            "experiments.epw_restart_roundtrip",
            "case",
            "--native-root",
            str(Path(freeze["native_root"]).resolve()),
            "--temperature-kelvin",
            "20",
            "--chemical-potential-offset-ev",
            repr(freeze["offset"]),
            "--output",
            str(case),
            "--freeze",
            str(freeze_path),
            "--repeat",
            str(i),
        ]
        require(
            len(cmd) == len(expected)
            and str(Path(cmd[0]).resolve()) == expected[0]
            and cmd[1:] == expected[1:]
            and child["exit_code"] == 0
            and sha(log) == child["log_sha256"]
            and sha(receipt) == child["receipt_sha256"],
            "child invocation/exit/hash",
        )
        r = read(receipt)
        require(
            r["schema"] == "kernelgecp.epw-restart-roundtrip-case.v1"
            and r["status"] == "passed"
            and r["repeat"] == i
            and r["freeze_sha256"] == sha(freeze_path)
            and r["native_acceptance"] == "not_tested"
            and r["corrected_values_are_fit_predictions"] is False
            and r["pure_model_heldout_passed"] is True
            and r["lossless_native_bits_equal"] is True
            and r["printed_native_fields_equal"] is True
            and r["fit_rows"] == train,
            "case receipt",
        )
        files = r["files_sha256"]
        expected_files = {
            "model.json",
            "model.npz",
            "bundle.json",
            "bundle.npz",
            "raw-replay.dat",
            "restored-native.dat",
            "scientific-arrays.npz",
        }
        require(
            set(files) == expected_files
            and all(sha(case / k) == v for k, v in files.items()),
            "case files/hash",
        )
        require(
            (case / "raw-replay.dat").read_bytes() == data["path"].read_bytes(),
            "raw replay",
        )
        a, meta, cond = bundle(case / "bundle.json", data)
        require(meta["diagnostics"] == r["diagnostics"], "diagnostics binding")
        model_record = read(case / "model.json")
        require(
            digest_json({k: v for k, v in model_record.items() if k != "record_sha256"})
            == model_record.get("record_sha256")
            and model_record.get("archive")
            == {"file": "model.npz", "sha256": sha(case / "model.npz")}
            and model_record.get("source_context") == data["context"]
            and model_record.get("units") == "eV",
            "model manifest",
        )
        with np.load(case / "model.npz", allow_pickle=False) as z:
            require(
                set(z.files) == (KEYS - {"corrections"}) | {"__metadata__"},
                "model inventory",
            )
            for k in KEYS - {"corrections"}:
                require(np.array_equal(z[k], a[k]), "model bundle " + k)
        with np.load(case / "scientific-arrays.npz", allow_pickle=False) as z:
            require(len(z.files) == len(set(z.files)), "duplicate scientific member")
            samples = {k: z[k] for k in z.files}
        expected_samples = {
            prefix + "_" + k
            for prefix in ("heldout", "training", "full")
            for k in FIELDS
        } | {
            "coefficients",
            "frequencies_ev",
            "state_energies_ev",
            "restored_z",
            "restored_gap",
            "restored_chi",
        }
        require(
            set(samples) == expected_samples
            and all(np.all(np.isfinite(v)) for v in samples.values())
            and all(
                np.array_equal(samples[k], a[k])
                for k in ("coefficients", "frequencies_ev", "state_energies_ev")
            ),
            "scientific inventory",
        )
        if first is not None:
            require(
                all(np.array_equal(first[k], samples[k]) for k in first)
                and all(np.array_equal(first_bundle[k], a[k]) for k in a),
                "repeat arrays",
            )
        else:
            first, first_bundle = samples, a
        # Primary oracle uses direct-product cardinals and independent Nambu algebra.
        detail = {}
        for name, rows in (
            ("heldout", HOLDOUT),
            ("training", train),
            ("full", list(range(n))),
        ):
            m, d = metrics(
                data, a, rows, samples, name + "_", r["pure_model_metrics"][name]
            )
            detail[name] = {"metrics": m, **d}
            if name == "heldout":
                require(passed(m, 1), "independent heldout budget")
        target = np.array([data["fields"][j] for j in range(n)])
        freq, energy, restored = parse_native_file(
            case / "restored-native.dat", n, len(data["energy"]), range(n)
        )
        require(
            np.array_equal(bits(freq), bits(data["frequency"]))
            and np.array_equal(bits(energy), bits(data["energy"]))
            and np.array_equal(
                bits(np.array([restored[j] for j in range(n)])), bits(target)
            ),
            "independent exported native bits",
        )
        for j, k in enumerate(("restored_z", "restored_gap", "restored_chi")):
            require(
                np.array_equal(bits(samples[k]), bits(target[..., j])),
                "saved restored fields",
            )
        # Public decoder is a consistency check, never the independent oracle.
        from kernelgecp.integrations.epw_restart_archive import EPWRestartArchive

        restored_data = EPWRestartArchive.load(
            case / "bundle.json"
        ).reconstruct_native()
        require(
            np.array_equal(
                bits(
                    np.stack(
                        (restored_data.z, restored_data.delta_ev, restored_data.chi_ev),
                        axis=-1,
                    )
                ),
                bits(target),
            ),
            "public decoder consistency",
        )
        require(
            r["bytes"]["raw_ascii"] == data["path"].stat().st_size
            and r["bytes"]["raw_binary_arrays"]
            == data["frequency"].nbytes + data["energy"].nbytes + target.nbytes
            and r["bytes"]["model_only"]
            == sum((case / x).stat().st_size for x in ("model.json", "model.npz"))
            and r["bytes"]["model_plus_xor_bundle"]
            == sum((case / x).stat().st_size for x in ("bundle.json", "bundle.npz")),
            "byte counts",
        )
        results.append(
            {
                "repeat": i,
                "results": detail,
                "training_design_condition_max": cond,
                "xor_correction_bytes": a["corrections"].nbytes,
                "xor_correction_values": a["corrections"].size,
                "bytes": r["bytes"],
            }
        )
    return {
        "schema": "kernelgecp.epw-restart-roundtrip-audit.v1",
        "status": "passed",
        "native_acceptance": "not_tested",
        "scope": (
            "finite 20 K fit-only holdout and all-field truth-assisted exact "
            "reconstruction; native solver comparator required"
        ),
        "public_decoder": "consistency only",
        "cases": results,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit_campaign(args.manifest)
        code = 0
    except Exception as e:
        result = {
            "status": "failed",
            "native_acceptance": "not_tested",
            "error": f"{type(e).__name__}: {e}",
        }
        code = 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n"
    )
    return code


if __name__ == "__main__":
    raise SystemExit(main())
