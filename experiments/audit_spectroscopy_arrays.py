import argparse
import collections
import glob
import hashlib
import json
import os
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser(
    description=(
        "Read-only independent audit of extracted EH evidence; "
        "no sampling or provider replay."
    )
)
parser.add_argument("evidence", help="extracted evidence root")
parser.add_argument("--output", required=True, help="JSON result path")
args = parser.parse_args()
P = os.path.abspath(args.evidence)
out = {"checks": {}, "groups": {}, "blockers": [], "limitations": []}


def read_json(path):
    return json.loads(Path(path).read_text())


def read_npz(path):
    with np.load(path) as archive:
        return {k: archive[k] for k in archive.files}


def sha(f):
    return hashlib.sha256(Path(f).read_bytes()).hexdigest()


def record(ok, msg):
    if not ok:
        out["blockers"].append(msg)


def examine(d, n):
    files = []
    for f in glob.glob(P + "/" + d + "/*.json"):
        j = read_json(f)
        if "array_sha256" in j and "method" in j and os.path.isfile(f[:-5] + ".npz"):
            files.append((f, j))
    record(len(files) == n, d + " case count " + str(len(files)))
    groups = {}
    ids = set()
    maxcov = 0
    maxrms = 0
    maxmass = 0
    allrank = []
    cond = []
    boot = []
    nonconv = 0
    shapes = 0
    for f, j in files:
        name = os.path.basename(f)
        a = read_npz(f[:-5] + ".npz")
        record(sha(f[:-5] + ".npz") == j["array_sha256"], name + " array hash")
        key = (j["method"], j["config"]["default_model"])
        groups.setdefault("/".join(key), []).append(j["whitened_rms"])
        ids.add(
            (j.get("case", "measured"), *key, j.get("seed", name.rsplit("-r", 1)[-1]))
        )
        record(
            j["status"] == "passed"
            and all(j["gates"].values())
            and j["whitened_rms"] <= 2,
            name + " gates",
        )
        record(
            len(j["converged"]) == 33
            and len(j["iterations"]) == 33
            and j["converged"][j["selected_index"]],
            name + " trials",
        )
        nonconv += sum(not x for x in j["converged"])
        shapes += int(
            a["trial_spectra"].shape == (33, 401) and a["chi2"].shape == (33,)
        )
        mass = np.trapezoid(a["density"], a["omega"])
        maxmass = max(maxmass, abs(mass - 1))
        record(np.all(a["density"] >= 0) and abs(mass - 1) <= 1e-6, name + " density")
        c = a["covariance"]
        pp = a["projection"]
        cp = pp @ c @ pp.T
        maxcov = max(maxcov, float(np.max(np.abs(cp - a["projected_covariance"]))))
        record(
            np.allclose(cp, a["projected_covariance"], rtol=1e-8, atol=1e-16),
            name + " projected covariance",
        )
        ck = (
            "validation_covariance"
            if "validation_covariance" in a
            else "query_covariance"
        )
        ref = "query_mean" if "query_mean" in a else "reference"
        C = a[ck]
        delta = a["prediction"] - a[ref]
        ev, U = np.linalg.eigh((C + C.T) / 2)
        if "validation_covariance_modes" in j:
            threshold = j["validation_covariance_modes"]["threshold"]
        else:
            threshold = 0
        keep = ev > threshold
        rms = np.sqrt(np.sum((U[:, keep].T @ delta) ** 2 / ev[keep]) / sum(keep))
        maxrms = max(maxrms, abs(rms - j["whitened_rms"]))
        record(
            abs(rms - j["whitened_rms"]) < 1e-6,
            name + " RMS " + str(rms) + " " + str(j["whitened_rms"]),
        )
        allrank.append(int(keep.sum()))
        cond.append(float(ev[keep].max() / ev[keep].min()))
        count = j["bootstrap_count"]
        record(
            a["bootstrap_spectra"].shape == ((count, 401) if count else (0,))
            and len(j["bootstrap_records"]) == count,
            name + " bootstrap retention",
        )
        if count:
            std = np.std(a["bootstrap_weights"], axis=0, ddof=1)
            m = j.get(
                "bootstrap_window_standard_deviation",
                j.get("metrics", {}).get("bootstrap_window_standard_deviation"),
            )
            record(np.allclose(std, m, atol=1e-12, rtol=1e-8), name + " bootstrap std")
            boot.extend(std.tolist())
            for rec in j["bootstrap_records"]:
                record(
                    len(rec["converged"]) == 33
                    and len(rec["iterations"]) == 33
                    and rec["converged"][rec["selected_index"]],
                    name + " bootstrap selected",
                )
    out["groups"][d] = {
        "count": len(files),
        "unique_ids": len(ids),
        "max_rms_by_method_default": {k: max(v) for k, v in groups.items()},
        "max_rms_discrepancy": maxrms,
        "projection_covariance_max_abs_error": maxcov,
        "max_mass_error": maxmass,
        "validation_rank_range": [min(allrank), max(allrank)],
        "validation_condition_range": [min(cond), max(cond)],
        "bootstrap_window_std_range": [min(boot), max(boot)] if boot else None,
        "nonconverged_main_trials": nonconv,
        "complete_trial_arrays": shapes,
    }


for d, n in [
    ("calibration", 24),
    ("campaign", 72),
    ("native-calibration", 8),
    ("native-final", 24),
]:
    examine(d, n)
m = read_json(P + "/fixed-bath/manifest.json")
r = read_json(P + "/fixed-bath/result.json")
seeds = m["protocol"]["seeds"]
chains = r["chains"]
record(
    len(seeds) == len(set(seeds)) == len(chains) == 256
    and all(s % 2 and str(s).startswith("321") for s in seeds),
    "seed identities",
)
record(
    [c["seed"] for c in chains] == seeds
    and [c["role"] for c in chains] == ["train"] * 128 + ["heldout"] * 128,
    "chain ordering",
)
for c in chains:
    stem = P + f"/fixed-bath/chains/chain-{c['index']:03d}"
    for suf, k in [
        (".json", "receipt_sha256"),
        (".npz", "arrays_sha256"),
        (".timing.json", "timing_sha256"),
    ]:
        record(sha(stem + suf) == c[k], stem + suf + " hash")
out["checks"]["fixed_bath"] = {
    "chain_count": len(chains),
    "unique_seeds": len(set(seeds)),
    "odd_321_seeds": all(s % 2 and str(s).startswith("321") for s in seeds),
    "roles": dict(collections.Counter(c["role"] for c in chains)),
    "manifest_sha256_matches": sha(P + "/fixed-bath/manifest.json")
    == r["manifest_sha256"],
    "calibration_audit_sha256_matches": sha(P + "/native-calibration-audit.json")
    == m["calibration"]["audit_sha256"],
}
for x, v in out["checks"]["fixed_bath"].items():
    if x.endswith("_matches"):
        record(v, x)
out["checks"]["native_final_binding"] = {
    "binding_file_sha256": sha(P + "/fixed-bath-binding.json"),
    "receipt_binding_sha256": read_json(P + "/native-final/final.json")[
        "binding_sha256"
    ],
    "all_sources_match": all(
        read_json(f)["source"]["binding_sha256"] == sha(P + "/fixed-bath-binding.json")
        for f in glob.glob(P + "/native-final/*-r?.json")
        if "source" in read_json(f)
    ),
}
record(out["checks"]["native_final_binding"]["all_sources_match"], "final binding")
out["limitations"] = [
    "Finite-grid held-out chain-mean diagnostics only; no continuous proof, "
    "true native spectrum, unique spectrum, nonlinear predictive coverage, "
    "or method-wide dominance.",
    "New cohort does not erase earlier 6.495 and 2.023 failures "
    "or establish mu0 as their sole cause.",
    "Independent chain means and covariance, stored-prediction RMS, P C P^T, "
    "chain bath bytes and runtime receipts checked; independent forward "
    "quadrature, runtime authenticity, and bootstrap resampling are not recomputed.",
    "Synthetic final 72 receipts reuse 24 method/default/case identities across "
    "repeated runs; these are not 72 independent cases.",
]


# Recompute fixed-bath means and sample-mean covariance independently.
m = read_json(P + "/fixed-bath/manifest.json")
x = []
bathok = True
maximag = 0
times = []
receipts = True
for i in range(256):
    stem = f"{P}/fixed-bath/chains/chain-{i:03d}"
    z = read_npz(stem + ".npz")
    bathok &= (
        hashlib.sha256(z["bath_iw"].tobytes()).hexdigest()
        == m["parent"]["bath_bytes_sha256"]
    )
    g = z["g_tau"]
    maximag = max(maximag, float(abs(g.imag).max()))
    x.append(g.real.mean(axis=1))
    t = read_json(stem + ".timing.json")
    r = read_json(stem + ".json")
    times.append(t["elapsed_seconds"])
    receipts &= (
        t["seed"] == r["seed"]
        and t["solve_status"] == r["solve_status"] == 0
        and t["elapsed_seconds"] > 0
        and t["peak_rss_bytes"] <= m["protocol"]["memory_limit_bytes"]
    )
x = np.array(x)
o = x[:128, ::48]
q = x[128:, 24 + 192 * np.arange(16)]
om = o.mean(axis=0)
qm = q.mean(axis=0)
oc = np.cov(o, rowvar=False, ddof=1) / 128
qc = np.cov(q, rowvar=False, ddof=1) / 128
maxdiff = {
    k: 0.0
    for k in ("observation_values", "query_mean", "covariance", "validation_covariance")
}
for f in glob.glob(P + "/native-final/*-r?.npz"):
    a = read_npz(f)
    for k, v in [
        ("observation_values", om),
        ("query_mean", qm),
        ("covariance", oc),
        ("validation_covariance", qc),
    ]:
        maxdiff[k] = max(maxdiff[k], float(np.max(np.abs(a[k] - v))))
check = {
    "all_bath_bytes_match": bool(bathok),
    "max_imaginary": maximag,
    "all_chain_runtime_receipts_valid": bool(receipts),
    "runtime_seconds_range": [min(times), max(times)],
    "chain_recomputed_max_abs_errors": maxdiff,
}
out["checks"]["chain_recomputation"] = check
if not bathok or not receipts or maximag > 1e-10 or max(maxdiff.values()) > 1e-12:
    out["blockers"].append("chain recomputation failed")

Path(args.output).write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
print("audit complete; blockers:", len(out["blockers"]))
raise SystemExit(bool(out["blockers"]))
