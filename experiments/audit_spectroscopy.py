"""Independent, fail-closed audit of frozen synthetic MaxEnt campaigns.

Checks saved solver outputs without running MaxEnt; no unique-peak or coverage claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad_vec
from scipy.special import expit, ndtr

ROOT = Path(__file__).resolve().parents[1]
CASES = ("offgrid_atoms_independent", "offgrid_atoms_correlated", "smooth_correlated")
METHODS = ("raw", "kernelgecp", "pydlr", "sparse_ir")
DEFAULTS = ("flat", "broad_gaussian")
THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)
ARRAYS = {
    "omega",
    "density",
    "trial_spectra",
    "tau",
    "values",
    "covariance",
    "projected_values",
    "projected_covariance",
    "projection",
    "observation_design",
    "query",
    "query_covariance",
    "prediction",
    "reference",
    "alphas",
    "chi2",
    "bootstrap_weights",
    "bootstrap_spectra",
    "bootstrap_values",
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    value = json.loads(Path(path).read_text())
    require(isinstance(value, dict) and bool(value), f"empty metadata: {path}")
    return value


def close(a, b, message, *, rtol=1e-8, atol=1e-11):
    require(
        np.shape(a) == np.shape(b) and np.allclose(a, b, rtol=rtol, atol=atol), message
    )


def modes(c, rtol=1e-10, atol=1e-18):
    c = np.asarray(c)
    require(
        c.ndim == 2 and c.shape[0] == c.shape[1] and c.size and np.all(np.isfinite(c)),
        "invalid covariance",
    )
    scale = np.max(np.abs(c))
    slack = 128 * np.finfo(float).eps * len(c) * scale
    require(np.max(np.abs(c - c.T)) <= slack, "asymmetric covariance")
    e, v = np.linalg.eigh((c + c.T) / 2)
    require(
        e[0] >= -slack and 0 < rtol < 1 and atol > 0,
        "indefinite covariance or invalid cutoffs",
    )
    keep = e > max(atol, rtol * e[-1])
    require(np.any(keep), "no retained covariance modes")
    return e[keep], v[:, keep]


def rms(residual, e, v):
    return float(np.linalg.norm((v.T @ residual) / np.sqrt(e)) / np.sqrt(len(e)))


def kernel(t, w):
    t = np.asarray(t)
    w = np.asarray(w)
    return -np.exp(
        np.where(w >= 0, -t[..., None] * w, (10 - t[..., None]) * w)
    ) * expit(10 * np.abs(w))


def truth(case, t):
    if case.startswith("offgrid_atoms"):
        return kernel(t, np.array([-1.371, 0.413, 2.117])) @ np.array([0.2, 0.55, 0.25])

    def integrand(w):
        return kernel(t, np.array([w])).reshape(-1) * smooth(w)

    return quad_vec(integrand, -4, 4, epsabs=2e-13, epsrel=2e-13)[0]


def smooth(w):
    out = 0.0
    for mu, sigma, mass in zip(
        (-1.1, 0.55, 1.9), (0.32, 0.45, 0.25), (0.2, 0.55, 0.25), strict=True
    ):
        norm = ndtr((4 - mu) / sigma) - ndtr((-4 - mu) / sigma)
        out += (
            mass
            * np.exp(-0.5 * ((w - mu) / sigma) ** 2)
            / (np.sqrt(2 * np.pi) * sigma * norm)
        )
    return out


def true_cdf(case, x):
    if case.startswith("offgrid_atoms"):
        return sum(
            m * (x >= p)
            for p, m in zip((-1.371, 0.413, 2.117), (0.2, 0.55, 0.25), strict=True)
        )
    return sum(
        m
        * (ndtr((x - mu) / s) - ndtr((-4 - mu) / s))
        / (ndtr((4 - mu) / s) - ndtr((-4 - mu) / s))
        for mu, s, m in zip(
            (-1.1, 0.55, 1.9), (0.32, 0.45, 0.25), (0.2, 0.55, 0.25), strict=True
        )
    )


def spectrum_metrics(case, omega, density):
    cdf = np.r_[0.0, np.cumsum(np.diff(omega) * (density[:-1] + density[1:]) / 2)]
    actual = true_cdf(case, omega)
    boundaries = np.array([-4.0, -1.0, 1.0, 4.0])
    windows = np.diff(np.interp(boundaries, omega, cdf))
    return (
        windows,
        windows - np.diff(true_cdf(case, boundaries)),
        float(np.max(np.abs(cdf - actual))),
        float(np.trapezoid(np.abs(cdf - actual), omega)),
    )


def fingerprint(
    beta, rtol, atol, tau, values, covariance, projection=None, source_tau=None
):
    h = hashlib.sha256()
    arrays = [np.array([beta, rtol, atol]), tau, values, covariance]
    if projection is not None:
        arrays.extend([projection, source_tau])
    for a in arrays:
        a = np.asarray(a)
        h.update(str(a.shape).encode())
        h.update(np.ascontiguousarray(a, dtype="<f8").tobytes())
    return h.hexdigest()


def check_provenance(p):
    require(
        p.get("tracked_clean") is True
        and all(p.get("threads", {}).get(k) == "1" for k in THREADS),
        "unclean source or threads",
    )
    require(
        all(
            p.get(k)
            for k in (
                "source_revision",
                "source_hashes",
                "python",
                "platform",
                "versions",
                "python_executable",
            )
        ),
        "missing runtime provenance",
    )
    rev = p["source_revision"]
    require(
        subprocess.run(
            ["git", "merge-base", "--is-ancestor", rev, "HEAD"],
            cwd=ROOT,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).returncode
        == 0,
        "source revision is not an ancestor",
    )
    required = {
        "experiments/spectroscopy.py",
        "experiments/spectroscopy_campaign.py",
        "experiments/audit_spectroscopy.py",
        "experiments/dlr_adapters.py",
        "kernelgecp/integrations/maxent.py",
        "kernelgecp/basis.py",
        "kernelgecp/fitting.py",
        "environments/locks/triqs-linux-64-conda-explicit.txt",
        "environments/locks/triqs-maxent-pip.txt",
        "environments/locks/dlr-physics-pip.txt",
    }
    require(
        isinstance(p["source_hashes"], dict) and required.issubset(p["source_hashes"]),
        "missing source inventory",
    )
    for name, expected in p["source_hashes"].items():
        blob = subprocess.check_output(["git", "show", f"{rev}:{name}"], cwd=ROOT)
        require(
            hashlib.sha256(blob).hexdigest() == expected, f"Git blob mismatch: {name}"
        )


def jobs(calibration):
    for i, case in enumerate(CASES):
        for method in METHODS:
            for model in DEFAULTS:
                for repeat in range(1 if calibration else 3):
                    job = {
                        "case": case,
                        "method": method,
                        "default_model": model,
                        "seed": (101 if calibration else 1001) + i,
                        "eps": 1e-8,
                        "n_omega": 401,
                        "n_alpha": 33,
                        "query_order": 24 if calibration else 36,
                        "bootstrap": 0 if calibration or repeat else 8,
                        "pilot": calibration,
                        "repeat": repeat,
                    }
                    yield f"{case}-{method}-{model}-r{repeat}", job


def dataset_check(path, a, result):
    d = read(path.with_suffix(".json"))
    require(
        d.get("schema") == "kernelgecp.spectral-dataset.v1"
        and d.get("observation_domain") == "imaginary_time"
        and d.get("beta") == 10.0,
        "dataset schema",
    )
    require(
        d.get("units")
        == {"tau": "inverse_energy", "omega": "energy", "G_tau": "dimensionless"}
        and d.get("state_ordering") == ["scalar"],
        "dataset units/order",
    )
    require(
        d.get("provider")
        == {"name": result["method"], "continuation": "TRIQS/MaxEnt4.0.1"}
        and d.get("calculation_settings") == result["config"],
        "dataset provider/settings",
    )
    require(sha(path) == d.get("file_sha256"), "dataset file hash")
    with np.load(path, allow_pickle=False) as z:
        arrays = {k: z[k] for k in z.files}
    expected = {"tau", "values", "covariance"} | (
        set() if result["method"] == "raw" else {"projection", "projection_source_tau"}
    )
    require(
        set(arrays) == expected and set(d.get("arrays", {})) == expected,
        "dataset array inventory",
    )
    for k, value in arrays.items():
        meta = d["arrays"][k]
        require(
            list(value.shape) == meta["shape"]
            and hashlib.sha256(
                np.ascontiguousarray(value, dtype="<f8").tobytes()
            ).hexdigest()
            == meta["sha256"],
            f"dataset array hash {k}",
        )
    close(arrays["tau"], a["tau"], "dataset tau")
    close(arrays["values"], a["projected_values"], "dataset values")
    close(
        arrays["covariance"],
        a["projected_covariance"],
        "dataset covariance",
        atol=1e-18,
    )
    rtol = d["covariance_rtol"]
    atol = d["covariance_atol"]
    require(rtol == 1e-10 and atol == 1e-18, "unexpected covariance cutoffs")
    e, v = modes(arrays["covariance"], rtol, atol)
    require(
        len(e) == d["effective_rank"] == result["effective_rank"], "covariance rank"
    )
    projection = arrays.get("projection")
    source_tau = arrays.get("projection_source_tau")
    if projection is not None:
        close(projection, a["projection"], "dataset projection")
        close(source_tau, a["tau"], "source tau")
        sourcefp = fingerprint(10.0, rtol, atol, a["tau"], a["values"], a["covariance"])
        require(d["source_fingerprint"] == sourcefp, "source fingerprint")
    else:
        require(d["source_fingerprint"] is None, "raw source fingerprint")
    require(
        d["fingerprint"]
        == fingerprint(
            10.0,
            rtol,
            atol,
            arrays["tau"],
            arrays["values"],
            arrays["covariance"],
            projection,
            source_tau,
        ),
        "dataset fingerprint",
    )
    return e, v


def audit_case(path, job=None, full=True):
    path = Path(path)
    result = read(path)
    npz = path.with_suffix(".npz")
    require(
        result.get("schema") == "kernelgecp.spectroscopy-case.v1"
        and result.get("status") == "passed",
        "case status/schema",
    )
    if job:
        for k in ("case", "method", "seed", "eps", "pilot", "query_order"):
            require(result.get(k) == job[k], f"case {k}")
        require(result.get("bootstrap_count") == job["bootstrap"], "bootstrap count")
        for k in ("default_model", "n_omega", "n_alpha"):
            require(result.get("config", {}).get(k) == job[k], f"config {k}")
    require(sha(npz) == result.get("array_sha256"), "case archive hash")
    with np.load(npz, allow_pickle=False) as z:
        a = {k: z[k] for k in z.files}
    require(
        set(a) == ARRAYS and all(np.all(np.isfinite(x)) for x in a.values()),
        "case arrays missing or nonfinite",
    )
    fit = result.get("coefficient_fit_archive")
    if result["method"] == "kernelgecp":
        require(
            isinstance(fit, dict)
            and fit.get("path") == path.stem + "-gecp-fit.npz"
            and sha(path.parent / fit["path"]) == fit.get("sha256"),
            "GECP archive",
        )
        require(
            result.get("fit_projection_max_difference") is not None
            and result["fit_projection_max_difference"] <= 1e-10,
            "GECP projection fit",
        )
    else:
        require(
            fit is None and result.get("fit_projection_max_difference") is None,
            "unexpected fit archive",
        )
    e, v = dataset_check(path.with_name(path.stem + "-dataset.npz"), a, result)
    if not full:
        return a
    omega = a["omega"]
    t = a["tau"]
    p = a["projection"]
    design = a["observation_design"]
    values = a["values"]
    c = a["covariance"]
    require(
        result["config"]
        == {
            "omega_min": -4.0,
            "omega_max": 4.0,
            "n_omega": 401,
            "alpha_min": 1e-4,
            "alpha_max": 1e4,
            "n_alpha": 33,
            "maxiter": 1000,
            "default_model": result["config"]["default_model"],
            "normalize_mass": 1.0,
        },
        "unfrozen solver configuration",
    )
    if fit:
        from kernelgecp.basis import FermionicFit

        restored = FermionicFit.load(path.parent / fit["path"])
        difference = float(
            np.max(np.abs(restored.evaluate_tau(t) - a["projected_values"]))
        )
        require(difference <= 1e-10, "GECP fit archive projection")
        close(
            result["fit_projection_max_difference"],
            difference,
            "GECP fit difference",
            atol=1e-10,
        )
    require(
        omega.shape == (401,)
        and np.all(np.diff(omega) > 0)
        and t.shape == (65,)
        and p.shape == (65, 65)
        and design.ndim == 2
        and design.shape[0] == 65,
        "mesh/design shapes",
    )
    close(omega, np.linspace(-4, 4, 401), "omega mesh")
    close(t, np.linspace(0, 10, 65), "tau mesh")
    close(p, p.T, "projection symmetry", atol=2e-10)
    close(p @ p, p, "projection idempotence", atol=2e-10)
    if result["method"] == "raw":
        close(p, np.eye(65), "raw projection")
    else:
        sv = np.linalg.svd(design, compute_uv=False)
        diag = result["projection_svd"]
        cutoff = diag["rcond"] * sv[0]
        require(
            diag["effective_rank"] == int(np.sum(sv > cutoff))
            and abs(diag["rank_threshold"] - cutoff) <= 1e-12,
            "design SVD rank",
        )
        require(
            np.linalg.norm(p @ design - design, 2) <= max(2e-9, 2 * cutoff),
            "projection truncation bound",
        )
    close(a["projected_values"], p @ values, "projected data")
    close(a["projected_covariance"], p @ c @ p.T, "projected covariance", atol=2e-18)
    case = result["case"]
    seed = result["seed"]
    corr = case.endswith("correlated")

    def noise_cov(x):
        return 1e-8 * (
            (0.25 * np.eye(len(x)) + 0.75 * np.exp(-np.abs(x[:, None] - x) / 0.6))
            if corr
            else np.eye(len(x))
        )

    close(c, noise_cov(t), "input covariance", atol=1e-18)
    close(
        values,
        truth(case, t)
        + np.linalg.cholesky(c) @ np.random.default_rng(seed).standard_normal(65),
        "synthetic values",
        atol=2e-12,
    )
    query = (np.polynomial.legendre.leggauss(result["query_order"])[0] + 1) * 5
    close(a["query"], query, "heldout query")
    require(
        not np.any(np.isclose(query[:, None], t[None, :], atol=1e-12, rtol=0)),
        "heldout overlaps observations",
    )
    qc = noise_cov(query)
    close(a["query_covariance"], qc, "query covariance", atol=1e-18)
    qe, qv = modes(qc)
    reference = truth(case, query)
    close(a["reference"], reference, "analytic truth", atol=2e-12)
    trials = a["trial_spectra"]
    alphas = a["alphas"]
    chi = a["chi2"]
    conv = np.asarray(result["converged"])
    it = result["iterations"]
    require(
        trials.shape == (33, 401)
        and alphas.shape == (33,)
        and chi.shape == (33,)
        and conv.shape == (33,)
        and len(it) == 33
        and all(type(x) is bool for x in result["converged"])
        and all(type(x) is int and x >= 0 for x in it),
        "trial inventory",
    )
    close(alphas, np.geomspace(1e4, 1e-4, 33), "alpha grid")
    widths = np.r_[
        np.diff(omega)[0] / 2, (omega[2:] - omega[:-2]) / 2, np.diff(omega)[-1] / 2
    ]
    forwards = p @ (kernel(t, omega) * widths) @ trials.T
    measured = np.sum(
        ((v.T @ (forwards - a["projected_values"][:, None])) / np.sqrt(e)[:, None])
        ** 2,
        axis=0,
    )
    close(chi, measured, "trial chi2", rtol=3e-5, atol=2e-5)
    eligible = conv & np.isfinite(chi) & (chi >= 0) & np.all(trials >= 0, axis=1)
    require(np.any(eligible), "no converged positive trial")
    scores = np.full(33, np.inf)
    scores[eligible] = np.abs(
        np.log(np.maximum(chi[eligible], np.finfo(float).tiny) / len(e))
    )
    idx = int(np.argmin(scores))
    require(result["selected_index"] == idx, "alpha selection")
    raw = float(np.trapezoid(trials[idx], omega))
    require(raw > 0, "raw mass")
    close(result["raw_mass"], raw, "raw mass")
    density = trials[idx] / raw
    close(a["density"], density, "normalized selected spectrum")
    require(
        np.all(density >= 0) and abs(np.trapezoid(density, omega) - 1) <= 1e-6,
        "positive normalized spectrum",
    )
    close(result["mass"], np.trapezoid(density, omega), "recorded mass")

    def forward(times):
        return quad_vec(
            lambda w: (
                kernel(times, np.array([w])).reshape(-1) * np.interp(w, omega, density)
            ),
            -4,
            4,
            points=omega[1:-1],
            epsabs=1e-11,
            epsrel=1e-11,
        )[0]

    train = rms(p @ forward(t) - a["projected_values"], e, v)
    close(
        result["normalization_training_whitened_rms"], train, "training RMS", rtol=2e-5
    )
    # Integrate piecewise-linear density independently of the native trapezoid forward.
    prediction = forward(query)
    close(a["prediction"], prediction, "saved prediction", atol=2e-12)
    heldout = rms(prediction - reference, qe, qv)
    require(heldout <= 2, f"heldout whitened RMS {heldout}")
    saved_rms = rms(a["prediction"] - reference, qe, qv)
    close(result["whitened_rms"], saved_rms, "saved heldout RMS")
    wm, err, cmax, cint = spectrum_metrics(case, omega, density)
    metric = result["metrics"]
    for k, x in [
        ("window_weights", wm),
        ("window_weight_errors", err),
        ("cdf_max_error", cmax),
        ("cdf_integrated_absolute_error", cint),
    ]:
        close(metric[k], x, k)
    n = result["bootstrap_count"]
    require(
        result["bootstrap_seed"] == seed + 100000
        and len(result["bootstrap_records"]) == n,
        "bootstrap metadata",
    )
    require(
        a["bootstrap_values"].shape == ((n, 65) if n else (0,))
        and a["bootstrap_spectra"].shape == ((n, 401) if n else (0,))
        and a["bootstrap_weights"].shape == ((n, 3) if n else (0,)),
        "bootstrap shapes",
    )
    rng = np.random.default_rng(seed + 100000)
    bs = []
    errors = []
    for i in range(n):
        close(
            a["bootstrap_values"][i],
            values + np.linalg.cholesky(c) @ rng.standard_normal(65),
            "bootstrap sample",
        )
        d = a["bootstrap_spectra"][i]
        rec = result["bootstrap_records"][i]
        require(
            np.all(d >= 0)
            and abs(np.trapezoid(d, omega) - 1) <= 1e-6
            and rec["raw_mass"] > 0
            and abs(rec["mass"] - 1) <= 1e-6
            and rec["effective_rank"] == len(e),
            "bootstrap spectrum",
        )
        require(
            len(rec["converged"]) == 33
            and len(rec["iterations"]) == 33
            and 0 <= rec["selected_index"] < 33
            and rec["converged"][rec["selected_index"]] is True,
            "bootstrap convergence",
        )
        w, _, cm, _ = spectrum_metrics(case, omega, d)
        bs.append(w)
        errors.append(cm)
    if n:
        close(a["bootstrap_weights"], bs, "bootstrap windows")
        close(
            metric["bootstrap_window_standard_deviation"],
            np.std(bs, axis=0, ddof=1),
            "bootstrap variation",
        )
        close(
            metric["bootstrap_cdf_error_range"],
            [min(errors), max(errors)],
            "bootstrap cdf",
        )
    return a


def audit_campaign(manifest_path, freeze_path):
    manifest_path = Path(manifest_path)
    freeze_path = Path(freeze_path)
    m = read(manifest_path)
    f = read(freeze_path)
    kind = m.get("kind")
    require(
        kind in ("calibration", "campaign")
        and m.get("schema") == "kernelgecp.spectroscopy-campaign.v1"
        and m.get("status") == "passed",
        "campaign schema/status",
    )
    require(
        f.get("schema") == "kernelgecp.spectroscopy-freeze.v1"
        and f.get("status") == "passed",
        "freeze status/schema",
    )
    controls = {
        "cases": list(CASES),
        "methods": list(METHODS),
        "defaults": list(DEFAULTS),
        "eps": 1e-8,
        "n_omega": 401,
        "n_alpha": 33,
        "calibration_seeds": [101, 102, 103],
        "calibration_query_order": 24,
        "final_seeds": [1001, 1002, 1003],
        "final_query_order": 36,
        "bootstrap_first_repeat": 8,
        "repeats": 3,
    }
    require(
        m.get("controls") == controls == f.get("controls"), "frozen controls mismatch"
    )
    p = m.get("provenance", {})
    check_provenance(p)

    def stable(x):
        return {
            k: v for k, v in x.items() if k not in ("source_revision", "tracked_clean")
        }

    require(
        stable(p) == stable(f.get("provenance", {})), "freeze runtime/source mismatch"
    )
    check_provenance(f["provenance"])
    if kind == "campaign":
        require(m.get("freeze_sha256") == sha(freeze_path), "freeze hash")
    else:
        require(
            f.get("calibration_manifest_sha256") == sha(manifest_path),
            "calibration manifest hash",
        )
    expected = dict(jobs(kind == "calibration"))
    children = m.get("children")
    require(
        isinstance(children, dict) and set(children) == set(expected),
        "incomplete or extra campaign children",
    )
    directory = manifest_path.parent
    prior = {}
    audited = 0
    for name, job in expected.items():
        record = children[name]
        path = directory / (name + ".json")
        require(
            record.get("job") == job
            and record.get("status") == "passed"
            and record.get("exit_code") == 0,
            f"{name}: job/exit",
        )
        require(
            sha(directory / (name + ".log")) == record.get("log_sha256"),
            f"{name}: log hash",
        )
        meta = read(path)
        require(
            meta.get("provenance")
            == {k: v for k, v in p.items() if k != "python_executable"},
            f"{name}: provenance",
        )
        hashes = record.get("hashes", {})
        fit = meta.get("coefficient_fit_archive")
        expected_hashes = {
            "json_sha256": sha(path),
            "npz_sha256": sha(path.with_suffix(".npz")),
            "dataset_sha256": sha(directory / (name + "-dataset.npz")),
            "dataset_json_sha256": sha(directory / (name + "-dataset.json")),
            "fit_sha256": sha(directory / fit["path"]) if fit else None,
        }
        require(hashes == expected_hashes, f"{name}: artifact hashes")
        key = (job["case"], job["method"], job["default_model"])
        a = audit_case(path, job, full=job["repeat"] == 0)
        if job["repeat"] == 0:
            prior[key] = a
            audited += 1
        else:
            first = prior[key]
            for k in ARRAYS:
                if not k.startswith("bootstrap_"):
                    require(
                        np.array_equal(a[k], first[k]), f"{name}: repeat mismatch {k}"
                    )
            # Metadata affecting selection and reported metrics must also reproduce.
            base = read(directory / (name.rsplit("-r", 1)[0] + "-r0.json"))
            for k in (
                "selected_index",
                "converged",
                "iterations",
                "raw_mass",
                "mass",
                "effective_rank",
                "projection_svd",
                "fit_projection_max_difference",
                "normalization_training_whitened_rms",
                "whitened_rms",
            ):
                require(meta.get(k) == base.get(k), f"{name}: repeat metadata {k}")
    return {
        "status": "passed",
        "kind": kind,
        "cases": len(expected),
        "independently_audited_cases": audited,
        "manifest_sha256": sha(manifest_path),
        "freeze_sha256": sha(freeze_path),
        "claim_scope": (
            "finite synthetic cases; bootstrap variation is not "
            "coverage or unique peak recovery"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--freeze", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit_campaign(args.manifest, args.freeze)
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
    sys.exit(main())
