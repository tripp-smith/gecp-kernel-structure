"""Independent audit of frozen CT-HYB spectroscopy; never runs MaxEnt."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from scipy.integrate import quad_vec
from scipy.special import expit

OBS = np.arange(0, 3073, 48)
QUERY = 24 + 192 * np.arange(16)
METHODS = ("raw", "kernelgecp", "pydlr", "sparse_ir")
MODELS = ("flat", "broad_gaussian")
REV = "9d28053060329cd20647dc03e038a383c351adf1"
BINDING = "bdc256787c761d95d00ecd59b66ba3763a85fd8da4571e64306ae8c9f5ecef88"
ROOT = Path(__file__).resolve().parents[1]


def require(x, msg):
    if not x:
        raise ValueError(msg)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    x = json.loads(Path(p).read_text())
    require(isinstance(x, dict) and bool(x), "invalid JSON")
    return x


def close(a, b, msg, rtol=1e-8, atol=1e-11):
    a, b = np.asarray(a), np.asarray(b)
    require(
        a.shape == b.shape
        and np.isfinite(a).all()
        and np.isfinite(b).all()
        and np.allclose(a, b, rtol=rtol, atol=atol),
        msg,
    )


def moments(a):
    a = np.asarray(a, float)
    require(a.ndim == 2 and len(a) > 1 and np.isfinite(a).all(), "sample matrix")
    # Separate two-pass accumulation, without calling producer covariance code.
    mu = np.sum(a, axis=0) / len(a)
    d = a - mu
    return mu, np.einsum("ni,nj->ij", d, d) / (len(a) * (len(a) - 1))


def modes(c):
    c = np.asarray(c)
    require(
        c.ndim == 2 and c.shape[0] == c.shape[1] and np.isfinite(c).all(), "covariance"
    )
    scale = np.max(np.abs(c))
    slack = 128 * np.finfo(float).eps * len(c) * scale
    require(np.max(np.abs(c - c.T)) <= slack, "covariance asymmetry")
    e, v = np.linalg.eigh((c + c.T) / 2)
    require(e[0] >= -slack, "negative covariance")
    keep = e > max(1e-18, 1e-4 * e[-1])
    require(np.any(keep), "no covariance modes")
    return e[keep], v[:, keep]


def weights(w, d):
    cum = np.r_[0, np.cumsum(np.diff(w) * (d[1:] + d[:-1]) / 2)]
    return np.diff(np.interp([-4, -1, 1, 4], w, cum))


def forward(t, w, d):
    t = np.asarray(t)

    def integrand(x):
        k = -np.exp(np.where(x >= 0, -t * x, (10 - t) * x)) * expit(10 * abs(x))
        return k * np.interp(x, w, d)

    return quad_vec(integrand, -4, 4, points=w[1:-1], epsabs=1e-11, epsrel=1e-11)[0]


def raw_cohort(root, audit, binding, stage):
    root = Path(root)
    require(sha(binding) == BINDING, "EE binding hash")
    b = read(binding)
    a = read(audit)
    manifest = read(root / "manifest.json")
    result = read(root / "result.json")
    require(
        b.get("status") == "passed"
        and b.get("generation_revision") == REV
        and b.get("audit_receipt_sha256") == sha(audit),
        "EE binding",
    )
    require(
        a.get("status") == "passed"
        and a.get("generation_revision") == REV
        and a.get("gates")
        and all(v is True for v in a["gates"].values()),
        "EE independent audit",
    )
    audit_blob = subprocess.check_output(
        ["git", "show", a["audit_revision"] + ":experiments/audit_cthyb_dmft.py"],
        cwd=ROOT,
    )
    require(
        hashlib.sha256(audit_blob).hexdigest() == a.get("audit_sha256"),
        "EE auditor source",
    )
    for key, name in [
        ("manifest_sha256", "manifest.json"),
        ("result_sha256", "result.json"),
        ("freeze_sha256", "freeze.json"),
    ]:
        require(b.get(key) == sha(root / name), "EE " + name)
    require(
        manifest.get("source_revision") == REV
        and result.get("identity") == manifest
        and result.get("status") == "passed",
        "EE provenance",
    )
    cfg = manifest["config"]
    require(
        cfg.get("seed_mapping") == "odd_stride_2_v1"
        and cfg.get("n_tau") == 3073
        and cfg.get("final_replicas") == 64,
        "EE configuration",
    )
    role, offset = (
        ("kernelgecp-final", 40000)
        if stage == "calibration"
        else ("dense-final", 30000)
    )
    samples = []
    fingerprints = set()
    bath = None
    for i in range(64):
        p = root / role / f"replica-{i:02d}.json"
        rel = p.relative_to(root).as_posix()
        d = read(p)
        require(
            b["measurement_archives"].get(rel)
            == sha(p)
            == result["measurement_archives"].get(rel),
            "chain receipt hash",
        )
        require(
            d.get("seed") == cfg["seed_base"] + offset + 1 + 2 * i
            and d.get("cycles") == 1600000
            and d.get("solve_status") == 0
            and d.get("identity") == manifest
            and d.get("average_sign", 0) > 0.99,
            "chain seed/cycles/sign",
        )
        require(
            d.get("arrays_sha256") == sha(p.with_suffix(".npz")), "chain array hash"
        )
        with np.load(p.with_suffix(".npz"), allow_pickle=False) as z:
            x = z["g_tau"]
            t = z["bath_iw"]
            require(
                x.shape == (3073, 2)
                and np.isfinite(x).all()
                and np.max(np.abs(x.imag)) <= 1e-10,
                "spin data",
            )
            if bath is None:
                bath = t.copy()
            else:
                require(np.array_equal(bath, t), "cohort bath")
            h = hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()
            require(h not in fingerprints, "duplicate chain data")
            fingerprints.add(h)
            samples.append((x[:, 0].real + x[:, 1].real) / 2)
    samples = np.array(samples)
    om, oc = moments(samples[:32, OBS])
    qm, qc = moments(samples[32:, QUERY])
    return samples, om, oc, qm, qc


def audit_case(path, method, model, repeat, stats):
    d = read(path)
    npz = path.with_suffix(".npz")
    require(
        d.get("schema") == "kernelgecp.cthyb-spectroscopy.v1"
        and d.get("status") == "passed"
        and d.get("method") == method
        and d.get("config", {}).get("default_model") == model
        and d.get("array_sha256") == sha(npz),
        "case identity",
    )
    require(
        d.get("covariance_rtol") == 1e-4
        and d.get("observation_indices") == OBS.tolist()
        and d.get("query_indices") == QUERY.tolist(),
        "partition",
    )
    require(
        d.get("eps") == 1e-8
        and d.get("config")
        == {
            "omega_min": -4.0,
            "omega_max": 4.0,
            "n_omega": 401,
            "alpha_min": 1e-4,
            "alpha_max": 1e4,
            "n_alpha": 33,
            "maxiter": 1000,
            "default_model": model,
            "normalize_mass": 1.0,
        },
        "solver controls",
    )
    count = 8 if repeat == 0 else 0
    require(
        d.get("bootstrap_count") == count and d.get("bootstrap_seed") == 193,
        "bootstrap count/seed",
    )
    with np.load(npz, allow_pickle=False) as z:
        a = {k: z[k] for k in z.files}
    required = {
        "omega",
        "density",
        "observation_tau",
        "observation_values",
        "covariance",
        "projected_values",
        "projected_covariance",
        "projection",
        "observation_design",
        "native_kernel",
        "projected_kernel",
        "query",
        "query_mean",
        "validation_covariance",
        "prediction",
        "alphas",
        "chi2",
        "trial_spectra",
        "bootstrap_weights",
        "bootstrap_spectra",
        "bootstrap_indices",
    }
    require(
        set(a) == required and all(np.isfinite(x).all() for x in a.values()),
        "array inventory",
    )
    samples, om, oc, qm, qc = stats
    t = np.linspace(0, 10, 3073)
    close(a["observation_tau"], t[OBS], "obs tau")
    close(a["query"], t[QUERY], "query tau")
    require(not np.intersect1d(OBS, QUERY).size, "overlap")
    close(a["observation_values"], om, "train mean")
    close(a["covariance"], oc, "train covariance", atol=1e-18)
    close(a["query_mean"], qm, "heldout mean")
    close(a["validation_covariance"], qc, "heldout covariance", atol=1e-18)
    p = a["projection"]
    require(p.shape == (65, 65), "projection shape")
    close(p, p.T, "projection symmetry", atol=2e-10)
    close(p @ p, p, "projection idempotence", atol=2e-10)
    if method == "raw":
        close(p, np.eye(65), "raw projection")
    design = a["observation_design"]
    require(design.ndim == 2 and design.shape[0] == 65, "design")
    if method != "raw":
        u, s, _ = np.linalg.svd(design, full_matrices=False)
        keep = s > 1e-10 * s[0]
        close(p, u[:, keep] @ u[:, keep].T, "independent projection", atol=2e-8)
    native_times = t[OBS]
    w0 = np.linspace(-4, 4, 401)
    widths = np.r_[np.diff(w0)[0] / 2, (w0[2:] - w0[:-2]) / 2, np.diff(w0)[-1] / 2]
    native = (
        -np.exp(
            np.where(
                w0[None, :] >= 0,
                -native_times[:, None] * w0[None, :],
                (10 - native_times[:, None]) * w0[None, :],
            )
        )
        * expit(10 * np.abs(w0))[None, :]
    )
    close(a["native_kernel"], native, "native kernel", rtol=1e-10, atol=1e-12)
    close(a["projected_values"], p @ om, "projected values")
    close(a["projected_covariance"], p @ oc @ p.T, "projected covariance", atol=2e-18)
    close(a["projected_kernel"], p @ a["native_kernel"], "projected kernel")
    for key, cov in [
        ("training_covariance_modes", oc),
        ("projected_covariance_modes", p @ oc @ p.T),
        ("validation_covariance_modes", qc),
    ]:
        eig = np.linalg.eigvalsh((cov + cov.T) / 2)
        record = d[key]
        close(
            record["eigenvalues"],
            eig,
            "covariance eigenspectrum",
            rtol=1e-6,
            atol=1e-18,
        )
        close(
            record["threshold"],
            max(1e-18, 1e-4 * eig[-1]),
            "covariance threshold",
            rtol=1e-6,
            atol=1e-18,
        )
        require(
            record["relative_cutoff"] == 1e-4
            and record["absolute_cutoff"] == 1e-18
            and record["effective_rank"] == len(modes(cov)[0]),
            "covariance mode record",
        )
    e, v = modes(p @ oc @ p.T)
    qe, qv = modes(qc)
    require(
        d["training_rank"] == len(modes(oc)[0])
        and d["projected_rank"] == len(e)
        and d["validation_rank"] == len(qe),
        "covariance rank",
    )
    w = a["omega"]
    den = a["density"]
    close(w, np.linspace(-4, 4, 401), "omega")
    require(
        den.shape == (401,)
        and np.min(den) >= -1e-6
        and abs(np.trapezoid(den, w) - 1) <= 1e-6,
        "spectrum positivity or mass",
    )
    trials = a["trial_spectra"]
    chi = a["chi2"]
    alph = a["alphas"]
    conv = d["converged"]
    require(
        trials.shape == (33, 401)
        and chi.shape == (33,)
        and len(conv) == len(d["iterations"]) == 33
        and all(type(x) is bool for x in conv),
        "trials",
    )
    close(alph, np.geomspace(1e4, 1e-4, 33), "alphas")
    calc = np.sum(
        (
            (
                v.T
                @ (
                    p @ a["native_kernel"] @ (widths[:, None] * trials.T)
                    - a["projected_values"][:, None]
                )
            )
            / np.sqrt(e)[:, None]
        )
        ** 2,
        axis=0,
    )
    close(chi, calc, "native chi2", rtol=3e-5, atol=2e-5)
    eligible = np.array(conv) & (chi >= 0) & np.all(trials >= 0, axis=1)
    require(eligible.any(), "convergence")
    score = np.full(33, np.inf)
    score[eligible] = np.abs(
        np.log(np.maximum(chi[eligible], np.finfo(float).tiny) / len(e))
    )
    idx = int(np.argmin(score))
    require(d["selected_index"] == idx, "alpha selection")
    raw = np.trapezoid(trials[idx], w)
    require(raw > 0, "raw mass")
    close(den, trials[idx] / raw, "selected spectrum")
    prediction = forward(a["query"], w, den)
    close(prediction, a["prediction"], "independent quadrature", rtol=2e-6, atol=2e-7)
    rms = float(
        np.linalg.norm((qv.T @ (prediction - qm)) / np.sqrt(qe)) / np.sqrt(len(qe))
    )
    require(rms <= 2, "heldout whitened RMS")
    close(d["whitened_rms"], rms, "reported RMS", rtol=2e-4, atol=2e-4)
    close(d["window_weights"], weights(w, den), "window weights")
    require(
        a["bootstrap_spectra"].shape == ((8, 401) if count else (0,))
        and a["bootstrap_indices"].shape == ((8, 32) if count else (0,)),
        "bootstrap shape",
    )
    if count:
        ix = np.random.default_rng(d["bootstrap_seed"]).integers(0, 32, size=(8, 32))
        require(np.array_equal(ix, a["bootstrap_indices"]), "bootstrap indices")
        bw = []
        for x in a["bootstrap_spectra"]:
            require(
                np.min(x) >= -1e-6 and abs(np.trapezoid(x, w) - 1) <= 1e-6,
                "bootstrap spectrum",
            )
            bw.append(weights(w, x))
        close(a["bootstrap_weights"], bw, "bootstrap weights")
        close(
            d["bootstrap_window_standard_deviation"],
            np.std(bw, axis=0, ddof=1),
            "bootstrap variation",
        )
        require(
            len(d["bootstrap_records"]) == 8
            and all(
                x["converged"][x["selected_index"]] for x in d["bootstrap_records"]
            ),
            "bootstrap convergence",
        )
    return d, a


def check_snapshot(snapshot):
    require(
        isinstance(snapshot, dict)
        and snapshot.get("revision")
        and snapshot.get("sources")
        and snapshot.get("provider_modules")
        and snapshot.get("packages", {}).get("triqs_maxent") == "4.0.1",
        "runtime snapshot",
    )
    require(
        all(v == "1" for v in snapshot.get("threads", {}).values())
        and len(snapshot.get("threads", {})) >= 5,
        "thread snapshot",
    )
    rev = snapshot["revision"]
    require(
        subprocess.run(
            ["git", "merge-base", "--is-ancestor", rev, "HEAD"],
            cwd=ROOT,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).returncode
        == 0,
        "snapshot revision",
    )
    tracked = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", rev, "kernelgecp"], cwd=ROOT, text=True
    ).splitlines()
    required = {p for p in tracked if p.endswith(".py")}
    required |= {
        "experiments/cthyb_spectroscopy.py",
        "experiments/audit_cthyb_spectroscopy.py",
        "experiments/cthyb_dmft.py",
        "experiments/spectroscopy_campaign.py",
        "experiments/audit_spectroscopy.py",
        "experiments/spectroscopy.py",
        "experiments/dlr_adapters.py",
        "environments/locks/triqs-linux-64-conda-explicit.txt",
        "environments/locks/triqs-maxent-pip.txt",
        "environments/locks/dlr-physics-pip.txt",
    }
    require(
        required.issubset(snapshot["sources"]), "incomplete scientific source inventory"
    )
    for name, h in snapshot["sources"].items():
        require(
            not Path(name).is_absolute() and ".." not in Path(name).parts, "source path"
        )
        blob = subprocess.check_output(["git", "show", rev + ":" + name], cwd=ROOT)
        require(hashlib.sha256(blob).hexdigest() == h, "source hash " + name)
    require(
        "experiments/cthyb_spectroscopy.py" in snapshot["sources"]
        and "experiments/audit_cthyb_spectroscopy.py" in snapshot["sources"],
        "source inventory",
    )


def check_workers(
    manifest, campaign, root, ee_audit, binding, stage, calibration, calibration_audit
):
    """Verify the complete pinned child-job and process receipt inventory."""
    workers = manifest.get("workers")
    require(
        isinstance(workers, dict) and set(workers) == set(manifest["cases"]),
        "worker inventory",
    )
    require(
        manifest.get("global_timeout_seconds") == 19800
        and manifest.get("case_timeout_seconds") == 3600,
        "worker deadline controls",
    )
    elapsed_total = manifest.get("elapsed_seconds")
    require(
        type(elapsed_total) in (int, float)
        and np.isfinite(elapsed_total)
        and 0 < elapsed_total <= 19800,
        "campaign duration",
    )
    seen = set()
    parent_pid = manifest.get("parent_pid")
    require(type(parent_pid) is int and parent_pid > 0, "campaign parent PID")
    campaign_dir = Path(campaign).resolve().parent
    for name, record in workers.items():
        method, model, repeat = name[:-5].rsplit("-", 2)
        require(
            record.get("status") == "passed" and record.get("exit_code") == 0,
            "worker exit",
        )
        pid = record.get("pid")
        require(type(pid) is int and pid > 0 and pid not in seen, "worker PID")
        seen.add(pid)
        elapsed = record.get("elapsed_seconds")
        require(
            type(elapsed) in (int, float)
            and np.isfinite(elapsed)
            and 0 < elapsed <= 3600,
            "worker duration",
        )

        def local(key, current_record=record):
            p = Path(current_record[key])
            p = p if p.is_absolute() else campaign_dir / p
            p = p.resolve()
            require(p.is_relative_to(campaign_dir) and p.is_file(), "worker file path")
            return p

        job_path, log_path, worker_path = (
            local(k) for k in ("job_file", "log_file", "worker_receipt_file")
        )
        require(
            sha(job_path) == record.get("job_sha256")
            and sha(log_path) == record.get("log_sha256")
            and sha(worker_path) == record.get("worker_receipt_sha256"),
            "worker artifact hash",
        )
        job = read(job_path)
        parent = job.get("parent_pid")
        require(
            type(parent) is int and parent > 0 and parent != pid, "worker parent PID"
        )
        require(parent == parent_pid, "worker parent mismatch")
        require(
            job.get("schema") == "kernelgecp.cthyb-spectroscopy-worker-job.v1"
            and job.get("stage") == stage
            and job.get("method") == method
            and job.get("model") == model
            and job.get("repeat") == int(repeat[1:])
            and job.get("protocol") == manifest["protocol"]
            and job.get("snapshot") == manifest["snapshot"],
            "worker job controls",
        )
        require(
            Path(job["folder"]).resolve() == Path(root).resolve()
            and Path(job["audit_receipt"]).resolve() == Path(ee_audit).resolve()
            and Path(job["audit_binding"]).resolve() == Path(binding).resolve()
            and job.get("binding_sha256") == BINDING
            and Path(job["output"]).resolve() == campaign_dir / name,
            "worker input paths",
        )
        if stage == "final":
            require(
                Path(job["calibration"]).resolve() == Path(calibration).resolve()
                and job.get("calibration_sha256") == sha(calibration)
                and Path(job["calibration_audit"]).resolve()
                == Path(calibration_audit).resolve()
                and job.get("calibration_audit_sha256") == sha(calibration_audit),
                "worker calibration binding",
            )
        else:
            require(
                all(
                    job.get(k) is None
                    for k in (
                        "calibration",
                        "calibration_sha256",
                        "calibration_audit",
                        "calibration_audit_sha256",
                    )
                ),
                "calibration worker controls",
            )
        command = record.get("command")
        require(
            isinstance(command, list)
            and len(command) == 7
            and Path(command[0]).resolve()
            == Path(manifest["snapshot"]["executable"]).resolve()
            and command[1:4] == ["-m", "experiments.cthyb_spectroscopy", "--worker"]
            and Path(command[4]).resolve() == job_path
            and command[5:] == ["--job-sha256", sha(job_path)],
            "worker command",
        )
        receipt = read(worker_path)
        require(
            receipt.get("schema") == "kernelgecp.cthyb-spectroscopy-worker.v1"
            and receipt.get("status") == "passed"
            and receipt.get("pid") == pid
            and receipt.get("job_sha256") == sha(job_path)
            and receipt.get("result_sha256") == manifest["cases"][name]["result_sha256"]
            and receipt.get("array_sha256") == manifest["cases"][name]["array_sha256"],
            "worker receipt",
        )
    require(len(seen) == len(workers), "worker PID inventory")


def campaign_audit(
    campaign, root, ee_audit, binding, calibration=None, calibration_audit=None
):
    campaign = Path(campaign)
    m = read(campaign)
    stage = (
        "calibration"
        if m.get("schema") == "kernelgecp.cthyb-spectroscopy-calibration.v1"
        else "final"
    )
    require(
        m.get("schema") == f"kernelgecp.cthyb-spectroscopy-{stage}.v1"
        and m.get("status") == "passed"
        and m.get("binding_sha256") == BINDING,
        "campaign schema",
    )
    protocol = m["protocol"]
    require(
        protocol
        == {
            "methods": list(METHODS),
            "default_models": list(MODELS),
            "eps": 1e-8,
            "n_omega": 401,
            "n_alpha": 33,
            "covariance_rtol": 1e-4,
            "bootstrap_count": 8,
            "bootstrap_repeat": 0,
            "bootstrap_seed": 193,
            "calibration_repeats": 1,
            "final_repeats": 3,
            "validation": "independent_heldout_chain_mean_conditional_on_training_fit",
        },
        "protocol",
    )
    check_snapshot(m["snapshot"])
    if stage == "final":
        require(
            calibration is not None and calibration_audit is not None,
            "independent calibration audit required",
        )
        ca = read(calibration_audit)
        cal = read(calibration)
        require(
            m.get("calibration_sha256") == sha(calibration)
            and ca.get("schema") == "kernelgecp.cthyb-spectroscopy-independent-audit.v1"
            and ca.get("status") == "passed"
            and ca.get("kind") == "calibration"
            and ca.get("calibration_sha256") == sha(calibration)
            and ca.get("binding_sha256") == BINDING
            and ca.get("gates")
            and all(v is True for v in ca["gates"].values()),
            "calibration audit binding",
        )
        require(
            cal.get("protocol") == protocol and cal.get("snapshot") == m["snapshot"],
            "calibration freeze",
        )
    samples = raw_cohort(root, ee_audit, binding, stage)
    repeats = 1 if stage == "calibration" else 3
    expected = {
        f"{method}-{model}-r{repeat}.json"
        for method in METHODS
        for model in MODELS
        for repeat in range(repeats)
    }
    require(set(m["cases"]) == expected, "child inventory")
    check_workers(
        m, campaign, root, ee_audit, binding, stage, calibration, calibration_audit
    )
    records = {}
    arrays = {}
    for name in sorted(expected):
        method, model, r = name[:-5].rsplit("-", 2)
        repeat = int(r[1:])
        p = campaign.parent / name
        h = m["cases"][name]
        require(
            sha(p) == h["result_sha256"]
            and sha(p.with_suffix(".npz")) == h["array_sha256"],
            "child hashes",
        )
        records[name], arrays[name] = audit_case(p, method, model, repeat, samples)
    if stage == "final":
        for method in METHODS:
            for model in MODELS:
                first = arrays[f"{method}-{model}-r0.json"]
                for rep in (1, 2):
                    other = arrays[f"{method}-{model}-r{rep}.json"]
                    for k in first:
                        if not k.startswith("bootstrap_"):
                            require(
                                np.array_equal(first[k], other[k]),
                                "repeated scientific arrays " + k,
                            )
    sensitivity = {}
    for repeat in range(repeats):
        for method in METHODS:
            a = records[f"{method}-flat-r{repeat}.json"]
            b = records[f"{method}-broad_gaussian-r{repeat}.json"]
            delta = np.array(b["window_weights"]) - a["window_weights"]
            key = f"{method}-r{repeat}"
            close(
                m["sensitivity"][key]["broad_minus_flat_window_weights"],
                delta,
                "default sensitivity",
            )
            sensitivity[key] = delta.tolist()
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    source = ROOT / "experiments/audit_cthyb_spectroscopy.py"
    blob = subprocess.check_output(
        ["git", "show", revision + ":experiments/audit_cthyb_spectroscopy.py"], cwd=ROOT
    )
    require(
        sha(source) == hashlib.sha256(blob).hexdigest(),
        "auditor source must be committed",
    )
    return {
        "schema": "kernelgecp.cthyb-spectroscopy-independent-audit.v1",
        "status": "passed",
        "kind": stage,
        "calibration_sha256": sha(campaign)
        if stage == "calibration"
        else sha(calibration),
        "campaign_sha256": sha(campaign),
        "binding_sha256": BINDING,
        "auditor_revision": revision,
        "auditor_sha256": sha(source),
        "gates": {
            "ee_source_and_raw_chains": True,
            "frozen_protocol_and_hashes": True,
            "covariance_and_projection": True,
            "native_trials_and_quadrature": True,
            "heldout_rms": True,
            "bootstrap_and_defaults": True,
            "repeat_arrays": True,
            "fresh_process_receipts": True,
        },
        "cases": len(records),
        "default_sensitivity": sensitivity,
        "scope": (
            "heldout chain mean only; finite native model, no truth CDF, "
            "nonlinear coverage or unique spectrum claim"
        ),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("campaign", type=Path)
    p.add_argument("--ee-root", type=Path, required=True)
    p.add_argument("--ee-audit", type=Path, required=True)
    p.add_argument("--ee-binding", type=Path, required=True)
    p.add_argument("--calibration", type=Path)
    p.add_argument("--calibration-audit", type=Path)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    try:
        r = campaign_audit(
            a.campaign,
            a.ee_root,
            a.ee_audit,
            a.ee_binding,
            a.calibration,
            a.calibration_audit,
        )
        code = 0
    except Exception as e:
        r = {
            "schema": "kernelgecp.cthyb-spectroscopy-independent-audit.v1",
            "status": "failed",
            "error": f"{type(e).__name__}: {e}",
        }
        code = 1
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(r, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
