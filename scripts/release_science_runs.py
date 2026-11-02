"""Select authenticated workflow metadata, never scientific eligibility."""

import re

REPOSITORY = "tripp-smith/gecp-kernel-structure"
SCOPE = (
    "run/artifact metadata only, numerical eligibility requires "
    "independent payload audits"
)
TRIQS_STEPS = [
    "Execute actual TRIQS, CT-HYB and MaxEnt calculations",
    "Test actual installed TRIQS interoperability outside the checkout",
    "Run the frozen interacting three-provider campaign and independent audit",
]
ROLES = {
    "providers": (
        "provider-runtimes.yml",
        {
            "triqs (ubuntu-24.04, linux-64)": TRIQS_STEPS,
            "triqs (macos-15, osx-arm64)": TRIQS_STEPS,
            "gp (ubuntu-24.04)": [
                "Execute real operator solves and Gaussian probability checks"
            ],
            "gp (macos-15)": [
                "Execute real operator solves and Gaussian probability checks"
            ],
            "native (green)": [
                "Execute native numerical calculation and validate receipt"
            ],
            "native (epw)": [
                "Execute native numerical calculation and validate receipt"
            ],
            "payload-audit": [
                "Audit provider payloads and independently replay ED campaigns"
            ],
        },
        [
            "provider-triqs-ubuntu-24.04",
            "provider-triqs-macos-15",
            "provider-gp-ubuntu-24.04",
            "provider-gp-macos-15",
            "provider-green-linux",
            "provider-epw-linux",
            "provider-payload-audit",
        ],
    ),
    "molecular": (
        "molecular-green.yml",
        {
            "actual-native": [
                "Execute both molecular geometries on both native grids",
                "Check actual HDF5 interoperability and independent FCI references",
                "Generate the frozen scientific comparison",
                "Independently audit the comparison",
            ]
        },
        ["molecular-green-native-linux"],
    ),
    "epw": (
        "epw-restart-archive.yml",
        {
            "native": [
                "Generate the pinned native 20 K source and provenance",
                "Run the fixed pure-model holdout and explicit correction archive",
                "Independently audit all three archive cases",
                "Compare actual unmodified native raw and restored restarts",
            ]
        },
        ["epw-explicit-restart-archive-linux"],
    ),
    "covariance": (
        "covariance.yml",
        {
            "real-data": [
                "Execute and independently audit a small real-data case",
                "Execute the full frozen EI branch or explicitly dispatched campaign",
            ]
        },
        ["real-data-covariance"],
    ),
    "cthyb": (
        "cthyb-evidence.yml",
        {"replay": ["Recompute independent audit from committed native data"]},
        ["cthyb-evidence-replay-{id}"],
    ),
    "spectroscopy": (
        "spectroscopy-evidence.yml",
        {"replay": ["Replay independent scientific audits"]},
        ["spectroscopy-evidence-{id}-{attempt}"],
    ),
}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _positive(value):
    return type(value) is int and value > 0


def select_science_runs(repository, revision, run_items, jobs_by_run, artifacts_by_run):
    """Select exact-main runs from authenticated, fully paginated API data.

    Runs are keyed by role; jobs and artifacts are keyed by integer run ID.
    Caller authenticates API responses. Artifact contents are not examined.
    """
    _require(repository == REPOSITORY, "unexpected repository")
    _require(
        isinstance(revision, str) and re.fullmatch(r"[0-9a-f]{40}", revision),
        "invalid revision",
    )
    selected = {}
    for role, (workflow, required_jobs, templates) in ROLES.items():
        candidates = [
            r
            for r in run_items.get(role, [])
            if r.get("head_sha") == revision
            and r.get("head_branch") == "main"
            and r.get("event") in ("push", "workflow_dispatch")
            and r.get("path") == f".github/workflows/{workflow}"
            and (r.get("repository") or {}).get("full_name") == repository
            and (r.get("head_repository") or {}).get("full_name") == repository
            and _positive(r.get("id"))
            and _positive(r.get("run_attempt"))
        ]
        _require(candidates, f"missing exact-main {role} run")
        run = max(candidates, key=lambda r: (r["id"], r["run_attempt"]))
        run_id, attempt = run["id"], run["run_attempt"]
        _require(
            run.get("status") == "completed" and run.get("conclusion") == "success",
            f"latest {role} run failed",
        )
        jobs = jobs_by_run.get(run_id, [])
        _require(
            len(jobs) == len(required_jobs)
            and {j.get("name") for j in jobs} == set(required_jobs),
            f"{role} jobs incomplete or duplicated",
        )
        job_records = []
        for job in jobs:
            _require(
                _positive(job.get("id"))
                and job.get("run_id") == run_id
                and job.get("run_attempt") == attempt
                and job.get("head_sha") == revision
                and job.get("status") == "completed"
                and job.get("conclusion") == "success",
                f"invalid {role} job",
            )
            for name in required_jobs[job["name"]]:
                found = [s for s in job.get("steps", []) if s.get("name") == name]
                _require(
                    len(found) == 1
                    and found[0].get("status") == "completed"
                    and found[0].get("conclusion") == "success",
                    f"{role} required step failed or skipped: {name}",
                )
            job_records.append({"id": job["id"], "name": job["name"]})
        expected = {name.format(id=run_id, attempt=attempt) for name in templates}
        artifacts = artifacts_by_run.get(run_id, [])
        names = [a.get("name") for a in artifacts]
        ids = [a.get("id") for a in artifacts]
        auxiliary = [name for name in names if name not in expected]
        limit = {"providers": 2, "molecular": 1, "epw": 1}.get(role, 0)
        _require(
            len(names) == len(set(names))
            and len(ids) == len(set(ids))
            and expected.issubset(names)
            and len(auxiliary) <= limit
            and all(
                isinstance(name, str)
                and re.fullmatch(
                    r"tripp-smith~gecp-kernel-structure~[A-Za-z0-9_-]+\.dockerbuild",
                    name,
                )
                for name in auxiliary
            ),
            f"{role} artifact inventory mismatch",
        )
        records = []
        auxiliary_records = []
        for artifact in artifacts:
            binding = artifact.get("workflow_run") or {}
            digest = artifact.get("digest")
            artifact_attempt = binding.get("run_attempt")
            _require(
                ("run_attempt" not in binding and attempt == 1)
                or (type(artifact_attempt) is int and artifact_attempt == attempt),
                f"ambiguous or invalid {role} artifact run attempt",
            )
            _require(
                _positive(artifact.get("id"))
                and artifact.get("expired") is False
                and _positive(artifact.get("size_in_bytes"))
                and isinstance(digest, str)
                and re.fullmatch(r"sha256:[0-9a-f]{64}", digest)
                and binding.get("id") == run_id
                and binding.get("head_sha") == revision
                and binding.get("head_branch") == "main",
                f"invalid {role} artifact",
            )
            record = {k: artifact[k] for k in ("id", "name", "digest", "size_in_bytes")}
            (records if artifact["name"] in expected else auxiliary_records).append(
                record
            )
        selected[role] = {
            "workflow": workflow,
            "id": run_id,
            "attempt": attempt,
            "jobs": sorted(job_records, key=lambda x: x["name"]),
            "artifacts": sorted(records, key=lambda x: x["name"]),
            "auxiliary_artifacts": sorted(auxiliary_records, key=lambda x: x["name"]),
        }
    return {
        "scope": SCOPE,
        "repository": repository,
        "revision": revision,
        "runs": selected,
    }
