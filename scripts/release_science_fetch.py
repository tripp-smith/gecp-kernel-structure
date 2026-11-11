"""Collect bounded authenticated workflow metadata, not scientific evidence."""

import argparse
import json
import re
import subprocess
from pathlib import Path

try:
    from scripts.release_science_runs import (
        REPOSITORY,
        profile_roles,
        select_science_runs,
        select_spectroscopy_origin,
    )
except ModuleNotFoundError:
    from release_science_runs import (
        REPOSITORY,
        profile_roles,
        select_science_runs,
        select_spectroscopy_origin,
    )

PAGE_SIZE = 100
MAX_PAGES = 100


def _gh(path, page):
    try:
        result = subprocess.run(
            [
                "gh",
                "api",
                "--method",
                "GET",
                path,
                "-f",
                f"per_page={PAGE_SIZE}",
                "-f",
                f"page={page}",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as exc:
        raise ValueError("GitHub metadata request failed") from exc
    if result.returncode:
        raise ValueError("GitHub metadata request failed")
    try:
        return json.loads(result.stdout)
    except (ValueError, TypeError) as exc:
        raise ValueError("invalid GitHub metadata JSON") from exc


def _pages(path, key, get):
    items, ids, total = [], set(), None
    for page in range(1, MAX_PAGES + 1):
        response = get(path, page)
        if not isinstance(response, dict):
            raise ValueError("invalid metadata response")
        count, batch = response.get("total_count"), response.get(key)
        if type(count) is not int or count < 0 or not isinstance(batch, list):
            raise ValueError("invalid metadata page")
        if total is None:
            total = count
            if total > PAGE_SIZE * MAX_PAGES:
                raise ValueError("metadata page budget exceeded")
        if count != total:
            raise ValueError("metadata total changed during pagination")
        if len(batch) != min(PAGE_SIZE, total - len(items)):
            raise ValueError("missing or excessive metadata page items")
        for item in batch:
            if (
                not isinstance(item, dict)
                or type(item.get("id")) is not int
                or item["id"] <= 0
            ):
                raise ValueError("invalid metadata item ID")
            if item["id"] in ids:
                raise ValueError("duplicate metadata item ID")
            ids.add(item["id"])
        items.extend(batch)
        if len(items) == total:
            return items
    raise ValueError("metadata page budget exceeded")


def fetch_science_runs(repository, revision, get=None, profile=None):
    """Fetch all pages and select the latest exact-main run in each role."""
    if (
        repository != REPOSITORY
        or not isinstance(revision, str)
        or not re.fullmatch(r"[0-9a-f]{40}", revision)
    ):
        raise ValueError("invalid repository or revision")
    get = _gh if get is None else get
    prefix = f"repos/{repository}/actions"
    runs, jobs, artifacts = {}, {}, {}
    for role, (workflow, _, _) in profile_roles(profile).items():
        runs[role] = _pages(
            f"{prefix}/workflows/{workflow}/runs?branch=main&head_sha={revision}",
            "workflow_runs",
            get,
        )
        candidates = [
            r
            for r in runs[role]
            if r.get("head_sha") == revision
            and r.get("head_branch") == "main"
            and r.get("event") in ("push", "workflow_dispatch")
            and r.get("path") == f".github/workflows/{workflow}"
            and isinstance(r.get("repository"), dict)
            and r["repository"].get("full_name") == repository
            and isinstance(r.get("head_repository"), dict)
            and r["head_repository"].get("full_name") == repository
            and type(r.get("run_attempt")) is int
            and r["run_attempt"] > 0
        ]
        if not candidates:
            raise ValueError(f"missing exact-main {role} run")
        run = max(candidates, key=lambda r: (r["id"], r["run_attempt"]))
        if run.get("status") != "completed" or run.get("conclusion") != "success":
            raise ValueError(f"latest {role} run failed")
        run_id, attempt = run["id"], run["run_attempt"]
        if run_id in jobs:
            raise ValueError("duplicate selected run ID")
        jobs[run_id] = _pages(
            f"{prefix}/runs/{run_id}/attempts/{attempt}/jobs", "jobs", get
        )
        artifacts[run_id] = _pages(
            f"{prefix}/runs/{run_id}/artifacts", "artifacts", get
        )
    result = select_science_runs(
        repository, revision, runs, jobs, artifacts, profile=profile
    )
    result["provenance"] = {
        "source": "authenticated gh api GET",
        "page_size": PAGE_SIZE,
        "max_pages_per_endpoint": MAX_PAGES,
        "payloads_downloaded": False,
    }
    return result


def fetch_spectroscopy_origin(repository, revision, run_id, get=None):
    if (
        repository != REPOSITORY
        or not isinstance(revision, str)
        or not re.fullmatch(r"[0-9a-f]{40}", revision)
        or type(run_id) is not int
        or run_id <= 0
    ):
        raise ValueError("invalid origin identity")
    get = _gh if get is None else get
    prefix = f"repos/{repository}/actions/runs/{run_id}"
    run = get(prefix, 1)
    # Query current attempt first: a rerun must not be hidden by historical attempt 1.
    if (
        not isinstance(run, dict)
        or type(run.get("run_attempt")) is not int
        or run["run_attempt"] != 1
    ):
        raise ValueError("canonical origin requires current attempt 1")
    jobs = _pages(prefix + "/attempts/1/jobs", "jobs", get)
    artifacts = _pages(prefix + "/artifacts", "artifacts", get)
    result = select_spectroscopy_origin(
        repository, revision, run_id, run, jobs, artifacts
    )
    result["provenance"] = {
        "source": "authenticated gh api GET",
        "page_size": PAGE_SIZE,
        "max_pages_per_endpoint": MAX_PAGES,
        "payloads_downloaded": False,
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["fetch", "fetch-spectroscopy-origin"])
    parser.add_argument("--repository", required=True)
    parser.add_argument("--revision", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--profile", choices=["production-components"])
    parser.add_argument("--run-id", type=int)
    args = parser.parse_args()
    try:
        if args.action == "fetch":
            if args.run_id is not None:
                raise ValueError("run-id is only for spectroscopy origin")
            result = fetch_science_runs(
                args.repository, args.revision, profile=args.profile
            )
        else:
            if args.profile is not None or args.run_id is None:
                raise ValueError("origin requires run-id and no profile")
            result = fetch_spectroscopy_origin(
                args.repository, args.revision, args.run_id
            )
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    except (ValueError, OSError) as exc:
        parser.exit(1, f"science metadata fetch failed: {exc}\n")


if __name__ == "__main__":
    main()
