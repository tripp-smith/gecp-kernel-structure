import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from test_release_science_payloads import fake, selection, zip_bytes
from test_release_spectroscopy_origin import origin, select

from scripts import release_science_gate as gate
from scripts import release_science_payloads as payloads


def fixtures(tmp_path):
    blob = zip_bytes([("uci-501.zip", "test")])
    selected = selection(blob)
    selected["profile"] = "production-components"
    del selected["runs"]["spectroscopy"]
    o = select(origin())
    o["artifact"]["size_in_bytes"] = len(blob)
    import hashlib

    o["artifact"]["digest"] = "sha256:" + hashlib.sha256(blob).hexdigest()
    sf, of = tmp_path / "selected.json", tmp_path / "origin.json"
    sf.write_text(json.dumps(selected))
    of.write_text(json.dumps(o))
    return blob, selected, o, sf, of


def test_single_role_validates_all_roles(tmp_path):
    blob, selected, _, _, _ = fixtures(tmp_path)
    result = payloads.materialize(
        selected,
        tmp_path / "payloads",
        tmp_path / "materialized",
        fake(blob),
        role="cthyb",
    )
    assert {r["role"] for r in result["artifacts"]} == {"cthyb"}
    selected["runs"]["epw"]["artifacts"].pop()
    with pytest.raises(ValueError, match="artifact names"):
        payloads.materialize(
            selected,
            tmp_path / "bad",
            tmp_path / "bad-receipt",
            fake(blob),
            role="cthyb",
        )


@pytest.mark.parametrize("role", gate.ROLES)
def test_role_mapping_and_fresh_receipt(tmp_path, monkeypatch, role):
    blob, selected, o, sf, of = fixtures(tmp_path)
    root = tmp_path / "repo"
    root.mkdir()
    monkeypatch.setattr(
        gate.subprocess,
        "check_output",
        lambda args, **kw: selected["revision"] if "rev-parse" in args else "",
    )

    def wrapper(args, **kw):
        assert args[2].endswith(
            "spectroscopy_v4"
            if role == "spectroscopy"
            else "provider_audits"
            if role == "providers"
            else role + "_audits"
        )
        out = Path(args[args.index("--output") + 1])
        r = dict(
            schema=gate.SCHEMAS[role],
            status="passed",
            source_revision=o["revision"]
            if role == "spectroscopy"
            else selected["revision"],
        )
        r.update({k: {"present": True} for k in gate.FIELDS[role]})
        if role == "providers":
            from scripts.release_provider_audits import LANES

            r["lanes"] = {
                lane: dict(
                    receipt_valid=True,
                    receipt_sha256="x",
                    bound_payload_sha256={},
                    independent_dmft="audited.json",
                )
                for lane in LANES
            }
        if role in ("cthyb", "epw"):
            r["payload_unchanged"] = True
        if role == "covariance":
            r["inputs_unchanged"] = True
        if role == "spectroscopy":
            assert (
                Path(args[args.index("--archive") + 1]).name
                == o["artifact"]["name"] + ".zip"
            )
            r.update(
                counts={
                    "synthetic_calibration": 24,
                    "synthetic_final": 72,
                    "measured_calibration": 8,
                    "measured_final": 24,
                    "raw_chains": 256,
                },
                receipt_sha256={
                    name: "a" * 64
                    for name in (
                        "synthetic-cal",
                        "synthetic-final",
                        "measured-cal",
                        "raw",
                        "measured-final",
                    )
                },
                archive_sha256=o["artifact"]["digest"][7:],
                archive_size=len(blob),
                run_id=str(o["id"]),
            )
        if role != "providers":
            out.mkdir()
            out = out / "receipt.json"
        out.write_text(json.dumps(r))
        return SimpleNamespace(returncode=0)

    result = gate.run_role(
        role,
        sf,
        of,
        "/usr/bin/python3",
        tmp_path / "out",
        downloader=fake(blob),
        runner=wrapper,
        root=root,
    )
    assert result["status"] == "passed" and result["role"] == role
    assert result["wrapper_receipt_sha256"] == gate.sha(
        tmp_path / "out/wrapper-receipt.json"
    )
    assert (tmp_path / "out/role-receipt.json").is_file()


def test_failure_no_success_receipt(tmp_path, monkeypatch):
    blob, selected, _, sf, of = fixtures(tmp_path)
    root = tmp_path / "repo"
    root.mkdir()
    monkeypatch.setattr(
        gate.subprocess,
        "check_output",
        lambda args, **kw: selected["revision"] if "rev-parse" in args else "",
    )
    with pytest.raises(ValueError, match="wrapper command failed"):
        gate.run_role(
            "cthyb",
            sf,
            of,
            "/usr/bin/python3",
            tmp_path / "out",
            downloader=fake(blob),
            runner=lambda *a, **k: SimpleNamespace(returncode=1),
            root=root,
        )
    assert not (tmp_path / "out/role-receipt.json").exists()
    assert (tmp_path / "out/failure.log").exists()


def test_rejects_wrong_receipts():
    for role in gate.ROLES:
        with pytest.raises(ValueError):
            gate.validate_wrapper(
                role,
                {
                    "schema": gate.SCHEMAS[role],
                    "status": "passed",
                    "source_revision": "a",
                },
                "a",
            )
    with pytest.raises(ValueError):
        gate.origin_check({"attempt": 2})
