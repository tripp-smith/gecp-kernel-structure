"""An outside source copy or editable package cannot impersonate a wheel."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.providers import installed_wheel

FILES = {"basis": "kernelgecp/basis.py", "green": "kernelgecp/integrations/green.py"}


class Distribution:
    version = "1.1.0"

    def __init__(self):
        self.files = list(FILES.values())
        self.record = "\n".join(path + ",sha256=example,100" for path in FILES.values())
        self.direct = None
        self.prefix = Path("/opt/genuine-wheel")

    def read_text(self, name):
        return self.record if name == "RECORD" else self.direct

    def locate_file(self, name):
        return self.prefix / name


@pytest.fixture
def distribution(monkeypatch):
    value = Distribution()
    monkeypatch.setattr(installed_wheel.metadata, "distribution", lambda _: value)
    return value


def paths(distribution):
    return {
        name: distribution.locate_file(relative) for name, relative in FILES.items()
    }


def test_actual_distribution_membership_and_record_are_required(distribution):
    receipt = installed_wheel.installed_wheel_identity(paths(distribution), FILES)
    assert receipt["module_count"] == 2
    assert receipt["editable"] is False
    assert len(receipt["record_sha256"]) == 64
    assert "url" not in receipt and "prefix" not in receipt


def test_copy_outside_checkout_is_not_enough(distribution):
    copied = {
        name: Path("/tmp/copied-source") / relative for name, relative in FILES.items()
    }
    with pytest.raises(RuntimeError, match="does not belong"):
        installed_wheel.installed_wheel_identity(copied, FILES)


@pytest.mark.parametrize(
    "attack",
    [
        "missing_files",
        "missing_record",
        "bad_record",
        "duplicate",
        "unlisted",
        "editable",
        "bad_direct",
        "wrong_editable_type",
    ],
)
def test_missing_malformed_or_editable_metadata_is_rejected(distribution, attack):
    if attack == "missing_files":
        distribution.files = None
    elif attack == "missing_record":
        distribution.record = None
    elif attack == "bad_record":
        distribution.record = "one,two"
    elif attack == "duplicate":
        distribution.record += "\nkernelgecp/basis.py,hash,1"
    elif attack == "unlisted":
        distribution.files.pop()
    elif attack == "editable":
        distribution.direct = json.dumps({"dir_info": {"editable": True}})
    elif attack == "bad_direct":
        distribution.direct = "invalid json"
    elif attack == "wrong_editable_type":
        distribution.direct = json.dumps({"dir_info": {"editable": 1}})
    with pytest.raises(RuntimeError):
        installed_wheel.installed_wheel_identity(paths(distribution), FILES)
