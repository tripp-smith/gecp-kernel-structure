"""The executable same-workload example is exercised, not just its protocol."""

import importlib.util

import pytest

from examples.dlr_interoperability import run_demo


def test_public_gecp_physical_interoperability_example():
    record = run_demo(("kernelgecp",))["records"][0]
    assert record["method"] == "kernelgecp"
    assert record["accuracy_passed"]
    assert record["solver_converged"]
    assert record["solver_iterations"] > 1


@pytest.mark.parametrize("method", ["pydlr", "sparse_ir"])
def test_actual_optional_physical_interoperability_example(method):
    if importlib.util.find_spec(method) is None:
        pytest.skip(f"optional {method} is not installed")
    record = run_demo((method,))["records"][0]
    assert record["method"] == method
    assert record["accuracy_passed"]
    assert record["solver_converged"]


@pytest.mark.parametrize("methods", [(), ("unknown",), ("kernelgecp", "kernelgecp")])
def test_example_rejects_invalid_provider_selection(methods):
    with pytest.raises(ValueError):
        run_demo(methods)
