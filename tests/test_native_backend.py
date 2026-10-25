"""Input, arithmetic-state and parity checks for the optional native backend."""

from __future__ import annotations

import importlib.util
import os
import re
import tracemalloc
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from kernelgecp import __version__
from kernelgecp._basis_numerics import BasisNumerics

ROOT = Path(__file__).resolve().parents[1]
MODULE = os.environ.get("KERNELGECP_NATIVE_MODULE")
if MODULE:
    spec = importlib.util.spec_from_file_location("_native", Path(MODULE))
    assert spec is not None and spec.loader is not None
    native = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(native)
else:
    native = pytest.importorskip("kernelgecp._native")


def args(rank: int = 2) -> list[object]:
    identity = [[str(int(i == j)) for j in range(rank)] for i in range(rank)]
    return [
        128,
        list(map(float, range(rank))),
        list(range(rank)),
        identity,
        [r[:] for r in identity],
    ]


def handle():
    return native.ReferenceCardinal128(*args())


def test_metadata_and_empty_little_endian_output():
    info = native.build_info()
    assert native.NATIVE_PROTOCOL_VERSION == info["protocol_version"] == 1
    assert native.WORKING_PRECISION_BITS == info["working_precision_bits"] == 128
    assert info["mpfr_tls"] is True
    assert info["package_version"] == __version__
    assert isinstance(info["build_target"], str) and info["build_target"]
    assert info["source_revision"] == "unrecorded" or re.fullmatch(
        r"[0-9a-f]{40}", info["source_revision"]
    )
    assert handle().rank == 2
    assert handle().design_tau([]) == b""
    scalar = native.ReferenceCardinal128(128, [0.0], [0], [["1"]], [["1"]])
    assert scalar.design_tau([0, 0.5, 1]) == np.full(3, 0.5, dtype="<f8").tobytes()


@pytest.mark.parametrize(
    "bits", [True, False, 53, 127, 129, 4096, 128.0, "128", -(2**128), 2**100]
)
def test_precision_rejected(bits):
    values = args()
    values[0] = bits
    with pytest.raises((ValueError, TypeError)):
        native.ReferenceCardinal128(*values)


@pytest.mark.parametrize(
    "omega",
    [
        [],
        list(range(129)),
        [True, 1],
        [0.0, -0.0],
        [1.0, 1.0],
        [0.0, float("inf")],
        [0.0, float("nan")],
        [0.0, 16385.0],
        [0.0, -16385.0],
    ],
)
def test_frequency_validation(omega):
    values = args()
    values[1] = omega
    with pytest.raises((ValueError, TypeError)):
        native.ReferenceCardinal128(*values)


@pytest.mark.parametrize(
    "order", [[0], [0, 0], [0, 2], [-1, 0], [True, 1], [0.0, 1], [0, 2**100], "01"]
)
def test_permutation_validation(order):
    values = args()
    values[2] = order
    with pytest.raises((ValueError, TypeError)):
        native.ReferenceCardinal128(*values)


@pytest.mark.parametrize(
    "value", [[], [["1"]], [["1", "0"], ["0"]], "identity", [[1, 0], [0, 1]]]
)
@pytest.mark.parametrize("slot", [3, 4])
def test_shape_and_token_type_validation(value, slot):
    values = args()
    values[slot] = value
    with pytest.raises((ValueError, TypeError)):
        native.ReferenceCardinal128(*values)


@pytest.mark.parametrize(
    "token",
    [
        "nan",
        "inf",
        "-Inf",
        "1_000",
        " ",
        "",
        "e4",
        "1e1e1",
        "1e99999999999999999999999999999",
        "1e65536",
        "1e-65536",
        "9" * 97,
        ".",
        "+",
        "１２",
    ],
)
def test_factor_token_validation(token):
    values = args()
    values[4][0][0] = token
    with pytest.raises((ValueError, TypeError)):
        native.ReferenceCardinal128(*values)


@pytest.mark.parametrize(
    "slot,i,j,value",
    [
        (3, 0, 0, "2"),
        (4, 0, 0, "0"),
        (4, 0, 0, "-0.0"),
        (3, 0, 1, "0.1"),
        (4, 1, 0, "0.1"),
    ],
)
def test_lu_structure_validation(slot, i, j, value):
    values = args()
    values[slot][i][j] = value
    with pytest.raises(ValueError):
        native.ReferenceCardinal128(*values)


@pytest.mark.parametrize(
    "times",
    [
        [True],
        [-1e-20],
        [1.0000000000000002],
        [float("nan")],
        [float("inf")],
        [1j],
        "0.5",
        np.array([0.5]),
        [0.5] * 16385,
    ],
)
def test_query_validation(times):
    with pytest.raises((ValueError, TypeError)):
        handle().design_tau(times)


def test_work_budget_checked_before_evaluation():
    h = native.ReferenceCardinal128(*args(128))
    with pytest.raises(ValueError, match="budget"):
        h.design_tau([0.5] * 8193)


def test_immutable_handle_and_concurrent_readers():
    h = handle()
    with pytest.raises((AttributeError, TypeError)):
        h.rank = 3
    times = np.linspace(0, 1, 257).tolist()
    expected = h.design_tau(times)
    with ThreadPoolExecutor(max_workers=4) as pool:
        assert all(x == expected for x in pool.map(h.design_tau, [times] * 8))


def test_nonfinite_binary64_export_fails():
    h = native.ReferenceCardinal128(128, [0.0], [0], [["1"]], [["1e-1000"]])
    with pytest.raises(RuntimeError, match="export"):
        h.design_tau([0.3])


def test_final_binary64_underflow_has_explicit_existing_export_semantics():
    h = native.ReferenceCardinal128(128, [0.0], [0], [["1"]], [["1e1000"]])
    assert h.design_tau([0.3]) == np.zeros(1, dtype="<f8").tobytes()


def test_caller_mpmath_numpy_contexts_unchanged():
    saved_mp, saved_iv = mp.mp.prec, mp.iv.prec
    saved_np = np.seterr(all="raise")
    try:
        mp.mp.prec, mp.iv.prec = 91, 81
        handle().design_tau([0, 0.123456, 1])
        assert (mp.mp.prec, mp.iv.prec) == (91, 81)
        assert np.geterr() == dict.fromkeys(saved_np, "raise")
    finally:
        mp.mp.prec, mp.iv.prec = saved_mp, saved_iv
        np.seterr(**saved_np)


def test_128_bit_lu_parity_and_independent_256_bit_qr():
    c = mp.mp.clone()
    c.prec = 128
    nodes = [0.0, 0.2, 0.8, 1.0]
    omegas = [-6.0, -1.0, 2.0, 7.0]

    def kernel(context, t, w):
        return (
            context.exp(-t * w) / (1 + context.exp(-w))
            if w >= 0
            else context.exp((1 - t) * w) / (1 + context.exp(w))
        )

    core = c.matrix([[kernel(c, c.mpf(t), c.mpf(w)) for w in omegas] for t in nodes])
    p, lo, up = c.lu(core.T)
    order = [next(j for j in range(4) if p[i, j]) for i in range(4)]
    h = native.ReferenceCardinal128(
        128,
        omegas,
        order,
        [[c.nstr(lo[i, j], 60) for j in range(4)] for i in range(4)],
        [[c.nstr(up[i, j], 60) for j in range(4)] for i in range(4)],
    )
    times = [
        0.0,
        float(np.nextafter(0.0, 1.0)),
        0.1234,
        0.500001,
        float(np.nextafter(1.0, 0.0)),
        1.0,
    ]
    actual = np.frombuffer(h.design_tau(times), dtype="<f8").reshape(-1, 4)
    high = mp.mp.clone()
    high.prec = 256
    hi_core = high.matrix(
        [[kernel(high, high.mpf(t), high.mpf(w)) for w in omegas] for t in nodes]
    )
    reference = BasisNumerics(1.0, 7.0, nodes, omegas, precision_bits=128)
    baseline = reference.design_tau(times)
    assert np.array_equal(actual, baseline)
    for k, t in enumerate(times):
        truth = high.qr_solve(
            hi_core.T,
            high.matrix([kernel(high, high.mpf(t), high.mpf(w)) for w in omegas]),
        )[0]
        for j in range(4):
            assert abs(high.mpf(float(actual[k, j])) - truth[j]) <= high.mpf("2e-15")


def test_mutating_float_hook_cannot_shrink_constructor_coordinates():
    values = args()
    coordinates = values[1]
    calls = []

    class MutatingFloat:
        def __float__(self):
            calls.append("called")
            coordinates.pop()
            return 0.0

    coordinates[0] = MutatingFloat()
    values[2] = [1, 0]
    with pytest.raises(TypeError):
        native.ReferenceCardinal128(*values)
    assert calls == []
    assert len(coordinates) == 2


def test_mutating_float_hook_cannot_shrink_query_coordinates():
    times = [None, 0.5]
    calls = []

    class MutatingFloat:
        def __float__(self):
            calls.append("called")
            times.pop()
            return 0.0

    times[0] = MutatingFloat()
    with pytest.raises(TypeError):
        handle().design_tau(times)
    assert calls == []
    assert len(times) == 2


def test_integer_subclasses_are_rejected_at_precision_and_order_boundaries():
    class IntegerSubclass(int):
        pass

    values = args()
    values[0] = IntegerSubclass(128)
    with pytest.raises((TypeError, ValueError)):
        native.ReferenceCardinal128(*values)
    values = args()
    values[2][0] = IntegerSubclass(0)
    with pytest.raises(TypeError):
        native.ReferenceCardinal128(*values)


@pytest.mark.parametrize("slot", [1, 2, 3, 4])
def test_list_subclasses_are_rejected_at_constructor_boundaries(slot):
    class ListSubclass(list):
        pass

    values = args()
    values[slot] = ListSubclass(values[slot])
    with pytest.raises(TypeError):
        native.ReferenceCardinal128(*values)


@pytest.mark.parametrize("slot", [3, 4])
def test_list_subclasses_are_rejected_for_matrix_rows(slot):
    class ListSubclass(list):
        pass

    values = args()
    values[slot][0] = ListSubclass(values[slot][0])
    with pytest.raises(TypeError):
        native.ReferenceCardinal128(*values)


def test_list_and_real_subclasses_are_rejected_for_query_values():
    class ListSubclass(list):
        pass

    class FloatSubclass(float):
        pass

    class IntegerSubclass(int):
        pass

    h = handle()
    for times in (ListSubclass([0.5]), [FloatSubclass(0.5)], [IntegerSubclass(0)]):
        with pytest.raises(TypeError):
            h.design_tau(times)


def test_string_subclasses_are_rejected_before_factor_conversion():
    class StringSubclass(str):
        pass

    values = args()
    values[3][0][0] = StringSubclass("1")
    with pytest.raises(TypeError):
        native.ReferenceCardinal128(*values)


def test_oversized_unicode_token_rejected_before_large_utf8_allocation():
    if tracemalloc.is_tracing():
        pytest.skip("memory-boundary test requires an isolated tracer")
    token = "é" * 1_000_000
    tracemalloc.start()
    try:
        with pytest.raises(ValueError):
            native.ReferenceCardinal128(128, [0.0], [0], [[token]], [["1"]])
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert peak < 512 * 1024
