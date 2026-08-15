from fractions import Fraction

from kernelgecp.surrogates import (
    exact_gecp,
    exact_gecp_pivot_sign_coherence,
    fraction_determinant,
    geometric_surrogate,
    inspect_surrogate,
)


def test_fraction_determinant_and_exact_pivots() -> None:
    matrix = geometric_surrogate(4, Fraction(2, 3))
    assert (
        fraction_determinant([[Fraction(2), Fraction(1)], [Fraction(1), Fraction(3)]])
        == 5
    )
    rows, columns, pivots = exact_gecp(matrix)
    assert len(rows) == len(columns) == len(pivots) == 4
    assert all(pivot != 0 for pivot in pivots)


def test_geometric_minors_have_expected_sign() -> None:
    for size in range(2, 6):
        record = inspect_surrogate(size, Fraction(3, 4))
        assert not any(record.minor_sign_mismatches.values())
        assert not any(record.zero_minors.values())


def test_geometric_surrogate_selected_crosses_are_sign_coherent() -> None:
    for size in range(2, 9):
        for q in (Fraction(1, 2), Fraction(2, 3), Fraction(3, 4)):
            checks = exact_gecp_pivot_sign_coherence(geometric_surrogate(size, q))
            assert len(checks) == size
            assert all(checks)


def test_minimized_sign_regular_pivot_obstruction() -> None:
    left = [[Fraction(4), Fraction(3)], [Fraction(2), Fraction(1)]]
    right = [[Fraction(1), Fraction(2)], [Fraction(3), Fraction(4)]]
    assert fraction_determinant(left) == fraction_determinant(right) == -2
    left_rows, left_columns, _ = exact_gecp(left)
    right_rows, right_columns, _ = exact_gecp(right)
    assert (left_rows[0], left_columns[0]) == (0, 0)
    assert (right_rows[0], right_columns[0]) == (1, 1)
