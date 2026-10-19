"""Exact finite audits of adjacent Taylor-cutoff transfer certificates."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
import subprocess
from fractions import Fraction

RationalMatrix = list[list[Fraction]]


def determinant_two(matrix: RationalMatrix) -> Fraction:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def conditioned_factors(
    rows: tuple[Fraction, ...],
    columns: tuple[Fraction, ...],
    features: int,
    depth: int,
    pivot_terms: int,
) -> tuple[RationalMatrix, RationalMatrix, RationalMatrix, list[int], list[int]]:
    matrix = [
        [
            sum(
                (
                    (x * y) ** degree / math.factorial(degree)
                    for degree in range(pivot_terms)
                ),
                Fraction(),
            )
            for y in columns
        ]
        for x in rows
    ]
    left = [
        [x**degree / math.factorial(degree) for degree in range(features)] for x in rows
    ]
    right = [[y**degree for y in columns] for degree in range(features)]
    selected_rows: list[int] = []
    selected_columns: list[int] = []
    for _ in range(depth):
        _, _, _, row, column = max(
            (abs(matrix[i][j]), -i, -j, i, j)
            for i in range(len(rows))
            for j in range(len(columns))
            if i not in selected_rows and j not in selected_columns
        )
        pivot = matrix[row][column]
        if not pivot:
            raise ArithmeticError("Exact pivot kernel is rank deficient")
        left = [
            [
                value - matrix[i][column] * left[row][k] / pivot
                for k, value in enumerate(line)
            ]
            for i, line in enumerate(left)
        ]
        right = [
            [
                value - right[k][column] * matrix[row][j] / pivot
                for j, value in enumerate(line)
            ]
            for k, line in enumerate(right)
        ]
        matrix = [
            [
                value - matrix[i][column] * matrix[row][j] / pivot
                for j, value in enumerate(line)
            ]
            for i, line in enumerate(matrix)
        ]
        selected_rows.append(row)
        selected_columns.append(column)
    return left, right, matrix, selected_rows, selected_columns


def matched_coordinates(
    left: RationalMatrix,
    right: RationalMatrix,
    fresh_rows: tuple[int, int],
    fresh_columns: tuple[int, int],
    orientation: int,
) -> dict[tuple[int, int], Fraction]:
    return {
        (i, j): orientation
        * (
            left[fresh_rows[0]][i] * left[fresh_rows[1]][j]
            - left[fresh_rows[0]][j] * left[fresh_rows[1]][i]
        )
        * (
            right[i][fresh_columns[0]] * right[j][fresh_columns[1]]
            - right[i][fresh_columns[1]] * right[j][fresh_columns[0]]
        )
        for i, j in itertools.combinations(range(len(right)), 2)
    }


def cutoff_coefficients(
    coordinates: dict[tuple[int, int], Fraction], features: int, cutoff: int
) -> tuple[Fraction, Fraction, Fraction]:
    return (
        sum(
            (coordinates[i, j] for i, j in itertools.combinations(range(cutoff), 2)),
            Fraction(),
        ),
        sum(
            (coordinates[i, j] for i in range(cutoff) for j in range(cutoff, features)),
            Fraction(),
        ),
        sum(
            (
                coordinates[i, j]
                for i, j in itertools.combinations(range(cutoff, features), 2)
            ),
            Fraction(),
        ),
    )


def transfer_numerator(
    coordinates: dict[tuple[int, int], Fraction], features: int, moved: int
) -> tuple[Fraction, Fraction, Fraction]:
    constant, mixed, quadratic = cutoff_coefficients(coordinates, features, moved)
    left = sum((coordinates[i, moved] for i in range(moved)), Fraction())
    right = sum((coordinates[moved, j] for j in range(moved + 1, features)), Fraction())
    numerator = constant * right**2 - mixed * left * right + quadratic * left**2
    return numerator, left, right


def audit(cases: int = 100, seed: int = 0xC212C0DE) -> dict[str, object]:
    if cases < 1:
        raise ValueError("cases must be positive")
    rng = random.Random(seed)
    transfers = 0
    covered = 0
    failures: list[dict[str, object]] = []
    for case in range(cases):
        size = rng.randrange(4, 10)
        features = rng.randrange(5, 15)
        depth = rng.randrange(1, size - 1)
        row_denominator = rng.randrange(size + 1, 31)
        column_denominator = rng.randrange(size + 1, 33)
        rows = tuple(
            Fraction(v, row_denominator)
            for v in sorted(rng.sample(range(1, row_denominator), size))
        )
        columns = tuple(
            Fraction(v, column_denominator)
            for v in sorted(rng.sample(range(1, column_denominator), size))
        )
        left, right, _, selected_rows, selected_columns = conditioned_factors(
            rows, columns, features, depth, features + 5
        )
        remaining_rows = [i for i in range(size) if i not in selected_rows]
        remaining_columns = [j for j in range(size) if j not in selected_columns]
        for row_pair in itertools.combinations(remaining_rows, 2):
            for column_pair in itertools.combinations(remaining_columns, 2):
                inversions = sum(i > j for i in selected_rows for j in row_pair) + sum(
                    i > j for i in selected_columns for j in column_pair
                )
                orientation = (-1) ** inversions
                coordinates = matched_coordinates(
                    left, right, row_pair, column_pair, orientation
                )
                nodes = {}
                for moved in range(1, features - 1):
                    constant, mixed, quadratic = cutoff_coefficients(
                        coordinates, features, moved
                    )
                    numerator, lsum, rsum = transfer_numerator(
                        coordinates, features, moved
                    )
                    if min(constant, mixed, quadratic) < 0:
                        raise AssertionError("Negative sampled coefficient")
                    if lsum <= 0 or rsum <= 0:
                        raise AssertionError("Nonpositive sampled transfer sum")
                    if mixed**2 < 4 * constant * quadratic:
                        raise AssertionError("Negative sampled discriminant")
                    nodes[moved] = numerator
                    transfers += 1
                    if numerator > 0:
                        failures.append(
                            {
                                "case": case,
                                "rows": list(map(str, rows)),
                                "columns": list(map(str, columns)),
                                "features": features,
                                "pivot_terms": features + 5,
                                "selected_rows": selected_rows,
                                "selected_columns": selected_columns,
                                "fresh_rows": row_pair,
                                "fresh_columns": column_pair,
                                "moved_degree": moved,
                            }
                        )
                for cutoff in range(2, features - 1):
                    if nodes[cutoff - 1] > 0 and nodes[cutoff] > 0:
                        raise AssertionError("Sampled adjacent-node coverage failed")
                    covered += 1
    return {
        "arithmetic": "exact rational",
        "precision": "unbounded rational arithmetic",
        "conditioning": "nonzero pivots checked exactly; no floating-point estimate",
        "cases": cases,
        "seed": seed,
        "pivot_rule": "complete pivot with lexicographically smallest ties",
        "pivot_certification": "exact finite-sample comparisons",
        "near_ties": "exact comparisons only; no tolerance is used",
        "transfers": transfers,
        "interior_cutoffs_covered": covered,
        "positive_node_numerators": len(failures),
        "failed_fixed_nodes": failures,
        "coverage_holds_on_samples": True,
        "certification_scope": (
            "finite samples only; no universal or infinite-tail claim"
        ),
        "stop_reason": "completed",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=int, default=100)
    parser.add_argument("--seed", type=int, default=0xC212C0DE)
    args = parser.parse_args()
    report = audit(args.cases, args.seed)
    report["revision"] = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], text=True
    ).strip()
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
