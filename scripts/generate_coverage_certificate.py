"""Reproduce the integer data checked by the Phase DB Lean certificate.

Emit Lean source to stdout, or use --check to compare the committed data.
This generator is not trusted by Lean: the proof checks its arithmetic.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
from fractions import Fraction
from pathlib import Path

ROWS = (7, 11, 12, 18, 19, 21, 22, 23, 27)
COLUMNS = (8, 9, 11, 13, 24, 33, 39, 57, 60)
FEATURES = 13
SIZE = 9


def _vector(values: list[int] | tuple[int, ...]) -> str:
    return "![" + ", ".join(str(value) for value in values) + "]"


def _matrix(name: str, values: list[list[int]]) -> str:
    return (
        f"def {name} : Matrix (Fin {len(values)}) (Fin {len(values[0])}) ℤ :=\n"
        "  Matrix.of ![\n    "
        + ",\n    ".join(_vector(line) for line in values)
        + "\n  ]\n"
    )


def generate_certificate() -> tuple[str, dict[str, object]]:
    degree = FEATURES - 1
    kernel_scale = (29 * 63) ** degree * math.factorial(degree)
    left = [
        [
            x**k
            * (29 * 63) ** (degree - k)
            * (math.factorial(degree) // math.factorial(k))
            for k in range(FEATURES)
        ]
        for x in ROWS
    ]
    right = [[y**k for y in COLUMNS] for k in range(FEATURES)]
    current = [
        [sum(left[i][k] * right[k][j] for k in range(FEATURES)) for j in range(SIZE)]
        for i in range(SIZE)
    ]
    previous = 1
    states = [current]
    divisors = [previous]
    selected_rows: list[int] = []
    selected_columns: list[int] = []
    row_projector = [[Fraction(i == j) for j in range(SIZE)] for i in range(SIZE)]
    column_projector = [[Fraction(i == j) for j in range(SIZE)] for i in range(SIZE)]
    for _ in range(7):
        _, _, _, row, column = max(
            (abs(current[i][j]), -i, -j, i, j)
            for i in range(SIZE)
            for j in range(SIZE)
            if i not in selected_rows and j not in selected_columns
        )
        pivot = current[row][column]
        assert pivot != 0
        row_projector = [
            [
                value - Fraction(current[i][column], pivot) * row_projector[row][j]
                for j, value in enumerate(line)
            ]
            for i, line in enumerate(row_projector)
        ]
        column_projector = [
            [
                value - column_projector[i][column] * Fraction(current[row][j], pivot)
                for j, value in enumerate(line)
            ]
            for i, line in enumerate(column_projector)
        ]
        numerators = [
            [
                pivot * current[i][j] - current[i][column] * current[row][j]
                for j in range(SIZE)
            ]
            for i in range(SIZE)
        ]
        assert all(value % previous == 0 for line in numerators for value in line)
        current = [[value // previous for value in line] for line in numerators]
        previous = pivot
        selected_rows.append(row)
        selected_columns.append(column)
        states.append(current)
        divisors.append(previous)

    # The reversed fresh-row orientation makes the determinant coefficients positive.
    fresh_rows, fresh_columns = (6, 4), (1, 2)
    scale = divisors[-1]
    el = [[row_projector[i][j] * scale for j in range(SIZE)] for i in fresh_rows]
    gr = [[column_projector[i][j] * scale for j in fresh_columns] for i in range(SIZE)]
    ln = [
        [sum(el[i][j] * left[j][k] for j in range(SIZE)) for k in range(FEATURES)]
        for i in range(2)
    ]
    rn = [
        [sum(right[k][j] * gr[j][i] for j in range(SIZE)) for i in range(2)]
        for k in range(FEATURES)
    ]
    assert all(value.denominator == 1 for line in el + gr + ln + rn for value in line)

    source = """import Mathlib.LinearAlgebra.Matrix.Notation

/-! Integer data for a coherent 13-term Taylor coverage obstruction.
Every table is checked by the companion proof module, not trusted as an axiom.
The common kernel scale is `(29*63)^12 * 12!`. -/

namespace GECPKernelStructure.Fermionic

"""
    source += "def phaseDBRowNumerators : Fin 9 → ℤ := " + _vector(ROWS) + "\n"
    source += "def phaseDBColumnNumerators : Fin 9 → ℤ := " + _vector(COLUMNS) + "\n"
    source += "def phaseDBPivotRows : Fin 7 → Fin 9 := " + _vector(selected_rows) + "\n"
    source += (
        "def phaseDBPivotColumns : Fin 7 → Fin 9 := " + _vector(selected_columns) + "\n"
    )
    source += "def phaseDBDivisors : Fin 8 → ℤ := " + _vector(divisors) + "\n"
    source += f"def phaseDBKernelScale : ℤ := {kernel_scale}\n\n"
    for step, state in enumerate(states):
        source += _matrix(f"phaseDBState{step}", state) + "\n"
    source += "def phaseDBState : Fin 8 → Matrix (Fin 9) (Fin 9) ℤ :=\n  !["
    source += ", ".join(f"phaseDBState{k}" for k in range(8)) + "]\n\n"
    for name, values in (
        ("phaseDBLeftProjectorNumerators", el),
        ("phaseDBRightProjectorNumerators", gr),
        ("phaseDBLeftNumerators", ln),
        ("phaseDBRightNumerators", rn),
    ):
        source += (
            _matrix(name, [[int(value) for value in line] for line in values]) + "\n"
        )
    source += "end GECPKernelStructure.Fermionic\n"
    return source, {
        "arithmetic": "exact integer and rational",
        "precision": "unbounded integer and rational arithmetic",
        "conditioning": "all pivots checked nonzero exactly",
        "pivot_rule": "complete pivot, lexicographically smallest exact ties",
        "near_ties": "exact comparisons, no tolerance",
        "selected_rows": selected_rows,
        "selected_columns": selected_columns,
        "features": FEATURES,
        "feature_scale_decimal_digits": len(str(abs(scale))),
        "certification_scope": (
            "fixed finite Taylor witness; Lean checks generated data"
        ),
        "stop_reason": "completed",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    source, report = generate_certificate()
    if args.check:
        root = Path(__file__).resolve().parents[1]
        target = (
            root / "GECPKernelStructure/Fermionic/AdaptiveCoverageObstructionData.lean"
        )
        if target.read_text() != source:
            raise SystemExit("Committed Lean data differs from the exact generator")
        report["revision"] = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True
        ).strip()
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(source, end="")


if __name__ == "__main__":
    main()
