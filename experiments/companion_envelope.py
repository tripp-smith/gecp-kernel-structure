"""Finite-grid stress evidence for the still-open conditioned companion bound.

Pivots always come from the original exponential kernel. Decimal agreement
is not interval certification or a continuous-domain complete-pivot proof.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from fractions import Fraction

import mpmath as mp

WEIGHTS = tuple(range(56, 129, 9))
CASES = (
    (0, 1),
    (2, 16),
    (8, 16),
    (16, 16),
    (32, 16),
    (8, 256),
    (16, 256),
    (32, 256),
    (64, 256),
)


def experiment(depth: int, cutoff: int, digits: int) -> dict[str, object]:
    if depth < 0 or cutoff < 1 or digits < 100:
        raise ValueError("Require depth >= 0, cutoff >= 1, and digits >= 100")
    size = depth + 18
    rows = [Fraction(i + 1, size) for i in range(size)]
    local = [Fraction(j + 1, size) for j in range(size)]
    outer = (
        [Fraction(1) + (cutoff - 1) * Fraction(j + 1, depth) ** 2 for j in range(depth)]
        if cutoff > 1
        else []
    )
    with mp.workdps(digits):
        x = [mp.mpf(v.numerator) / v.denominator for v in rows]
        y = [mp.mpf(v.numerator) / v.denominator for v in local + outer]
        denominator = [1 + mp.exp(v) for v in y]
        matrix = [[mp.exp(a * b) / denominator[j] for j, b in enumerate(y)] for a in x]
        left = [[a**k / mp.factorial(k) for k in range(8)] for a in x]
        right = [[b**k / denominator[j] for j, b in enumerate(y)] for k in range(8)]
        selected_rows: list[int] = []
        selected_columns: list[int] = []
        gaps, pivots = [], []
        for _ in range(depth):
            candidates = sorted(
                (abs(matrix[i][j]), -i, -j, i, j)
                for i in range(size)
                for j in range(len(y))
                if i not in selected_rows and j not in selected_columns
            )
            peak, _, _, row, column = candidates[-1]
            if not peak:
                raise ArithmeticError("Uncertified zero pivot")
            gaps.append((peak - candidates[-2][0]) / peak)
            pivot = matrix[row][column]
            pivots.append(abs(pivot))
            left = [
                [
                    v - matrix[i][column] * left[row][k] / pivot
                    for k, v in enumerate(line)
                ]
                for i, line in enumerate(left)
            ]
            right = [
                [
                    v - right[k][column] * matrix[row][j] / pivot
                    for j, v in enumerate(line)
                ]
                for k, line in enumerate(right)
            ]
            matrix = [
                [
                    v - matrix[i][column] * matrix[row][j] / pivot
                    for j, v in enumerate(line)
                ]
                for i, line in enumerate(matrix)
            ]
            selected_rows.append(row)
            selected_columns.append(column)
        fresh_rows = [i for i in range(size) if i not in selected_rows]
        fresh_columns = [j for j in range(size) if j not in selected_columns]
        sample = [
            fresh_columns[round(k * (len(fresh_columns) - 1) / 17)] for k in range(18)
        ]
        residual = mp.matrix([[matrix[i][j] for j in sample] for i in fresh_rows])
        low = mp.matrix(
            [
                [sum(left[i][k] * right[k][j] for k in range(8)) for j in sample]
                for i in fresh_rows
            ]
        )
        ratios: list[str] = []
        envelope_holds = True
        peak = mp.mpf(0)
        for weight in WEIGHTS:
            companion = weight * residual - (weight - 1) * low
            peak = max(peak, max(abs(v) for v in companion))
            determinant = mp.det(companion)
            if not determinant:
                raise ArithmeticError("Uncertified zero determinant")
            ratio = abs(determinant) / 2**18
            envelope_holds &= ratio <= 1
            ratios.append(mp.nstr(ratio, 45))
        return {
            "depth": depth,
            "cutoff": cutoff,
            "digits": digits,
            "selected_rows": selected_rows,
            "selected_columns": selected_columns,
            "sample_columns": sample,
            "weights": list(WEIGHTS),
            "max_entry": mp.nstr(peak, 30),
            "max_low": mp.nstr(max(abs(v) for v in low), 30),
            "max_residual": mp.nstr(max(abs(v) for v in residual), 30),
            "determinant_ratios": ratios,
            "sampled_envelope_holds": envelope_holds,
            "relative_pivot_gap": mp.nstr(min(gaps), 20) if gaps else None,
            "smallest_pivot": mp.nstr(min(pivots), 20) if pivots else None,
            "conditioning": (
                "Not estimated; inspect precision agreement and pivot diagnostics"
            ),
            "tie_policy": "Largest magnitude, then smallest row and column indices",
            "scope": "Finite-grid numerical evidence only",
            "certification": "Not interval-certified or a theorem",
            "stop_reason": "completed_fixed_depth",
        }


def audit_case(depth: int, cutoff: int, digits: int = 320) -> dict[str, object]:
    lower, upper = (
        experiment(depth, cutoff, precision) for precision in (digits, 2 * digits)
    )
    keys = ("selected_rows", "selected_columns", "sample_columns", "determinant_ratios")
    agreement = all(lower[key] == upper[key] for key in keys)
    envelope_holds = all(run["sampled_envelope_holds"] for run in (lower, upper))
    stop_reason = (
        "uncertified_precision_disagreement"
        if not agreement
        else "precision_agreement"
        if envelope_holds
        else "sampled_envelope_violation"
    )
    return {
        "configuration": {"depth": depth, "cutoff": cutoff},
        "precision_agreement_45_significant_digits": agreement,
        "stop_reason": stop_reason,
        "runs": [lower, upper],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", type=int)
    parser.add_argument("--cutoff", type=int, default=256)
    parser.add_argument("--digits", type=int, default=320)
    args = parser.parse_args()
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], text=True))
    cases = CASES if args.depth is None else ((args.depth, args.cutoff),)
    failed = False
    for depth, cutoff in cases:
        report = audit_case(depth, cutoff, args.digits)
        report["revision"] = revision
        report["working_tree_dirty"] = dirty
        print(json.dumps(report, sort_keys=True), flush=True)
        failed |= report["stop_reason"] != "precision_agreement"
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
