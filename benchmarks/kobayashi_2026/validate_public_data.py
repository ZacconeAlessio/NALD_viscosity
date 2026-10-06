#!/usr/bin/env python3
"""Validate the public collaborator-supplied Kobayashi benchmark subset."""

from __future__ import annotations

import csv
from pathlib import Path

DATA = Path(__file__).with_name("data")
SUMMARY = DATA / "benchmark_summary.csv"
SAMPLED = DATA / "bulk_corr_sampled.csv"

EXPECTED = {
    "s_1.00": {
        "rows": 8001,
        "c0": 0.6729e9,
        "zeta": 0.000380497267485,
        "published": 0.000379,
    },
    "s_0.00": {
        "rows": 6001,
        "c0": 0.3908e9,
        "zeta": 0.000035109834144,
        "published": 0.0000344,
    },
}


def close(a: float, b: float, rtol: float = 1.0e-10) -> bool:
    return abs(a - b) <= rtol * max(abs(a), abs(b), 1.0)


def main() -> None:
    with SUMMARY.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if {r["model"] for r in rows} != set(EXPECTED):
        raise SystemExit("FAIL: unexpected models in benchmark_summary.csv")

    for row in rows:
        state = row["model"]
        exp = EXPECTED[state]
        if int(row["raw_corr_rows"]) != exp["rows"]:
            raise SystemExit(f"FAIL: {state} row count")
        if not close(float(row["corr_t0_Pa"]), exp["c0"]):
            raise SystemExit(f"FAIL: {state} C(0)")
        if not close(float(row["integrated_zeta_Pa_s"]), exp["zeta"]):
            raise SystemExit(f"FAIL: {state} integrated zeta")
        pub = float(row["published_zeta_Pa_s"])
        if abs(exp["zeta"] / pub - 1.0) > 0.03:
            raise SystemExit(f"FAIL: {state} supplied integral vs published zeta")

    last = {}
    first = {}
    with SAMPLED.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            state = row["model"]
            first.setdefault(state, row)
            last[state] = row

    for state, exp in EXPECTED.items():
        if not close(float(first[state]["bulk_corr_Pa"]), exp["c0"]):
            raise SystemExit(f"FAIL: {state} sampled C(0)")
        if not close(float(last[state]["cumulative_zeta_Pa_s"]), exp["zeta"]):
            raise SystemExit(f"FAIL: {state} sampled endpoint cumulative zeta")

    print("PASS: public Kobayashi benchmark subset is internally consistent.")
    for state, exp in EXPECTED.items():
        print(f"{state}: zeta={exp['zeta']:.12g} Pa s")


if __name__ == "__main__":
    main()
