#!/usr/bin/env python3
"""Summarize supplied Kobayashi pressure-relaxation correlations."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from benchmarks.kobayashi_2026.pressure_relaxation import relaxation_summary


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("data_dir", type=Path)
    p.add_argument("--out", type=Path, default=Path("pressure_relaxation_summary.csv"))
    args = p.parse_args()

    rows = []
    for model in ("s_1.00", "s_0.00"):
        rows.append({
            "model": model,
            **relaxation_summary(args.data_dir / model / "bulk_corr.dat"),
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
