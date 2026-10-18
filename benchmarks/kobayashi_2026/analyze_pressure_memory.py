#!/usr/bin/env python3
"""Analyze the second-order pressure memory from supplied bulk correlations."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from benchmarks.kobayashi_2026.pressure_memory import (
    pressure_memory_summary,
    second_order_memory_laplace,
)


RATES_PS = (0.0, 0.05, 0.10, 0.20, 0.50, 1.0, 2.0, 5.0, 10.0, 20.0)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("data_dir", type=Path)
    p.add_argument("--outdir", type=Path, default=Path("pressure_memory_results"))
    p.add_argument("--curvature-window-ps", type=float, default=0.04)
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    scan_rows = []

    for model in ("s_1.00", "s_0.00"):
        path = args.data_dir / model / "bulk_corr.dat"
        data = np.loadtxt(path, dtype=float)
        window_s = args.curvature_window_ps * 1.0e-12
        summary = pressure_memory_summary(path, curvature_window_s=window_s)
        summary_rows.append({"model": model, **summary})

        omega0_sq = summary["omega0_sq_s-2"]
        for rate_ps in RATES_PS:
            nu = second_order_memory_laplace(
                data[:, 0],
                data[:, 1],
                rate_ps * 1.0e12,
                omega0_sq_s2=omega0_sq,
            )
            scan_rows.append(
                {
                    "model": model,
                    "laplace_rate_ps-1": rate_ps,
                    "memory_laplace_ps-1": nu / 1.0e12,
                }
            )

    for name, rows in (
        ("pressure_memory_summary.csv", summary_rows),
        ("pressure_memory_laplace.csv", scan_rows),
    ):
        with (args.outdir / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
