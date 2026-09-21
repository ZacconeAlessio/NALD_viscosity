#!/usr/bin/env python3
"""Check zeta=(K_inf-K_0)*tau_zeta for the transcribed benchmark."""

from __future__ import annotations

import csv
from pathlib import Path


def main() -> None:
    path = Path(__file__).with_name("reference.csv")
    worst = 0.0
    checked = 0
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["model"] == "experiment" or not row["Kinf_minus_K0_GPa"]:
                continue
            zeta_table = float(row["zeta_1e-1_mPa_s"])
            kdiff = float(row["Kinf_minus_K0_GPa"])
            tau = float(row["tau_zeta_ps"])
            # GPa * ps = 1e-3 Pa s = 1 mPa s.
            reconstructed_table_units = 10.0 * kdiff * tau
            rel = abs(reconstructed_table_units - zeta_table) / zeta_table
            worst = max(worst, rel)
            checked += 1
    print(f"checked {checked} MD rows")
    print(f"worst relative discrepancy: {worst:.4%}")
    if worst > 0.05:
        raise SystemExit("FAIL: relation differs by more than 5%")
    print("PASS")


if __name__ == "__main__":
    main()
