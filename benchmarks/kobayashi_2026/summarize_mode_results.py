#!/usr/bin/env python3
"""Summarize frame_metrics.csv from the Kobayashi mode ensemble.

This produces the mean structural-prefactor cutoff scan and a cross-model
comparison against the published 124 K Green-Kubo ratios.  It intentionally
does not select a preferred cutoff or claim a fitted viscosity prediction.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


MODELS = ("s_1.00", "s_0.00")
PUBLISHED = {
    "s_1.00": {"zeta_Pa_s": 3.79e-4, "eta_Pa_s": 0.958e-4},
    "s_0.00": {"zeta_Pa_s": 0.344e-4, "eta_Pa_s": 0.790e-4},
}


def _cutoff_columns(fieldnames: list[str]) -> list[tuple[float, str]]:
    result = []
    for name in fieldnames:
        if not name.startswith("A_bulk_") or not name.endswith("_kg_m"):
            continue
        key = name[len("A_bulk_") : -len("_kg_m")]
        cutoff = float(key.replace("p", "."))
        result.append((cutoff, key))
    return sorted(result)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("frame_metrics", type=Path)
    p.add_argument("--outdir", type=Path, default=Path("summary"))
    args = p.parse_args()

    with args.frame_metrics.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        if not rows:
            raise SystemExit("frame_metrics.csv is empty")
        cutoffs = _cutoff_columns(list(rows[0].keys()))

    args.outdir.mkdir(parents=True, exist_ok=True)

    scan = []
    for model in MODELS:
        model_rows = [r for r in rows if r["model"] == model]
        if not model_rows:
            raise SystemExit(f"missing rows for {model}")

        for cutoff, key in cutoffs:
            bulk = np.array([float(r[f"A_bulk_{key}_kg_m"]) for r in model_rows])
            shear = np.array([float(r[f"A_shear_{key}_kg_m"]) for r in model_rows])
            scan.append(
                {
                    "model": model,
                    "cutoff_THz": cutoff,
                    "A_bulk_mean_kg_m": float(bulk.mean()),
                    "A_bulk_std_kg_m": float(bulk.std(ddof=1)),
                    "A_shear_mean_kg_m": float(shear.mean()),
                    "A_shear_std_kg_m": float(shear.std(ddof=1)),
                    "A_bulk_over_A_shear": float(bulk.mean() / shear.mean()),
                    "published_zeta_over_eta": (
                        PUBLISHED[model]["zeta_Pa_s"]
                        / PUBLISHED[model]["eta_Pa_s"]
                    ),
                }
            )

    with (args.outdir / "channel_prefactor_scan.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(scan[0].keys()))
        writer.writeheader()
        writer.writerows(scan)

    cross = []
    for cutoff, _ in cutoffs:
        by_model = {
            model: next(
                row
                for row in scan
                if row["model"] == model and row["cutoff_THz"] == cutoff
            )
            for model in MODELS
        }
        cross.append(
            {
                "cutoff_THz": cutoff,
                "A_bulk_s1_over_s0": (
                    by_model["s_1.00"]["A_bulk_mean_kg_m"]
                    / by_model["s_0.00"]["A_bulk_mean_kg_m"]
                ),
                "published_zeta_s1_over_s0": (
                    PUBLISHED["s_1.00"]["zeta_Pa_s"]
                    / PUBLISHED["s_0.00"]["zeta_Pa_s"]
                ),
                "A_shear_s1_over_s0": (
                    by_model["s_1.00"]["A_shear_mean_kg_m"]
                    / by_model["s_0.00"]["A_shear_mean_kg_m"]
                ),
                "published_eta_s1_over_s0": (
                    PUBLISHED["s_1.00"]["eta_Pa_s"]
                    / PUBLISHED["s_0.00"]["eta_Pa_s"]
                ),
            }
        )

    with (args.outdir / "cross_model_ratios.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cross[0].keys()))
        writer.writeheader()
        writer.writerows(cross)


if __name__ == "__main__":
    main()
