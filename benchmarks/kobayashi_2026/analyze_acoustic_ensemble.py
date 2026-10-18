#!/usr/bin/env python3
"""Evaluate independent acoustic q_min scales on the Kobayashi snapshots."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from benchmarks.kobayashi_2026.acoustic_finite_size import plane_wave_rayleigh
from benchmarks.kobayashi_2026.static_nald import load_rxyz_frames, source_box_length_m


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("data_dir", type=Path)
    p.add_argument("--out", type=Path, default=Path("acoustic_scan.csv"))
    p.add_argument("--max-harmonic", type=int, default=3)
    args = p.parse_args()

    rows = []
    for model, scale in (("s_1.00", 1.0), ("s_0.00", 0.0)):
        frames = load_rxyz_frames(args.data_dir / model / "rxyz.dat")
        length = source_box_length_m(natoms=len(frames[0][1]))
        for frame_index, (step, xyz) in enumerate(frames):
            for n in range(1, args.max_harmonic + 1):
                rows.append({
                    "model": model,
                    "frame_index": frame_index,
                    "step": step,
                    **plane_wave_rayleigh(xyz, length, scale, n),
                })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
