#!/usr/bin/env python3
"""
average_ko_counts.py

Produces a celltype_info output directory suitable for the Fig4B-C pipeline
where HEK293T-KO counts_sum = average of (HEK293T-KO-T3-8 + HEK293T-KO-T3-17).

Inputs:
  --split-dir   celltype_info_split_outputs/ (get_celltype_info run with
                libmeta_ko_split_template — has HEK293T-KO-T3-8 and
                HEK293T-KO-T3-17 columns instead of HEK293T-KO)
  --orig-dir    celltype_info_outputs/ (standard run — used as-is for
                perc_avg and norm_perc files)
  --out-dir     destination directory for the averaged outputs

For *_counts_sum* files: reads split-dir, replaces HEK293T-KO-T3-8 and
  HEK293T-KO-T3-17 with a single HEK293T-KO = (T3-8 + T3-17) / 2 column.
  Column order matches the original counts_sum column layout.

For *_perc_avg* and *_norm_perc* files: copied unchanged from orig-dir.
"""

import argparse
import glob
import os
import shutil

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--split-dir", required=True,
                   help="celltype_info dir generated with the split KO libmeta")
    p.add_argument("--orig-dir", required=True,
                   help="standard celltype_info dir (for perc_avg/norm_perc)")
    p.add_argument("--out-dir", required=True,
                   help="output directory for averaged files")
    return p.parse_args()


def average_counts_sum(split_path: str, out_path: str) -> None:
    df = pd.read_csv(split_path, sep="\t")

    t38_col  = "HEK293T-KO-T3-8"
    t317_col = "HEK293T-KO-T3-17"

    if t38_col not in df.columns or t317_col not in df.columns:
        raise ValueError(
            f"Expected columns '{t38_col}' and '{t317_col}' in {split_path}.\n"
            f"Found: {list(df.columns)}"
        )

    avg = (df[t38_col] + df[t317_col]) / 2.0

    # Insert HEK293T-KO at the position of the first clone column, drop both clones.
    insert_pos = df.columns.get_loc(t38_col)
    df.drop(columns=[t38_col, t317_col], inplace=True)
    df.insert(insert_pos, "HEK293T-KO", avg)

    df.to_csv(out_path, sep="\t", index=False)
    print(f"  averaged → {os.path.basename(out_path)}")


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # ── counts_sum: average from split dir ───────────────────────────────────
    split_files = glob.glob(os.path.join(args.split_dir, "*_counts_sum*.tsv"))
    if not split_files:
        print(f"WARNING: no *_counts_sum* files found in {args.split_dir}")
    for src in split_files:
        fname = os.path.basename(src)
        # Rename: remove any _split suffix if present, keep standard naming
        out_fname = fname.replace("_split", "")
        average_counts_sum(src, os.path.join(args.out_dir, out_fname))

    # ── perc_avg and norm_perc: copy unchanged from orig dir ─────────────────
    for pattern in ("*_perc_avg*.tsv", "*_norm_perc*.tsv"):
        for src in glob.glob(os.path.join(args.orig_dir, pattern)):
            dst = os.path.join(args.out_dir, os.path.basename(src))
            shutil.copy2(src, dst)
            print(f"  copied   → {os.path.basename(dst)}")


if __name__ == "__main__":
    main()
