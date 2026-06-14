#!/usr/bin/env python3
"""
Unified merge script that joins RNA-seq TPM expression data into rNMP output
matrices. Replaces the four separate scripts:
  - merge_exp.py            (celltype_info, TSS, no strand, perc_avg/counts_sum/norm_perc)
  - merge_exp_perc.py       (celltype_info, TSS, strand-aware, perc_avg)
  - merge_exp_EF.py         (EF_avg, TSS, strand-aware, EF)
  - merge_exp_EF_genebody.py(EF_avg, gene body, strand-aware, EF)

Join path:
  RNA-seq (Ensembl ID) -> hg38_gene_mapping.tsv (Ensembl->gene name) -> rNMP data (id = gene name)

Usage:
  python3 merge_exp_all.py                          # run all jobs, default TPM files
  python3 merge_exp_all.py --jobs 4 5               # run only jobs 4 and 5
  python3 merge_exp_all.py --tpm /path/to/file.tsv  # use a specific TPM file
  python3 merge_exp_all.py --tpm f1.tsv f2.tsv      # merge multiple TPM files
  python3 merge_exp_all.py --out-root /path/to/dir  # redirect all outputs under dir
  python3 merge_exp_all.py --ef-dir /path/to/ef1 /path/to/ef2  # ad-hoc EF regions
  python3 merge_exp_all.py --list                    # list available jobs
"""

import argparse
import pandas as pd
import os
import sys

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE        = "/storage/home/hcoda1/8/twarner31/r-fstorici3-0"
GENE_MAP_F  = f"{BASE}/RNA-seq-pipeline/scripts/hg38_gene_mapping.tsv"
FC_DIR      = f"{BASE}/G4_analysis/feature_correlations"

# Default TPM files (used when --tpm is not specified)
DEFAULT_TPM_FILES = [
    f"{BASE}/RNA-seq-pipeline/output/WT-KO_03-05-26/WT-KO.tpms_by_condition.tsv",
    f"{BASE}/RNA-seq-pipeline/output/Run_02-19-26/siTOP.tpms_by_condition.tsv",
    f"{BASE}/RNA-seq-pipeline/output/CD4T_04-07-26_hisat2/CD4T.tpms_by_condition.tsv",
    f"{BASE}/RNA-seq-pipeline/output/H9_04-07-26_hisat2/H9.tpms_by_condition.tsv",
]

CELLTYPE_DIR = f"{FC_DIR}/celltype_info_outputs"
EF_TSS_DIR      = f"{BASE}/psoralem_20kbupstream-check/rNMP_EF/annotate_outputs/hg38_refGene_TSS_downstream_0_1kb_1000bp"
EF_GB_DIR       = f"{BASE}/psoralem_20kbupstream-check/rNMP_EF/annotate_outputs/hg38_refGene_collapsed"
EF_TSS_FLANK_DIR = f"{BASE}/psoralem_20kbupstream-check/rNMP_EF/annotate_outputs/hg38_refGene_TSS_flanks_1kb_each_side_2kb_total"

TSS_DATASET       = "hg38_refGene_TSS_downstream_0_1kb"
GB_DATASET        = "hg38_refGene_collapsed"
TSS_FLANK_DATASET = "hg38_refGene_TSS_flanks_1kb_each_side_2kb_total"

NUCS_4 = ["rA", "rC", "rG", "rU"]
NUCS_5 = ["rA", "rC", "rG", "rN", "rU"]

STRANDS_ALL = ["both", "same", "opp"]
# Strand infix for celltype_info filenames (no infix for "both")
STRAND_INFIX = {"both": "", "same": "_same", "opp": "_opp"}

META_COLS = ("Geneid", "Chr", "Start", "End", "Strand", "Length", "gene_name")

# ── Job definitions ───────────────────────────────────────────────────────────
# Each job is a dict with:
#   name        : human-readable label
#   source      : "celltype_info" or "ef_avg"
#   source_dir  : directory containing input files
#   dataset     : dataset name string
#   nucs        : list of nucleotides
#   dtypes      : list of data type suffixes (celltype_info) or ["EF_avg"] (EF)
#   strands     : list of strand types, or [None] for no strand dimension
#   out_base    : output base directory
#   ko_avg_keep : if True, average KO clones AND keep individual clone columns

JOBS = [
    {
        "name": "celltype_info TSS (no strand): perc_avg, counts_sum, norm_perc",
        "source": "celltype_info",
        "source_dir": CELLTYPE_DIR,
        "dataset": TSS_DATASET,
        "nucs": NUCS_4,
        "dtypes": ["perc_avg", "counts_sum", "norm_perc"],
        "strands": [None],
        "out_base": f"{FC_DIR}/merged_exp_outputs",
        "ko_avg_keep": False,
    },
    {
        "name": "celltype_info TSS (strand-aware): perc_avg",
        "source": "celltype_info",
        "source_dir": CELLTYPE_DIR,
        "dataset": TSS_DATASET,
        "nucs": NUCS_4,
        "dtypes": ["perc_avg", "counts_sum", "norm_perc"],
        "strands": STRANDS_ALL,
        "out_base": f"{FC_DIR}/merged_exp_outputs",
        "ko_avg_keep": False,
    },
    {
        "name": "EF TSS (strand-aware)",
        "source": "ef_avg",
        "source_dir": EF_TSS_DIR,
        "dataset": TSS_DATASET,
        "nucs": NUCS_5,
        "dtypes": ["EF_avg"],
        "strands": STRANDS_ALL,
        "out_base": f"{FC_DIR}/merged_exp_outputs",
        "ko_avg_keep": True,
    },
    {
        "name": "EF gene body (strand-aware)",
        "source": "ef_avg",
        "source_dir": EF_GB_DIR,
        "dataset": GB_DATASET,
        "nucs": NUCS_5,
        "dtypes": ["EF_avg"],
        "strands": STRANDS_ALL,
        "out_base": f"{FC_DIR}/merged_exp_outputs_genebody",
        "ko_avg_keep": True,
    },
    {
        "name": "EF TSS flanks 2kb (strand-aware)",
        "source": "ef_avg",
        "source_dir": EF_TSS_FLANK_DIR,
        "dataset": TSS_FLANK_DATASET,
        "nucs": NUCS_5,
        "dtypes": ["EF_avg"],
        "strands": STRANDS_ALL,
        "out_base": f"{FC_DIR}/merged_exp_outputs_tss_flanks",
        "ko_avg_keep": True,
    },
    {
        "name": "celltype_info TSS upstream 2-3kb (strand-aware): perc_avg, counts_sum, norm_perc",
        "source": "celltype_info",
        "source_dir": CELLTYPE_DIR,
        "dataset": "hg38_refGene_TSS_upstream_2_3kb",
        "nucs": NUCS_4,
        "dtypes": ["perc_avg", "counts_sum", "norm_perc"],
        "strands": STRANDS_ALL,
        "out_base": f"{FC_DIR}/merged_exp_outputs",
        "ko_avg_keep": False,
    },
    {
        "name": "celltype_info TTS downstream 4-5kb (strand-aware): perc_avg, counts_sum, norm_perc",
        "source": "celltype_info",
        "source_dir": CELLTYPE_DIR,
        "dataset": "hg38_refGene_TTS_downstream_4_5kb",
        "nucs": NUCS_4,
        "dtypes": ["perc_avg", "counts_sum", "norm_perc"],
        "strands": STRANDS_ALL,
        "out_base": f"{FC_DIR}/merged_exp_outputs",
        "ko_avg_keep": False,
    },
]


def load_expression(gene_map_f, tpm_files):
    """Load gene mapping and one or more TPM files, returning a join-ready DataFrame.

    Multiple TPM files are outer-joined on Geneid so all expression columns
    are available in a single table. Genes present in any TPM file are kept;
    expression values default to NaN where a gene is absent from a file.
    """
    gene_map = pd.read_csv(gene_map_f, sep="\t", header=None,
                           usecols=[0, 1], names=["Geneid", "gene_name"])

    combined = None
    for tpm_f in tpm_files:
        tpm = pd.read_csv(tpm_f, sep="\t")
        tpm = tpm.rename(columns=lambda c: c.replace("+", "p"))

        # Extract expression columns (everything except featureCounts metadata)
        exp_cols = [c for c in tpm.columns if c not in META_COLS and c != "Geneid"]
        subset = tpm[["Geneid"] + exp_cols]

        if combined is None:
            combined = subset
        else:
            combined = combined.merge(subset, on="Geneid", how="outer")

    # Add gene names, drop rows without a mapping or on chrX/Y
    combined = combined.merge(gene_map, on="Geneid", how="left")
    combined = combined.dropna(subset=["gene_name"])

    # Need Chr info to filter X/Y — get it from the first TPM file
    chr_info = pd.read_csv(tpm_files[0], sep="\t", usecols=["Geneid", "Chr"])
    combined = combined.merge(chr_info, on="Geneid", how="left")
    combined = combined[~combined["Chr"].astype(str).isin(["X", "Y"])]
    combined = combined.drop(columns=["Geneid", "Chr"])

    return combined.drop_duplicates(subset="gene_name")


def build_input_path(job, nuc, strand, dtype):
    """Construct the input file path based on source type."""
    if job["source"] == "celltype_info":
        infix = STRAND_INFIX.get(strand, "") if strand else ""
        return os.path.join(
            job["source_dir"],
            f"{nuc}_{job['dataset']}{infix}_regions_{dtype}.tsv"
        )
    else:  # ef_avg
        ds = job["dataset"]
        # EF source directory has an extra suffix for TSS downstream
        if ds == TSS_DATASET:
            fname = f"{ds}_1000bp_{nuc}_{strand}_EF_regions_EF_avg.tsv"
        else:
            fname = f"{ds}_{nuc}_{strand}_EF_regions_EF_avg.tsv"
        return os.path.join(job["source_dir"], nuc, "out", fname)


def build_output_path(job, nuc, strand, dtype):
    """Construct the output file path."""
    ds = job["dataset"]
    if job["source"] == "celltype_info":
        if strand:
            return f"{nuc}_{ds}_regions_{dtype}_{strand}_withexp.tsv"
        else:
            return f"{nuc}_{ds}_regions_{dtype}_withexp.tsv"
    else:  # ef_avg
        return f"{nuc}_{ds}_regions_EF_{strand}_withexp.tsv"


def process_job(job, tpm_join):
    """Run a single merge job."""
    print(f"\n{'='*70}")
    print(f"Job: {job['name']}")
    print(f"{'='*70}")

    for strand in job["strands"]:
        if strand:
            out_dir = os.path.join(job["out_base"], strand)
        else:
            out_dir = job["out_base"]
        os.makedirs(out_dir, exist_ok=True)

        if strand:
            print(f"\n  --- strand: {strand} ---")

        for dtype in job["dtypes"]:
            for nuc in job["nucs"]:
                in_f = build_input_path(job, nuc, strand, dtype)
                df = pd.read_csv(in_f, sep="\t")

                # Replace + with p in column names
                df = df.rename(columns=lambda c: c.replace("+", "p"))

                # For EF sources: average KO clones and keep individual columns
                if job["ko_avg_keep"]:
                    ko_t38  = "HEK293T-RNASEH2A-KO-T3-8"
                    ko_t317 = "HEK293T-RNASEH2A-KO-T3-17"
                    df.insert(df.columns.get_loc(ko_t38), "HEK293T-KO",
                              df[[ko_t38, ko_t317]].mean(axis=1))
                    df = df.rename(columns={ko_t38: "KO-T3-8", ko_t317: "KO-T3-17"})

                # Rename any expression columns that collide with data columns
                data_cols = set(df.columns) - {"id"}
                exp_rename = {c: f"{c}_exp" for c in tpm_join.columns
                              if c in data_cols and c != "gene_name"}
                tpm_to_merge = tpm_join.rename(columns=exp_rename) if exp_rename else tpm_join
                if exp_rename:
                    print(f"      Renamed expression cols to avoid collision: {exp_rename}")

                # Merge with expression data
                merged = df.merge(tpm_to_merge, left_on="id", right_on="gene_name", how="inner")
                merged = merged.drop(columns=["gene_name"])
                merged = merged[~merged["chr"].isin(["chrX", "chrY"])]

                out_fname = build_output_path(job, nuc, strand, dtype)
                out_f = os.path.join(out_dir, out_fname)
                merged.to_csv(out_f, sep="\t", index=False)
                label = f"{nuc} {dtype}" if job["source"] == "celltype_info" else nuc
                print(f"    {label}: {len(merged)} genes -> {out_f}")

    # Print column map for one sample file
    sample_strand = job["strands"][0]
    sample_dtype = job["dtypes"][0]
    sample_nuc = job["nucs"][0]
    sample_fname = build_output_path(job, sample_nuc, sample_strand, sample_dtype)
    if sample_strand:
        sample_dir = os.path.join(job["out_base"], sample_strand)
    else:
        sample_dir = job["out_base"]
    sample_f = os.path.join(sample_dir, sample_fname)
    if os.path.exists(sample_f):
        cols = pd.read_csv(sample_f, sep="\t", nrows=0).columns.tolist()
        print(f"\n  Column positions in {os.path.basename(sample_f)} (1-indexed):")
        for i, c in enumerate(cols, 1):
            print(f"    {i:2d}  {c}")


def main():
    parser = argparse.ArgumentParser(
        description="Unified rNMP + expression merge script")
    parser.add_argument("--jobs", nargs="+", type=int,
                        help="Job numbers to run (1-indexed). Default: all.")
    parser.add_argument("--tpm", nargs="+", metavar="FILE",
                        help="TPM file(s) to merge. Default: Ribome + CD4T + H9.")
    parser.add_argument("--out-root", metavar="DIR",
                        help="Redirect all job outputs under this directory. "
                             "Each job's out_base dirname is remapped relative to this root.")
    parser.add_argument("--ef-dir", nargs="+", metavar="DIR",
                        help="Ad-hoc EF region directories. Each creates a strand-aware "
                             "EF merge job. Dataset name = directory basename. "
                             "Ignores --jobs; can be combined with --out-root.")
    parser.add_argument("--celltype-dir", metavar="DIR",
                        help="Override the celltype_info source directory for all "
                             "celltype_info jobs (default: celltype_info_outputs/).")
    parser.add_argument("--list", action="store_true",
                        help="List available jobs and exit.")
    args = parser.parse_args()

    if args.list:
        print("Available jobs:")
        for i, job in enumerate(JOBS, 1):
            print(f"  {i}: {job['name']}")
        sys.exit(0)

    tpm_files = args.tpm if args.tpm else DEFAULT_TPM_FILES

    # Load expression data once
    print("Loading expression data...")
    for f in tpm_files:
        print(f"  TPM: {f}")
    tpm_join = load_expression(GENE_MAP_F, tpm_files)
    exp_cols = [c for c in tpm_join.columns if c != "gene_name"]
    print(f"  {len(tpm_join)} unique genes, {len(exp_cols)} expression columns: {exp_cols}")

    # Build job list
    if args.ef_dir:
        # Ad-hoc EF jobs from directory paths
        out_root = os.path.abspath(args.out_root) if args.out_root else FC_DIR
        selected = []
        for ef_path in args.ef_dir:
            ef_path = os.path.abspath(ef_path)
            dataset = os.path.basename(ef_path)
            selected.append({
                "name": f"EF {dataset} (strand-aware)",
                "source": "ef_avg",
                "source_dir": ef_path,
                "dataset": dataset,
                "nucs": NUCS_5,
                "dtypes": ["EF_avg"],
                "strands": STRANDS_ALL,
                "out_base": os.path.join(out_root, f"merged_exp_{dataset}"),
                "ko_avg_keep": True,
            })
    else:
        # Select from predefined jobs
        if args.jobs:
            selected = [JOBS[j - 1] for j in args.jobs]
        else:
            selected = JOBS

        # Override output directories if --out-root is given
        if args.out_root:
            out_root = os.path.abspath(args.out_root)
            for job in selected:
                # Keep only the last component of the original out_base path
                job["out_base"] = os.path.join(out_root, os.path.basename(job["out_base"]))

        # Override celltype_info source directory if requested
        if args.celltype_dir:
            celltype_dir = os.path.abspath(args.celltype_dir)
            for job in selected:
                if job.get("source") == "celltype_info":
                    job["source_dir"] = celltype_dir

    for job in selected:
        process_job(job, tpm_join)

    print("\nDone.")


if __name__ == "__main__":
    main()
