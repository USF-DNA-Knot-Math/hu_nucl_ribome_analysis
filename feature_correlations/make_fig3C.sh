#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# make_fig3C.sh
#
# Minimal pipeline to reproduce Fig 3C: overlays the rNMP EF (rN, all
# nucleotides combined) vs. expression trend for two genomic regions —
# TSS +/- 1 kb flanks (red) and TSS-to-TTS gene body (blue) — for
# HEK293T-WT, across all three strand relationships (both/same/opp).
#
# Edit the paths below before running.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"

# ── Inputs ───────────────────────────────────────────────────────────────────
# rNMP_EF/annotate_outputs directories for each region (see rNMP_EF/README.md)
ef_dir_genebody="/path/to/rNMP_EF/annotate_outputs/hg38_refGene_collapsed"
ef_dir_tss_flanks="/path/to/rNMP_EF/annotate_outputs/hg38_refGene_TSS_flanks_1kb_each_side_2kb_total"

# RNA-seq TPM file(s) (outer-joined on Geneid by merge_exp_all.py)
tpm_files=(
  "/path/to/RNA-seq/WT-KO.tpms_by_condition.tsv"
)

out_dir="${script_dir}/output_fig3C"

# ── Column layout of the merged EF tables (see merge_exp_all.py output) ───────
# Cols 7-16: EF cell types (col 9 = HEK293T-WT)
# Cols 17+:  RNA-seq TPM columns appended in the order TPM files are joined
#            (with the WT-KO file above, col 19 = HEK293T-WT TPM)
# Check the header of a *_withexp.tsv file and adjust these if your TPM
# file set differs.
ef_ycol=9    # HEK293T-WT EF
ef_xcol=19   # HEK293T-WT RNA-seq TPM
ymax=4
tick=1

# ── Stage 1: merge RNA-seq expression into the EF matrices ────────────────────
# Jobs 4 and 5 = EF gene body and EF TSS flanks (strand-aware)
python3 "${script_dir}/merge_exp_all.py" --jobs 4 5 \
  --tpm "${tpm_files[@]}" \
  --out-root "$out_dir"

genebody_merged="${out_dir}/merged_exp_outputs_genebody"
tss_flanks_merged="${out_dir}/merged_exp_outputs_tss_flanks"
genebody_dataset="hg38_refGene_collapsed"
tss_flanks_dataset="hg38_refGene_TSS_flanks_1kb_each_side_2kb_total"

# ── Stage 2: bin rN EF by expression for each region/strand ────────────────────
bins_genebody="${out_dir}/pair_bins_genebody"
bins_tss="${out_dir}/pair_bins_tss"

for strand in both same opp; do
  for region in genebody tss_flanks; do
    case "$region" in
      genebody)
        merged_dir="$genebody_merged"; dataset="$genebody_dataset"; out="$bins_genebody" ;;
      tss_flanks)
        merged_dir="$tss_flanks_merged"; dataset="$tss_flanks_dataset"; out="$bins_tss" ;;
    esac

    out_strand="${out}/EF_${strand}"
    mkdir -p "$out_strand"
    input="${merged_dir}/${strand}/rN_${dataset}_regions_EF_${strand}_withexp.tsv"
    run_name="rN_${dataset}_EF_${strand}_HEK293T-WT"

    Rscript "${script_dir}/pair_bins.R" \
      -m "$input" \
      -n "$run_name" \
      -c "#777777" \
      -p \
      -e "$ef_xcol" \
      -r "$ef_ycol" \
      -y "$ymax" \
      -o "$out_strand"
  done
done

# ── Stage 3: overlay the two regions on a single trend plot per strand ────────
plots_dir="${out_dir}/two_region_plots"

for strand in both same opp; do
  tss_stats="${bins_tss}/EF_${strand}/rN_${tss_flanks_dataset}_EF_${strand}_HEK293T-WT_stats.tsv"
  gb_stats="${bins_genebody}/EF_${strand}/rN_${genebody_dataset}_EF_${strand}_HEK293T-WT_stats.tsv"
  out_strand="${plots_dir}/EF_${strand}"
  mkdir -p "$out_strand"
  run_name="${tss_flanks_dataset}_EF_${strand}_HEK293T-WT"

  # Labeled version
  Rscript "${script_dir}/two_region_trends.R" \
    -a "$tss_stats" \
    -b "$gb_stats" \
    -n "$run_name" \
    -l "rNMP in TSS +/- 1 kb" \
    -k "rNMP in TSS to TTS" \
    -y "$ymax" \
    -t "$tick" \
    -o "$out_strand"

  # No-label version (for figure assembly)
  Rscript "${script_dir}/two_region_trends.R" \
    -a "$tss_stats" \
    -b "$gb_stats" \
    -n "$run_name" \
    -y "$ymax" \
    -t "$tick" \
    -q \
    -o "$out_strand"
done

echo "Done. Outputs under: ${out_dir}/"
