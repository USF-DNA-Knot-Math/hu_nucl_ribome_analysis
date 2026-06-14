#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# make_fig4BC.sh
#
# Minimal pipeline to reproduce Fig 4B-C: per-nucleotide (rA/rC/rG/rU) rNMP
# enrichment factor (EF) and nucleotide-composition (perc_avg, norm_perc,
# counts_sum) trends vs. expression, overlaid on one plot per data type, for
# HEK293T-WT on the TSS 0-1 kb downstream region across all three strand
# relationships (both/same/opp).
#
# Edit the paths below before running.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"

# ── Inputs ───────────────────────────────────────────────────────────────────
# Per-nucleotide withnuc files: bedtools annotate output (rNMP counts per
# library) merged with bedtools nuc output (nucleotide composition), for the
# TSS 0-1 kb downstream region. See rNMP_EF/README.md for how these are made.
withnuc_dir="/path/to/rNMP_EF/withnuc_outputs/hg38_refGene_TSS_downstream_0_1kb_1000bp"
rN_withnuc="${withnuc_dir}/rN_annotated_counts_withnuc.tsv"

# Per-nucleotide counts files (library filename -> total rNMP count)
counts_dir="/path/to/rNMP_EF/annotate_outputs/hg38_refGene_TSS_downstream_0_1kb_1000bp"

# RNA-seq TPM file(s) (outer-joined on Geneid by merge_exp_all.py)
tpm_files=(
  "/path/to/RNA-seq/WT-KO.tpms_by_condition.tsv"
)

out_dir="${script_dir}/output_fig4BC"
celltype_dir="${out_dir}/celltype_info_outputs"

# ── Stage 1: per-nucleotide composition + counts (get_celltype_info.R) ────────
mkdir -p "$celltype_dir"
declare -A nuc_letter=( [rA]=A [rC]=C [rG]=G [rU]=U )
declare -A nuc_bg=( [rA]=7 [rU]=7 [rC]=8 [rG]=8 )  # background composition column

for nuc in rA rC rG rU; do
  letter="${nuc_letter[$nuc]}"
  sed "s/NUC/${letter}/g" "${script_dir}/libmeta_template.tsv" > "${celltype_dir}/libmeta_${nuc}.tsv"

  Rscript "${script_dir}/get_celltype_info.R" \
    -a "${withnuc_dir}/${nuc}_annotated_counts_withnuc.tsv" \
    -t "$rN_withnuc" \
    -c "${counts_dir}/${nuc}/out/all_counts.tsv" \
    -b "${nuc_bg[$nuc]}" \
    -f "${celltype_dir}/libmeta_${nuc}.tsv" \
    -o "${celltype_dir}/${nuc}_hg38_refGene_TSS_downstream_0_1kb"
done

# ── Stage 2: merge RNA-seq expression into EF + celltype_info matrices ────────
# Job 2 = celltype_info TSS (strand-aware): perc_avg/counts_sum/norm_perc
# Job 3 = EF TSS (strand-aware)
python3 "${script_dir}/merge_exp_all.py" --jobs 2 3 \
  --tpm "${tpm_files[@]}" \
  --celltype-dir "$celltype_dir" \
  --out-root "$out_dir"

ef_merged="${out_dir}/merged_exp_outputs"
ci_merged="${out_dir}/merged_exp_outputs"
ef_dataset="hg38_refGene_TSS_downstream_0_1kb"
ci_dataset="hg38_refGene_TSS_downstream_0_1kb"

# ── Column layout of the merged tables (see merge_exp_all.py output) ──────────
# EF:           cols 7-16 = EF cell types (col 9 = HEK293T-WT)
#               cols 17+  = RNA-seq TPM (col 19 = HEK293T-WT TPM with the
#                           WT-KO TPM file above)
# celltype_info: cols 16-23 = rNMP cell types (col 18 = HEK293T-WT)
#               cols 24+   = RNA-seq TPM (col 26 = HEK293T-WT TPM)
# Check the header of a *_withexp.tsv file and adjust these if your TPM
# file set differs.
ef_xcol=19;  ef_ycol=9
ci_xcol=26;  ci_ycol=18

colors=( [rA]="#D55E00" [rC]="#0072B2" [rG]="#F0C742" [rU]="#009E73" )

# data type -> bins ymax | trends ymax | tick interval
declare -A dtype_params=(
  [EF]="6|6|1"
  [perc_avg]="55|55|10"
  [norm_perc]="5|2|0.5"
  [counts_sum]="4|4|1"
)

pair_bins_out="${out_dir}/pair_bins_outputs"

for strand in both same opp; do
  for dtype in EF perc_avg norm_perc counts_sum; do
    IFS='|' read -r bins_ymax trends_ymax tick <<< "${dtype_params[$dtype]}"
    dtype_out="${pair_bins_out}/${dtype}_${strand}"
    mkdir -p "$dtype_out"

    for nuc in rA rC rG rU; do
      if [[ "$dtype" == "EF" ]]; then
        input="${ef_merged}/${strand}/${nuc}_${ef_dataset}_regions_EF_${strand}_withexp.tsv"
        xcol="$ef_xcol"; ycol="$ef_ycol"
      else
        input="${ci_merged}/${strand}/${nuc}_${ci_dataset}_regions_${dtype}_${strand}_withexp.tsv"
        xcol="$ci_xcol"; ycol="$ci_ycol"
      fi
      run_name="${nuc}_${ef_dataset}_${dtype}_${strand}_HEK293T-WT"

      Rscript "${script_dir}/pair_bins.R" \
        -m "$input" \
        -n "$run_name" \
        -c "${colors[$nuc]}" \
        -p \
        -e "$xcol" \
        -r "$ycol" \
        -y "$bins_ymax" \
        -o "$dtype_out"
    done

    # ── Stage 3: overlay rA/rC/rG/rU trends on one plot per data type ─────────
    trends_dir="${dtype_out}/trends"
    mkdir -p "$trends_dir"
    pattern="${ef_dataset}_${dtype}_${strand}_HEK293T-WT"

    Rscript "${script_dir}/trends_merge.R" \
      -i "$dtype_out" \
      -f "$pattern" \
      -y "$trends_ymax" \
      -b "$tick" \
      -o "$trends_dir"

    Rscript "${script_dir}/trends_merge.R" \
      -i "$dtype_out" \
      -f "$pattern" \
      -y "$trends_ymax" \
      -b "$tick" \
      -o "$trends_dir" \
      -q
  done
done

echo "Done. Outputs under: ${out_dir}/"
