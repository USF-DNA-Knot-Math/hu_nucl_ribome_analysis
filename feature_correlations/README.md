# Feature Correlations

Scripts for correlating rNMP enrichment factor (EF) with RNA-seq expression across genomic regions. These analyses produced **Fig 3C** (two-region rNMP EF vs. expression trend comparison) and **Fig 4B–C** (per-nucleotide rNMP EF and composition vs. expression across multiple regions).

---

## Overview of the pipeline

The pipeline has four conceptual stages:

```
[1] get_celltype_info.R    ← per-nucleotide rNMP counts/perc per cell type
         ↓
[2] average_ko_counts.py   ← average T3-8 and T3-17 KO clones into a single KO column
         ↓
[3] merge_exp_all.py       ← join RNA-seq TPM expression into rNMP matrices
         ↓
[4] pair_bins.R             ← bin genes by expression, compute rNMP EF mean per bin
    trends_merge.R          ← overlay rA/rC/rG/rU trends on a single plot (Fig 4B-C)
    two_region_trends.R     ← overlay two genomic regions on a single trend plot (Fig 3C)
```

EF inputs (from `rNMP_EF/`) are fed directly into Stage 3; they bypass Stage 1–2.  
Stages 1–2 are only needed for the nucleotide-composition data types (`counts_sum`, `perc_avg`, `norm_perc`).

For a runnable end-to-end example of each figure, see [Reproducing Fig 3C and Fig 4B-C](#reproducing-fig-3c-and-fig-4b-c) below.

---

## Scripts

### Stage 1 — `get_celltype_info.R`

**Original author:** Deepali Kundnani  
**Updated by:** TPW — fixed unit mismatch in `norm_perc` normalization: background nucleotide proportion (0–1) is now scaled by 50 before dividing the rNMP percentage (0–100), giving an EF-like quantity near 1 at background composition.

Reads a `bedtools annotate`-style withnuc TSV (rNMP counts + nucleotide composition from `bedtools nuc`) and outputs three files per nucleotide per region:

| Output | Description |
|---|---|
| `*_regions_counts_sum.tsv` | Raw rNMP count summed across replicates per cell type |
| `*_regions_perc_avg.tsv` | rNMP percentage per region, averaged across replicates |
| `*_regions_norm_perc.tsv` | `perc_avg` normalized by background nucleotide composition |

```bash
Rscript get_celltype_info.R \
  -a rA_annotated_counts_withnuc.tsv \   # withnuc annotated file for this nucleotide
  -t rN_annotated_counts_withnuc.tsv \   # total (rN) withnuc file for normalization
  -c all_counts.tsv \                    # per-library total rNMP counts
  -b 7 \                                 # background composition column (e.g. 7=pct_AT, 8=pct_GC)
  -f libmeta.tsv \                       # library metadata (order file)
  -o rA_hg38_refGene_TSS_downstream_0_1kb
```

**libmeta format** (`libmeta_template.tsv` / `libmeta_ko_split_template.tsv`):

| Col | Field | Description |
|---|---|---|
| 1 | num | Sort order |
| 2 | Filename | BED filename (without path) |
| 3 | Library | Library ID |
| 4 | Cellline | Cell type label used for aggregation |
| 5 | Rep | Replicate label |

Use `libmeta_ko_split_template.tsv` when T3-8 and T3-17 KO clones should remain separate (needed before Stage 2).

---

### Stage 2 — `average_ko_counts.py`

Averages the T3-8 and T3-17 KO clone `counts_sum` files produced by the split-libmeta `get_celltype_info.R` run into a single `HEK293T-KO` column. `perc_avg` and `norm_perc` are copied unchanged from the standard (non-split) run.

```bash
python3 average_ko_counts.py \
  --split-dir celltype_info_split_outputs/ \
  --orig-dir  celltype_info_outputs/ \
  --out-dir   celltype_info_avg_outputs/
```

---

### Stage 3 — `merge_exp_all.py`

Unified script that joins RNA-seq TPM expression data into rNMP EF or celltype_info matrices. Replaces four earlier per-job scripts (`merge_exp.py`, `merge_exp_perc.py`, `merge_exp_EF.py`, `merge_exp_EF_genebody.py`).

Join path:  
`RNA-seq TPM (Ensembl ID)` → `hg38_gene_mapping.tsv` → `rNMP matrix (id = gene name)`

Hardcoded default paths at the top of the script (lines 29–44) point to cluster filesystem locations and **must be updated** before running on a different system.

```bash
# List available predefined jobs
python3 merge_exp_all.py --list

# Run specific jobs (e.g. jobs 4 and 5 for EF gene body + TSS flanks)
python3 merge_exp_all.py --jobs 4 5 \
  --tpm WT-KO.tpms_by_condition.tsv CD4T.tpms_by_condition.tsv H9.tpms_by_condition.tsv \
  --out-root /path/to/output/

# Ad-hoc EF merge for any annotate_outputs directory
python3 merge_exp_all.py \
  --ef-dir /path/to/rNMP_EF/annotate_outputs/hg38_refGene_TSS_downstream_0_1kb_1000bp \
  --tpm WT-KO.tpms_by_condition.tsv \
  --out-root /path/to/output/
```

**Output column layout** (merged EF files, used by `pair_bins.R`):

| Cols | Content |
|---|---|
| 1–6 | chr, start, stop, id, length, Subtype |
| 7 | CD4T |
| 8 | hESC-H9 |
| 9 | HEK293T-WT |
| 10 | HEK293T-KO (avg of T3-8 and T3-17) |
| 11 | KO-T3-8 |
| 12 | KO-T3-17 |
| 13–16 | siRNA cell types |
| 17–23 | RNA-seq TPM columns (expression, used as X axis in pair_bins) |

---

### Stage 4a — `pair_bins.R`

**Original author:** Deepali Kundnani  
**Updated by:** TPW  
**Used for:** Fig 3C (rN only) and Fig 4B–C (rA, rC, rG, rU separately)

Bins genes into 10 percentile groups by expression (X axis column), computes the mean rNMP EF in each bin, fits a linear model, and reports Pearson R and p. Outputs:

| Output | Description |
|---|---|
| `*_stats.tsv` | Per-bin mean EF, SE, N — input for `trends_merge.R` and `two_region_trends.R` |
| `*_corr.tsv` | Pearson R, p, Spearman rho, slope |
| `*bins.svg` | Violin + boxplot per expression bin |
| `*_scatter.png` | Per-gene scatter with regression line |

```bash
Rscript pair_bins.R \
  -m rN_hg38_refGene_TSS_downstream_0_1kb_regions_EF_opp_withexp.tsv \
  -n rN_TSS_WTexp_EF_opp_HEK293T-WT \
  -c "#777777" \         # point/line color (hex)
  -p \                   # use percentile binning (recommended)
  -e 19 \                # X axis column (e.g. col 19 = WT RNA-seq TPM)
  -r 9 \                 # Y axis column (e.g. col 9 = HEK293T-WT EF)
  -y 4 \                 # Y axis maximum
  -o pair_bins_outputs/EF_opp/
```

---

### Stage 4b — `trends_merge.R`

**Original author:** Deepali Kundnani  
**Updated by:** TPW  
**Used for:** Fig 4B–C multi-nucleotide trend overlays

Reads the `*_stats.tsv` files for rA, rC, rG, rU from a single directory and overlays them on one plot. Produces both a labeled and a no-label (`_nolab`) SVG for figure assembly, plus a legend SVG.

```bash
Rscript trends_merge.R \
  -i pair_bins_outputs/EF_opp/ \       # directory containing *_stats.tsv files
  -f "TSS_WTexp_EF_opp_HEK293T-WT" \  # filename pattern (shared across rA/rC/rG/rU)
  -y 6 \                               # Y axis maximum
  -b 1 \                               # Y axis tick interval
  -o pair_bins_outputs/EF_opp/trends/

# No-label version (suppress axis text for figure assembly)
Rscript trends_merge.R ... -q
```

---

### Stage 4c — `two_region_trends.R`

**Updated by:** TPW  
**Used for:** Fig 3C two-region overlay (TSS ± 1 kb flanks vs. TSS-to-TTS gene body)

Takes the `*_stats.tsv` and `*_corr.tsv` output from two separate `pair_bins.R` runs (one per region) and overlays them on a single plot with Pearson R and p annotated per line. Region 1 is plotted in red, region 2 in blue.

```bash
Rscript two_region_trends.R \
  -a pair_bins_tss/EF_opp/rN_TSS_flanks_WTexp_EF_opp_HEK293T-WT_stats.tsv \
  -b pair_bins_genebody/EF_opp/rN_genebody_WTexp_EF_opp_HEK293T-WT_stats.tsv \
  -n "TSS_flanks_WTexp_EF_opp_HEK293T-WT" \
  -l "rNMP in TSS +/- 1 kb" \   # legend label for region 1 (red)
  -k "rNMP in TSS to TTS" \      # legend label for region 2 (blue)
  -y 4 \
  -t 1 \
  -o two_region_plots/EF_opp/

# No-label version
Rscript two_region_trends.R ... -q
```

---

## Nucleotide colors

All trend scripts use a consistent per-nucleotide color scheme:

| Nucleotide | Color | Hex |
|---|---|---|
| rA | Orange | `#D55E00` |
| rC | Blue | `#0072B2` |
| rG | Yellow | `#F0C742` |
| rU | Green | `#009E73` |

---

## Genomic regions used

| Region | BED annotation | Used in |
|---|---|---|
| TSS 0–1 kb downstream | `hg38_refGene_TSS_downstream_0_1kb_1000bp` | Fig 4B–C |
| TSS ± 1 kb flanks (2 kb total) | `hg38_refGene_TSS_flanks_1kb_each_side_2kb_total` | Fig 3C |
| TSS-to-TTS gene body | `hg38_refGene_collapsed` | Fig 3C |
| TSS 2–3 kb upstream | `hg38_refGene_TSS_upstream_2_3kb` | Fig 4B–C |
| TTS 4–5 kb downstream | `hg38_refGene_TTS_downstream_4_5kb` | Fig 4B–C |

EF inputs come from `rNMP_EF/annotate_outputs/`. Nucleotide-composition inputs require running `get_celltype_info.R` on `bedtools nuc`-annotated files first.

---

## Reproducing Fig 3C and Fig 4B-C

Two self-contained wrapper scripts demonstrate the full pipeline end-to-end for HEK293T-WT (the patterns generalize to the other cell types and expression profiles used in the paper by swapping the relevant column indices and TPM files).

### `make_fig3C.sh`

Reproduces Fig 3C: overlays the rN EF vs. expression trend for **TSS ± 1 kb flanks** (red) and **TSS-to-TTS gene body** (blue), for HEK293T-WT, across `both`/`same`/`opp` strands.

Stages: `merge_exp_all.py` (jobs 4 & 5) → `pair_bins.R` (rN, per region/strand) → `two_region_trends.R` (overlay).

```bash
bash make_fig3C.sh
```

### `make_fig4BC.sh`

Reproduces Fig 4B-C: per-nucleotide (rA/rC/rG/rU) EF, `perc_avg`, `norm_perc`, and `counts_sum` vs. expression, overlaid into one trend plot per data type, for HEK293T-WT on the TSS 0-1 kb downstream region across `both`/`same`/`opp` strands.

Stages: `get_celltype_info.R` (per nucleotide) → `merge_exp_all.py` (jobs 2 & 3) → `pair_bins.R` (per nucleotide/data type/strand) → `trends_merge.R` (overlay).

```bash
bash make_fig4BC.sh
```

Both scripts use placeholder input paths (`/path/to/...`) at the top — update these to point at your `rNMP_EF` outputs and RNA-seq TPM files before running. Outputs are written to `output_fig3C/` and `output_fig4BC/` respectively.
