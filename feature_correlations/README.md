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

---

## Final figure scripts: Fig 3E, Fig 4D-E, Fig 4F

These three standalone scripts generate the final versions of Fig 3E, Fig 4D-E, and Fig 4F directly from EF and TPM tables (independent of the `pair_bins.R`/`merge_exp_all.py` pipeline above). Each is self-contained: edit the input paths/constants near the top of the script, then run with `Rscript`.

### `make_fig3E.R`

**Used for:** Fig 3E — strand-bias (opposite-strand minus same-strand rN EF) boxplots, stratified by expression level (Low/Mod/High), for 7 cell types across two genomic regions (Promoter: TSS 0-1kb downstream; GeneBody: TSS-to-TTS).

Inputs:

| Path | Description |
|---|---|
| `Tyler_rNMPs/hg38_refGene_TSS_downstream_0_1kb_1000bp/rN/out/*_rN_{opp,same}_EF_regions_EF.tsv` | Promoter region rN EF tables |
| `Tyler_rNMPs/hg38_refGene_collapsed/rN/out/*_rN_{opp,same}_EF_regions_EF.tsv` | Gene body region rN EF tables |
| `Expression-level-lists/bowtie/{CD4T,H9,WT,KO-T3-8,KO-T3-17,T3-8siRNA-NC_5_5pmol,T3-8siRNA-TOP1_5_5pmol}{Low,Mod,High}_pc.bed` | Per-cell-type, per-expression-level gene ID lists (col 4 = gene ID) |

For each library, the script computes `bias = opp - same` EF, groups by cell type/expression level/region, and runs Wilcoxon tests (Low vs Mod, Low vs High) to annotate significance brackets.

```bash
Rscript make_fig3E.R
```

Outputs (under `Fig3E/050526/`):

| Output | Description |
|---|---|
| `Fig3E_0-1kb-downstream-of-TSS.svg` / `_nolabs.svg` | Promoter region boxplots |
| `Fig3E_TSS-to-TTS.svg` / `_nolabs.svg` | Gene body region boxplots |
| `Fig3E_pvalues_summary.csv` | Wilcoxon p-values for all cell type / region / expression-group comparisons |
| `Fig3E_legend_vertical.svg` / `_notext.svg` | Standalone vertical legend |

### `make_fig4D-E.R`

**Used for:** Fig 4D-E — rG EF difference (asinh-scaled) vs. RNA log2TPM expression difference, quadrant analysis, for four comparisons (WTvKO with low/high rNMP bias, NCvTOP1, TOP1vWT). Genes are binned into Low/Medium/High expression by WT RNA quartiles before any bias threshold is applied.

Inputs:

| Path | Description |
|---|---|
| `Expression_TPMs/Bowtie/Ribome-data.tpms_by_condition.tsv` | RNA-seq TPMs (WT, KO-T3-8, KO-T3-17, siRNA-NC, siRNA-TOP1 columns) |
| `Tyler_rNMPs/hg38_refGene_TSS_downstream_0_1kb_1000bp/rG/out/*_rG_both_EF_regions_EF_avg.tsv` | rG EF (both strands) |
| `Tyler_rNMPs/hg38_refGene_TSS_downstream_0_1kb_1000bp/rN/out/*_rN_{same,opp}_EF_regions_EF_avg.tsv` | rN EF (same/opp strand), used to compute rNMP strand bias for KO and TOP1 |
| `Static_annotations/ribome_ensembl_to_symbol_map_pc.tsv` | Static Ensembl-to-symbol map (cols `ENSEMBL`, `SYMBOL`), restricted to protein-coding genes |

TOP1vWT RNA values are batch-corrected with `limma::removeBatchEffect` (WT/KO vs. NC/TOP1 as batches) before computing differences.

```bash
Rscript make_fig4D-E.R
```

Outputs (under `Fig4D/05-14-26_final/static_gene_mapping_pc/dynamic_pre_filter_with_batchnorm_keep/`):

| Output | Description |
|---|---|
| `plots/*_{Low,Medium,High}.svg` / `_nolabs.svg` | Scatter plots per comparison/expression bin |
| `annotated/*_annotated.svg` | Scatter plots annotated with quadrant counts and Q2/Q3, Q1/Q4 ratios |
| `dot_plots/high_expression_*_dot_plot*.svg` | High-expression quadrant-ratio summary dot plots and legend |
| `*_quadrant_tally.csv`, `*_stage_summary.csv`, `*_removed_rG_both0_genes.csv` | Per-comparison QC tables |
| `all_comparisons_*.csv`, `all_modes_ratio_summary.csv`, `mapping_mode_gene_counts.csv` | Combined summary tables across all bias-threshold cutoffs |

### `make_fig4F.R`

**Used for:** Fig 4F — Expression Ratio (Q3/Q1 of log2TPM, expressed genes only) per genotype, comparing splice-aware (joined WT-KO + siTOP TPM tables) and splice-unaware (single Ribome-data TPM table) datasets, with `limma::removeBatchEffect` batch correction.

Inputs:

| Path | Description |
|---|---|
| `WT-KO.tpms_by_condition.tsv`, `siTOP.tpms_by_condition.tsv` | Splice-aware dataset (inner-joined on `Geneid`) |
| `Ribome-data.tpms_by_condition.tsv` | Splice-unaware dataset |

Each TPM table is expected to have `Chr`, `Start`, `End`, `Strand`, `Length`, `Geneid` columns plus one column per sample; `Chr` in `{M, X, Y}` is excluded.

```bash
Rscript make_fig4F.R
```

Outputs (under `Fig4F/<MM-DD-YY>/`, dated by run date):

| Output | Description |
|---|---|
| `fig4F_ribome_splice_aware.svg` / `_noy.svg` | Splice-aware Expression Ratio dot plot |
| `fig4F_ribome_splice_unaware.svg` / `_noy.svg` | Splice-unaware Expression Ratio dot plot |
| `fig4F_legend_full.svg`, `fig4F_legend_dots_only.svg` | Shared legends |

The script also runs a series of diagnostic comparisons (density plots, PCA before/after batch correction, median-centering vs. limma correction) and prints summary tables to the console; these are exploratory and not required to reproduce the final figure.
