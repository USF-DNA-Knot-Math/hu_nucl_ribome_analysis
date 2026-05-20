# DNA_rez

This directory contains scripts and data for visualizing REZ (Ribo-Enriched Zone) tracks derived from DNA-sequencing (input/control) libraries.

---

## Files

| File | Description |
|---|---|
| `DNAseq_REZ_plot.R` | Karyotype plot of DNA-seq enrichment factor across the genome |
| `DNA_rez_500K_counts.tsv` | 500 kb bin rNMP counts for DNA-seq samples (dsF, RE1, RE2, RE3) |

---

## `DNAseq_REZ_plot.R`

### What it does

Produces genome-wide karyotype plots of enrichment factor (EF) calculated from DNA-seq library read counts in 500 kb bins, overlaid with CpG island annotations. Outputs SVG files for the full genome (all autosomes) and a chr19 zoom, plus standalone x-axis strip SVGs.

This script is the DNA-seq counterpart to `../plot_REZ.R`. It uses four DNA-seq library labels (dsF, RE1, RE2, RE3) and the same color scheme.

### Input files

| Variable | Description |
|---|---|
| `anno_file` | CpG island BED annotation (hg38) |
| `ef_file` | 500 kb bin counts TSV for DNA-seq samples (`DNA_rez_500K_counts.tsv`) |
| `genome_fai` | `.fai` index of the reference genome (hg38 autosomes, no alts) |

### Outputs

| File | Description |
|---|---|
| `out_all_chr` | Whole-genome karyotype SVG |
| `out_chr19` | chr19 zoom SVG |
| `out_legend_horizontal` | x-axis strip SVG (chr1, full width) |
| `out_legend_chr19` | x-axis strip SVG (chr19 width) |

### TPW changes (from original)

**Specific changes made:**
- Paths generalized to configurable variables at the top of the script rather than hardcoded values
- Genome chromosome lengths now loaded from a `.fai` index file via an explicit path variable
- Added `stopifnot(file.exists(...))` checks before reading any input file
- Removed unused library imports (`bedr`, `regioneR`, `stringr`, `zoo`, `httpgd`)
- Changed output from PNG to SVG
- Replaced `numticks=2` with `tick.pos = c(0, 4)` to display only the 0 and 4 axis ticks
- Increased `data1height`/`data2height` to 600 and `leftmargin` to 0.15 to match the rNMP-seq figure layout

### How to run

```r
# Set all path variables at the top of the script, then:
source("DNAseq_REZ_plot.R")
```
