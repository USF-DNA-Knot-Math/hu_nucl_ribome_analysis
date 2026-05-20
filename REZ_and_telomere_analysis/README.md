# REZ and Telomere Analysis

This directory contains scripts for Ribo-enriched Zone (REZ) identification and visualization across the human nuclear genome (hg38, autosomes only, no alts).

---

## Directory Structure

```
REZ_and_telomere_analysis/
├── plot_REZ.R                  # Whole-genome REZ karyotype plot (rNMP-seq data)
├── plot_freq.R                 # rNMP frequency plots
├── hg38_cpg_islands.bed        # CpG island annotation (hg38)
├── genome_telomere_unit.tsv    # Telomere repeat unit lengths
├── telomere.tsv                # Telomere coordinates
├── DNA_rez/                    # REZ analysis for DNA-seq control data
│   ├── DNAseq_REZ_plot.R       # Karyotype plot for DNA-seq REZ tracks
│   └── DNA_rez_500K_counts.tsv # 500 kb bin counts for DNA-seq
├── RLregions/                  # R-loop region (RLseq, S96) data
│   ├── RLregions_rez.R         # REZ overlap with RLseq peaks
│   └── RLregions_ucsc          # UCSC track file for RLseq regions
├── Rian-seq/                   # RIAN-seq R-loop peak data
│   ├── Rianseq_rez.R           # REZ overlap with RIAN-seq peaks
│   ├── bgcall_overlaps.bed     # Background-called peak overlaps
│   └── get_peaks.sh            # Shell script to extract peaks from bigwig
├── callG4s/                    # G-quadruplex peak processing scripts
├── rez_analysis_output/        # REZ calling outputs at various thresholds
└── pxu/                        # Penghao Xu pipeline (rNMP binning, telomere, exon/intron)
```

---

## `plot_REZ.R`

### What it does

Produces genome-wide karyotype plots of rNMP enrichment factor (EF) relative to REZ (Ribo-Enriched Zone) positions. EF is plotted per strand (forward = top panel, reverse = bottom panel) as colored lines for each cell type. Annotation features (CpG islands, G4 structures, R-loops) are overlaid as shaded bands below the zero line.

Outputs a pair of SVG files per annotation:
- `<name>_22_chr_large_horizontal.svg` — all autosomes tiled horizontally
- `<name>_chr19.svg` — chr19 zoom

Also outputs standalone x-axis strip SVGs for figure assembly:
- `xaxis_horizontal.svg`
- `xaxis_chr19.svg`

### Input files

| Variable | Description |
|---|---|
| `ef_file` | TSV of rNMP EF values per 500 kb bin per strand per cell type (from `pxu/rez/`) |
| `rez_file` | TSV of common REZ coordinates across cell types (from `pxu/rez/`) |
| `genome_fai` | `.fai` index of the reference genome (hg38 autosomes, no alts) |
| `anno_specs` | List of annotation BED files; each entry specifies name, path, color, and coordinate system |

### TPW changes (from original)

The original script (`plot_REZ.R` in the inherited codebase) used `setwd()` to a Windows local path and loaded EF data with `read.table()` pointing to hardcoded Windows file paths. It also used `numticks=2` for axis labels, which produced a spurious tick at EF=2. The plot function relied on global variables set earlier in the script rather than taking arguments, making it difficult to re-run for different annotations.

**Specific changes made:**
- Removed `setwd()` and replaced all paths with configurable variables at the top of the script
- Added genome loading from a `.fai` index file (rather than the original `filtered_hg38-nucleus-noXY.fa.fai` with no path)
- Replaced `numticks=2` with `tick.pos = c(0, 4)` to show only the 0 and 4 tick marks
- Wrapped plotting logic in `plot_one_annotation()` function accepting a GRanges object, allowing iteration over multiple annotations
- Added `read_anno_3col()` helper with defensive loading: strips UCSC track/browser header lines, coerces coordinates to integer, clamps to valid BED range, harmonizes seqlevels against the plotted genome, and calls `trim()` to clip any out-of-bounds intervals
- Changed output format from PNG to SVG for lossless vector output
- Removed unused library imports (`bedr`, `regioneR`, `stringr`, `zoo`, `httpgd`)
- Added `plot_axis_only_svgs()` to produce standalone x-axis strips for multi-panel figure assembly

### How to run

```r
# Set genome_fai, ef_file, rez_file, and anno_specs paths at the top of the script, then:
source("plot_REZ.R")
```

The script skips any annotation whose file path does not exist (with a warning) rather than erroring out, so partial annotation sets work fine.
