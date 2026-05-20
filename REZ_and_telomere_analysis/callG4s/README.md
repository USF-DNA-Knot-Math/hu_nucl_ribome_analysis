# callG4s

Scripts for generating a consensus G-quadruplex (G4) peak set from published ChIP-seq data (GSE133379, 293T G4P replicate experiments, hg19) and lifting over coordinates to hg38.

The resulting `intersect_hg38.nochrMXY.bed` file is used as the G4 experimental annotation track in `../plot_REZ.R`.

---

## Scripts

### `callG4sFromNarrowPeaks.sh` *(used for final figures)*

Downloads pre-called narrowPeak files for two ChIP-seq replicates, intersects them to keep only reproducible peaks, and lifts over to hg38. This mirrors the approach used for the RIAN-seq data (intersect rather than union of replicates).

**Output:** `intersect_hg38.nochrMXY.bed`

**Prerequisites:** `bedtools`, `liftOver`, `hg19ToHg38.over.chain`

**Steps:**
1. Download `rep1.narrowPeak` and `rep2.narrowPeak` from GEO (wget lines are commented out; run once)
2. Sort each narrowPeak file
3. `bedtools intersect` to keep peaks present in both replicates
4. Trim to first 3 columns (full narrowPeak format causes liftOver issues)
5. `liftOver` hg19 → hg38
6. Remove chrM, chrX, chrY

---

### `callG4sFrombw.sh` *(exploratory, not used in final figures)*

Alternative approach starting from BigWig signal tracks rather than pre-called peaks. Merges the two replicate BigWig signals and re-calls peaks with MACS2 before lifting over. Retained here as documentation of the approach that was explored.

**Prerequisites:** `bigWigToBedGraph` (UCSC tools), `bedtools`, `MACS2`, `liftOver`

---

## Intermediate files

The following intermediate and output files produced by `callG4sFromNarrowPeaks.sh` are checked in for reference:

| File | Description |
|---|---|
| `rep1.narrowPeak` / `rep2.narrowPeak` | Downloaded replicate peak files (hg19) |
| `rep1.sorted.narrowPeak` / `rep2.sorted.narrowPeak` | Sorted replicate files |
| `intersect.narrowPeak` | Bedtools intersect of both replicates (hg19) |
| `intersect.bed` | First 3 columns of intersect (for liftOver input) |
| `intersect_hg38.bed` | LiftOver output (hg38) |
| `intersect_unmapped.bed` | Regions that failed liftOver |
| `intersect_hg38.nochrMXY.bed` | Final file used in analysis (hg38, no chrM/X/Y) |
| `hg19ToHg38.over.chain` | UCSC chain file used for liftOver |
