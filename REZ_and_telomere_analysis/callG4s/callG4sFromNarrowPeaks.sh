#!/bin/bash
# callG4sFromNarrowPeaks.sh
# Author: Tyler P. Warner (TPW)
#
# Generates a consensus G-quadruplex (G4) peak set from two ChIP-seq replicates
# in narrowPeak format (hg19), lifts over coordinates to hg38, and removes
# chrM/chrX/chrY.
#
# Data source: GSE133379 (293T G4P ChIP-seq, hg19)
# Approach: intersect replicate peaks (same approach used for the RIAN-seq data)
# rather than union, to keep only reproducible peaks.
#
# Prerequisites: bedtools, liftOver, UCSC chain file (hg19ToHg38.over.chain)

# Step 1: Download replicate narrowPeak files
#wget -O rep1.narrowPeak.gz ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133379/suppl/GSE133379_293T-G4P-hg19-rep1.narrowPeak.gz
#wget -O rep2.narrowPeak.gz ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133379/suppl/GSE133379_293T-G4P-hg19-rep2.narrowPeak.gz

# Step 2: Unzip
#gunzip -f rep1.narrowPeak.gz
#gunzip -f rep2.narrowPeak.gz

# Step 3: Sort the peak files
#sort -k1,1 -k2,2n rep1.narrowPeak > rep1.sorted.narrowPeak
#sort -k1,1 -k2,2n rep2.narrowPeak > rep2.sorted.narrowPeak

# Step 4: Intersect peaks (bedtools)
bedtools intersect -a rep1.sorted.narrowPeak -b rep2.sorted.narrowPeak > intersect.narrowPeak

# Step 5: Download UCSC chain file (hg19->hg38)
#wget http://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz
#gunzip hg19ToHg38.over.chain.gz

# Step 6: LiftOver to hg38
# The full narrowPeak format caused issues with liftOver; trim to first 3 columns first.
cut -f1-3 intersect.narrowPeak > intersect.bed
liftOver intersect.bed hg19ToHg38.over.chain intersect_hg38.bed intersect_unmapped.bed

# Step 7: Remove chrM, chrX, chrY
grep -v -E "^(chrM|chrX|chrY)" intersect_hg38.bed > intersect_hg38.nochrMXY.bed

echo "Done. Final hg38 peaks (no chrM/X/Y): intersect_hg38.nochrMXY.bed"
