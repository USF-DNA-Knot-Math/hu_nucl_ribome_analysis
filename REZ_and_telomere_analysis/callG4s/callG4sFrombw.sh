#!/bin/bash
# callG4sFrombw.sh
# Author: Tyler P. Warner (TPW)
#
# Alternative G4 peak calling starting from BigWig signal tracks rather than
# pre-called narrowPeak files. Downloads hg19 ratio BigWigs for two 293T
# replicates (GSE133379), merges the signal, calls peaks with MACS2, and
# lifts over to hg38.
#
# This approach was explored as an alternative to callG4sFromNarrowPeaks.sh.
# The narrowPeak intersection method (callG4sFromNarrowPeaks.sh) was used
# for the final paper figures.
#
# Prerequisites: bigWigToBedGraph (UCSC tools), bedtools, MACS2, liftOver

GSE="GSE133379"
DIR="${GSE}_data"
mkdir -p "$DIR" && cd "$DIR"

# Download BigWig signal tracks (hg19 ratios) for 293T replicates
#wget -O rep1.bw "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/${GSE}/suppl/${GSE}_293T-G4P-hg19-ratio-rep1.bigWig"
#wget -O rep2.bw "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/${GSE}/suppl/${GSE}_293T-G4P-hg19-ratio-rep2.bigWig"

# Convert BigWig to bedGraph
bigWigToBedGraph rep1.bw rep1.bedGraph
bigWigToBedGraph rep2.bw rep2.bedGraph

# Sort bedGraphs
sort -k1,1 -k2,2n rep1.bedGraph > rep1.sorted.bedGraph
sort -k1,1 -k2,2n rep2.bedGraph > rep2.sorted.bedGraph

# Merge signals: sum scores across the two replicates
bedtools unionbedg -i rep1.sorted.bedGraph rep2.sorted.bedGraph > union.bedGraph
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4+$5}' union.bedGraph > merged_signal.bedGraph

# Call peaks on merged signal with MACS2
macs2 callpeak \
    -t merged_signal.bedGraph \
    -f BED \
    -n 293T_G4P \
    -g hs \
    --outdir macs2_out

# Download UCSC chain file for hg19->hg38 if not already present
#wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz
#gunzip hg19ToHg38.over.chain.gz

# LiftOver peaks to hg38
liftOver macs2_out/293T_G4P_peaks.narrowPeak \
    hg19ToHg38.over.chain \
    293T_G4P_peaks.hg38.bed \
    293T_G4P_peaks.unmapped.bed

# Extract first 3 columns if only coordinates are needed
cut -f1-3 293T_G4P_peaks.hg38.bed > 293T_G4P_peaks.hg38.coords.bed
