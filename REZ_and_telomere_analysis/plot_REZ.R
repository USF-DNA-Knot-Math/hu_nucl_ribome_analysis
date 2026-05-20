# plot_REZ.R
# Original author: Deepali Kundnani
# Updated by: Tyler P. Warner (TPW)
#
# TPW changes:
#   - Load genome chromosome lengths from a .fai index file (no setwd required)
#   - Use configurable input paths instead of hardcoded working-directory assumptions
#   - Support multiple annotation tracks via anno_specs list (CpG islands, G4, R-loops)
#   - Defensive annotation loading: strip UCSC track/browser header lines, coerce
#     numeric coords, clamp to chromosome bounds, harmonize seqinfo, and trim
#   - Explicit tick.pos = c(0, 4) to suppress spurious intermediate axis labels
#   - Output SVG (vector) rather than PNG

library(karyoploteR)
library(GenomicRanges)
library(dplyr)
library(data.table)

# -----------------------
# Set input paths
# Replace these with the actual paths on your system.
# -----------------------
ef_file    <- "pxu/rez/REZ_output_ef.tsv"
rez_file   <- "pxu/rez/REZ_output_common_rezs.tsv"
genome_fai <- "/path/to/hg38-autosomes-noalts.fa.fai"

# Annotation inputs.
# is_1based_inclusive:
#   TRUE  => input start/stop are 1-based inclusive; script converts to 0-based BED
#   FALSE => input is already 0-based BED half-open (no conversion needed)
anno_specs <- list(
  list(
    name  = "cpg_islands",
    path  = "hg38_cpg_islands.bed",
    color = transparent("#9E1717", amount = 0.9),
    is_1based_inclusive = FALSE
  ),
  list(
    name  = "g4_experimental",
    path  = "callG4s/intersect_hg38.nochrMXY.bed",
    color = transparent("#2CA02C", amount = 0.9),
    is_1based_inclusive = FALSE
  ),
  list(
    name  = "g4_comp",
    path  = "/path/to/G4_computational.noMXY.reordered.bed",
    color = transparent("#F2E300", amount = 0.9),
    is_1based_inclusive = FALSE
  ),
  list(
    name  = "rLoop-47",
    path  = "RLregions/RLregions_ucsc.nochrMXY.bed",
    color = transparent("#D9731E", amount = 0.9),
    is_1based_inclusive = FALSE
  ),
  list(
    name  = "rLoop-48",
    path  = "Rian-seq/bgcall_overlaps.bed",
    color = transparent("#A63EDE", amount = 0.9),
    is_1based_inclusive = FALSE
  )
)

# -----------------------
# Genome build (autosomes only, no alts)
# -----------------------
genome <- read.table(genome_fai, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
genome_lengths <- setNames(as.integer(genome[[2]]), genome[[1]])

custom.genome <- toGRanges(
  data.frame(chr   = names(genome_lengths),
             start = rep(1L, length(genome_lengths)),
             end   = as.integer(genome_lengths))
)

# Set seqlengths explicitly so trim/seqinfo works correctly downstream
seqlengths(custom.genome) <- genome_lengths

# -----------------------
# Colors for each cell line
# -----------------------
colors <- c(
  "CD4T"           = "#FF7F0E",
  "hESC"           = "#2CA02C",
  "HEK293T"        = "#D62728",
  "RNH2A-KO T3-8"  = "#E377C2",
  "RNH2A-KO T3-17" = "#9467BD"
)

rez_col <- transparent("#86ECFA", amount = 0.4)

# -----------------------
# Load REZ and EF data
# -----------------------
EF  <- fread(ef_file)
REZ <- fread(rez_file) %>%
  rename(chr = Chromosome, start = Start, stop = End, strand = Strand)

# Average EF across replicates per bin/strand/celltype
EF <- EF %>%
  group_by(Chromosome, Start, End, Strand, Celltype) %>%
  summarise(EF = mean(EF, na.rm = TRUE), .groups = "drop") %>%
  rename(chr = Chromosome, start = Start, stop = End, strand = Strand)

EF.range  <- makeGRangesFromDataFrame(EF,  keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
REZ.range <- makeGRangesFromDataFrame(REZ, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
REZ.split <- sort(split(REZ.range, strand(REZ.range)))

# -----------------------
# Helper: read annotation using only first 3 columns (chr/start/stop).
# Applies optional 1-based inclusive -> BED conversion.
# -----------------------
read_anno_3col <- function(path, is_1based_inclusive = FALSE) {
  dt <- fread(path, header = FALSE, sep = "\t", fill = TRUE, quote = "", data.table = TRUE)

  if (ncol(dt) < 3) stop("Annotation file has <3 columns: ", path)

  dt <- dt[, 1:3]
  setnames(dt, c("chr", "start", "stop"))

  dt <- dt[!is.na(chr)]
  dt <- dt[!grepl("^(track|browser)$", chr)]

  dt[, start := suppressWarnings(as.integer(start))]
  dt[, stop  := suppressWarnings(as.integer(stop))]

  if (anyNA(dt$start) || anyNA(dt$stop)) {
    bad <- dt[is.na(start) | is.na(stop)][1:10]
    stop("Annotation has non-numeric start/stop in: ", path, "\nFirst offending rows:\n",
         paste(capture.output(print(bad)), collapse = "\n"))
  }

  if (is_1based_inclusive) dt[, start := start - 1L]

  dt[start < 0, start := 0L]
  dt <- dt[start < stop]

  gr <- makeGRangesFromDataFrame(
    dt,
    seqnames.field = "chr",
    start.field    = "start",
    end.field      = "stop",
    starts.in.df.are.0based = TRUE
  )

  common_levels  <- intersect(seqlevels(gr), seqlevels(custom.genome))
  dropped_levels <- setdiff(seqlevels(gr), seqlevels(custom.genome))
  if (length(dropped_levels) > 0) {
    message("Dropping annotation seqlevels not in custom.genome for ", basename(path), ": ",
            paste(dropped_levels, collapse = ", "))
  }

  gr <- keepSeqlevels(gr, common_levels, pruning.mode = "coarse")
  seqinfo(gr) <- seqinfo(custom.genome)[seqlevels(gr)]
  gr <- trim(gr)

  oob <- end(gr) > seqlengths(gr)[as.character(seqnames(gr))]
  if (any(oob, na.rm = TRUE)) {
    bad <- gr[which(oob)][1:10]
    stop("Still have out-of-bounds ranges AFTER trim() for ", path, ". Example offenders:\n",
         paste(capture.output(show(bad)), collapse = "\n"))
  }

  gr
}

# -----------------------
# Plotting function for one annotation GRanges.
# Produces a pair of SVGs: whole-genome and chr19 zoom.
# -----------------------
plot_one_annotation <- function(annot_gr, annot_color, out_prefix) {

  pp <- getDefaultPlotParams(plot.type = 2)
  pp$ideogramheight <- 0
  pp$data1height    <- 600
  pp$data2height    <- 600
  pp$leftmargin     <- 0.15
  pp$data1inmargin  <- 0
  pp$data2inmargin  <- 0

  y_ticks <- c(0, 4)

  # Whole genome horizontal plot
  svg(paste0(out_prefix, "_22_chr_large_horizontal.svg"), height = 15, width = 11)
  par(mar = c(2, 2, 0, 0))

  kp1 <- plotKaryotype(genome = custom.genome, plot.type = 2,
                       plot.params = pp, labels.plotter = NULL)
  kpAddChromosomeNames(kp1, cex = 1.2, xoffset = 0.03, yoffset = 0, srt = 0)

  kpAxis(kp1, data.panel = 1, r0 = 0.2, r1 = 0.85, cex = 0.8,
         ymin = 0, ymax = 4, tick.pos = y_ticks,
         tick.len = 10e5, label.margin = -20e5)
  kpAxis(kp1, data.panel = 2, r0 = 0.2, r1 = 0.85, cex = 0.8,
         ymin = 0, ymax = 4, tick.pos = y_ticks,
         tick.len = 10e5, label.margin = -20e5)

  kpPlotRegions(kp1, data = REZ.split$`+`, r0 = 0.2, r1 = 0.85,
                col = rez_col, border = NA, data.panel = 1)
  kpPlotRegions(kp1, data = REZ.split$`-`, r0 = 0.2, r1 = 0.85,
                col = rez_col, border = NA, data.panel = 2)

  kpPlotRegions(kp1, data = annot_gr, r0 = -0.17, r1 = 0.17,
                col = annot_color, avoid.overlapping = FALSE, data.panel = 2)

  for (n in names(colors)) {
    EF.split <- EF.range[elementMetadata(EF.range)[, "Celltype"] == n]
    gr_pos <- sort(split(EF.split, strand(EF.split))$`+`)
    gr_neg <- sort(split(EF.split, strand(EF.split))$`-`)
    kpLines(kp1, data = gr_pos, y = gr_pos$EF, data.panel = 1, col = colors[n],
            ymin = 0, ymax = 4, r0 = 0.2, r1 = 0.85, lwd = 1)
    kpLines(kp1, data = gr_neg, y = gr_neg$EF, data.panel = 2, col = colors[n],
            ymin = 0, ymax = 4, r0 = 0.2, r1 = 0.85, lwd = 1)
  }
  dev.off()

  # chr19 zoom plot
  svg(paste0(out_prefix, "_chr19.svg"), height = 1.5, width = 4)
  par(mfrow = c(1, 1), mar = c(5, 5, 0, 1))

  kp2 <- plotKaryotype(genome = custom.genome, chromosomes = "chr19",
                       plot.type = 2, plot.params = pp, labels.plotter = NULL)

  kpAxis(kp2, data.panel = 1, r0 = 0.2, r1 = 0.85, cex = 1.2,
         ymin = 0, ymax = 4, tick.pos = y_ticks,
         tick.len = 10e5, label.margin = -10e5)
  kpAxis(kp2, data.panel = 2, r0 = 0.2, r1 = 0.85, cex = 1.2,
         ymin = 0, ymax = 4, tick.pos = y_ticks,
         tick.len = 10e5, label.margin = -10e5)

  rez_pos_chr19 <- REZ.split$`+`[seqnames(REZ.split$`+`) == "chr19"]
  rez_neg_chr19 <- REZ.split$`-`[seqnames(REZ.split$`-`) == "chr19"]

  kpPlotRegions(kp2, data = rez_pos_chr19, r0 = 0.2, r1 = 0.85,
                col = rez_col, border = NA, data.panel = 1)
  kpPlotRegions(kp2, data = rez_neg_chr19, r0 = 0.2, r1 = 0.85,
                col = rez_col, border = NA, data.panel = 2)

  kpPlotRegions(kp2, annot_gr[seqnames(annot_gr) == "chr19"],
                r0 = -0.17, r1 = 0.17,
                col = annot_color, avoid.overlapping = FALSE, data.panel = 2)

  for (n in names(colors)) {
    EF.split <- EF.range[elementMetadata(EF.range)[, "Celltype"] == n]
    gr_pos <- sort(split(EF.split, strand(EF.split))$`+`)
    gr_neg <- sort(split(EF.split, strand(EF.split))$`-`)

    gr_pos_chr19 <- gr_pos[seqnames(gr_pos) == "chr19"]
    gr_neg_chr19 <- gr_neg[seqnames(gr_neg) == "chr19"]

    kpLines(kp2, data = gr_pos_chr19, y = gr_pos_chr19$EF, data.panel = 1, col = colors[n],
            ymin = 0, ymax = 4, r0 = 0.2, r1 = 0.85, lwd = 1)
    kpLines(kp2, data = gr_neg_chr19, y = gr_neg_chr19$EF, data.panel = 2, col = colors[n],
            ymin = 0, ymax = 4, r0 = 0.2, r1 = 0.85, lwd = 1)
  }
  dev.off()

  invisible(TRUE)
}

# -----------------------
# Standalone x-axis SVGs (for figure assembly / legend strips)
# -----------------------
plot_axis_only_svgs <- function(out_prefix = "xaxis") {
  pp <- getDefaultPlotParams(plot.type = 2)
  pp$ideogramheight <- 0
  pp$data1height    <- 500
  pp$data2height    <- 500
  pp$leftmargin     <- 0.075
  pp$data1inmargin  <- 0
  pp$data2inmargin  <- 0

  svg(paste0(out_prefix, "_horizontal.svg"), height = 1, width = 11)
  par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
  kp3 <- plotKaryotype(genome = custom.genome, chromosomes = "chr1",
                       plot.type = 2, plot.params = pp, labels.plotter = NULL)
  kpAddBaseNumbers(kp3,
                   tick.dist = 5e7, minor.tick.dist = 1e7,
                   tick.len = 50, minor.tick.len = 25,
                   tick.col = "black", minor.tick.col = "black",
                   cex = 1, units = "Mb")
  dev.off()

  svg(paste0(out_prefix, "_chr19.svg"), height = 0.6, width = 4)
  par(mfrow = c(1, 1), mar = c(0.2, 0, 0, 0))
  kp4 <- plotKaryotype(genome = custom.genome, chromosomes = "chr19",
                       plot.type = 2, plot.params = pp, labels.plotter = NULL)
  kpAddBaseNumbers(kp4,
                   tick.dist = 1e7,
                   tick.len  = 120,
                   tick.col  = "black",
                   cex = 0.9, units = "Mb")
  dev.off()

  invisible(TRUE)
}

# -----------------------
# Run all annotation specs.
# Each produces: <name>_22_chr_large_horizontal.svg and <name>_chr19.svg
# -----------------------
for (spec in anno_specs) {
  if (is.null(spec$name) || is.null(spec$path)) {
    stop("Each anno_specs entry must include at least 'name' and 'path'.")
  }
  if (!file.exists(spec$path)) {
    warning("Annotation file not found (skipping): ", spec$path)
    next
  }

  annot_col <- spec$color
  if (is.null(annot_col)) annot_col <- transparent("gray50", amount = 0.7)

  is1 <- isTRUE(spec$is_1based_inclusive)

  message("Plotting annotation: ", spec$name,
          " | file: ", spec$path,
          " | is_1based_inclusive: ", is1)

  annot_gr <- read_anno_3col(spec$path, is_1based_inclusive = is1)
  plot_one_annotation(annot_gr, annot_col, out_prefix = spec$name)
}

plot_axis_only_svgs(out_prefix = "xaxis")
