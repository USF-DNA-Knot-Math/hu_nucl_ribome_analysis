# DNAseq_REZ_plot.R
# Original author: Deepali Kundnani
# Updated by: Tyler P. Warner (TPW)
#
# TPW changes:
#   - Replaced setwd() with explicit configurable path variables
#   - Load genome chromosome lengths from a .fai index file
#   - Added file existence checks (stopifnot) before reading inputs
#   - Removed unused library imports (bedr, regioneR, stringr, zoo, httpgd)
#   - Changed output to SVG for vector-quality figures
#   - Explicit tick.pos = c(0, 4) to suppress spurious intermediate axis labels

library(karyoploteR)
library(dplyr)
library(data.table)
library(GenomicRanges)

# -----------------------
# Set input/output paths
# Replace these with the actual paths on your system.
# -----------------------
anno_file  <- "/path/to/hg38_cpg_islands.bed"
ef_file    <- "/path/to/DNA_rez_500K_counts.tsv"
genome_fai <- "/path/to/hg38-autosomes-noalts.fa.fai"

out_all_chr           <- "/path/to/output/DNA_22_rez_large.svg"
out_chr19             <- "/path/to/output/DNA_chr19.svg"
out_legend_horizontal <- "/path/to/output/legend_horizontal.svg"
out_legend_chr19      <- "/path/to/output/legend_chr19.svg"

colors <- c(
  "dsF" = "#C55A11",
  "RE1" = "#FFC000",
  "RE2" = "#2E75B6",
  "RE3" = "#00B050"
)

annot_col <- transparent("#9E1717", amount = 0.90)
y_ticks   <- c(0, 4)

plot_rez <- function(kp, annot.range, EF.range, axis_cex = 0.8) {
  kpAxis(kp, data.panel = 1, r0 = 0.2, r1 = 0.85, cex = axis_cex,
         ymin = 0, ymax = 4, tick.pos = y_ticks,
         tick.len = 10e5, label.margin = -20e5)
  kpAxis(kp, data.panel = 2, r0 = 0.2, r1 = 0.85, cex = axis_cex,
         ymin = 0, ymax = 4, tick.pos = y_ticks,
         tick.len = 10e5, label.margin = -20e5)

  kpPlotRegions(kp, annot.range, r0 = -0.17, r1 = 0.17,
                col = annot_col, avoid.overlapping = FALSE, data.panel = 2)

  for (n in names(colors)) {
    EF.split <- EF.range[elementMetadata(EF.range)[, "Sample"] == n]
    gr_pos <- sort(split(EF.split, strand(EF.split))$`+`)
    gr_neg <- sort(split(EF.split, strand(EF.split))$`-`)

    kpLines(kp, data = gr_pos, y = gr_pos$EF, data.panel = 1, col = unname(colors[n]),
            ymin = 0, ymax = 4, r0 = 0.2, r1 = 0.85, lwd = 1)
    kpLines(kp, data = gr_neg, y = gr_neg$EF, data.panel = 2, col = unname(colors[n]),
            ymin = 0, ymax = 4, r0 = 0.2, r1 = 0.85, lwd = 1)
  }
}

stopifnot(file.exists(anno_file))
stopifnot(file.exists(ef_file))
stopifnot(file.exists(genome_fai))

genome <- read.table(genome_fai, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
genome_lengths <- setNames(as.integer(genome[[2]]), genome[[1]])

custom.genome <- toGRanges(
  data.frame(
    chr   = names(genome_lengths),
    start = rep(1L, length(genome_lengths)),
    end   = as.integer(genome_lengths)
  )
)
seqlengths(custom.genome) <- genome_lengths

EF <- fread(ef_file)

EF$length <- EF$End - EF$Start
EF$EF <- (EF$Count / EF$length) / (sum(EF$Count) / sum(EF$length))

EF.range <- makeGRangesFromDataFrame(EF, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)

annot <- fread(anno_file, header = FALSE, sep = "\t", fill = TRUE)
annot <- annot[, 1:7]
setnames(annot, c("chr", "start", "stop", "id", "length", "strand", "label"))
annot.range <- makeGRangesFromDataFrame(annot, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
strand(annot.range) <- "*"

pp <- getDefaultPlotParams(plot.type = 2)
pp$ideogramheight <- 0
pp$data1height    <- 600
pp$data2height    <- 600
pp$leftmargin     <- 0.15
pp$data1inmargin  <- 0
pp$data2inmargin  <- 0

svg(out_all_chr, height = 15, width = 11)
par(mar = c(2, 2, 0, 0))
kp1 <- plotKaryotype(genome = custom.genome, plot.type = 2, plot.params = pp, labels.plotter = NULL)
kpAddChromosomeNames(kp1, cex = 1.2, xoffset = 0.03, yoffset = 0, srt = 0)
plot_rez(kp1, annot.range, EF.range, axis_cex = 0.8)
dev.off()

svg(out_chr19, height = 1.5, width = 4)
par(mfrow = c(1, 1), mar = c(5, 5, 0, 1))
kp2 <- plotKaryotype(genome = custom.genome, chromosomes = "chr19",
                     plot.type = 2, plot.params = pp, labels.plotter = NULL)

kpAxis(kp2, data.panel = 1, r0 = 0.2, r1 = 0.85, cex = 1.2,
       ymin = 0, ymax = 4, tick.pos = y_ticks,
       tick.len = 10e5, label.margin = -10e5)
kpAxis(kp2, data.panel = 2, r0 = 0.2, r1 = 0.85, cex = 1.2,
       ymin = 0, ymax = 4, tick.pos = y_ticks,
       tick.len = 10e5, label.margin = -10e5)

kpPlotRegions(kp2, annot.range[seqnames(annot.range) == "chr19"],
              r0 = -0.17, r1 = 0.17,
              col = annot_col, avoid.overlapping = FALSE, data.panel = 2)

for (n in names(colors)) {
  EF.split <- EF.range[elementMetadata(EF.range)[, "Sample"] == n]
  gr_pos <- sort(split(EF.split, strand(EF.split))$`+`)
  gr_neg <- sort(split(EF.split, strand(EF.split))$`-`)

  gr_pos_chr19 <- gr_pos[seqnames(gr_pos) == "chr19"]
  gr_neg_chr19 <- gr_neg[seqnames(gr_neg) == "chr19"]

  kpLines(kp2, data = gr_pos_chr19, y = gr_pos_chr19$EF, data.panel = 1,
          col = unname(colors[n]), ymin = 0, ymax = 4, r0 = 0.2, r1 = 0.85, lwd = 1)
  kpLines(kp2, data = gr_neg_chr19, y = gr_neg_chr19$EF, data.panel = 2,
          col = unname(colors[n]), ymin = 0, ymax = 4, r0 = 0.2, r1 = 0.85, lwd = 1)
}
dev.off()

svg(out_legend_horizontal, height = 1, width = 11)
par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
kp3 <- plotKaryotype(genome = custom.genome, chromosomes = "chr1",
                     plot.type = 2, plot.params = pp, labels.plotter = NULL)
kpAddBaseNumbers(kp3,
                 tick.dist = 5e7, minor.tick.dist = 1e7,
                 tick.len = 50, minor.tick.len = 25,
                 tick.col = "black", minor.tick.col = "black",
                 cex = 1, units = "Mb")
dev.off()

svg(out_legend_chr19, height = 0.6, width = 4)
par(mfrow = c(1, 1), mar = c(0.2, 0, 0, 0))
kp4 <- plotKaryotype(genome = custom.genome, chromosomes = "chr19",
                     plot.type = 2, plot.params = pp, labels.plotter = NULL)
kpAddBaseNumbers(kp4,
                 tick.dist = 1e7,
                 tick.len  = 120,
                 tick.col  = "black",
                 cex = 0.9, units = "Mb")
dev.off()
