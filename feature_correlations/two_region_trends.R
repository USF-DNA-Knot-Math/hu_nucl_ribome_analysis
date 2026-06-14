#!/usr/bin/env Rscript
# Overlays rNMP EF trends for two genomic regions (e.g. TSS window vs gene body)
# on a single plot per cell type, with Pearson R and p annotated per line.
#
# Inputs: _stats.tsv and _corr.tsv files produced by pair_bins.R for rN.
# Output: *_two_region_trend.svg per cell type.

library("optparse")

option_list = list(
  make_option(c("-a", "--stats1"), type="character", default=NULL,
              help="Path to _stats.tsv for region 1 (TSS window, plotted red)",
              metavar="file"),
  make_option(c("-b", "--stats2"), type="character", default=NULL,
              help="Path to _stats.tsv for region 2 (gene body, plotted blue)",
              metavar="file"),
  make_option(c("-n", "--name"), type="character", default=NULL,
              help="Output file base name (no extension)", metavar="character"),
  make_option(c("-l", "--label1"), type="character", default="rNMP in TSS +/- 1 kb",
              help="Legend label for region 1 [default '%default']", metavar="character"),
  make_option(c("-k", "--label2"), type="character", default="rNMP in TSS to TTS",
              help="Legend label for region 2 [default '%default']", metavar="character"),
  make_option(c("-y", "--ymax"), type="numeric", default=6,
              help="Y axis maximum [default %default]", metavar="numeric"),
  make_option(c("-s", "--ymin"), type="numeric", default=0,
              help="Y axis minimum [default %default]", metavar="numeric"),
  make_option(c("-t", "--tick_interval"), type="numeric", default=1,
              help="Major Y axis tick interval [default %default]", metavar="numeric"),
  make_option(c("-o", "--outdir"), type="character", default=".",
              help="Output directory [default '%default']", metavar="character"),
  make_option(c("-q", "--nolab"), action="store_true", default=FALSE,
              help="Suppress axis labels, legend, and R/p annotations; append _nolab to filename [default %default]")
)

opt_parser = OptionParser(option_list = option_list)
opt        = parse_args(opt_parser)

if (is.null(opt$stats1) || is.null(opt$stats2) || is.null(opt$name)) {
  print_help(opt_parser)
  stop("--stats1, --stats2, and --name are required.", call. = FALSE)
}

suppressMessages({
  library(ggplot2)
  library(dplyr)
})

# ── Helper: read corr file and return Pearson R + p ───────────────────────────
read_corr <- function(stats_path) {
  corr_path <- sub("_stats\\.tsv$", "_corr.tsv", stats_path)
  if (!file.exists(corr_path)) {
    warning("corr file not found: ", corr_path)
    return(list(R = NA_real_, p = NA_real_))
  }
  corr <- read.table(corr_path, sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE, fill = TRUE, quote = "")
  # The file has an unnamed leading index column; the stats names live in
  # the column literally called "stats" (or the second column if unnamed).
  key_col <- if ("stats" %in% names(corr)) "stats" else names(corr)[2]
  val_col <- if ("value" %in% names(corr)) "value" else names(corr)[3]
  get_val <- function(name) {
    idx <- which(corr[[key_col]] == name)
    if (length(idx) == 0) return(NA_real_)
    as.numeric(corr[[val_col]][idx[1]])
  }
  list(
    R = get_val("Pearson'sR"),
    p = get_val("Pearson'sPvalue")
  )
}

# ── Format p-value the same way as the reference figure ──────────────────────
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) formatC(p, format = "e", digits = 2) else round(p, 3)
}

# ── Read data ─────────────────────────────────────────────────────────────────
s1 <- read.table(opt$stats1, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
s2 <- read.table(opt$stats2, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

c1 <- read_corr(opt$stats1)
c2 <- read_corr(opt$stats2)

line_colors <- c("#FF0000", "#0000FF")
fill_colors <- c("#FF9999", "#9999FF")
names(line_colors) <- c(opt$label1, opt$label2)
names(fill_colors) <- c(opt$label1, opt$label2)

s1$region <- opt$label1
s2$region <- opt$label2
s1$ylo <- s1$ymean - s1$se
s1$yhi <- s1$ymean + s1$se
s2$ylo <- s2$ymean - s2$se
s2$yhi <- s2$ymean + s2$se

dat <- bind_rows(s1, s2)
dat$region <- factor(dat$region, levels = c(opt$label1, opt$label2))

# ── Annotation positions ──────────────────────────────────────────────────────
# Place R/p text at x = midpoint, y near each line's mean
x_mid   <- median(dat$groups)
y1_mean <- mean(s1$ymean, na.rm = TRUE)
y2_mean <- mean(s2$ymean, na.rm = TRUE)
y_range <- opt$ymax - opt$ymin
offset  <- y_range * 0.06   # small vertical nudge

ann1 <- sprintf("italic(R) == %.3f*','~italic(p) == %s", c1$R, fmt_p(c1$p))
ann2 <- sprintf("italic(R) == %.3f*','~italic(p) == %s", c2$R, fmt_p(c2$p))

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(dat, aes(x = groups, y = ymean, colour = region, fill = region,
                     group = region)) +
  geom_ribbon(aes(ymin = ylo, ymax = yhi), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.4)

if (!opt$nolab) {
  p <- p +
    annotate("text", x = x_mid, y = y1_mean + offset,
             label = ann1, parse = TRUE,
             colour = line_colors[opt$label1], size = 2.5, hjust = 0.5) +
    annotate("text", x = x_mid, y = y2_mean - offset,
             label = ann2, parse = TRUE,
             colour = line_colors[opt$label2], size = 2.5, hjust = 0.5)
}

p <- p +
  scale_colour_manual(values = line_colors,
                      guide  = guide_legend(title = NULL,
                                            override.aes = list(linewidth = 0.8))) +
  scale_fill_manual(values = fill_colors, guide = "none") +
  scale_x_continuous(expand = c(0, 0),
                     breaks = s1$groups,
                     labels = as.character(s1$groups),
                     limits = c(min(s1$groups), max(s1$groups) + 0.5)) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(opt$ymin, opt$ymax + opt$ymax * 0.004),
                     breaks = seq(opt$ymin, opt$ymax, by = opt$tick_interval),
                     name   = "rNMP EF") +
  theme_classic(base_size = 20) +
  theme(
    panel.border         = element_blank(),
    axis.title           = element_blank(),
    axis.ticks           = element_line(linewidth = 0.28, colour = "black"),
    axis.ticks.length    = unit(0.075, "cm"),
    axis.line            = element_line(linewidth = 0.28, colour = "black"),
    panel.background     = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor     = element_blank(),
    panel.grid.major     = element_blank(),
    plot.background      = element_rect(fill = "transparent", colour = NA),
    strip.background     = element_blank(),
    plot.margin          = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )

if (opt$nolab) {
  p <- p + theme(
    legend.position      = "none",
    axis.text            = element_blank()
  )
} else {
  p <- p + theme(
    legend.position      = "top",
    legend.key.width     = unit(0.4, "cm"),
    legend.text          = element_text(size = 7),
    legend.margin        = margin(0, 0, 0, 0),
    legend.box.margin    = margin(-4, 0, -4, 0),
    axis.text            = element_text(colour = "black"),
    axis.text.x          = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y          = element_text(size = 9)
  )
}

# ── Save ──────────────────────────────────────────────────────────────────────
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
suffix <- if (opt$nolab) "_nolab" else ""
out_path <- file.path(opt$outdir, paste0(opt$name, "_two_region_trend", suffix, ".svg"))
svg(out_path, width = 2.5, height = 2.2)
print(p)
dev.off()

writeLines(paste("Written:", out_path))
