# =============================================================================
# Expression Ratio (75th/25th percentile of log2TPM) per genotype (Fig 4F)
#
# This refactor keeps each dataset in its own analysis scope so we can generate
# figures for multiple inputs without overwriting objects from a previous run.
# =============================================================================

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(limma)
library(scales)
library(svglite)
library(ggrepel)
library(cowplot)

output_dir <- file.path("Fig4F", format(Sys.Date(), "%m-%d-%y"))

cell_type_levels <- c("HEK293T", "RNH2A KO-T3-8", "RNH2A KO-T3-17", "NC", "TOP1")

# The first four shapes follow the high-expression dot plot order in
# make_fig4D-E.R.
point_shapes <- c(
  "HEK293T" = 16,
  "RNH2A KO-T3-8" = 17,
  "RNH2A KO-T3-17" = 15,
  "NC" = 18,
  "TOP1" = 8
)

point_sizes <- c(
  "HEK293T" = 6,
  "RNH2A KO-T3-8" = 6,
  "RNH2A KO-T3-17" = 6,
  "NC" = 8,
  "TOP1" = 5
)

filter_tpm_table <- function(tbl) {
  tbl %>%
    filter(!Chr %in% c("M", "X", "Y")) %>%
    dplyr::select(-c(Chr, Start, End, Strand, Length))
}

load_joined_tpm <- function(
  wt_ko_file = "WT-KO.tpms_by_condition.tsv",
  sitop_file = "siTOP.tpms_by_condition.tsv"
) {
  hek <- read_tsv(wt_ko_file, show_col_types = FALSE) %>%
    filter_tpm_table()

  sitop <- read_tsv(sitop_file, show_col_types = FALSE) %>%
    filter_tpm_table()

  inner_join(hek, sitop, by = "Geneid") %>%
    column_to_rownames("Geneid")
}

load_single_tpm <- function(input_file) {
  tpm <- read_tsv(input_file, show_col_types = FALSE) %>%
    filter_tpm_table() %>%
    as.data.frame()

  rownames(tpm) <- tpm$Geneid
  tpm$Geneid <- NULL
  tpm
}

pca_plot <- function(pca, title, batch) {
  df <- as.data.frame(pca$x[, 1:2])
  df$batch <- batch
  df$sample <- rownames(df)

  ggplot(df, aes(PC1, PC2, color = batch, label = sample)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3) +
    ggtitle(title) +
    theme_classic()
}

save_svg_bundle <- function(results, prefix) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(results$splice_aware)) {
    ggsave(
      filename = file.path(output_dir, paste0(prefix, ".svg")),
      plot = results$splice_aware,
      device = svglite,
      width = 2.8,
      height = 3.2,
      units = "in"
    )
  }

  if (!is.null(results$splice_unaware)) {
    ggsave(
      filename = file.path(output_dir, paste0(prefix, ".svg")),
      plot = results$splice_unaware,
      device = svglite,
      width = 2.8,
      height = 3.2,
      units = "in"
    )
  }

  if (!is.null(results$splice_aware_noy)) {
    ggsave(
      filename = file.path(output_dir, paste0(prefix, "_noy.svg")),
      plot = results$splice_aware_noy,
      device = svglite,
      width = 2.8,
      height = 3.2,
      units = "in"
    )
  }

  if (!is.null(results$splice_unaware_noy)) {
    ggsave(
      filename = file.path(output_dir, paste0(prefix, "_noy.svg")),
      plot = results$splice_unaware_noy,
      device = svglite,
      width = 2.8,
      height = 3.2,
      units = "in"
    )
  }
}

save_shared_legends <- function(results) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ggsave(
    filename = file.path(output_dir, "fig4F_legend_full.svg"),
    plot = results$legend_full_plot,
    device = svglite,
    width = 2.5,
    height = 2.35,
    units = "in"
  )

  ggsave(
    filename = file.path(output_dir, "fig4F_legend_dots_only.svg"),
    plot = results$legend_dots_only_plot,
    device = svglite,
    width = 0.5,
    height = 2.35,
    units = "in"
  )
}

run_fig4d_analysis <- function(
  tpm,
  dataset_name,
  batch1_cols = c("WT", "KO-T3-8", "KO-T3-17"),
  batch2_cols = c("T3-8siRNA-NC_5+5pmol", "T3-8siRNA-TOP1_5+5pmol")
) {
  cat("\n==============================\n")
  cat("Running dataset:", dataset_name, "\n")
  cat("==============================\n")
  cat("Dimensions:", nrow(tpm), "genes x", ncol(tpm), "samples\n")
  cat("Columns:", paste(colnames(tpm), collapse = ", "), "\n")

  log2tpm <- log2(tpm + 1)

  batch <- ifelse(colnames(log2tpm) %in% batch1_cols, "batch1", "batch2")
  batch <- factor(batch)

  cat("\nBatch assignments:\n")
  print(data.frame(sample = colnames(log2tpm), batch = batch))

  genotype <- factor(case_when(
    colnames(log2tpm) %in% c("WT") ~ "WT",
    colnames(log2tpm) %in% c("KO-T3-8", "T3-8siRNA-NC_5+5pmol") ~ "KO-T3-8",
    colnames(log2tpm) %in% c("KO-T3-17") ~ "KO-T3-17",
    colnames(log2tpm) %in% c("T3-8siRNA-TOP1_5+5pmol") ~ "TOP1_KD",
    TRUE ~ "OTHER"
  ))

  log2tpm_corrected <- removeBatchEffect(as.matrix(log2tpm), batch = batch)

  cat("\nPre-correction column means:\n")
  print(round(colMeans(log2tpm), 3))
  cat("\nPost-correction column means:\n")
  print(round(colMeans(log2tpm_corrected), 3))

  expressed_mask <- log2tpm > 0

  compute_ratio_masked <- function(col_name) {
    x <- log2tpm_corrected[expressed_mask[, col_name], col_name]
    quantile(x, 0.75, na.rm = TRUE) / quantile(x, 0.25, na.rm = TRUE)
  }

  ratios <- sapply(colnames(log2tpm_corrected), compute_ratio_masked)
  names(ratios) <- colnames(log2tpm_corrected)

  cat("\nExpression ratios (Q3/Q1 of log2(TPM+1), expressed genes only):\n")
  print(round(ratios, 4))

  label_map <- c(
    "HEK293T-WT" = "HEK293T",
    "WT" = "HEK293T",
    "HEK293T-KO-T3-8" = "RNH2A KO-T3-8",
    "KO-T3-8" = "RNH2A KO-T3-8",
    "HEK293T-KO-T3-17" = "RNH2A KO-T3-17",
    "KO-T3-17" = "RNH2A KO-T3-17",
    "T3-8siRNA-NC_5+5pmol" = "NC",
    "T3-8siRNA-TOP1_5+5pmol" = "TOP1"
  )

  df <- data.frame(
    cell_type = factor(
      label_map[names(ratios)],
      levels = cell_type_levels
    ),
    ratio = as.numeric(ratios)
  )

  cat("\nPlotting data:\n")
  print(df)

  fill_colors <- c(
    "HEK293T" = "#E67D7E",
    "RNH2A KO-T3-8" = "#EEADDA",
    "RNH2A KO-T3-17" = "#BFA4D7",
    "TOP1" = "#79ADD2",
    "NC" = "#A6CEE6"
  )

  border_colors <- c(
    "HEK293T" = "#D62728",
    "RNH2A KO-T3-8" = "#E377C2",
    "RNH2A KO-T3-17" = "#9467BD",
    "TOP1" = "#1F77B4",
    "NC" = "#6BAED6"
  )

  p <- ggplot(df, aes(x = cell_type, y = ratio, color = cell_type, shape = cell_type, size = cell_type)) +
    geom_point(stroke = 1, alpha = 0.75) +
    scale_y_continuous(
      limits = c(3.1, 3.7),
      breaks = seq(3.1, 3.7, by = 0.1),
      expand = c(0, 0)
    ) +
    scale_color_manual(values = border_colors, name = NULL) +
    scale_shape_manual(values = point_shapes, name = NULL) +
    scale_size_manual(values = point_sizes, name = NULL) +
    labs(x = NULL, y = "Expression Ratio\n(75%/25%)") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.text = element_text(size = 15),
      legend.key.height = unit(0.34, "in"),
      legend.key.width = unit(0.28, "in"),
      legend.spacing.y = unit(0.04, "in"),
      legend.background = element_rect(fill = NA, color = NA),
      panel.grid.major.y = element_blank()
    )

  p_no_legend <- p + theme(legend.position = "none")
  legend_full <- get_legend(p)

  p_dots_only_legend <- p +
    theme(
      legend.text = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )

  legend_full_plot <- as_ggplot(legend_full) +
    theme(plot.margin = margin(0, 0, 0, 0))

  legend_dots_only_plot <- as_ggplot(get_legend(p_dots_only_legend)) +
    theme(plot.margin = margin(6, 0, 0, 0, unit = "pt"))

  splice_aware <- p_no_legend +
    scale_y_continuous(
      limits = c(3.8, 5),
      breaks = seq(3.8, 5.0, by = 0.2),
      expand = c(0, 0)
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      legend.position = "none"
    )

  splice_unaware <- p_no_legend +
    scale_y_continuous(
      limits = c(3.1, 3.7),
      breaks = seq(3.1, 3.7, by = 0.1),
      expand = c(0, 0)
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      legend.position = "none"
    )

  splice_aware_noy <- splice_aware +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
    )

  splice_unaware_noy <- splice_unaware +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
    )

  log2tpm_mat <- as.matrix(log2tpm)
  keep <- apply(log2tpm_mat, 1, var) > 0

  cat("Genes removed (zero variance):", sum(!keep), "\n")
  cat("Genes kept:", sum(keep), "\n")

  pca_before <- prcomp(t(log2tpm_mat[keep, ]), scale. = TRUE)
  pca_after <- prcomp(t(log2tpm_corrected[keep, ]), scale. = TRUE)

  list(
    dataset_name = dataset_name,
    tpm = tpm,
    log2tpm = log2tpm,
    log2tpm_corrected = log2tpm_corrected,
    batch = batch,
    genotype = genotype,
    ratios = ratios,
    df = df,
    p = p,
    p_no_legend = p_no_legend,
    legend_full = legend_full,
    legend_full_plot = legend_full_plot,
    legend_dots_only_plot = legend_dots_only_plot,
    splice_aware = splice_aware,
    splice_unaware = splice_unaware,
    splice_aware_noy = splice_aware_noy,
    splice_unaware_noy = splice_unaware_noy,
    pca_before = pca_plot(pca_before, paste("Before batch correction -", dataset_name), batch),
    pca_after = pca_plot(pca_after, paste("After batch correction -", dataset_name), batch)
  )
}

# -----------------------------------------------------------------------------
# Run both datasets without overwriting objects
# -----------------------------------------------------------------------------

ribome_splice_aware_results <- run_fig4d_analysis(
  tpm = load_joined_tpm(),
  dataset_name = "ribome_splice_aware"
)
ribome_splice_aware_results$splice_unaware <- NULL
ribome_splice_aware_results$splice_unaware_noy <- NULL

ribome_splice_unaware_results <- run_fig4d_analysis(
  tpm = load_single_tpm("Ribome-data.tpms_by_condition.tsv"),
  dataset_name = "ribome_splice_unaware"
)
ribome_splice_unaware_results$splice_aware <- NULL
ribome_splice_unaware_results$splice_aware_noy <- NULL

# -----------------------------------------------------------------------------
# Save dataset-specific SVG outputs
# -----------------------------------------------------------------------------

save_shared_legends(ribome_splice_aware_results)
save_svg_bundle(ribome_splice_aware_results, "fig4F_ribome_splice_aware")
save_svg_bundle(ribome_splice_unaware_results, "fig4F_ribome_splice_unaware")

message("Done. Dataset-specific SVGs saved without overwriting prior outputs.")



# -----------------------------------------------------------------------------
# Diagnostics — pull from results objects instead of global vars
# -----------------------------------------------------------------------------

make_density_df <- function(results) {
  as.data.frame(results$log2tpm) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "log2tpm") %>%
    filter(log2tpm > 0) %>%
    mutate(
      batch    = as.character(results$batch[match(sample, colnames(results$log2tpm))]),
      dataset  = results$dataset_name,
      corrected = FALSE
    )
}

make_density_df_corrected <- function(results) {
  as.data.frame(results$log2tpm_corrected) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "log2tpm") %>%
    filter(log2tpm > 0) %>%
    mutate(
      batch    = as.character(results$batch[match(sample, colnames(results$log2tpm_corrected))]),
      dataset  = results$dataset_name,
      corrected = TRUE
    )
}

# Run for both datasets
for (res in list(ribome_splice_aware_results, ribome_splice_unaware_results)) {

  dens_df <- bind_rows(
    make_density_df(res),
    make_density_df_corrected(res)
  ) %>%
    mutate(corrected = ifelse(corrected, "After correction", "Before correction"))

  p_dens <- ggplot(dens_df, aes(x = log2tpm, color = sample, linetype = batch)) +
    geom_density() +
    facet_wrap(~ corrected) +
    labs(title = res$dataset_name) +
    theme_classic()

  print(p_dens)
  print(res$pca_before)
  print(res$pca_after)
}


# Run compute_ratio_masked on uncorrected data for comparison
compute_ratio_raw <- function(col_name, log2tpm) {
  x <- log2tpm[log2tpm[, col_name] > 0, col_name]
  quantile(x, 0.75, na.rm = TRUE) / quantile(x, 0.25, na.rm = TRUE)
}

ratios_uncorrected <- sapply(
  colnames(ribome_splice_aware_results$log2tpm),
  compute_ratio_raw,
  log2tpm = as.matrix(ribome_splice_aware_results$log2tpm)
)

comparison <- data.frame(
  sample      = names(ratios_uncorrected),
  batch       = as.character(ribome_splice_aware_results$batch),
  uncorrected = round(ratios_uncorrected, 4),
  corrected   = round(ribome_splice_aware_results$ratios, 4),
  delta       = round(ribome_splice_aware_results$ratios - ratios_uncorrected, 4)
)
print(comparison)

# -----------------------------------------------------------------------------
# One-off: Median-centering diagnostic
# -----------------------------------------------------------------------------

run_medcenter_comparison <- function(results) {
  log2tpm_mat <- as.matrix(results$log2tpm)

  # Median-center each sample
  log2tpm_medcentered <- sweep(
    log2tpm_mat,
    2,
    apply(log2tpm_mat, 2, median),
    FUN = "-"
  )

  # Reuse existing expressed-gene mask
  expressed_mask <- log2tpm_mat > 0

  compute_ratio_from_mat <- function(mat) {
    sapply(colnames(mat), function(col) {
      x <- mat[expressed_mask[, col], col]
      quantile(x, 0.75, na.rm = TRUE) / quantile(x, 0.25, na.rm = TRUE)
    })
  }

  data.frame(
    sample      = colnames(log2tpm_mat),
    batch       = as.character(results$batch),
    uncorrected = round(compute_ratio_from_mat(log2tpm_mat), 4),
    limma       = round(results$ratios, 4),
    medcentered = round(compute_ratio_from_mat(log2tpm_medcentered), 4)
  ) %>%
    mutate(
      delta_limma  = round(limma - uncorrected, 4),
      delta_median = round(medcentered - uncorrected, 4)
    )
}



# -----------------------------------------------------------------------------
# One-off: Compare Q3/Q1 dot plots across correction methods
# -----------------------------------------------------------------------------

make_ratio_plot <- function(ratios, title, fill_colors, border_colors, label_map) {
  df <- data.frame(
    cell_type = factor(
      label_map[names(ratios)],
      levels = cell_type_levels
    ),
    ratio = as.numeric(ratios)
  )

  ggplot(df, aes(x = cell_type, y = ratio, color = cell_type, shape = cell_type, size = cell_type)) +
    geom_point(stroke = 1, alpha = 0.75) +
    scale_color_manual(values = border_colors, name = NULL) +
    scale_shape_manual(values = point_shapes, name = NULL) +
    scale_size_manual(values = point_sizes, name = NULL) +
    labs(x = NULL, y = "Expression Ratio\n(75%/25%)", title = title) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x    = element_blank(),
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.text    = element_text(size = 11),
      legend.key.size = unit(0.6, "lines"),
      legend.background = element_rect(fill = NA, color = NA)
    )
}

run_correction_comparison_plots <- function(results) {
  log2tpm_mat   <- as.matrix(results$log2tpm)
  expressed_mask <- log2tpm_mat > 0

  # Median-centered matrix
  log2tpm_medcentered <- sweep(
    log2tpm_mat, 2,
    apply(log2tpm_mat, 2, median),
    FUN = "-"
  )

  compute_ratios <- function(mat) {
    sapply(colnames(mat), function(col) {
      x <- mat[expressed_mask[, col], col]
      unname(quantile(x, 0.75, na.rm = TRUE)) /
        unname(quantile(x, 0.25, na.rm = TRUE))
    })
  }

  label_map <- c(
    "WT"                        = "HEK293T",
    "KO-T3-8"                   = "RNH2A KO-T3-8",
    "KO-T3-17"                  = "RNH2A KO-T3-17",
    "T3-8siRNA-NC_5+5pmol"      = "NC",
    "T3-8siRNA-TOP1_5+5pmol"    = "TOP1"
  )

  fill_colors <- c(
    "HEK293T"        = "#E67D7E",
    "RNH2A KO-T3-8"  = "#EEADDA",
    "RNH2A KO-T3-17" = "#BFA4D7",
    "TOP1"           = "#79ADD2",
    "NC"             = "#A6CEE6"
  )

  border_colors <- c(
    "HEK293T"        = "#D62728",
    "RNH2A KO-T3-8"  = "#E377C2",
    "RNH2A KO-T3-17" = "#9467BD",
    "TOP1"           = "#1F77B4",
    "NC"             = "#6BAED6"
  )

  p_uncorrected <- make_ratio_plot(
    compute_ratios(log2tpm_mat),
    "Uncorrected", fill_colors, border_colors, label_map
  )

  p_limma <- make_ratio_plot(
    results$ratios,
    "limma", fill_colors, border_colors, label_map
  )

  p_median <- make_ratio_plot(
    compute_ratios(log2tpm_medcentered),
    "Median-centered", fill_colors, border_colors, label_map
  )

  # Shared y-axis range across all three for fair visual comparison
  all_ratios <- c(
    compute_ratios(log2tpm_mat),
    results$ratios,
    compute_ratios(log2tpm_medcentered)
  )
  y_min <- floor(min(all_ratios) * 10) / 10 - 0.1
  y_max <- ceiling(max(all_ratios) * 10) / 10 + 0.1

  fix_y <- function(p) {
    p + scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = seq(y_min, y_max, by = 0.1)
    )
  }

  plot_grid(
    fix_y(p_uncorrected),
    fix_y(p_limma) + theme(axis.title.y = element_blank(),
                           axis.text.y  = element_blank()),
    nrow       = 1,
    rel_widths = c(1.3, 1),
    labels     = results$dataset_name
  )
}

print(run_correction_comparison_plots(ribome_splice_aware_results))
print(run_correction_comparison_plots(ribome_splice_unaware_results))
cat("\n--- Splice-aware ---\n")
print(run_medcenter_comparison(ribome_splice_aware_results))

cat("\n--- Splice-unaware ---\n")
print(run_medcenter_comparison(ribome_splice_unaware_results))


# -----------------------------------------------------------------------------
# Standalone: Median-centering as batch correction for Q3/Q1 comparison
# Shifts batch2 samples by the offset between their median and the
# batch1 reference median — preserves absolute scale, only corrects location
# -----------------------------------------------------------------------------

run_median_correction_comparison <- function(results) {

  log2tpm_mat    <- as.matrix(results$log2tpm)
  expressed_mask <- log2tpm_mat > 0

  # Compute per-sample medians on expressed genes only
  sample_medians <- sapply(colnames(log2tpm_mat), function(col) {
    median(log2tpm_mat[expressed_mask[, col], col], na.rm = TRUE)
  })

  cat("\nPer-sample medians (expressed genes):\n")
  print(round(sample_medians, 4))

  # Reference = mean of batch1 sample medians
  batch1_mask    <- as.character(results$batch) == "batch1"
  batch1_ref     <- mean(sample_medians[batch1_mask])

  cat("\nBatch1 reference median:", round(batch1_ref, 4), "\n")

  # Compute per-sample offsets — batch1 samples get offset of 0
  offsets <- ifelse(
    batch1_mask,
    0,
    batch1_ref - sample_medians
  )

  cat("\nCorrection offsets applied:\n")
  print(data.frame(
    sample  = colnames(log2tpm_mat),
    batch   = as.character(results$batch),
    median  = round(sample_medians, 4),
    offset  = round(offsets, 4)
  ))

  # Apply offsets column-wise (only batch2 samples shift)
  log2tpm_medcorrected <- sweep(log2tpm_mat, 2, offsets, FUN = "+")

  # Compute Q3/Q1 ratios
  compute_ratios <- function(mat) {
    sapply(colnames(mat), function(col) {
      x <- mat[expressed_mask[, col], col]
      unname(quantile(x, 0.75, na.rm = TRUE)) /
        unname(quantile(x, 0.25, na.rm = TRUE))
    })
  }

  ratios_uncorrected  <- compute_ratios(log2tpm_mat)
  ratios_medcorrected <- compute_ratios(log2tpm_medcorrected)
  ratios_limma        <- results$ratios

  # Summary table
  comparison <- data.frame(
    sample      = colnames(log2tpm_mat),
    batch       = as.character(results$batch),
    uncorrected = round(ratios_uncorrected,  4),
    limma       = round(ratios_limma,        4),
    median_corr = round(ratios_medcorrected, 4),
    delta_limma = round(ratios_limma        - ratios_uncorrected, 4),
    delta_median= round(ratios_medcorrected - ratios_uncorrected, 4)
  )

  cat("\nRatio comparison:\n")
  print(comparison)

  comparison
}

cat("\n--- Splice-aware ---\n")
comp_aware   <- run_median_correction_comparison(ribome_splice_aware_results)

cat("\n--- Splice-unaware ---\n")
comp_unaware <- run_median_correction_comparison(ribome_splice_unaware_results)

# -----------------------------------------------------------------------------
# Plot comparison: Uncorrected vs Limma vs Median-corrected
# -----------------------------------------------------------------------------

label_map <- c(
  "WT"                        = "HEK293T",
  "KO-T3-8"                   = "RNH2A KO-T3-8",
  "KO-T3-17"                  = "RNH2A KO-T3-17",
  "T3-8siRNA-NC_5+5pmol"      = "NC",
  "T3-8siRNA-TOP1_5+5pmol"    = "TOP1"
)

fill_colors <- c(
  "HEK293T"        = "#E67D7E",
  "RNH2A KO-T3-8"  = "#EEADDA",
  "RNH2A KO-T3-17" = "#BFA4D7",
  "NC"             = "#A6CEE6",
  "TOP1"           = "#79ADD2"
)

border_colors <- c(
  "HEK293T"        = "#D62728",
  "RNH2A KO-T3-8"  = "#E377C2",
  "RNH2A KO-T3-17" = "#9467BD",
  "NC"             = "#6BAED6",
  "TOP1"           = "#1F77B4"
)

make_ratio_plot <- function(ratios, title, y_limits) {
  df <- data.frame(
    cell_type = factor(
      label_map[names(ratios)],
      levels = cell_type_levels
    ),
    ratio = as.numeric(ratios)
  )

  ggplot(df, aes(x = cell_type, y = ratio, color = cell_type, shape = cell_type, size = cell_type)) +
    geom_point(stroke = 1, alpha = 0.75) +
    scale_y_continuous(
      limits = y_limits,
      breaks = seq(y_limits[1], y_limits[2], by = 0.1),
      expand = c(0, 0)
    ) +
    scale_color_manual(values = border_colors, name = NULL) +
    scale_shape_manual(values = point_shapes, name = NULL) +
    scale_size_manual(values = point_sizes, name = NULL) +
    labs(x = NULL, y = "Expression Ratio\n(75%/25%)", title = title) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x          = element_blank(),
      axis.title.y         = element_text(size = 13),
      axis.text.y          = element_text(size = 13),
      legend.position      = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.text          = element_text(size = 11),
      legend.key.size      = unit(0.6, "lines"),
      legend.background    = element_rect(fill = NA, color = NA),
      plot.title           = element_text(size = 12, hjust = 0.5)
    )
}

plot_correction_comparison <- function(results, comp) {

  log2tpm_mat    <- as.matrix(results$log2tpm)
  expressed_mask <- log2tpm_mat > 0

  # Recompute median-corrected matrix to get ratios with correct names
  sample_medians <- sapply(colnames(log2tpm_mat), function(col) {
    median(log2tpm_mat[expressed_mask[, col], col], na.rm = TRUE)
  })
  batch1_mask <- as.character(results$batch) == "batch1"
  batch1_ref  <- mean(sample_medians[batch1_mask])
  offsets     <- ifelse(batch1_mask, 0, batch1_ref - sample_medians)
  log2tpm_medcorrected <- sweep(log2tpm_mat, 2, offsets, FUN = "+")

  compute_ratios <- function(mat) {
    sapply(colnames(mat), function(col) {
      x <- mat[expressed_mask[, col], col]
      unname(quantile(x, 0.75, na.rm = TRUE)) /
        unname(quantile(x, 0.25, na.rm = TRUE))
    })
  }

  r_uncorr  <- compute_ratios(log2tpm_mat)
  r_limma   <- results$ratios
  r_median  <- compute_ratios(log2tpm_medcorrected)

  # Shared y-axis across all three panels
  all_vals <- c(r_uncorr, r_limma, r_median)
  y_limits <- c(
    floor(min(all_vals) * 10) / 10 - 0.1,
    ceiling(max(all_vals) * 10) / 10 + 0.1
  )

  p_uncorr <- make_ratio_plot(r_uncorr, "Uncorrected", y_limits)
  p_limma  <- make_ratio_plot(r_limma,  "limma",       y_limits) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank())
  p_median <- make_ratio_plot(r_median, "Median-corrected", y_limits) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank())

  plot_grid(
    p_uncorr, p_limma, p_median,
    nrow       = 1,
    rel_widths = c(1.4, 1, 1),
    labels     = results$dataset_name,
    label_size = 10
  )
}

print(plot_correction_comparison(ribome_splice_aware_results,   comp_aware))
print(plot_correction_comparison(ribome_splice_unaware_results, comp_unaware))
