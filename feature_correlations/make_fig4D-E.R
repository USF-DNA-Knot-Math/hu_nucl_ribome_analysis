suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
})

# -----------------------------
# Settings
# -----------------------------
# rG values are asinh-scaled for plotting so large rG differences do not
# dominate the x-axis. RNA differences remain on their original log2 scale.
ASINH_SCALE <- 5
Y_LIM_FIXED <- 5

# Cutoffs used for the comprehensive ratio summary. The figure uses only a
# subset, but the CSV output includes every cutoff below.
LE_CUTOFFS <- c(0.10, 0.15, 0.25, 0.30, 0.40, 0.50, 0.60, 0.75, 0.80, 0.90, 1.00)
GT_CUTOFFS <- c(1.00, 1.50, 2.00, 2.50, 3.00, 4.00)

# All outputs from this final script are written under this folder.
OUT_DIR <- "Fig4D/05-14-26_final"

# Use the static PC Ensembl-to-symbol map. This avoids changing results when
# Bioconductor annotation packages differ between machines.
SYMBOL_MAP_FILE_PC <- "Static_annotations/ribome_ensembl_to_symbol_map_pc.tsv"

high <- "#A31621"
med  <- "#EF5B5B"
low  <- "#F3A712"

facet_levels <- c("Low", "Medium", "High")
facet_colors <- c(Low = low, Medium = med, High = high)

RNA_FILE <- "Expression_TPMs/Bowtie/Ribome-data.tpms_by_condition.tsv"
RG_FILE  <- "Tyler_rNMPs/hg38_refGene_TSS_downstream_0_1kb_1000bp/rG/out/hg38_refGene_TSS_downstream_0_1kb_1000bp_rG_both_EF_regions_EF_avg.tsv"
RN_SAME  <- "Tyler_rNMPs/hg38_refGene_TSS_downstream_0_1kb_1000bp/rN/out/hg38_refGene_TSS_downstream_0_1kb_1000bp_rN_same_EF_regions_EF_avg.tsv"
RN_OPP   <- "Tyler_rNMPs/hg38_refGene_TSS_downstream_0_1kb_1000bp/rN/out/hg38_refGene_TSS_downstream_0_1kb_1000bp_rN_opp_EF_regions_EF_avg.tsv"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Helpers
# -----------------------------
avg_cols <- function(df, pattern, label = pattern) {
  cols <- grep(pattern, names(df), value = TRUE)
  if (length(cols) == 0) stop("No columns matched pattern for ", label, ": ", pattern)
  rowMeans(df[, cols, drop = FALSE], na.rm = TRUE)
}

# Normalize IDs before joins: strip Ensembl version suffixes and standardize
# symbols to uppercase.
clean_ens <- function(x) sub("\\.\\d+$", "", trimws(as.character(x)))
clean_sym <- function(x) toupper(trimws(as.character(x)))
threshold_tag <- function(x) gsub("\\.", "_", sprintf("%.2f", x))

# Load a precomputed static mapping table and keep only RNA-input Ensembl IDs.
build_symbol_map_static <- function(ensembl_ids, map_file) {
  if (!file.exists(map_file)) stop("Static Ensembl-to-symbol map is missing: ", map_file)

  map_tbl <- readr::read_tsv(
    map_file,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )

  required_cols <- c("ENSEMBL", "SYMBOL")
  missing_cols <- setdiff(required_cols, names(map_tbl))
  if (length(missing_cols) > 0) {
    stop("Static symbol map is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  map_tbl %>%
    dplyr::transmute(
      ENSEMBL = clean_ens(ENSEMBL),
      SYMBOL = clean_sym(SYMBOL)
    ) %>%
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
    dplyr::filter(ENSEMBL %in% unique(clean_ens(ensembl_ids))) %>%
    dplyr::distinct(ENSEMBL, SYMBOL)
}

# Assign expression bins from WT RNA quartiles. This is intentionally called
# before bias thresholds are applied, so Low/Medium/High bins are not dependent
# on any particular cutoff.
assign_dynamic_prefilter_bins <- function(df, expr_ref_col = "expr_ref") {
  qs <- stats::quantile(df[[expr_ref_col]], probs = c(0.25, 0.75), na.rm = TRUE)
  df %>%
    dplyr::mutate(
      Exp.level = dplyr::case_when(
        .data[[expr_ref_col]] < qs[1] ~ "Low",
        .data[[expr_ref_col]] >= qs[1] & .data[[expr_ref_col]] <= qs[2] ~ "Medium",
        TRUE ~ "High"
      ),
      Exp.level = factor(Exp.level, levels = facet_levels)
    )
}

# Count plotted points in the four quadrants. Points on either zero axis are not
# included in Q1-Q4 ratios, but are tracked in n_axis/check_sum.
tally_quadrants <- function(df, facet_col = "Exp.level") {
  df %>%
    dplyr::mutate(
      on_axis = (x_plot == 0) | (exp.diff.wt == 0),
      Q1 = (!on_axis) & (x_plot > 0) & (exp.diff.wt > 0),
      Q2 = (!on_axis) & (x_plot < 0) & (exp.diff.wt > 0),
      Q3 = (!on_axis) & (x_plot < 0) & (exp.diff.wt < 0),
      Q4 = (!on_axis) & (x_plot > 0) & (exp.diff.wt < 0)
    ) %>%
    dplyr::group_by(.data[[facet_col]]) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_axis = sum(on_axis),
      Q1 = sum(Q1),
      Q2 = sum(Q2),
      Q3 = sum(Q3),
      Q4 = sum(Q4),
      .groups = "drop"
    ) %>%
    dplyr::rename(Exp.level = 1) %>%
    dplyr::mutate(
      Exp.level = factor(Exp.level, levels = facet_levels),
      check_sum = Q1 + Q2 + Q3 + Q4 + n_axis
    ) %>%
    dplyr::right_join(
      tibble::tibble(Exp.level = factor(facet_levels, levels = facet_levels)),
      by = "Exp.level"
    ) %>%
    dplyr::mutate(
      dplyr::across(c(n_total, n_axis, Q1, Q2, Q3, Q4, check_sum), ~ dplyr::coalesce(., 0L))
    )
}

# Avoid infinite ratios when a quadrant denominator is zero.
safe_ratio_num <- function(num, den) {
  dplyr::if_else(is.na(den) | den == 0, NA_real_, num / den)
}

# Human-readable two-decimal ratio labels for annotated SVGs.
format_ratio <- function(num, den) {
  if (is.na(den) || den == 0) return("NA")
  sprintf("%.2f", num / den)
}

# Base scatter plot for one expression bin. x_plot is the asinh-scaled rG
# difference; exp.diff.wt is the RNA expression difference.
build_single_plot <- function(df_lvl, lvl, x_max, show_axis_text = TRUE, extra_bottom_margin = 18) {
  ggplot2::ggplot(df_lvl, ggplot2::aes(x = x_plot, y = exp.diff.wt)) +
    ggplot2::geom_point(
      color = facet_colors[[lvl]],
      size = 1.75,
      alpha = 0.2,
      shape = 19,
      stroke = 0
    ) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
    ggplot2::scale_x_continuous(limits = c(-x_max, x_max), expand = c(0.02, 0.02)) +
    ggplot2::scale_y_continuous(
      limits = c(-Y_LIM_FIXED, Y_LIM_FIXED),
      breaks = c(-5, -2.5, 0, 2.5, 5),
      minor_breaks = c(-3.75, -1.25, 1.25, 3.75),
      expand = c(0.02, 0.02)
    ) +
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.9, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.9, color = "black"),
      axis.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      legend.position = "none",
      axis.text = if (show_axis_text) ggplot2::element_text(size = 22) else ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 10, r = 30, b = extra_bottom_margin, l = 30, unit = "pt")
    )
}

# Add quadrant counts, Q2/Q3 and Q1/Q4 ratios, and plotted/removed counts to a
# single scatter panel.
build_annotated_plot <- function(df_lvl, lvl, x_max, tally_row) {
  base_plot <- build_single_plot(df_lvl, lvl, x_max, show_axis_text = TRUE, extra_bottom_margin = 42)

  q_text <- tibble::tibble(
    x = c(0.42 * x_max, -0.42 * x_max, -0.42 * x_max, 0.42 * x_max),
    y = c(0.68 * Y_LIM_FIXED, 0.68 * Y_LIM_FIXED, -0.68 * Y_LIM_FIXED, -0.68 * Y_LIM_FIXED),
    hjust = c(0, 1, 1, 0),
    label = c(
      paste0("Q1: ", tally_row$Q1),
      paste0("Q2: ", tally_row$Q2),
      paste0("Q3: ", tally_row$Q3),
      paste0("Q4: ", tally_row$Q4)
    )
  )

  ratio_text <- tibble::tibble(
    x = c(-0.58 * x_max, 0.58 * x_max),
    y = c(0.45, 0.45),
    label = c(
      paste0("Q2/Q3: ", format_ratio(tally_row$Q2, tally_row$Q3)),
      paste0("Q1/Q4: ", format_ratio(tally_row$Q1, tally_row$Q4))
    )
  )

  summary_label <- paste0("Plotted: ", tally_row$n_plotted, " | Removed: ", tally_row$n_removed_rG_both0)

  base_plot +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::geom_text(
      data = q_text,
      ggplot2::aes(x = x, y = y, label = label, hjust = hjust),
      inherit.aes = FALSE,
      size = 4.2,
      fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = ratio_text,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 4.2
    ) +
    ggplot2::labs(caption = summary_label) +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(size = 12, hjust = 0.5, margin = ggplot2::margin(t = 12))
    )
}

# Save one labeled, one no-label, and one annotated SVG for each expression bin.
save_scatter_set <- function(df_plot, quad_tally, comparison_tag, plot_dir, annotated_dir, x_max) {
  if (nrow(df_plot) == 0) return(invisible(NULL))

  for (lvl in facet_levels) {
    df_lvl <- df_plot %>% dplyr::filter(Exp.level == lvl)
    tally_row <- quad_tally %>% dplyr::filter(Exp.level == lvl)

    p <- build_single_plot(df_lvl, lvl, x_max, show_axis_text = TRUE)
    p_nolabs <- build_single_plot(df_lvl, lvl, x_max, show_axis_text = FALSE)
    p_annotated <- build_annotated_plot(df_lvl, lvl, x_max, tally_row)

    ggplot2::ggsave(file.path(plot_dir, paste0(comparison_tag, "_", lvl, ".svg")), p, device = "svg", width = 4, height = 4, units = "in")
    ggplot2::ggsave(file.path(plot_dir, paste0(comparison_tag, "_", lvl, "_nolabs.svg")), p_nolabs, device = "svg", width = 4, height = 4, units = "in")
    ggplot2::ggsave(file.path(annotated_dir, paste0(comparison_tag, "_", lvl, "_annotated.svg")), p_annotated, device = "svg", width = 4, height = 4, units = "in")
  }
}

# Save high-expression summary dot plots. The default SVGs use the compact
# figure-layout size. A second no-label version removes titles and axis text for
# assembly. The point_size column implements the manual diamond size adjustment.
save_dot_plots <- function(high_ratio_df, dot_dir, comparison_shapes, comparison_point_sizes, mode_label, condition_label) {
  high_ratio_df <- high_ratio_df %>%
    dplyr::mutate(
      comparison_label = factor(
        comparison_label,
        levels = c("WTvKO <=0.1", "WTvKO >2", "NCvTOP1 >2", "TOP1vWT >2")
      )
    )

  ratio_long_df <- high_ratio_df %>%
    tidyr::pivot_longer(
      cols = c(Q2_Q3_ratio, Q1_Q4_ratio),
      names_to = "ratio_type",
      values_to = "ratio_value"
    ) %>%
    dplyr::mutate(
      ratio_label = dplyr::case_when(
        ratio_type == "Q2_Q3_ratio" ~ "Q2/Q3",
        ratio_type == "Q1_Q4_ratio" ~ "Q1/Q4",
        TRUE ~ ratio_type
      ),
      ratio_label = factor(ratio_label, levels = c("Q2/Q3", "Q1/Q4")),
      x_pos = dplyr::case_when(
        ratio_label == "Q2/Q3" ~ as.numeric(comparison_label) - 0.09,
        ratio_label == "Q1/Q4" ~ as.numeric(comparison_label) + 0.09,
        TRUE ~ as.numeric(comparison_label)
      ),
      point_size = comparison_point_sizes[as.character(comparison_label)]
    )

  y_pad_ratio <- max(0.08, 0.08 * diff(range(ratio_long_df$ratio_value, na.rm = TRUE)))
  y_limits_ratio <- c(min(ratio_long_df$ratio_value, na.rm = TRUE) - y_pad_ratio, max(ratio_long_df$ratio_value, na.rm = TRUE) + y_pad_ratio * 1.8)

  p_ratio <- ggplot2::ggplot(ratio_long_df, ggplot2::aes(x = x_pos, y = ratio_value, shape = comparison_label, size = point_size)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4, color = "gray50") +
    ggplot2::geom_point(color = high) +
    ggplot2::geom_text(ggplot2::aes(label = ratio_label), vjust = -0.9, size = 3.4, color = high, show.legend = FALSE) +
    ggplot2::scale_shape_manual(values = comparison_shapes) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_x_continuous(breaks = seq_along(levels(high_ratio_df$comparison_label)), labels = levels(high_ratio_df$comparison_label)) +
    ggplot2::coord_cartesian(ylim = y_limits_ratio, clip = "off") +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.9, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.9, color = "black"),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, color = "black"),
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 18, r = 20, b = 20, l = 20, unit = "pt")
    ) +
    ggplot2::labs(title = paste(mode_label, "|", condition_label), x = "Comparison", y = "Quadrant ratio")

  p_ratio_nolabs <- p_ratio +
    ggplot2::labs(title = NULL, x = NULL, y = NULL) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank()
    )

  y_pad_delta <- max(0.05, 0.10 * diff(range(high_ratio_df$Q1_Q4_minus_Q2_Q3, na.rm = TRUE)))
  y_limits_delta <- c(min(high_ratio_df$Q1_Q4_minus_Q2_Q3, na.rm = TRUE) - y_pad_delta, max(high_ratio_df$Q1_Q4_minus_Q2_Q3, na.rm = TRUE) + y_pad_delta)

  high_ratio_df <- high_ratio_df %>%
    dplyr::mutate(point_size = comparison_point_sizes[as.character(comparison_label)])

  p_delta <- ggplot2::ggplot(high_ratio_df, ggplot2::aes(x = comparison_label, y = Q1_Q4_minus_Q2_Q3, shape = comparison_label, size = point_size)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "gray50") +
    ggplot2::geom_point(color = high) +
    ggplot2::scale_shape_manual(values = comparison_shapes) +
    ggplot2::scale_size_identity() +
    ggplot2::coord_cartesian(ylim = y_limits_delta, clip = "off") +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.9, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.9, color = "black"),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, color = "black"),
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 16, r = 20, b = 20, l = 20, unit = "pt")
    ) +
    ggplot2::labs(title = paste(mode_label, "|", condition_label), x = "Comparison", y = "Q1/Q4 ratio minus Q2/Q3 ratio")

  p_delta_nolabs <- p_delta +
    ggplot2::labs(title = NULL, x = NULL, y = NULL) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank()
    )

  legend_df <- tibble::tibble(
    comparison_label = factor(names(comparison_shapes), levels = names(comparison_shapes)),
    x = 1,
    y = rev(seq_along(comparison_shapes)),
    point_size = comparison_point_sizes[names(comparison_shapes)]
  )

  p_legend <- ggplot2::ggplot(
    legend_df,
    ggplot2::aes(x = x, y = y, shape = comparison_label, size = point_size)
  ) +
    ggplot2::geom_point(color = high) +
    ggplot2::geom_text(
      ggplot2::aes(label = as.character(comparison_label)),
      x = 1.22,
      hjust = 0,
      size = 4,
      color = "black"
    ) +
    ggplot2::scale_shape_manual(values = comparison_shapes) +
    ggplot2::scale_size_identity() +
    ggplot2::coord_cartesian(xlim = c(0.75, 2.6), ylim = c(0.5, length(comparison_shapes) + 0.5), clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(4, 8, 4, 8, unit = "pt"))

  p_legend_nolabs <- ggplot2::ggplot(
    legend_df,
    ggplot2::aes(x = x, y = y, shape = comparison_label, size = point_size)
  ) +
    ggplot2::geom_point(color = high) +
    ggplot2::scale_shape_manual(values = comparison_shapes) +
    ggplot2::scale_size_identity() +
    ggplot2::coord_cartesian(xlim = c(0.75, 1.25), ylim = c(0.5, length(comparison_shapes) + 0.5), clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(4, 8, 4, 8, unit = "pt"))

  ggplot2::ggsave(file.path(dot_dir, "high_expression_quadrant_ratios_dot_plot.svg"), p_ratio, device = "svg", width = 4.4, height = 4.8, units = "in")
  ggplot2::ggsave(file.path(dot_dir, "high_expression_quadrant_ratios_dot_plot_nolabs.svg"), p_ratio_nolabs, device = "svg", width = 4.4, height = 4.8, units = "in")
  ggplot2::ggsave(file.path(dot_dir, "high_expression_ratio_difference_dot_plot.svg"), p_delta, device = "svg", width = 4.4, height = 4.8, units = "in")
  ggplot2::ggsave(file.path(dot_dir, "high_expression_ratio_difference_dot_plot_nolabs.svg"), p_delta_nolabs, device = "svg", width = 4.4, height = 4.8, units = "in")
  ggplot2::ggsave(file.path(dot_dir, "high_expression_dot_plot_legend.svg"), p_legend, device = "svg", width = 1.9, height = 1.6, units = "in")
  ggplot2::ggsave(file.path(dot_dir, "high_expression_dot_plot_legend_nolabs.svg"), p_legend_nolabs, device = "svg", width = 0.45, height = 1.6, units = "in")

  readr::write_csv(high_ratio_df, file.path(dot_dir, "high_expression_ratio_summary.csv"))
}

# Construct the filter functions for the summary table. These functions are
# applied after expression binning, so bin membership is cutoff-independent.
build_summary_filter_modes <- function() {
  c(
    list(list(filter_tag = "unfiltered", threshold_dir = "unfiltered", threshold_value = NA_real_, cutoff_direction = "none", bias_filter = "none", fn = function(df) df)),
    lapply(LE_CUTOFFS, function(threshold) {
      list(filter_tag = paste0("le_", threshold_tag(threshold)), threshold_dir = paste0("le_", sprintf("%.2f", threshold)), threshold_value = threshold, cutoff_direction = "<=", bias_filter = paste0("abs bias <= ", threshold), fn = function(df) dplyr::filter(df, abs_bias <= threshold))
    }),
    lapply(GT_CUTOFFS, function(threshold) {
      list(filter_tag = paste0("gt_", threshold_tag(threshold)), threshold_dir = paste0("gt_", sprintf("%.2f", threshold)), threshold_value = threshold, cutoff_direction = ">", bias_filter = paste0("abs bias > ", threshold), fn = function(df) dplyr::filter(df, abs_bias > threshold))
    })
  )
}

# Load all source tables and merge them into one gene-symbol table. This function
# is shared by every plot and summary output so all downstream files use the same
# static-PC gene universe.
read_shared_tables <- function(symbol_map_tbl) {
  df_rna_raw <- readr::read_tsv(RNA_FILE, show_col_types = FALSE)
  required_rna_cols <- c("Geneid", "WT", "KO-T3-8", "KO-T3-17", "T3-8siRNA-NC_5+5pmol", "T3-8siRNA-TOP1_5+5pmol")
  missing_rna_cols <- setdiff(required_rna_cols, names(df_rna_raw))
  if (length(missing_rna_cols) > 0) stop("RNA file is missing required columns: ", paste(missing_rna_cols, collapse = ", "))

  # Treat zero TPM as missing before taking log2.
  df_rna_raw <- df_rna_raw %>% dplyr::mutate(dplyr::across(dplyr::all_of(required_rna_cols[-1]), ~ dplyr::na_if(., 0)))
  df_rna_raw$HEK293T_KO_AVG <- rowMeans(df_rna_raw[, c("KO-T3-8", "KO-T3-17")], na.rm = TRUE)

  # Raw log2TPM RNA table used for WTvKO and NCvTOP1 RNA differences, and for
  # WT expression reference values used in dynamic binning.
  df_rna <- df_rna_raw %>%
    dplyr::transmute(
      ENSEMBL = clean_ens(Geneid),
      `HEK293T-WT.rna` = log2(WT),
      `HEK293T-KO.rna` = log2(HEK293T_KO_AVG),
      `HEK293T-NC.rna` = log2(`T3-8siRNA-NC_5+5pmol`),
      `HEK293T-TOP1.rna` = log2(`T3-8siRNA-TOP1_5+5pmol`)
    ) %>%
    dplyr::filter(dplyr::if_any(dplyr::ends_with(".rna"), is.finite)) %>%
    dplyr::inner_join(symbol_map_tbl, by = "ENSEMBL") %>%
    dplyr::select(SYMBOL, dplyr::ends_with(".rna")) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(., na.rm = TRUE)), .groups = "drop")

  # Batch-normalized RNA table used for TOP1vWT RNA differences. WT/KO and
  # NC/TOP1 samples are treated as separate batches while preserving group means.
  rna_batch_matrix <- df_rna_raw %>%
    dplyr::transmute(
      ENSEMBL = clean_ens(Geneid),
      WT = log2(WT),
      KO_T3_8 = log2(`KO-T3-8`),
      KO_T3_17 = log2(`KO-T3-17`),
      NC = log2(`T3-8siRNA-NC_5+5pmol`),
      TOP1 = log2(`T3-8siRNA-TOP1_5+5pmol`)
    ) %>%
    dplyr::inner_join(symbol_map_tbl, by = "ENSEMBL") %>%
    dplyr::select(SYMBOL, WT, KO_T3_8, KO_T3_17, NC, TOP1) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(., na.rm = TRUE)), .groups = "drop")

  rna_batch_matrix_values <- as.matrix(rna_batch_matrix[, c("WT", "KO_T3_8", "KO_T3_17", "NC", "TOP1")])
  rownames(rna_batch_matrix_values) <- rna_batch_matrix$SYMBOL
  batch <- factor(c("wt_ko_batch", "wt_ko_batch", "wt_ko_batch", "sirna_batch", "sirna_batch"))
  group <- factor(c("WT", "KO", "KO", "NC", "TOP1"), levels = c("WT", "KO", "NC", "TOP1"))
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- sub("^group", "", colnames(design))

  rna_batch_matrix_corrected <- limma::removeBatchEffect(x = rna_batch_matrix_values, batch = batch, design = design)

  # Keep only genes with finite batch-normalized WT and TOP1 RNA values for the
  # final analysis condition.
  df_rna_batchnorm <- tibble::tibble(
    SYMBOL = rownames(rna_batch_matrix_corrected),
    `HEK293T-WT.batchnorm.rna` = rna_batch_matrix_corrected[, "WT"],
    `HEK293T-NC.batchnorm.rna` = rna_batch_matrix_corrected[, "NC"],
    `HEK293T-TOP1.batchnorm.rna` = rna_batch_matrix_corrected[, "TOP1"]
  ) %>%
    dplyr::mutate(batchnorm_keep = is.finite(`HEK293T-WT.batchnorm.rna`) & is.finite(`HEK293T-TOP1.batchnorm.rna`))

  # rG signal table. Condition-specific columns are averaged to one value per
  # condition, then collapsed to one row per symbol.
  df_rg_raw <- readr::read_tsv(RG_FILE, show_col_types = FALSE)
  if (!("id" %in% names(df_rg_raw))) stop("rG file missing required column: id")
  df_rg <- df_rg_raw %>%
    dplyr::mutate(
      SYMBOL = clean_sym(id),
      `HEK293T-WT.rG` = avg_cols(., "HEK293T-WT", "rG WT libraries"),
      `HEK293T-KO.rG` = avg_cols(., "HEK293T-RNASEH2A-KO-T3-(8|17)", "rG KO libraries"),
      `HEK293T-NC.rG` = avg_cols(., "HEK293T-T3-8-siRNA-NC-5\\+5pmol", "rG NC libraries"),
      `HEK293T-TOP1.rG` = avg_cols(., "HEK293T-T3-8-siRNA-TOP1-5\\+5pmol", "rG TOP1 libraries")
    ) %>%
    dplyr::transmute(SYMBOL, `HEK293T-WT.rG`, `HEK293T-KO.rG`, `HEK293T-NC.rG`, `HEK293T-TOP1.rG`) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(., na.rm = TRUE)), .groups = "drop")

  # rNMP bias is computed as opposite-strand signal minus same-strand signal for
  # KO and TOP1. These bias values are what the cutoffs threshold.
  df_same_raw <- readr::read_tsv(RN_SAME, show_col_types = FALSE)
  df_opp_raw <- readr::read_tsv(RN_OPP, show_col_types = FALSE)
  if (!("id" %in% names(df_same_raw))) stop("rN SAME file missing required column: id")
  if (!("id" %in% names(df_opp_raw))) stop("rN OPP file missing required column: id")

  df_same <- df_same_raw %>%
    dplyr::mutate(SYMBOL = clean_sym(id)) %>%
    dplyr::transmute(
      SYMBOL,
      `HEK293T-KO.same` = rowMeans(dplyr::across(c(`HEK293T-RNASEH2A-KO-T3-8`, `HEK293T-RNASEH2A-KO-T3-17`)), na.rm = TRUE),
      `HEK293T-NC.same` = `HEK293T-T3-8-siRNA-NC-5+5pmol`,
      `HEK293T-TOP1.same` = `HEK293T-T3-8-siRNA-TOP1-5+5pmol`
    ) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(., na.rm = TRUE)), .groups = "drop")

  df_opp <- df_opp_raw %>%
    dplyr::mutate(SYMBOL = clean_sym(id)) %>%
    dplyr::transmute(
      SYMBOL,
      `HEK293T-KO.opp` = rowMeans(dplyr::across(c(`HEK293T-RNASEH2A-KO-T3-8`, `HEK293T-RNASEH2A-KO-T3-17`)), na.rm = TRUE),
      `HEK293T-NC.opp` = `HEK293T-T3-8-siRNA-NC-5+5pmol`,
      `HEK293T-TOP1.opp` = `HEK293T-T3-8-siRNA-TOP1-5+5pmol`
    ) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(., na.rm = TRUE)), .groups = "drop")

  df_rnmp <- dplyr::inner_join(df_opp, df_same, by = "SYMBOL") %>%
    dplyr::mutate(
      `HEK293T-KO.rNMP.bias` = `HEK293T-KO.opp` - `HEK293T-KO.same`,
      `HEK293T-TOP1.rNMP.bias` = `HEK293T-TOP1.opp` - `HEK293T-TOP1.same`
    ) %>%
    dplyr::select(SYMBOL, `HEK293T-KO.rNMP.bias`, `HEK293T-TOP1.rNMP.bias`)

  # Final joined table used downstream. Every comparison starts from this same
  # symbol-level data frame, then filters to batchnorm_keep before binning.
  merged_core <- df_rna %>%
    dplyr::inner_join(df_rna_batchnorm, by = "SYMBOL") %>%
    dplyr::inner_join(df_rg, by = "SYMBOL") %>%
    dplyr::inner_join(df_rnmp, by = "SYMBOL")

  list(
    merged_core = merged_core,
    batchnorm_keep_symbols = df_rna_batchnorm %>% dplyr::filter(batchnorm_keep) %>% dplyr::pull(SYMBOL),
    symbol_map_tbl = symbol_map_tbl,
    merged_symbols = merged_core$SYMBOL,
    map_stats = tibble::tibble(
      n_symbol_map_rows = nrow(symbol_map_tbl),
      n_merged_symbols = nrow(merged_core),
      n_batchnorm_keep_symbols = sum(df_rna_batchnorm$batchnorm_keep, na.rm = TRUE)
    )
  )
}

# Figure-panel comparisons. x_* columns define rG difference direction;
# y_* columns define RNA difference direction; bias_col defines the threshold.
comparisons_for_plots <- list(
  list(id = "wt_v_ko_le_0_1", comparison_label = "WTvKO <=0.1", base_comparison = "WTvKO", file_tag = "WTvKO_le0_1", x_high_col = "HEK293T-KO.rG", x_low_col = "HEK293T-WT.rG", y_high_col = "HEK293T-WT.rna", y_low_col = "HEK293T-KO.rna", expr_ref_col = "HEK293T-WT.rna", bias_col = "HEK293T-KO.rNMP.bias", bias_filter_label = "abs bias <= 0.1", filter_fn = function(df) dplyr::filter(df, abs_bias <= 0.1)),
  list(id = "wt_v_ko_gt_2", comparison_label = "WTvKO >2", base_comparison = "WTvKO", file_tag = "WTvKO_gt2", x_high_col = "HEK293T-KO.rG", x_low_col = "HEK293T-WT.rG", y_high_col = "HEK293T-WT.rna", y_low_col = "HEK293T-KO.rna", expr_ref_col = "HEK293T-WT.rna", bias_col = "HEK293T-KO.rNMP.bias", bias_filter_label = "abs bias > 2", filter_fn = function(df) dplyr::filter(df, abs_bias > 2)),
  list(id = "nc_v_top1_gt_2", comparison_label = "NCvTOP1 >2", base_comparison = "NCvTOP1", file_tag = "NCvTOP1_gt2", x_high_col = "HEK293T-TOP1.rG", x_low_col = "HEK293T-NC.rG", y_high_col = "HEK293T-NC.rna", y_low_col = "HEK293T-TOP1.rna", expr_ref_col = "HEK293T-WT.rna", bias_col = "HEK293T-TOP1.rNMP.bias", bias_filter_label = "abs bias > 2", filter_fn = function(df) dplyr::filter(df, abs_bias > 2)),
  list(id = "top1_v_wt_gt_2", comparison_label = "TOP1vWT >2", base_comparison = "TOP1vWT", file_tag = "TOP1vWT_gt2_batchnormalized", x_high_col = "HEK293T-TOP1.rG", x_low_col = "HEK293T-WT.rG", y_high_col = "HEK293T-WT.batchnorm.rna", y_low_col = "HEK293T-TOP1.batchnorm.rna", expr_ref_col = "HEK293T-WT.rna", bias_col = "HEK293T-TOP1.rNMP.bias", bias_filter_label = "abs bias > 2", filter_fn = function(df) dplyr::filter(df, abs_bias > 2))
)

# Summary comparisons use the same directions but are expanded across all
# thresholds from build_summary_filter_modes().
summary_comparisons <- list(
  list(id = "wt_v_ko", comparison_label = "WTvKO", x_high_col = "HEK293T-KO.rG", x_low_col = "HEK293T-WT.rG", y_high_col = "HEK293T-WT.rna", y_low_col = "HEK293T-KO.rna", expr_ref_col = "HEK293T-WT.rna", bias_col = "HEK293T-KO.rNMP.bias"),
  list(id = "nc_v_top1", comparison_label = "NCvTOP1", x_high_col = "HEK293T-TOP1.rG", x_low_col = "HEK293T-NC.rG", y_high_col = "HEK293T-NC.rna", y_low_col = "HEK293T-TOP1.rna", expr_ref_col = "HEK293T-WT.rna", bias_col = "HEK293T-TOP1.rNMP.bias"),
  list(id = "wt_v_top1", comparison_label = "WTvTOP1", x_high_col = "HEK293T-TOP1.rG", x_low_col = "HEK293T-WT.rG", y_high_col = "HEK293T-WT.batchnorm.rna", y_low_col = "HEK293T-TOP1.batchnorm.rna", expr_ref_col = "HEK293T-WT.rna", bias_col = "HEK293T-TOP1.rNMP.bias")
)

comparison_shapes <- c(`WTvKO <=0.1` = 16, `WTvKO >2` = 17, `NCvTOP1 >2` = 15, `TOP1vWT >2` = 18)

# Shape 18 is a diamond and appears visually smaller than the other point
# shapes, so TOP1vWT is drawn one size larger in dot plots.
# Dot plot symbol sizes. The diamond shape appears smaller than circles and
# triangles at the same ggplot size value, so TOP1vWT keeps a manual size boost.
# These are larger than the original final values of 6 and 7, while keeping
# the manually boosted diamond close in apparent area to the other symbols.
comparison_point_sizes <- c(`WTvKO <=0.1` = 12, `WTvKO >2` = 12, `NCvTOP1 >2` = 12, `TOP1vWT >2` = 14)
summary_filter_modes <- build_summary_filter_modes()

# Only the static PC mapping mode belongs in the final output.
mapping_modes <- list(
  list(
    mode_id = "static_gene_mapping_pc",
    mode_label = "static gene mapping pc",
    build_map = function(ensembl_ids) build_symbol_map_static(ensembl_ids, SYMBOL_MAP_FILE_PC)
  )
)

# Final condition: keep genes with finite batch-normalized WT/TOP1 RNA, then
# assign dynamic bins before threshold cutoffs.
conditions <- list(
  list(condition_id = "dynamic_pre_filter_with_batchnorm_keep", condition_label = "dynamic pre-filter quantile bins + batchnorm keep", apply_batchnorm_keep = TRUE)
)

all_mode_ratio_summaries <- list()
all_mode_map_stats <- list()
shared_by_mode <- list()

for (map_mode in mapping_modes) {
  # Build the static PC gene map once, then use it for every downstream table.
  symbol_map_tbl <- map_mode$build_map(readr::read_tsv(RNA_FILE, show_col_types = FALSE)$Geneid)
  shared <- read_shared_tables(symbol_map_tbl)
  all_mode_map_stats[[map_mode$mode_id]] <- shared$map_stats %>% dplyr::mutate(mode_id = map_mode$mode_id, mode_label = map_mode$mode_label)
  shared_by_mode[[map_mode$mode_id]] <- shared

  for (condition in conditions) {
    condition_dir <- file.path(OUT_DIR, map_mode$mode_id, condition$condition_id)
    plot_dir <- file.path(condition_dir, "plots")
    annotated_dir <- file.path(condition_dir, "annotated")
    dot_dir <- file.path(condition_dir, "dot_plots")
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(annotated_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(dot_dir, showWarnings = FALSE, recursive = TRUE)

    all_quad_tallies <- list()
    all_stage_summaries <- list()
    all_removed_gene_tables <- list()
    all_scatter_plot_data <- list()

    # Main figure-panel loop. Each comparison gets scatter plots, annotated
    # scatter plots, quadrant tallies, stage summaries, and removed-gene tables.
    for (cmp in comparisons_for_plots) {
      cmp_source_df <- shared$merged_core
      if (condition$apply_batchnorm_keep) {
        # Apply batchnorm keep before dynamic binning. This defines the analysis
        # universe used for Low/Medium/High quartiles.
        cmp_source_df <- cmp_source_df %>% dplyr::filter(SYMBOL %in% shared$batchnorm_keep_symbols)
      }

      # Compute comparison-specific differences and assign expression bins before
      # applying the comparison's bias threshold.
      cmp_df_base <- cmp_source_df %>%
        dplyr::mutate(
          rG.diff = .data[[cmp$x_high_col]] - .data[[cmp$x_low_col]],
          exp.diff.wt = .data[[cmp$y_high_col]] - .data[[cmp$y_low_col]],
          expr_ref = .data[[cmp$expr_ref_col]],
          abs_bias = abs(.data[[cmp$bias_col]])
        ) %>%
        dplyr::filter(is.finite(rG.diff), is.finite(exp.diff.wt), is.finite(expr_ref), is.finite(abs_bias)) %>%
        assign_dynamic_prefilter_bins("expr_ref")

      # Thresholding happens after binning. Therefore the bin labels are not
      # recalculated separately for <=0.1, >2, or any other cutoff.
      filtered_df <- cmp$filter_fn(cmp_df_base)
      removed_rg_both0 <- filtered_df %>%
        dplyr::filter(is.finite(.data[[cmp$x_high_col]]), is.finite(.data[[cmp$x_low_col]]), .data[[cmp$x_high_col]] == 0, .data[[cmp$x_low_col]] == 0)

      # Remove genes with zero rG in both compared conditions from the plot and
      # quadrant counts, while preserving a removal tally for reporting.
      df_plot <- filtered_df %>%
        dplyr::filter(!(.data[[cmp$x_high_col]] == 0 & .data[[cmp$x_low_col]] == 0)) %>%
        dplyr::mutate(x_plot = asinh(rG.diff / ASINH_SCALE)) %>%
        dplyr::filter(is.finite(x_plot), is.finite(exp.diff.wt))

      removed_tally <- removed_rg_both0 %>%
        dplyr::count(Exp.level, name = "n_removed_rG_both0") %>%
        dplyr::right_join(tibble::tibble(Exp.level = factor(facet_levels, levels = facet_levels)), by = "Exp.level") %>%
        dplyr::mutate(n_removed_rG_both0 = dplyr::coalesce(n_removed_rG_both0, 0L))

      quad_tally <- tally_quadrants(df_plot) %>%
        dplyr::left_join(removed_tally, by = "Exp.level") %>%
        dplyr::mutate(
          mode_id = map_mode$mode_id,
          mode_label = map_mode$mode_label,
          condition_id = condition$condition_id,
          condition_label = condition$condition_label,
          comparison = cmp$id,
          comparison_label = cmp$comparison_label,
          base_comparison = cmp$base_comparison,
          bias_filter = cmp$bias_filter_label,
          n_removed_rG_both0 = dplyr::coalesce(n_removed_rG_both0, 0L),
          n_plotted = n_total,
          n_plotted_check = check_sum,
          n_pre_rG_zero_removal = n_plotted + n_removed_rG_both0,
          Q2_Q3_ratio = safe_ratio_num(Q2, Q3),
          Q1_Q4_ratio = safe_ratio_num(Q1, Q4),
          Q1_Q4_minus_Q2_Q3 = Q1_Q4_ratio - Q2_Q3_ratio
        )

      stage_summary <- tibble::tibble(
        mode_id = map_mode$mode_id,
        mode_label = map_mode$mode_label,
        condition_id = condition$condition_id,
        condition_label = condition$condition_label,
        comparison = cmp$id,
        comparison_label = cmp$comparison_label,
        base_comparison = cmp$base_comparison,
        bias_filter = cmp$bias_filter_label,
        n_after_shared_filters = nrow(cmp_df_base),
        n_after_cutoff_filter = nrow(filtered_df),
        n_removed_rG_both0 = nrow(removed_rg_both0),
        n_plotted = nrow(df_plot)
      )

      removed_gene_table <- removed_rg_both0 %>%
        dplyr::transmute(
          mode_id = map_mode$mode_id,
          mode_label = map_mode$mode_label,
          condition_id = condition$condition_id,
          condition_label = condition$condition_label,
          comparison = cmp$id,
          comparison_label = cmp$comparison_label,
          SYMBOL,
          Exp.level,
          x_high_value = .data[[cmp$x_high_col]],
          x_low_value = .data[[cmp$x_low_col]],
          rG.diff,
          exp.diff.wt,
          expr_ref,
          abs_bias
        )

      all_quad_tallies[[cmp$id]] <- quad_tally
      all_stage_summaries[[cmp$id]] <- stage_summary
      all_removed_gene_tables[[cmp$id]] <- removed_gene_table
      all_scatter_plot_data[[cmp$id]] <- list(df_plot = df_plot, quad_tally = quad_tally, file_tag = cmp$file_tag)

      readr::write_csv(quad_tally, file.path(condition_dir, paste0(cmp$file_tag, "_quadrant_tally.csv")))
      readr::write_csv(stage_summary, file.path(condition_dir, paste0(cmp$file_tag, "_stage_summary.csv")))
      readr::write_csv(removed_gene_table, file.path(condition_dir, paste0(cmp$file_tag, "_removed_rG_both0_genes.csv")))
    }

    # Use a common symmetric x-axis limit for all scatter panels so visual scale
    # is comparable across comparisons and expression bins.
    global_x_max <- all_scatter_plot_data %>% purrr::map("df_plot") %>% dplyr::bind_rows() %>% dplyr::pull(x_plot) %>% abs() %>% max(na.rm = TRUE)
    if (!is.finite(global_x_max) || global_x_max == 0) global_x_max <- 1

    purrr::walk(all_scatter_plot_data, ~ save_scatter_set(.x$df_plot, .x$quad_tally, .x$file_tag, plot_dir, annotated_dir, global_x_max))

    all_quad_tallies_df <- dplyr::bind_rows(all_quad_tallies)
    all_stage_summaries_df <- dplyr::bind_rows(all_stage_summaries)
    all_removed_gene_tables_df <- dplyr::bind_rows(all_removed_gene_tables)

    readr::write_csv(all_quad_tallies_df, file.path(condition_dir, "all_comparisons_quadrant_tallies.csv"))
    readr::write_csv(all_stage_summaries_df, file.path(condition_dir, "all_comparisons_stage_summary.csv"))
    readr::write_csv(all_removed_gene_tables_df, file.path(condition_dir, "all_comparisons_removed_rG_both0_genes.csv"))

    # Dot plots summarize only the High-expression rows from the four main
    # figure-panel comparisons.
    high_ratio_df <- all_quad_tallies_df %>%
      dplyr::filter(Exp.level == "High") %>%
      dplyr::select(comparison_label, Q2_Q3_ratio, Q1_Q4_ratio, Q1_Q4_minus_Q2_Q3) %>%
      dplyr::distinct()
    save_dot_plots(high_ratio_df, dot_dir, comparison_shapes, comparison_point_sizes, map_mode$mode_label, condition$condition_label)

    # Comprehensive summary output. This repeats the same calculation pattern as
    # the plotted comparisons but expands each comparison over all thresholds.
    summary_ratio_rows <- list()
    for (cmp in summary_comparisons) {
      cmp_source_df <- shared$merged_core
      if (condition$apply_batchnorm_keep) {
        cmp_source_df <- cmp_source_df %>% dplyr::filter(SYMBOL %in% shared$batchnorm_keep_symbols)
      }

      # Dynamic bins are assigned once before the threshold loop. This is the
      # critical pre-filter behavior: changing threshold does not change bins.
      cmp_df_base <- cmp_source_df %>%
        dplyr::mutate(
          rG.diff = .data[[cmp$x_high_col]] - .data[[cmp$x_low_col]],
          exp.diff.wt = .data[[cmp$y_high_col]] - .data[[cmp$y_low_col]],
          expr_ref = .data[[cmp$expr_ref_col]],
          abs_bias = abs(.data[[cmp$bias_col]])
        ) %>%
        dplyr::filter(is.finite(rG.diff), is.finite(exp.diff.wt), is.finite(expr_ref), is.finite(abs_bias)) %>%
        assign_dynamic_prefilter_bins("expr_ref")

      for (filter_mode in summary_filter_modes) {
        # Apply one threshold to the already-binned comparison table.
        filtered_df <- filter_mode$fn(cmp_df_base)
        if (nrow(filtered_df) == 0) {
          summary_ratio_rows[[paste(cmp$id, filter_mode$filter_tag, sep = "__")]] <- tibble::tibble(
            mode_id = map_mode$mode_id,
            mode_label = map_mode$mode_label,
            condition_id = condition$condition_id,
            condition_label = condition$condition_label,
            comparison = cmp$id,
            comparison_label = cmp$comparison_label,
            filter_tag = filter_mode$filter_tag,
            threshold_dir = filter_mode$threshold_dir,
            threshold_value = filter_mode$threshold_value,
            cutoff_direction = filter_mode$cutoff_direction,
            bias_filter = filter_mode$bias_filter,
            Exp.level = factor(facet_levels, levels = facet_levels),
            n_plotted = 0L,
            n_removed_rG_both0 = 0L,
            Q1 = 0L, Q2 = 0L, Q3 = 0L, Q4 = 0L,
            Q2_Q3_ratio = NA_real_,
            Q1_Q4_ratio = NA_real_,
            Q1_Q4_minus_Q2_Q3 = NA_real_
          )
          next
        }

        # Mirror the plotted-panel handling of rG-both-zero genes.
        removed_rg_both0 <- filtered_df %>%
          dplyr::filter(is.finite(.data[[cmp$x_high_col]]), is.finite(.data[[cmp$x_low_col]]), .data[[cmp$x_high_col]] == 0, .data[[cmp$x_low_col]] == 0)
        df_plot <- filtered_df %>%
          dplyr::filter(!(.data[[cmp$x_high_col]] == 0 & .data[[cmp$x_low_col]] == 0)) %>%
          dplyr::mutate(x_plot = asinh(rG.diff / ASINH_SCALE)) %>%
          dplyr::filter(is.finite(x_plot), is.finite(exp.diff.wt))
        removed_tally <- removed_rg_both0 %>%
          dplyr::count(Exp.level, name = "n_removed_rG_both0") %>%
          dplyr::right_join(tibble::tibble(Exp.level = factor(facet_levels, levels = facet_levels)), by = "Exp.level") %>%
          dplyr::mutate(n_removed_rG_both0 = dplyr::coalesce(n_removed_rG_both0, 0L))

        quad_tally <- tally_quadrants(df_plot) %>%
          dplyr::left_join(removed_tally, by = "Exp.level") %>%
          dplyr::mutate(
            mode_id = map_mode$mode_id,
            mode_label = map_mode$mode_label,
            condition_id = condition$condition_id,
            condition_label = condition$condition_label,
            comparison = cmp$id,
            comparison_label = cmp$comparison_label,
            filter_tag = filter_mode$filter_tag,
            threshold_dir = filter_mode$threshold_dir,
            threshold_value = filter_mode$threshold_value,
            cutoff_direction = filter_mode$cutoff_direction,
            bias_filter = filter_mode$bias_filter,
            n_plotted = n_total,
            n_removed_rG_both0 = dplyr::coalesce(n_removed_rG_both0, 0L),
            Q2_Q3_ratio = safe_ratio_num(Q2, Q3),
            Q1_Q4_ratio = safe_ratio_num(Q1, Q4),
            Q1_Q4_minus_Q2_Q3 = Q1_Q4_ratio - Q2_Q3_ratio
          ) %>%
          dplyr::select(mode_id, mode_label, condition_id, condition_label, comparison, comparison_label, filter_tag, threshold_dir, threshold_value, cutoff_direction, bias_filter, Exp.level, n_plotted, n_removed_rG_both0, Q1, Q2, Q3, Q4, Q2_Q3_ratio, Q1_Q4_ratio, Q1_Q4_minus_Q2_Q3)

        summary_ratio_rows[[paste(cmp$id, filter_mode$filter_tag, sep = "__")]] <- quad_tally
      }
    }

    summary_ratio_df <- dplyr::bind_rows(summary_ratio_rows)
    readr::write_csv(summary_ratio_df, file.path(condition_dir, "all_comparisons_ratio_summary.csv"))
    all_mode_ratio_summaries[[paste(map_mode$mode_id, condition$condition_id, sep = "__")]] <- summary_ratio_df
  }
}

readr::write_csv(dplyr::bind_rows(all_mode_ratio_summaries), file.path(OUT_DIR, "all_modes_ratio_summary.csv"))
readr::write_csv(dplyr::bind_rows(all_mode_map_stats), file.path(OUT_DIR, "mapping_mode_gene_counts.csv"))
