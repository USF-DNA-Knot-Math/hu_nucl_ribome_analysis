library(tidyverse)
library(ggpubr)
library(parallel)
library(grid)

# -----------------------------
# User-adjustable plot settings
# -----------------------------
promoter_ylims <- c(-1.0, 1.4)
gene_body_ylims <- c(-0.5, 0.61)

output_dir <- "Fig3E/050526"


# -----------------------------
# Sample mapping
# -----------------------------
library_to_celltype <- function(lib_name) {
  case_when(
    grepl("CD4T",                 lib_name) ~ "CD4T",
    grepl("hESC-H9",              lib_name) ~ "hESC-H9",
    grepl("RNASEH2A-KO-T3-8",     lib_name) ~ "HEK293T-RNASEH2A-KO-T3-8",
    grepl("RNASEH2A-KO-T3-17",    lib_name) ~ "HEK293T-RNASEH2A-KO-T3-17",
    grepl("HEK293T-WT",           lib_name) ~ "HEK293T-WT",
    grepl("siRNA-NC-5\\+5pmol",   lib_name) ~ "HEK293T-T3-8-siRNA-NC-5+5pmol",
    grepl("siRNA-TOP1-5\\+5pmol", lib_name) ~ "HEK293T-T3-8-siRNA-TOP1-5+5pmol",
    TRUE ~ NA_character_
  )
}

celltype_order <- c(
  "CD4T",
  "hESC-H9",
  "HEK293T-WT",
  "HEK293T-RNASEH2A-KO-T3-8",
  "HEK293T-RNASEH2A-KO-T3-17",
  "HEK293T-T3-8-siRNA-NC-5+5pmol",
  "HEK293T-T3-8-siRNA-TOP1-5+5pmol"
)

celltype_labels <- c(
  "CD4+T",
  "hESC-H9",
  "HEK293T",
  "RNH2KO T3-8",
  "RNH2KO T3-17",
  "siRNA-NC 5+5pmol",
  "siRNA-TOP1 5+5pmol"
)

color_map <- c("Low" = "#F8CF7D", "Mod" = "#F6A5A5", "High" = "#CC7F85")
point_outline_map <- c("Low" = "#F0B02F", "Mod" = "#F27272", "High" = "#BD434D")
comparisons_list <- list(c("Low", "Mod"), c("Low", "High"))
dodge_offsets <- c("Low" = -0.267, "Mod" = 0, "High" = 0.267)
point_outline_alpha <- 1
point_fill_alpha <- 0.5
point_stroke <- 0.45
point_size <- 1.9

outline_from_fill <- function(fill_col) {
  matched_name <- names(color_map)[match(fill_col, color_map)]

  if (length(matched_name) == 0 || is.na(matched_name)) {
    return(fill_col)
  }

  point_outline_map[[matched_name]]
}

draw_key_boxplot_with_point <- function(data, params, size) {
  box_data <- data
  box_data$colour <- "black"
  box_data$fill <- dplyr::coalesce(data$fill, data$colour, "grey80")
  box_data$alpha <- 1
  box_grob <- ggplot2::draw_key_boxplot(box_data, params, size)

  point_fill <- dplyr::coalesce(data$fill, data$colour, "black")
  point_colour <- outline_from_fill(point_fill)

  fill_data <- data
  fill_data$shape <- 16
  fill_data$size <- point_size
  fill_data$stroke <- 0
  fill_data$colour <- point_fill
  fill_data$fill <- NA
  fill_data$alpha <- point_fill_alpha

  outline_data <- data
  outline_data$shape <- 21
  outline_data$size <- point_size
  outline_data$stroke <- point_stroke
  outline_data$colour <- point_colour
  outline_data$fill <- NA
  outline_data$alpha <- point_outline_alpha

  fill_grob <- ggplot2::draw_key_point(
    fill_data,
    params,
    size
  )

  outline_grob <- ggplot2::draw_key_point(
    outline_data,
    params,
    size
  )

  grid::grobTree(box_grob, fill_grob, outline_grob)
}

extract_legend <- function(plot_obj) {
  plot_grob <- ggplotGrob(plot_obj)
  guide_index <- which(vapply(plot_grob$grobs, function(x) x$name, character(1)) == "guide-box")

  if (length(guide_index) == 0) {
    stop("No legend found in plot.")
  }

  plot_grob$grobs[[guide_index[1]]]
}

save_vertical_legend <- function(output_file, show_text = TRUE) {
  legend_df <- tibble(
    expr_group = factor(c("Low", "Mod", "High"), levels = c("Low", "Mod", "High")),
    x = 1,
    y = c(3, 2, 1),
    point_fill_col = unname(color_map[c("Low", "Mod", "High")]),
    point_outline_col = unname(point_outline_map[c("Low", "Mod", "High")])
  )

  legend_plot <- ggplot(legend_df, aes(x = x, y = y, fill = expr_group)) +
    geom_boxplot(
      aes(group = expr_group),
      width = 0.35,
      color = "black",
      key_glyph = draw_key_boxplot_with_point,
      show.legend = TRUE
    ) +
    scale_fill_manual(
      values = color_map,
      aesthetics = c("fill")
    ) +
    guides(
      fill = guide_legend(ncol = 1, byrow = TRUE),
      color = "none"
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = if (show_text) element_text(size = 12) else element_blank(),
      legend.key.width = unit(0.9, "cm"),
      legend.key.height = unit(0.7, "cm")
    ) +
    labs(fill = NULL)

  legend_grob <- extract_legend(legend_plot)

  grDevices::svg(filename = output_file, width = 2.0, height = 2.8)
  grid::grid.newpage()
  grid::grid.draw(legend_grob)
  dev.off()

  invisible(output_file)
}

# -----------------------------
# IO helpers
# -----------------------------
read_ef <- function(file) {
  df <- readr::read_tsv(file, col_types = cols())

  if (!"id" %in% colnames(df)) {
    stop(sprintf("Missing 'id' column in %s", file))
  }

  df <- dplyr::rename(df, gene_id = id)
  library_cols <- df %>% dplyr::select(where(is.numeric)) %>% colnames()

  df %>%
    dplyr::select(gene_id, all_of(library_cols)) %>%
    pivot_longer(-gene_id, names_to = "library", values_to = "EF") %>%
    mutate(cell_type = library_to_celltype(library))
}

read_bed_ids <- function(file, gene_col = 4) {
  df <- readr::read_tsv(file, col_names = FALSE, col_types = cols())
  unique(df[[gene_col]])
}

# -----------------------------
# EF data
# -----------------------------
promoter_opp_file <- "Tyler_rNMPs/hg38_refGene_TSS_downstream_0_1kb_1000bp/rN/out/hg38_refGene_TSS_downstream_0_1kb_1000bp_rN_opp_EF_regions_EF.tsv"
promoter_same_file <- "Tyler_rNMPs/hg38_refGene_TSS_downstream_0_1kb_1000bp/rN/out/hg38_refGene_TSS_downstream_0_1kb_1000bp_rN_same_EF_regions_EF.tsv"
genebody_opp_file <- "Tyler_rNMPs/hg38_refGene_collapsed/rN/out/hg38_refGene_collapsed_rN_opp_EF_regions_EF.tsv"
genebody_same_file <- "Tyler_rNMPs/hg38_refGene_collapsed/rN/out/hg38_refGene_collapsed_rN_same_EF_regions_EF.tsv"

ef_long <- bind_rows(
  read_ef(promoter_opp_file) %>%
    mutate(strand = "opp", region = "Promoter"),
  read_ef(promoter_same_file) %>%
    mutate(strand = "same", region = "Promoter"),
  read_ef(genebody_opp_file) %>%
    mutate(strand = "opp", region = "GeneBody"),
  read_ef(genebody_same_file) %>%
    mutate(strand = "same", region = "GeneBody")
) %>%
  filter(cell_type %in% celltype_order)

# -----------------------------
# Expression groups
# -----------------------------
expr_lists <- list(
  "CD4T" = list(
    Low  = read_bed_ids("Expression-level-lists/bowtie/CD4TLow_pc.bed", gene_col = 4),
    Mod  = read_bed_ids("Expression-level-lists/bowtie/CD4TMod_pc.bed", gene_col = 4),
    High = read_bed_ids("Expression-level-lists/bowtie/CD4THigh_pc.bed", gene_col = 4)
  ),
  "hESC-H9" = list(
    Low  = read_bed_ids("Expression-level-lists/bowtie/H9Low_pc.bed", gene_col = 4),
    Mod  = read_bed_ids("Expression-level-lists/bowtie/H9Mod_pc.bed", gene_col = 4),
    High = read_bed_ids("Expression-level-lists/bowtie/H9High_pc.bed", gene_col = 4)
  ),
  "HEK293T-WT" = list(
    Low  = read_bed_ids("Expression-level-lists/bowtie/WTLow_pc.bed", gene_col = 4),
    Mod  = read_bed_ids("Expression-level-lists/bowtie/WTMod_pc.bed", gene_col = 4),
    High = read_bed_ids("Expression-level-lists/bowtie/WTHigh_pc.bed", gene_col = 4)
  ),
  "HEK293T-RNASEH2A-KO-T3-8" = list(
    Low  = read_bed_ids("Expression-level-lists/bowtie/KO-T3-8Low_pc.bed", gene_col = 4),
    Mod  = read_bed_ids("Expression-level-lists/bowtie/KO-T3-8Mod_pc.bed", gene_col = 4),
    High = read_bed_ids("Expression-level-lists/bowtie/KO-T3-8High_pc.bed", gene_col = 4)
  ),
  "HEK293T-RNASEH2A-KO-T3-17" = list(
    Low  = read_bed_ids("Expression-level-lists/bowtie/KO-T3-17Low_pc.bed", gene_col = 4),
    Mod  = read_bed_ids("Expression-level-lists/bowtie/KO-T3-17Mod_pc.bed", gene_col = 4),
    High = read_bed_ids("Expression-level-lists/bowtie/KO-T3-17High_pc.bed", gene_col = 4)
  ),
  "HEK293T-T3-8-siRNA-NC-5+5pmol" = list(
    Low  = read_bed_ids("Expression-level-lists/bowtie/T3-8siRNA-NC_5_5pmolLow_pc.bed", gene_col = 4),
    Mod  = read_bed_ids("Expression-level-lists/bowtie/T3-8siRNA-NC_5_5pmolMod_pc.bed", gene_col = 4),
    High = read_bed_ids("Expression-level-lists/bowtie/T3-8siRNA-NC_5_5pmolHigh_pc.bed", gene_col = 4)
  ),
  "HEK293T-T3-8-siRNA-TOP1-5+5pmol" = list(
    Low  = read_bed_ids("Expression-level-lists/bowtie/T3-8siRNA-TOP1_5_5pmolLow_pc.bed", gene_col = 4),
    Mod  = read_bed_ids("Expression-level-lists/bowtie/T3-8siRNA-TOP1_5_5pmolMod_pc.bed", gene_col = 4),
    High = read_bed_ids("Expression-level-lists/bowtie/T3-8siRNA-TOP1_5_5pmolHigh_pc.bed", gene_col = 4)
  )
)

lookup_df <- bind_rows(lapply(names(expr_lists), function(ct) {
  tibble(
    cell_type = ct,
    gene_id = c(expr_lists[[ct]]$Low, expr_lists[[ct]]$Mod, expr_lists[[ct]]$High),
    expr_group = c(
      rep("Low", length(expr_lists[[ct]]$Low)),
      rep("Mod", length(expr_lists[[ct]]$Mod)),
      rep("High", length(expr_lists[[ct]]$High))
    )
  )
}))

ef_annotated <- ef_long %>%
  left_join(lookup_df, by = c("cell_type", "gene_id")) %>%
  filter(!is.na(expr_group))

# -----------------------------
# Strand bias
# -----------------------------
bias_df <- ef_annotated %>%
  group_by(library, cell_type, expr_group, region, strand) %>%
  summarise(mean_EF = mean(EF, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = strand, values_from = mean_EF) %>%
  filter(!is.na(opp), !is.na(same)) %>%
  mutate(
    bias = opp - same,
    cell_type = factor(cell_type, levels = celltype_order, labels = celltype_labels),
    expr_group = factor(expr_group, levels = c("Low", "Mod", "High")),
    region = factor(region, levels = c("Promoter", "GeneBody"))
  )

compute_signif_df <- function(region_df) {
  region_df %>%
    group_by(cell_type) %>%
    group_modify(~ {
      purrr::map_dfr(comparisons_list, function(comp) {
        g1 <- filter(.x, expr_group == comp[1])$bias
        g2 <- filter(.x, expr_group == comp[2])$bias

        if (length(g1) < 2 || length(g2) < 2) {
          return(NULL)
        }

        tibble(
          group1 = comp[1],
          group2 = comp[2],
          p.value = wilcox.test(g1, g2)$p.value
        )
      })
    }) %>%
    ungroup() %>%
    mutate(
      p.signif = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    ) %>%
    filter(p.signif != "ns") %>%
    mutate(
      x_num = as.numeric(cell_type),
      xmin = x_num + dodge_offsets[group1],
      xmax = x_num + dodge_offsets[group2],
      comparison = paste(group1, group2, sep = "-")
    )
}

compute_pvalue_summary <- function(all_bias_df) {
  all_bias_df %>%
    group_by(cell_type, region) %>%
    group_modify(~ {
      purrr::map_dfr(comparisons_list, function(comp) {
        g1 <- filter(.x, expr_group == comp[1])$bias
        g2 <- filter(.x, expr_group == comp[2])$bias

        if (length(g1) < 2 || length(g2) < 2) {
          return(NULL)
        }

        tibble(
          group1 = comp[1],
          group2 = comp[2],
          n_group1 = length(g1),
          n_group2 = length(g2),
          p_value = wilcox.test(g1, g2)$p.value
        )
      })
    }) %>%
    ungroup() %>%
    mutate(
      p_signif = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    ) %>%
    arrange(region, cell_type, group1, group2)
}

make_breaks <- function(ylims) {
  pretty(ylims, n = 5) %>% .[. >= ylims[1] & . <= ylims[2]]
}

plot_and_save_region <- function(region_name, ylims, output_file) {
  region_df <- bias_df %>% filter(region == region_name)
  region_points_df <- region_df %>%
    mutate(
      point_fill_col = unname(color_map[as.character(expr_group)]),
      point_outline_col = unname(point_outline_map[as.character(expr_group)])
    )
  signif_df <- compute_signif_df(region_df)

  if (nrow(signif_df) > 0) {
    y_span <- diff(ylims)
    base_height <- ylims[2] - (0.12 * y_span)
    step_height <- 0.08 * y_span

    signif_df <- signif_df %>%
      group_by(cell_type) %>%
      arrange(comparison, .by_group = TRUE) %>%
      mutate(
        y_pos = base_height + (row_number() - 1) * step_height,
        tick_height = 0.02 * y_span,
        label_height = 0.025 * y_span,
        x_mid = (xmin + xmax) / 2
      ) %>%
      ungroup()
  }

  p_base <- ggplot(region_df, aes(x = cell_type, y = bias, fill = expr_group)) +
    geom_boxplot(
      outlier.shape = NA,
      position = position_dodge(0.8),
      color = "black",
      key_glyph = draw_key_boxplot_with_point,
      show.legend = TRUE
    ) +
    geom_point(
      data = region_points_df,
      aes(x = cell_type, y = bias, color = point_fill_col, group = expr_group),
      inherit.aes = FALSE,
      shape = 16,
      size = point_size,
      stroke = 0,
      alpha = point_fill_alpha,
      position = position_dodge(0.8),
      show.legend = FALSE
    ) +
    geom_point(
      data = region_points_df,
      aes(x = cell_type, y = bias, color = point_outline_col, group = expr_group),
      inherit.aes = FALSE,
      shape = 21,
      size = point_size,
      stroke = point_stroke,
      fill = NA,
      alpha = point_outline_alpha,
      position = position_dodge(0.8),
      show.legend = FALSE
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 0.4
    ) +
    scale_fill_manual(values = color_map) +
    scale_color_identity() +
    guides(
      fill = guide_legend(),
      color = "none"
    ) +
    scale_y_continuous(
      limits = ylims,
      breaks = make_breaks(ylims)
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "none",
      legend.title = element_blank(),
      legend.key.width = unit(0.9, "cm"),
      legend.key.height = unit(0.55, "cm"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black")
    ) +
    labs(
      x = NULL,
      y = "rNMP EF difference (template - non-template)",
      fill = NULL
    )

  p <- p_base

  if (nrow(signif_df) > 0) {
    p <- p +
      geom_segment(
        data = signif_df,
        aes(x = xmin, xend = xmax, y = y_pos, yend = y_pos),
        inherit.aes = FALSE,
        linewidth = 0.4
      ) +
      geom_segment(
        data = signif_df,
        aes(x = xmin, xend = xmin, y = y_pos - tick_height, yend = y_pos),
        inherit.aes = FALSE,
        linewidth = 0.4
      ) +
      geom_segment(
        data = signif_df,
        aes(x = xmax, xend = xmax, y = y_pos - tick_height, yend = y_pos),
        inherit.aes = FALSE,
        linewidth = 0.4
      ) +
      geom_text(
        data = signif_df,
        aes(x = x_mid, y = y_pos + label_height, label = p.signif),
        inherit.aes = FALSE,
        size = 3
      )
  }

  grDevices::svg(filename = output_file, width = 8.5, height = 5.5)
  print(p)
  dev.off()

  p_nolabs <- p_base +
    labs(x = NULL, y = NULL) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    )

  nolabs_file <- sub("\\.svg$", "_nolabs.svg", output_file)
  grDevices::svg(filename = nolabs_file, width = 8.5, height = 5.5)
  print(p_nolabs)
  dev.off()

  invisible(c(output_file, nolabs_file))
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

region_jobs <- list(
  list(
    region_name = "Promoter",
    ylims = promoter_ylims,
    output_file = file.path(output_dir, "Fig3E_0-1kb-downstream-of-TSS.svg")
  ),
  list(
    region_name = "GeneBody",
    ylims = gene_body_ylims,
    output_file = file.path(output_dir, "Fig3E_TSS-to-TTS.svg")
  )
)

cl <- makeCluster(length(region_jobs))
on.exit(stopCluster(cl), add = TRUE)

clusterEvalQ(cl, {
  library(tidyverse)
  library(ggpubr)
  NULL
})

clusterExport(
  cl,
  varlist = c(
    "bias_df",
    "comparisons_list",
    "dodge_offsets",
    "color_map",
    "point_outline_map",
    "point_outline_alpha",
    "point_fill_alpha",
    "point_stroke",
    "point_size",
    "outline_from_fill",
    "draw_key_boxplot_with_point",
    "compute_signif_df",
    "make_breaks",
    "plot_and_save_region"
  ),
  envir = environment()
)

all_pvalues <- compute_pvalue_summary(bias_df)
readr::write_csv(all_pvalues, file.path(output_dir, "Fig3E_pvalues_summary.csv"))

parallel::parLapply(
  cl,
  region_jobs,
  function(job) {
    do.call(plot_and_save_region, job)
  }
)

save_vertical_legend(file.path(output_dir, "Fig3E_legend_vertical.svg"), show_text = TRUE)
save_vertical_legend(file.path(output_dir, "Fig3E_legend_vertical_notext.svg"), show_text = FALSE)
