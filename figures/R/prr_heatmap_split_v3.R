options(bitmapType = "cairo")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readxl)
  library(tidyr)
  library(scales)
  library(patchwork)
})

base  <- "/users/fyp/fyp5/project/PRJNA945175/deseq2_clean"
nlr_d <- file.path(base, "nlr_prr")
res_d <- file.path(base, "results")

# -- 1. Load PRR landscape ---------------------------------------------------
prr_landscape <- read_excel(file.path(nlr_d, "nlr_prr_full_landscape.xlsx"),
                             sheet = "PRR_landscape")

# -- 2. Classification functions (TWO-TIER: ecto x backbone) ----------------
get_ecto <- function(rules) {
  if (!is.character(rules) || is.na(rules)) return("Other")
  tokens <- strsplit(rules, ",")[[1]][-1]
  ACCESSORY <- c("EGF","PAN","SDOM","NEW","RCC1","TNFR","GP_PDE","UNIDENTIFIED")
  primary <- tokens[!tokens %in% ACCESSORY][1]
  if (is.na(primary)) return("Other")
  if (primary == "LRR")                                    return("LRR")
  if (primary == "LysM")                                   return("LysM")
  if (primary == "WAK")                                    return("WAK")
  if (primary %in% c("MAL","SPARK"))                       return("Malectin/CrRLK")
  if (primary %in% c("BLEC","LLEC","CLEC","GLEC","GNK2")) return("Lectin")
  return("Other")
}

get_backbone <- function(rules) {
  if (!is.character(rules) || is.na(rules)) return("Other")
  strsplit(rules, ",")[[1]][1]
}

# Apply classification
prr_annot <- prr_landscape %>%
  mutate(
    ecto_group = sapply(matched_rules, get_ecto),
    backbone   = sapply(matched_rules, get_backbone),
    prr_class  = paste(ecto_group, backbone, sep = " - ")
  ) %>%
  select(gene_id, ecto_group, backbone, prr_class)

cat("PRR class distribution:\n")
print(sort(table(prr_annot$prr_class), decreasing = TRUE))
cat("Total PRR genes:", nrow(prr_annot), "\n\n")

# -- 3. Load all 9 DESeq2 results CSVs --------------------------------------
result_files <- list.files(res_d, pattern = "results_.*\\.csv", full.names = TRUE)

all_res <- bind_rows(lapply(result_files, function(f) {
  d     <- read.csv(f)
  parts <- strsplit(tools::file_path_sans_ext(basename(f)), "_")[[1]]
  virus <- parts[2]
  dpi   <- gsub("dpi", "", parts[3])
  d$contrast <- paste0(virus, "_", dpi, "dpi")
  d
}))

cat("Contrasts found:", paste(sort(unique(all_res$contrast)), collapse = ", "), "\n\n")

# -- 4. Filter PRR DEGs and join annotation ----------------------------------
prr_degs <- all_res %>%
  inner_join(prr_annot, by = "gene_id") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

n_unique <- length(unique(prr_degs$gene_id))
cat("Unique PRR DEG genes:", n_unique, "/ 666 expected\n\n")
if (n_unique != 666) warning("PRR DEG count mismatch! Expected 666, got ", n_unique)

# -- 5. Define row and column order ------------------------------------------
ecto_order     <- c("LRR", "Lectin", "WAK", "Malectin/CrRLK", "LysM", "Other")
backbone_order <- c("Secreted", "RLP", "RLK")  # RLK last = top of y-axis

col_order <- c(
  "HCRV_1dpi", "HCRV_7dpi", "HCRV_14dpi",
  "TSWV_1dpi", "TSWV_7dpi", "TSWV_14dpi",
  "TH_1dpi",   "TH_7dpi",   "TH_14dpi"
)

x_labels <- c(
  "HCRV_1dpi"  = "1dpi",  "HCRV_7dpi"  = "7dpi",  "HCRV_14dpi" = "14dpi",
  "TSWV_1dpi"  = "1dpi",  "TSWV_7dpi"  = "7dpi",  "TSWV_14dpi" = "14dpi",
  "TH_1dpi"    = "1dpi",  "TH_7dpi"    = "7dpi",  "TH_14dpi"   = "14dpi"
)

# -- 6. Split DEGs into up and down ------------------------------------------
up_degs   <- prr_degs %>% filter(log2FoldChange > 0)
down_degs <- prr_degs %>% filter(log2FoldChange < 0)

cat("Up DEG records:", nrow(up_degs), "\n")
cat("Down DEG records:", nrow(down_degs), "\n\n")

# -- 7. Compute summaries per ecto x backbone x contrast --------------------
summarise_panel <- function(degs, value_col = "mean_lfc") {
  degs %>%
    group_by(ecto_group, backbone, contrast) %>%
    summarise(
      n        = n(),
      mean_lfc = mean(log2FoldChange),
      .groups  = "drop"
    )
}

up_raw   <- summarise_panel(up_degs)
down_raw <- summarise_panel(down_degs)

# Build full grid and join
full_grid <- expand.grid(
  ecto_group = ecto_order,
  backbone   = backbone_order,
  contrast   = col_order,
  stringsAsFactors = FALSE
)

build_plot_df <- function(raw, full_grid) {
  full_grid %>%
    left_join(raw, by = c("ecto_group", "backbone", "contrast")) %>%
    mutate(
      ecto_group = factor(ecto_group, levels = ecto_order),
      backbone   = factor(backbone,   levels = backbone_order),
      contrast   = factor(contrast,   levels = col_order),
      n_label    = ifelse(!is.na(n) & n > 0, as.character(n), NA_character_)
    )
}

up_df   <- build_plot_df(up_raw,   full_grid)
down_df <- build_plot_df(down_raw, full_grid)

# -- 8. Shared theme helper --------------------------------------------------
shared_theme <- function(show_y_labels = TRUE) {
  t <- theme_minimal(base_family = "sans", base_size = 9) +
    theme(
      axis.text.x.top    = element_text(size = 8, colour = "#333333",
                                        angle = 0, hjust = 0.5),
      panel.grid         = element_blank(),
      legend.position    = "right",
      legend.title       = element_text(size = 8, colour = "#333333"),
      legend.text        = element_text(size = 7, colour = "#444444"),
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.background   = element_rect(fill = "white", colour = NA),
      strip.placement    = "outside",
      strip.background.y = element_rect(fill = "#F0F0F0", colour = NA),
      panel.spacing.y    = unit(0, "lines"),
      plot.margin        = margin(t = 4, r = 4, b = 8, l = 6)
    )
  if (show_y_labels) {
    t <- t + theme(
      axis.text.y        = element_text(size = 8, colour = "#333333"),
      strip.text.y.left  = element_text(angle = 0, hjust = 1, size = 8,
                                        face = "bold", colour = "#1A202C")
    )
  } else {
    t <- t + theme(
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      strip.text.y.left  = element_blank(),
      strip.background.y = element_blank()
    )
  }
  t
}

# -- 9. Build TOP panel (Upregulated) ----------------------------------------
p_up <- ggplot(up_df, aes(x = contrast, y = backbone)) +

  geom_tile(aes(fill = mean_lfc), colour = "white", linewidth = 0.4) +

  geom_text(
    data   = dplyr::filter(up_df, !is.na(n) & n > 0),
    aes(label = n),
    size   = 2.5,
    colour = "white",
    family = "sans"
  ) +

  scale_fill_gradient(
    name     = expression("Mean log"[2]*"FC (up)"),
    low      = "white",
    high     = "#6B2D8B",
    limits   = c(0, 4),
    oob      = scales::squish,
    na.value = "#CCCCCC",
    breaks   = c(0, 1, 2, 3, 4),
    labels   = c("0", "1", "2", "3", "4"),
    guide    = guide_colourbar(
      barwidth       = 0.8,
      barheight      = 5,
      ticks.colour   = "#555555",
      frame.colour   = "#AAAAAA",
      title.position = "top",
      label.theme    = element_text(size = 7, colour = "#333333")
    )
  ) +

  geom_vline(xintercept = c(3.5, 6.5), colour = "#888888",
             linewidth = 0.4, linetype = "solid") +

  scale_x_discrete(labels = x_labels, position = "top") +
  scale_y_discrete(drop = FALSE) +

  facet_grid(
    ecto_group ~ .,
    scales = "free_y",
    space  = "free_y",
    switch = "y"
  ) +

  coord_cartesian(clip = "off") +

  labs(x = NULL, y = NULL, title = "Upregulated PRR DEGs") +

  shared_theme(show_y_labels = TRUE) +
  theme(
    plot.title = element_text(size = 9, face = "bold", colour = "#1A202C",
                              hjust = 0.5, margin = margin(b = 4))
  )

# -- 10. Build BOTTOM panel (Downregulated) -----------------------------------
down_df <- down_df %>%
  mutate(abs_mean_lfc = abs(mean_lfc))

p_down <- ggplot(down_df, aes(x = contrast, y = backbone)) +

  geom_tile(aes(fill = abs_mean_lfc), colour = "white", linewidth = 0.4) +

  geom_text(
    data   = dplyr::filter(down_df, !is.na(n) & n > 0),
    aes(label = n),
    size   = 2.5,
    colour = "white",
    family = "sans"
  ) +

  scale_fill_gradient(
    name     = expression("Mean |log"[2]*"FC| (down)"),
    low      = "white",
    high     = "#0D7377",
    limits   = c(0, 4),
    oob      = scales::squish,
    na.value = "#CCCCCC",
    breaks   = c(0, 1, 2, 3, 4),
    labels   = c("0", "1", "2", "3", "4"),
    guide    = guide_colourbar(
      barwidth       = 0.8,
      barheight      = 5,
      ticks.colour   = "#555555",
      frame.colour   = "#AAAAAA",
      title.position = "top",
      label.theme    = element_text(size = 7, colour = "#333333")
    )
  ) +

  geom_vline(xintercept = c(3.5, 6.5), colour = "#888888",
             linewidth = 0.4, linetype = "solid") +

  # dpi labels at top — no virus group headers (those are on the virus strip above p_up only)
  scale_x_discrete(labels = x_labels, position = "top") +
  scale_y_discrete(drop = FALSE) +

  facet_grid(
    ecto_group ~ .,
    scales = "free_y",
    space  = "free_y",
    switch = "y"
  ) +

  coord_cartesian(clip = "off") +

  labs(x = NULL, y = NULL, title = "Downregulated PRR DEGs") +

  shared_theme(show_y_labels = TRUE) +
  theme(
    plot.title = element_text(size = 9, face = "bold", colour = "#1A202C",
                              hjust = 0.5, margin = margin(b = 4))
  )

# -- 11. Virus group label strip (above top panel only) ----------------------
make_virus_strip <- function() {
  virus_strip_df <- data.frame(
    label = c("HCRV", "TSWV", "TH (co-infection)"),
    x     = c(2, 5, 8),
    y     = 1
  )
  ggplot(virus_strip_df, aes(x = x, y = y, label = label)) +
    geom_text(size = 3.2, fontface = "bold", colour = "#1A202C", family = "sans") +
    scale_x_continuous(limits = c(0.5, 9.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.8, 1.2), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin     = margin(t = 2, r = 4, b = 0, l = 4)
    )
}

p_virus <- make_virus_strip()

# -- 12. Combine panels vertically --------------------------------------------
# Virus strip sits above the top (up) panel only.
# Stack: virus strip / up panel / down panel
# Give the virus strip a small fixed height relative to the heatmap panels.
# Up and down panels share equal height.
p_combined <- p_virus / p_up / p_down +
  plot_layout(ncol = 1, nrow = 3, heights = c(0.04, 0.48, 0.48))

# -- 13. Title / subtitle / caption ------------------------------------------
p_title <- ggplot() +
  labs(
    title    = "PRR ectodomain class response: upregulated vs downregulated",
    subtitle = expression(
      italic("N. benthamiana") ~
        "  |  666 PRR DEGs  |  padj < 0.05, |log"[2]*"FC| > 1"
    )
  ) +
  theme_void() +
  theme(
    plot.title      = element_text(size = 10, face = "bold", colour = "#1A202C",
                                   hjust = 0, margin = margin(b = 3)),
    plot.subtitle   = element_text(size = 8, colour = "#718096",
                                   hjust = 0, margin = margin(b = 2)),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin     = margin(t = 8, r = 10, b = 2, l = 10)
  )

p_caption <- ggplot() +
  labs(
    caption = paste(
      "Cell colour: mean log₂FC of upregulated (violet) or |log₂FC| of downregulated (teal) DEGs, capped at 4.",
      "Cell numbers = count of DEGs in that class × contrast.",
      "Grey cells = no DEGs.",
      sep = "  "
    )
  ) +
  theme_void() +
  theme(
    plot.caption    = element_text(size = 7, colour = "#999999",
                                   hjust = 0, margin = margin(t = 2)),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin     = margin(t = 0, r = 10, b = 6, l = 10)
  )

p_final <- p_title / p_combined / p_caption +
  plot_layout(heights = c(0.05, 0.91, 0.04)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

# -- 14. Save -----------------------------------------------------------------
out_pdf <- file.path(base, "figures/prr_heatmap_split_v3.pdf")
out_png <- file.path(base, "figures/prr_heatmap_split_v3.png")

ggsave(out_pdf, p_final, width = 8, height = 14, device = "pdf")
cat("Saved PDF:", out_pdf, "\n")

system(paste0(
  "gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 ",
  "-dGraphicsAlphaBits=4 -dTextAlphaBits=4 ",
  "-sOutputFile=", shQuote(out_png), " ", shQuote(out_pdf)
))
cat("Saved PNG:", out_png, "\n")
