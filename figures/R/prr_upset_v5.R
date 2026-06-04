suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readxl)
})
options(bitmapType = "cairo")

base  <- "/users/fyp/fyp5/project/PRJNA945175/deseq2_clean"
nlr_d <- file.path(base, "nlr_prr")
res_d <- file.path(base, "results")

# ── 1. PRR annotation — same ecto classification as Figures 1 & 2 ─────────
ACCESSORY <- c("EGF","PAN","SDOM","NEW","RCC1","TNFR","GP_PDE","UNIDENTIFIED")

get_ecto <- function(rules) {
  if (!is.character(rules) || is.na(rules)) return("Other")
  tokens <- strsplit(rules, ",")[[1]][-1]
  primary <- tokens[!tokens %in% ACCESSORY][1]
  if (is.na(primary))           return("Other")
  if (primary == "LRR")         return("LRR")
  if (primary == "LysM")        return("LysM")
  if (primary == "WAK")         return("WAK")
  if (primary %in% c("MAL","SPARK"))                    return("Malectin/CrRLK")
  if (primary %in% c("BLEC","LLEC","CLEC","GLEC","GNK2")) return("Lectin")
  return("Other")
}

prr_landscape <- read_excel(file.path(nlr_d, "nlr_prr_full_landscape.xlsx"),
                             sheet = "PRR_landscape")
prr_annot <- prr_landscape %>%
  mutate(ecto = sapply(matched_rules, get_ecto)) %>%
  select(gene_id, ecto)

prr_genes <- prr_annot$gene_id

# ── 2. Load all 9 DESeq2 results ──────────────────────────────────────────
all_res <- do.call(rbind, lapply(
  list.files(res_d, pattern = "results_.*\\.csv", full.names = TRUE),
  function(f) {
    d <- read.csv(f)
    p <- strsplit(tools::file_path_sans_ext(basename(f)), "_")[[1]]
    d$virus <- p[2]; d$dpi <- p[3]; d
  }
))

prr_degs <- all_res %>%
  filter(gene_id %in% prr_genes, !is.na(padj),
         padj < 0.05, abs(log2FoldChange) > 1)

hcrv_genes <- unique(prr_degs$gene_id[prr_degs$virus == "HCRV"])
tswv_genes <- unique(prr_degs$gene_id[prr_degs$virus == "TSWV"])
th_genes   <- unique(prr_degs$gene_id[prr_degs$virus == "TH"])
all_deg_genes <- unique(c(hcrv_genes, tswv_genes, th_genes))

cat(sprintf("Total PRR DEGs: %d  |  HCRV: %d  TSWV: %d  TH: %d\n",
            length(all_deg_genes),
            length(hcrv_genes), length(tswv_genes), length(th_genes)))

# ── 3. Intersection matrix ────────────────────────────────────────────────
mat <- data.frame(gene_id = all_deg_genes,
                  HCRV = as.integer(all_deg_genes %in% hcrv_genes),
                  TSWV = as.integer(all_deg_genes %in% tswv_genes),
                  TH   = as.integer(all_deg_genes %in% th_genes),
                  stringsAsFactors = FALSE) %>%
  left_join(prr_annot, by = "gene_id")
mat$ecto[is.na(mat$ecto)] <- "Other"

mat <- mat %>%
  mutate(intersection = case_when(
    HCRV==1 & TSWV==1 & TH==1 ~ "All three",
    HCRV==1 & TSWV==1 & TH==0 ~ "HCRV & TSWV",
    HCRV==0 & TSWV==1 & TH==1 ~ "TSWV & TH",
    HCRV==1 & TSWV==0 & TH==1 ~ "HCRV & TH",
    HCRV==1 & TSWV==0 & TH==0 ~ "HCRV only",
    HCRV==0 & TSWV==1 & TH==0 ~ "TSWV only",
    HCRV==0 & TSWV==0 & TH==1 ~ "TH only",
    TRUE ~ "None"
  ))

int_order <- names(sort(table(mat$intersection), decreasing = TRUE))
mat$intersection <- factor(mat$intersection, levels = int_order)

# ── 4. Colour palette — warm/cool tones distinct from NLR palette ─────────
# NLR uses Tableau oranges/greens/blues/purples; PRR uses reds/amber/teal/indigo
ecto_order <- c("LRR","Lectin","WAK","Malectin/CrRLK","LysM","Other")
ecto_cols  <- c(
  "LRR"           = "#5255A8",   # indigo
  "Lectin"        = "#6E8FD4",   # periwinkle
  "WAK"           = "#22AEAD",   # aqua
  "Malectin/CrRLK"= "#B8365A",   # raspberry pink
  "LysM"          = "#E07848",   # coral orange
  "Other"         = "#9A9A9A"    # cool grey
)

present_ectos <- ecto_order[ecto_order %in% mat$ecto]
mat$ecto <- factor(mat$ecto, levels = present_ectos)

col_hcrv <- "#4E7DAB"; col_tswv <- "#1A9BC9"; col_th <- "#E07B2E"

# ── 5. Intersection membership (for dot matrix) ───────────────────────────
int_members <- mat %>%
  group_by(intersection) %>%
  summarise(HCRV = max(HCRV), TSWV = max(TSWV), TH = max(TH), .groups = "drop") %>%
  mutate(intersection = as.character(intersection))

# ── 6. Layout: dots sit below y = 0 ──────────────────────────────────────
dot_y <- c("HCRV" = -22, "TSWV" = -42, "TH" = -62)
divider_y <- -10

bar_data   <- mat %>% count(intersection, ecto) %>%
  mutate(intersection = factor(intersection, levels = int_order),
         ecto = factor(ecto, levels = present_ectos))
bar_totals <- mat %>% count(intersection) %>%
  mutate(intersection = factor(intersection, levels = int_order))

dot_data <- expand.grid(virus = c("HCRV","TSWV","TH"),
                        intersection = int_order,
                        stringsAsFactors = FALSE) %>%
  left_join(int_members, by = "intersection") %>%
  mutate(
    active = case_when(virus=="HCRV"~HCRV, virus=="TSWV"~TSWV, TRUE~TH),
    y_pos  = dot_y[virus],
    intersection = factor(intersection, levels = int_order)
  )

seg_data <- int_members %>%
  filter(HCRV + TSWV + TH > 1) %>%
  mutate(
    y_top    = case_when(HCRV==1 ~ dot_y["HCRV"], TSWV==1 ~ dot_y["TSWV"],
                         TRUE ~ dot_y["TH"]),
    y_bottom = case_when(TH==1 ~ dot_y["TH"], TSWV==1 ~ dot_y["TSWV"],
                         TRUE ~ dot_y["HCRV"]),
    intersection = factor(intersection, levels = int_order)
  )

label_x <- -0.1

# ── 7. Plot ───────────────────────────────────────────────────────────────
p <- ggplot() +
  # Dot-section shaded background
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = min(dot_y) - 12, ymax = divider_y,
           fill = "#F6F6F6", colour = NA) +
  annotate("segment", x = -Inf, xend = Inf,
           y = divider_y, yend = divider_y,
           colour = "#BBBBBB", linewidth = 0.6) +
  # Subtle row guides in dot section
  annotate("segment", x = -Inf, xend = Inf,
           y = dot_y["HCRV"], yend = dot_y["HCRV"],
           colour = "#E0E0E0", linewidth = 0.3) +
  annotate("segment", x = -Inf, xend = Inf,
           y = dot_y["TSWV"], yend = dot_y["TSWV"],
           colour = "#E0E0E0", linewidth = 0.3) +
  annotate("segment", x = -Inf, xend = Inf,
           y = dot_y["TH"], yend = dot_y["TH"],
           colour = "#E0E0E0", linewidth = 0.3) +

  # Stacked bars
  geom_col(data = bar_data,
           aes(x = intersection, y = n, fill = ecto),
           width = 0.62, colour = "white", linewidth = 0.15) +
  geom_text(data = bar_totals,
            aes(x = intersection, y = n, label = n),
            inherit.aes = FALSE, vjust = -0.5, size = 3.2,
            fontface = "bold", colour = "#2D2D2D", family = "sans") +

  # Connecting segments
  geom_segment(data = seg_data,
               aes(x = intersection, xend = intersection,
                   y = y_top, yend = y_bottom),
               colour = "#2D2D2D", linewidth = 2.0, inherit.aes = FALSE) +
  # Inactive dots
  geom_point(data = filter(dot_data, active == 0),
             aes(x = intersection, y = y_pos),
             colour = "#D8D8D8", size = 5.5, shape = 19, inherit.aes = FALSE) +
  # Active dots
  geom_point(data = filter(dot_data, active == 1),
             aes(x = intersection, y = y_pos),
             colour = "#1A1A1A", size = 7.5, shape = 19, inherit.aes = FALSE) +

  # Virus labels left of dot matrix
  annotate("text", x = label_x, y = dot_y["HCRV"], hjust = 1, size = 3.3,
           label = paste0("HCRV  (n=", length(hcrv_genes), ")"),
           fontface = "bold", colour = col_hcrv, family = "sans") +
  annotate("text", x = label_x, y = dot_y["TSWV"], hjust = 1, size = 3.3,
           label = paste0("TSWV  (n=", length(tswv_genes), ")"),
           fontface = "bold", colour = col_tswv, family = "sans") +
  annotate("text", x = label_x, y = dot_y["TH"], hjust = 1, size = 3.3,
           label = paste0("TH  (n=", length(th_genes), ")"),
           fontface = "bold", colour = col_th, family = "sans") +

  # Scales
  scale_fill_manual(values = ecto_cols[present_ectos],
                    name   = "PRR ectodomain",
                    breaks = present_ectos,
                    guide  = guide_legend(ncol = 1, keywidth = 0.75, keyheight = 0.7,
                                          override.aes = list(colour = NA))) +
  scale_x_discrete(limits = int_order, expand = expansion(add = 0.55)) +
  scale_y_continuous(
    breaks = seq(0, 200, 50),
    expand = expansion(mult = c(0.02, 0.07)),
    name   = "Number of PRR DEGs"
  ) +
  coord_cartesian(clip = "off") +

  labs(
    title    = "PRR DEG overlap across HCRV, TSWV and TH infection",
    subtitle = expression(italic("N. benthamiana") ~
      "  |  666 PRR DEGs  |  padj < 0.05, |log2FC| > 1  |  union across 1, 7, 14 dpi")
  ) +
  theme_minimal(base_family = "sans", base_size = 10) +
  theme(
    plot.title         = element_text(size = 13, face = "bold", colour = "#1A202C"),
    plot.subtitle      = element_text(size = 9,  colour = "#718096"),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.title.x       = element_blank(),
    axis.text.y        = element_text(size = 9,  colour = "#555555"),
    axis.title.y       = element_text(size = 9,  colour = "#444444"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.35),
    legend.title       = element_text(size = 9,  face = "bold"),
    legend.text        = element_text(size = 8.5),
    legend.position    = "right",
    plot.margin        = margin(10, 10, 10, 95),
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.background   = element_rect(fill = "white", colour = NA)
  )

out_pdf <- file.path(base, "figures/prr_upset_v5.pdf")
out_png <- file.path(base, "figures/prr_upset_v5.png")
ggsave(out_pdf, p, width = 9.5, height = 7.0, device = "pdf")
system(paste0(
  "gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 ",
  "-dGraphicsAlphaBits=4 -dTextAlphaBits=4 ",
  "-sOutputFile=", shQuote(out_png), " ", shQuote(out_pdf)
))
cat("Saved prr_upset_v5.\n")
