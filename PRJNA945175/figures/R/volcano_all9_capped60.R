suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
})

CLEAN_DIR <- "/users/fyp/fyp5/project/PRJNA945175/deseq2_clean"
RESULTS   <- file.path(CLEAN_DIR, "results")
FIG_DIR   <- file.path(CLEAN_DIR, "figures")
dir.create(FIG_DIR, showWarnings=FALSE)

COL_UP   <- "#D62728"
COL_DOWN <- "#1F77B4"
COL_NS   <- "#CCCCCC"

PADJ_THRESH <- 0.05
LFC_THRESH  <- 1
Y_MAX       <- 60

load_contrast <- function(virus, tp) {
  f  <- file.path(RESULTS, sprintf("results_%s_%sdpi_vs_CK.csv", virus, tp))
  df <- read_csv(f, show_col_types=FALSE)
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  df$neg_log10_padj <- pmin(-log10(df$padj), Y_MAX)
  df$sig <- df$padj < PADJ_THRESH & abs(df$log2FoldChange) > LFC_THRESH
  df$direction <- "NS"
  df$direction[df$sig & df$log2FoldChange >  LFC_THRESH] <- "Up"
  df$direction[df$sig & df$log2FoldChange < -LFC_THRESH] <- "Down"
  df$direction <- factor(df$direction, levels=c("Up","Down","NS"))
  df
}

make_volcano <- function(df, title) {
  n_up   <- sum(df$direction == "Up")
  n_down <- sum(df$direction == "Down")

  ggplot(df, aes(x=log2FoldChange, y=neg_log10_padj, colour=direction)) +
    geom_point(data=df[df$direction=="NS",],
               size=0.4, alpha=0.25, stroke=0) +
    geom_point(data=df[df$direction!="NS",],
               size=0.6, alpha=0.65, stroke=0) +
    geom_vline(xintercept=c(-LFC_THRESH, LFC_THRESH),
               linetype="dashed", colour="grey40", linewidth=0.3) +
    geom_hline(yintercept=-log10(PADJ_THRESH),
               linetype="dashed", colour="grey40", linewidth=0.3) +
    annotate("text", x=9.5,  y=Y_MAX*0.97, hjust=1, vjust=1,
             label=paste0(n_up, " \u25b2"),
             colour=COL_UP, size=3, fontface="bold") +
    annotate("text", x=-9.5, y=Y_MAX*0.97, hjust=0, vjust=1,
             label=paste0("\u25bc ", n_down),
             colour=COL_DOWN, size=3, fontface="bold") +
    scale_colour_manual(values=c(Up=COL_UP, Down=COL_DOWN, NS=COL_NS),
                        guide="none") +
    scale_x_continuous(limits=c(-10, 10), oob=scales::squish,
                       breaks=c(-8,-4,0,4,8)) +
    scale_y_continuous(limits=c(0, Y_MAX),
                       breaks=seq(0, Y_MAX, 10),
                       expand=expansion(mult=c(0, 0.02))) +
    labs(title=title,
         x=expression(log[2]~"fold change"),
         y=expression(-log[10](p[adj]))) +
    theme_classic(base_size=10) +
    theme(
      plot.title       = element_text(face="bold", size=10),
      axis.text        = element_text(size=7, colour="grey20"),
      axis.title       = element_text(size=8, colour="grey20"),
      axis.line        = element_line(colour="grey60", linewidth=0.35),
      axis.ticks       = element_line(colour="grey60", linewidth=0.35),
      panel.grid.major = element_line(colour="grey94", linewidth=0.25),
      plot.margin      = margin(6, 8, 6, 6)
    )
}

cat("Building all 9 panels capped at 60...\n")

contrasts <- list(
  list(v="HCRV", tp=1,  lab="A   HCRV — 1 dpi"),
  list(v="HCRV", tp=7,  lab="B   HCRV — 7 dpi"),
  list(v="HCRV", tp=14, lab="C   HCRV — 14 dpi"),
  list(v="TSWV", tp=1,  lab="D   TSWV — 1 dpi"),
  list(v="TSWV", tp=7,  lab="E   TSWV — 7 dpi"),
  list(v="TSWV", tp=14, lab="F   TSWV — 14 dpi"),
  list(v="TH",   tp=1,  lab="G   TH (co-infection) — 1 dpi"),
  list(v="TH",   tp=7,  lab="H   TH (co-infection) — 7 dpi"),
  list(v="TH",   tp=14, lab="I   TH (co-infection) — 14 dpi")
)

panels <- lapply(contrasts, function(x)
  make_volcano(load_contrast(x$v, x$tp), title=x$lab))

fig <- wrap_plots(panels, ncol=3) +
  plot_annotation(
    title    = "Differential gene expression in virus-infected N. benthamiana",
    subtitle = "padj < 0.05, |log\u2082FC| > 1 | y-axis capped at 60; HCRV 7 dpi (max 86.6) and TSWV 7 dpi (max 90.4) exceed axis range | TH = TSWV + HCRV co-infection",
    theme = theme(
      plot.title    = element_text(size=12, face="bold"),
      plot.subtitle = element_text(size=8, colour="grey45")
    )
  )

ggsave(file.path(FIG_DIR, "volcano_all9_capped60.png"),
       fig, width=30, height=30, units="cm", dpi=300)
ggsave(file.path(FIG_DIR, "volcano_all9_capped60.pdf"),
       fig, width=30, height=30, units="cm")
cat("Saved: volcano_all9_capped60\n=== Done ===\n")
