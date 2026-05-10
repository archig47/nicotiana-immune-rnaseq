library(ggplot2)
library(dplyr)

# Load agriGO results
load_agrigo <- function(filepath, label) {
  df <- read.table(filepath, header=TRUE, sep="\t", quote="", 
                   fill=TRUE, comment.char="")
  df <- df[, 1:9]
  colnames(df) <- c("GO_acc", "term_type", "Term", "queryitem", 
                    "querytotal", "bgitem", "bgtotal", "pvalue", "FDR")
  df$subset <- label
  df$FDR <- as.numeric(df$FDR)
  df$queryitem <- as.numeric(df$queryitem)
  df$bgitem <- as.numeric(df$bgitem)
  df$querytotal <- as.numeric(df$querytotal)
  df$bgtotal <- as.numeric(df$bgtotal)
  df$GeneRatio <- df$queryitem / df$querytotal
  df <- df[df$FDR < 0.05, ]
  df
}

hcrv_up   <- load_agrigo("/users/fyp/fyp5/project/PRJNA945175/go_enrichment/agriGO_HCRV_up.txt", "HCRV_up")
tswv_spec <- load_agrigo("/users/fyp/fyp5/project/PRJNA945175/go_enrichment/agriGO_TSWV_specific_down.txt", "TSWV_specific_down")
tswv_down <- load_agrigo("/users/fyp/fyp5/project/PRJNA945175/go_enrichment/agriGO_TSWV_down.txt", "TSWV_down")

# Keep only Biological Process terms
hcrv_up   <- hcrv_up[hcrv_up$term_type == "P", ]
tswv_spec <- tswv_spec[tswv_spec$term_type == "P", ]
tswv_down <- tswv_down[tswv_down$term_type == "P", ]

# Select top 10 most significant terms per subset
top_terms <- function(df, n=10) {
  df[order(df$FDR), ][1:min(n, nrow(df)), ]
}

hcrv_top   <- top_terms(hcrv_up, 10)
tswv_spec_top <- top_terms(tswv_spec, 10)
tswv_down_top <- top_terms(tswv_down, 10)

# Combine
combined <- rbind(hcrv_top, tswv_spec_top, tswv_down_top)
combined$subset <- factor(combined$subset, 
                          levels=c("HCRV_up", "TSWV_specific_down", "TSWV_down"))
combined$Term <- factor(combined$Term, 
                        levels=rev(unique(combined$Term[order(combined$FDR)])))
combined$neglog10FDR <- -log10(combined$FDR)

# Plot
p <- ggplot(combined, aes(x=subset, y=Term, 
                           size=queryitem, 
                           color=neglog10FDR)) +
  geom_point() +
  scale_color_gradient(low="#fee8c8", high="#cc0000",
                       name="-log10(FDR)") +
  scale_size_continuous(name="Gene count", range=c(3, 10)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=30, hjust=1, size=11),
    axis.text.y = element_text(size=9),
    axis.title = element_blank(),
    panel.grid.major = element_line(colour="grey90"),
    legend.position = "right"
  ) +
  labs(title="GO Biological Process Enrichment — PRR DEG Subsets",
       subtitle="agriGO v2 SEA, NiBen v101 background, FDR < 0.05")

outdir <- "/users/fyp/fyp5/project/PRJNA945175/go_enrichment/"
ggsave(paste0(outdir, "GO_dotplot_PRR_subsets.pdf"), p, 
       width=12, height=8, units="in")
ggsave(paste0(outdir, "GO_dotplot_PRR_subsets.png"), p, 
       width=12, height=8, units="in", dpi=300)

cat("Dot plot saved\n")
print(p)
