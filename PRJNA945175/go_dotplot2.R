library(ggplot2)
library(dplyr)

load_agrigo <- function(filepath, label) {
  df <- read.table(filepath, header=TRUE, sep="\t", quote="",
                   fill=TRUE, comment.char="")
  df <- df[, 1:9]
  colnames(df) <- c("GO_acc","term_type","Term","queryitem",
                    "querytotal","bgitem","bgtotal","pvalue","FDR")
  df$subset <- label
  df$FDR <- as.numeric(df$FDR)
  df$queryitem <- as.numeric(df$queryitem)
  df$querytotal <- as.numeric(df$querytotal)
  df$GeneRatio <- df$queryitem / df$querytotal
  df <- df[df$FDR < 0.05 & df$term_type == "P", ]
  df
}

hcrv_up   <- load_agrigo("/users/fyp/fyp5/project/PRJNA945175/go_enrichment/agriGO_HCRV_up.txt", "HCRV_up")
tswv_spec <- load_agrigo("/users/fyp/fyp5/project/PRJNA945175/go_enrichment/agriGO_TSWV_specific_down.txt", "TSWV_specific_down")
tswv_down <- load_agrigo("/users/fyp/fyp5/project/PRJNA945175/go_enrichment/agriGO_TSWV_down.txt", "TSWV_down")

# Get terms unique to HCRV_up (not in TSWV subsets)
hcrv_unique_terms <- setdiff(hcrv_up$GO_acc, 
                              union(tswv_spec$GO_acc, tswv_down$GO_acc))
hcrv_unique <- hcrv_up[hcrv_up$GO_acc %in% hcrv_unique_terms, ]
hcrv_unique <- hcrv_unique[order(hcrv_unique$FDR), ][1:min(5, nrow(hcrv_unique)), ]

# Get top 5 shared terms by HCRV_up FDR
shared_terms <- intersect(hcrv_up$GO_acc, tswv_spec$GO_acc)
hcrv_shared <- hcrv_up[hcrv_up$GO_acc %in% shared_terms, ]
hcrv_shared <- hcrv_shared[order(hcrv_shared$FDR), ][1:min(5, nrow(hcrv_shared)), ]

# Selected terms to show
selected_terms <- unique(c(hcrv_unique$GO_acc, hcrv_shared$GO_acc))

# Build combined df with all three subsets for selected terms
all_data <- rbind(hcrv_up, tswv_spec, tswv_down)
plot_data <- all_data[all_data$GO_acc %in% selected_terms, ]

# Add missing combinations as NA rows so all subsets appear
all_combos <- expand.grid(GO_acc=selected_terms, 
                          subset=c("HCRV_up","TSWV_specific_down","TSWV_down"),
                          stringsAsFactors=FALSE)
plot_data <- merge(all_combos, 
                   plot_data[, c("GO_acc","Term","subset","queryitem",
                                 "querytotal","FDR","GeneRatio")],
                   by=c("GO_acc","subset"), all.x=TRUE)

# Fill Term for NAs
for (acc in selected_terms) {
  term_name <- all_data$Term[all_data$GO_acc == acc][1]
  plot_data$Term[plot_data$GO_acc == acc] <- term_name
}

plot_data$neglog10FDR <- -log10(plot_data$FDR)
plot_data$subset <- factor(plot_data$subset,
                           levels=c("HCRV_up","TSWV_specific_down","TSWV_down"))

# Label unique vs shared
plot_data$category <- ifelse(plot_data$GO_acc %in% hcrv_unique$GO_acc,
                             "HCRV_up unique", "Shared")

# Order terms: unique first, then shared
term_order <- c(hcrv_unique$Term, 
                hcrv_up$Term[hcrv_up$GO_acc %in% hcrv_shared$GO_acc])
term_order <- unique(term_order)
plot_data$Term <- factor(plot_data$Term, levels=rev(term_order))

p <- ggplot(plot_data, aes(x=subset, y=Term,
                            size=queryitem,
                            color=neglog10FDR)) +
  geom_point(na.rm=TRUE) +
  scale_color_gradient(low="#fee8c8", high="#cc0000",
                       name="-log10(FDR)", na.value="grey90") +
  scale_size_continuous(name="Gene count", range=c(2,10)) +
  facet_grid(category ~ ., scales="free_y", space="free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=30, hjust=1, size=11),
    axis.text.y = element_text(size=9),
    axis.title = element_blank(),
    strip.background = element_rect(fill="grey85"),
    strip.text = element_text(face="bold", size=10),
    panel.grid.major = element_line(colour="grey90"),
    legend.position = "right"
  ) +
  labs(title="GO Biological Process Enrichment - PRR DEG Subsets",
       subtitle="agriGO v2 SEA, NiBen v101 background, FDR < 0.05")

outdir <- "/users/fyp/fyp5/project/PRJNA945175/go_enrichment/"
ggsave(paste0(outdir, "GO_dotplot_v2.pdf"), p, width=12, height=8)
ggsave(paste0(outdir, "GO_dotplot_v2.png"), p, width=12, height=8, dpi=300)
cat("Done\n")

# Second plot — Cellular Component terms
hcrv_cc <- load_agrigo("/users/fyp/fyp5/project/PRJNA945175/go_enrichment/agriGO_HCRV_up.txt", "HCRV_up")
tswv_spec_cc <- load_agrigo("/users/fyp/fyp5/project/PRJNA945175/go_enrichment/agriGO_TSWV_specific_down.txt", "TSWV_specific_down")
tswv_down_cc <- load_agrigo("/users/fyp/fyp5/project/PRJNA945175/go_enrichment/agriGO_TSWV_down.txt", "TSWV_down")

# Override the BP filter — keep CC
load_agrigo_cc <- function(filepath, label) {
  df <- read.table(filepath, header=TRUE, sep="\t", quote="",
                   fill=TRUE, comment.char="")
  df <- df[, 1:9]
  colnames(df) <- c("GO_acc","term_type","Term","queryitem",
                    "querytotal","bgitem","bgtotal","pvalue","FDR")
  df$subset <- label
  df$FDR <- as.numeric(df$FDR)
  df$queryitem <- as.numeric(df$queryitem)
  df$querytotal <- as.numeric(df$querytotal)
  df$GeneRatio <- df$queryitem / df$querytotal
  df <- df[df$FDR < 0.05 & df$term_type == "C", ]
  df
}

hcrv_cc   <- load_agrigo_cc("/users/fyp/fyp5/project/PRJNA945175/go_enrichment/agriGO_HCRV_up.txt", "HCRV_up")
tswv_spec_cc <- load_agrigo_cc("/users/fyp/fyp5/project/PRJNA945175/go_enrichment/agriGO_TSWV_specific_down.txt", "TSWV_specific_down")
tswv_down_cc <- load_agrigo_cc("/users/fyp/fyp5/project/PRJNA945175/go_enrichment/agriGO_TSWV_down.txt", "TSWV_down")

cc_combined <- rbind(hcrv_cc, tswv_spec_cc, tswv_down_cc)
cc_combined$neglog10FDR <- -log10(cc_combined$FDR)
cc_combined$subset <- factor(cc_combined$subset,
                             levels=c("HCRV_up","TSWV_specific_down","TSWV_down"))
cc_combined$Term <- factor(cc_combined$Term,
                           levels=rev(unique(cc_combined$Term[order(cc_combined$FDR)])))

p_cc <- ggplot(cc_combined, aes(x=subset, y=Term,
                                 size=queryitem,
                                 color=neglog10FDR)) +
  geom_point(na.rm=TRUE) +
  scale_color_gradient(low="#deebf7", high="#084594",
                       name="-log10(FDR)") +
  scale_size_continuous(name="Gene count", range=c(3,10)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=30, hjust=1, size=11),
    axis.text.y = element_text(size=9),
    axis.title = element_blank(),
    panel.grid.major = element_line(colour="grey90"),
    legend.position = "right"
  ) +
  labs(title="GO Cellular Component Enrichment - PRR DEG Subsets",
       subtitle="agriGO v2 SEA, NiBen v101 background, FDR < 0.05")

ggsave(paste0(outdir, "GO_dotplot_CC.pdf"), p_cc, width=12, height=6)
ggsave(paste0(outdir, "GO_dotplot_CC.png"), p_cc, width=12, height=6, dpi=300)
cat("CC plot done\n")
