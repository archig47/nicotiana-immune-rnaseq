library(DESeq2)
library(ggplot2)

rld <- readRDS("/users/fyp/fyp5/project/PRJNA945175/deseq2/rld_group.rds")

pca_data <- plotPCA(rld, intgroup = c("condition", "timepoint"),
                    returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

virus_cols <- c(
  "CK"   = "#B4B2A9",
  "HCRV" = "#5DCAA5",
  "TSWV" = "#F0997B",
  "TH"   = "#AFA9EC"
)

tp_shapes <- c("1" = 16, "7" = 17, "14" = 15)

pca_data$virus     <- gsub("_.*", "", pca_data$condition)
pca_data$timepoint <- as.character(pca_data$timepoint)

p <- ggplot(pca_data, aes(x = PC1, y = PC2,
                           colour = virus,
                           fill   = virus,
                           shape  = timepoint)) +
  ## ellipses per condition (95% confidence, no outline)
  stat_ellipse(aes(group = virus),
               geom  = "polygon",
               alpha = 0.08,
               level = 0.95,
               show.legend = FALSE) +
  ## points on top
  geom_point(size = 3, stroke = 0.4, alpha = 0.95) +
  scale_colour_manual(values = virus_cols, name = "Condition") +
  scale_fill_manual(values = virus_cols,   name = "Condition") +
  scale_shape_manual(values = tp_shapes,
                     breaks = c("1", "7", "14"),
                     name   = "Timepoint (dpi)") +
  labs(
    x = paste0("PC1: ", pct_var[1], "% variance"),
    y = paste0("PC2: ", pct_var[2], "% variance")
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.line  = element_line(colour = "grey60", linewidth = 0.4),
    axis.ticks = element_line(colour = "grey60", linewidth = 0.4)
  )

ggsave("/users/fyp/fyp5/project/PRJNA945175/figures/pca_final.pdf",
       p, width = 14, height = 10, units = "cm")
ggsave("/users/fyp/fyp5/project/PRJNA945175/figures/pca_final.png",
       p, width = 14, height = 10, units = "cm", dpi = 300)

message("PCA final saved")
