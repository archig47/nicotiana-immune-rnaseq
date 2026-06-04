suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readxl)
})
options(bitmapType = "cairo")

base  <- "/users/fyp/fyp5/project/PRJNA945175/deseq2_clean"
nlr_d <- file.path(base, "nlr_prr")
res_d <- file.path(base, "results")

# ── Data ──────────────────────────────────────────────────────────────────
nlr_genes <- readLines(file.path(nlr_d, "nlr_gene_ids.txt"))
nlr_genes  <- nlr_genes[nlr_genes != ""]

phylo_raw <- read_excel(file.path(nlr_d, "NLR_phylogeneti_classification.xlsx"), sheet="Sheet2")
colnames(phylo_raw) <- c("protein_id","raw_class")
phylo_raw <- phylo_raw[phylo_raw$protein_id != "ID",]
phylo_class <- data.frame(
  gene_id = sub("-mRNA.*","", phylo_raw$protein_id),
  class   = dplyr::case_when(
    phylo_raw$raw_class == "CC-clade2-NRC-helper"    ~ "NRC-helper",
    phylo_raw$raw_class == "CC-clade2-NRC-sensor-01" ~ "NRC-sensor-01",
    phylo_raw$raw_class == "CC-clade2-NRC-sensor-02" ~ "NRC-sensor-02",
    phylo_raw$raw_class == "CC-clade2-ZAR1"          ~ "ZAR1",
    phylo_raw$raw_class == "CC-clade2"               ~ "CC-clade2",
    phylo_raw$raw_class == "CC-clade1"               ~ "CC-clade1",
    phylo_raw$raw_class == "TIR"                     ~ "TIR",
    phylo_raw$raw_class == "TNP"                     ~ "TNP",
    phylo_raw$raw_class == "RPW8-NRG1"               ~ "RPW8-NRG1",
    phylo_raw$raw_class == "RPW8-ADR1"               ~ "RPW8-ADR1",
    phylo_raw$raw_class == "CCG10"                   ~ "CCG10",
    TRUE ~ "Unclassified"),
  stringsAsFactors = FALSE
)

all_res <- do.call(rbind, lapply(
  list.files(res_d, pattern="results_.*\\.csv", full.names=TRUE),
  function(f) {
    d <- read.csv(f)
    p <- strsplit(tools::file_path_sans_ext(basename(f)),"_")[[1]]
    d$virus <- p[2]; d$dpi <- p[3]; d
  }
))

nlr_degs <- all_res[all_res$gene_id %in% nlr_genes & !is.na(all_res$padj) &
                    all_res$padj < 0.05 & abs(all_res$log2FoldChange) > 1,]
hcrv_genes <- unique(nlr_degs$gene_id[nlr_degs$virus=="HCRV"])
tswv_genes <- unique(nlr_degs$gene_id[nlr_degs$virus=="TSWV"])
th_genes   <- unique(nlr_degs$gene_id[nlr_degs$virus=="TH"])
all_genes  <- unique(c(hcrv_genes, tswv_genes, th_genes))

mat <- merge(
  data.frame(gene_id=all_genes,
             HCRV=as.integer(all_genes %in% hcrv_genes),
             TSWV=as.integer(all_genes %in% tswv_genes),
             TH  =as.integer(all_genes %in% th_genes)),
  phylo_class, by="gene_id", all.x=TRUE)
mat$class[is.na(mat$class)] <- "Unclassified"
mat$intersection <- dplyr::case_when(
  mat$HCRV==1&mat$TSWV==1&mat$TH==1 ~ "All three",
  mat$HCRV==1&mat$TSWV==1&mat$TH==0 ~ "HCRV & TSWV",
  mat$HCRV==0&mat$TSWV==1&mat$TH==1 ~ "TSWV & TH",
  mat$HCRV==1&mat$TSWV==0&mat$TH==1 ~ "HCRV & TH",
  mat$HCRV==1&mat$TSWV==0&mat$TH==0 ~ "HCRV only",
  mat$HCRV==0&mat$TSWV==1&mat$TH==0 ~ "TSWV only",
  mat$HCRV==0&mat$TSWV==0&mat$TH==1 ~ "TH only", TRUE ~ "None")

int_order <- names(sort(table(mat$intersection), decreasing=TRUE))
mat$intersection <- factor(mat$intersection, levels=int_order)

# ── Palettes ──────────────────────────────────────────────────────────────
class_order_all <- c("NRC-sensor-01","NRC-sensor-02","NRC-helper","CC-clade2",
                     "CC-clade1","ZAR1","CCG10","TIR","TNP","RPW8-NRG1",
                     "RPW8-ADR1","Unclassified")
# Tableau 20 palette — muted, harmonious pairs, easy on the eyes
class_cols <- c(
  "NRC-sensor-01" = "#F28E2B",  # orange         \
  "NRC-sensor-02" = "#FFBE7D",  # light orange    > NRC network
  "NRC-helper"    = "#B6992D",  # golden amber   /
  "CC-clade2"     = "#59A14F",  # green          \  other CC-NLRs
  "CC-clade1"     = "#8CD17D",  # light green    /
  "ZAR1"          = "#499894",  # teal           \  atypical CC singletons
  "CCG10"         = "#86BCB6",  # light teal     /
  "TIR"           = "#4E79A7",  # blue           \  TIR-NLRs
  "TNP"           = "#A0CBE8",  # light blue     /
  "RPW8-NRG1"     = "#B07AA1",  # purple         \  RPW8-NLRs
  "RPW8-ADR1"     = "#D4A6C8",  # light lavender /
  "Unclassified"  = "#BAB0AC"   # warm grey
)
present_classes <- class_order_all[class_order_all %in% mat$class]
mat$class <- factor(mat$class, levels=present_classes)

col_hcrv <- "#4E7DAB"; col_tswv <- "#1A9BC9"; col_th <- "#E07B2E"

# ── Intersection membership ───────────────────────────────────────────────
int_members <- mat %>%
  group_by(intersection) %>%
  summarise(HCRV=max(HCRV), TSWV=max(TSWV), TH=max(TH), .groups="drop") %>%
  mutate(intersection=as.character(intersection))

# ── y layout: bars = 0 to 25+;  dots spread wider for legibility ──────────
dot_y <- c("HCRV"=-2.6, "TSWV"=-5.0, "TH"=-7.4)

bar_data   <- mat %>% count(intersection, class) %>%
  mutate(intersection=factor(intersection,levels=int_order),
         class=factor(class,levels=present_classes))
bar_totals <- mat %>% count(intersection) %>%
  mutate(intersection=factor(intersection,levels=int_order))

dot_data <- expand.grid(virus=c("HCRV","TSWV","TH"), intersection=int_order,
                        stringsAsFactors=FALSE) %>%
  left_join(int_members, by="intersection") %>%
  mutate(active=dplyr::case_when(virus=="HCRV"~HCRV,virus=="TSWV"~TSWV,TRUE~TH),
         y_pos=dot_y[virus],
         intersection=factor(intersection,levels=int_order))

seg_data <- int_members %>%
  filter(HCRV+TSWV+TH > 1) %>%
  mutate(
    y_top    = dplyr::case_when(HCRV==1~dot_y["HCRV"],TSWV==1~dot_y["TSWV"],TRUE~dot_y["TH"]),
    y_bottom = dplyr::case_when(TH==1~dot_y["TH"],TSWV==1~dot_y["TSWV"],TRUE~dot_y["HCRV"]),
    intersection=factor(intersection,levels=int_order)
  )

# ── Virus label x position (left of all dots, outside panel) ─────────────
# x scale: levels 1..7, expand add=0.55 → leftmost dot ≈ 0.45
label_x <- -0.1   # will be clipped unless coord_cartesian(clip="off")

# ── Plot ──────────────────────────────────────────────────────────────────
p <- ggplot() +
  # Shaded dot section background
  annotate("rect", xmin=-Inf, xmax=Inf,
           ymin=min(dot_y)-1.2, ymax=-1.1,
           fill="#F6F6F6", colour=NA) +
  # Divider line
  annotate("segment", x=-Inf, xend=Inf, y=-1.1, yend=-1.1,
           colour="#BBBBBB", linewidth=0.6) +
  # Subtle row guide lines inside dot section
  annotate("segment", x=-Inf, xend=Inf, y=dot_y["HCRV"], yend=dot_y["HCRV"],
           colour="#E0E0E0", linewidth=0.3) +
  annotate("segment", x=-Inf, xend=Inf, y=dot_y["TSWV"], yend=dot_y["TSWV"],
           colour="#E0E0E0", linewidth=0.3) +
  annotate("segment", x=-Inf, xend=Inf, y=dot_y["TH"], yend=dot_y["TH"],
           colour="#E0E0E0", linewidth=0.3) +

  # ── Stacked bars ──────────────────────────────────────────────────────
  geom_col(data=bar_data,
           aes(x=intersection, y=n, fill=class),
           width=0.62, colour="white", linewidth=0.15) +
  geom_text(data=bar_totals,
            aes(x=intersection, y=n, label=n),
            inherit.aes=FALSE, vjust=-0.5, size=4.2,
            fontface="bold", colour="#2D2D2D", family="sans") +

  # ── Dot matrix ────────────────────────────────────────────────────────
  geom_segment(data=seg_data,
               aes(x=intersection, xend=intersection, y=y_top, yend=y_bottom),
               colour="#2D2D2D", linewidth=2.0, inherit.aes=FALSE) +
  geom_point(data=filter(dot_data, active==0),
             aes(x=intersection, y=y_pos),
             colour="#D8D8D8", size=5.5, shape=19, inherit.aes=FALSE) +
  geom_point(data=filter(dot_data, active==1),
             aes(x=intersection, y=y_pos),
             colour="#1A1A1A", size=7.5, shape=19, inherit.aes=FALSE) +

  # ── Virus name + count labels (left of dots) ──────────────────────────
  annotate("text", x=label_x, y=dot_y["HCRV"], hjust=1, size=4.2,
           label=paste0("HCRV  (n=",length(hcrv_genes),")"),
           fontface="bold", colour=col_hcrv, family="sans") +
  annotate("text", x=label_x, y=dot_y["TSWV"], hjust=1, size=4.2,
           label=paste0("TSWV  (n=",length(tswv_genes),")"),
           fontface="bold", colour=col_tswv, family="sans") +
  annotate("text", x=label_x, y=dot_y["TH"], hjust=1, size=4.2,
           label=paste0("TH  (n=",length(th_genes),")"),
           fontface="bold", colour=col_th, family="sans") +

  # ── Scales ────────────────────────────────────────────────────────────
  scale_fill_manual(values=class_cols[present_classes], name="NLR class",
                    breaks=present_classes,
                    guide=guide_legend(ncol=1, keywidth=0.75, keyheight=0.7,
                                       override.aes=list(colour=NA))) +
  scale_x_discrete(limits=int_order, expand=expansion(add=0.55)) +
  scale_y_continuous(
    breaks = c(0, 5, 10, 15, 20, 25),
    labels = c("0","5","10","15","20","25"),
    expand = expansion(mult=c(0.02, 0.07)),
    name   = "Number of NLR DEGs"
  ) +
  coord_cartesian(clip="off") +

  # ── Theme ─────────────────────────────────────────────────────────────
  theme_minimal(base_family="sans", base_size=12) +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.title.x       = element_blank(),
    axis.text.y        = element_text(size=11, colour="#555555"),
    axis.title.y       = element_text(size=12, colour="#444444"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour="#EEEEEE", linewidth=0.35),
    legend.title       = element_text(size=11, face="bold"),
    legend.text        = element_text(size=10.5),
    legend.position    = "right",
    plot.margin        = margin(10, 10, 10, 95),
    plot.background    = element_rect(fill="white", colour=NA),
    panel.background   = element_rect(fill="white", colour=NA)
  )

out_pdf <- file.path(base, "figures/nlr_upset_custom_v3.pdf")
out_png <- file.path(base, "figures/nlr_upset_custom_v3.png")
ggsave(out_pdf, p, width=9.5, height=7.0, device="pdf")
system(paste0(
  "gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 ",
  "-dGraphicsAlphaBits=4 -dTextAlphaBits=4 ",
  "-sOutputFile=", shQuote(out_png), " ", shQuote(out_pdf)
))
cat("Saved v3.\n")
