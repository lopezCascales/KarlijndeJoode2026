################################################################################
# Figure 6C — TCR gene expression in 4-1BB+ vs 4-1BB− CD4 T cells
#
# Rationale: MiXCR yielded 0 TCR clonotypes from bulk RNA-seq (expected:
# V(D)J sequences are ~0.01-0.1% of reads in non-enriched libraries).
# Instead, we quantify expression of TCR variable/joining/constant genes
# directly from HTSeq counts as a proxy for TCR repertoire activity.
#
# Samples:
#   4-1BB+ : F01counts, C01counts, E01counts  (n=3)
#   4-1BB− : G01counts, H01counts, D01counts  (n=3)
#
# Output plots:
#   1. Boxplot — total TCR gene expression per sample (with stats)
#   2. Heatmap — expression per TCR gene per sample (z-score, log2 CPM)
#   3. Barplot — mean expression per chain (TRA / TRB / TRD / TRG)
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(ggpubr)      # for stat_compare_means
library(pheatmap)    # for heatmap
library(patchwork)   # for combining plots

# ── 0. TCR gene reference table (ENSG → gene symbol, chain) ──────────────────
tcr_ref <- tribble(
  ~ENSG,                ~Symbol,         ~Chain,
  # Constant
  "ENSG00000277734",    "TRAC",          "TRA",
  "ENSG00000211751",    "TRBC1",         "TRB",
  "ENSG00000211772",    "TRBC2",         "TRB",
  "ENSG00000211666",    "TRDC",          "TRD",
  "ENSG00000211673",    "TRGC1",         "TRG",
  "ENSG00000211668",    "TRGC2",         "TRG",
  # Alpha variable
  "ENSG00000276308",    "TRAV1-1",       "TRA",
  "ENSG00000229134",    "TRAV1-2",       "TRA",
  "ENSG00000211801",    "TRAV2",         "TRA",
  "ENSG00000211800",    "TRAV3",         "TRA",
  "ENSG00000211799",    "TRAV4",         "TRA",
  "ENSG00000211798",    "TRAV5",         "TRA",
  "ENSG00000211797",    "TRAV6",         "TRA",
  "ENSG00000211796",    "TRAV7",         "TRA",
  "ENSG00000211795",    "TRAV8-1",       "TRA",
  "ENSG00000211794",    "TRAV8-2",       "TRA",
  "ENSG00000211793",    "TRAV8-3",       "TRA",
  "ENSG00000211792",    "TRAV8-4",       "TRA",
  "ENSG00000211791",    "TRAV9-1",       "TRA",
  "ENSG00000211790",    "TRAV9-2",       "TRA",
  "ENSG00000211789",    "TRAV10",        "TRA",
  "ENSG00000211788",    "TRAV12-1",      "TRA",
  "ENSG00000211787",    "TRAV12-2",      "TRA",
  "ENSG00000211786",    "TRAV12-3",      "TRA",
  "ENSG00000211785",    "TRAV13-1",      "TRA",
  "ENSG00000211784",    "TRAV13-2",      "TRA",
  "ENSG00000211783",    "TRAV14DV4",     "TRA",
  "ENSG00000211782",    "TRAV16",        "TRA",
  "ENSG00000211781",    "TRAV17",        "TRA",
  "ENSG00000211780",    "TRAV18",        "TRA",
  "ENSG00000211779",    "TRAV19",        "TRA",
  "ENSG00000211778",    "TRAV20",        "TRA",
  "ENSG00000211777",    "TRAV21",        "TRA",
  "ENSG00000211776",    "TRAV22",        "TRA",
  "ENSG00000211775",    "TRAV23DV6",     "TRA",
  "ENSG00000211774",    "TRAV24",        "TRA",
  "ENSG00000211773",    "TRAV25",        "TRA",
  "ENSG00000211771",    "TRAV26-1",      "TRA",
  "ENSG00000211770",    "TRAV26-2",      "TRA",
  "ENSG00000211769",    "TRAV27",        "TRA",
  "ENSG00000211768",    "TRAV29DV5",     "TRA",
  "ENSG00000211767",    "TRAV30",        "TRA",
  "ENSG00000211766",    "TRAV34",        "TRA",
  "ENSG00000211765",    "TRAV35",        "TRA",
  "ENSG00000211764",    "TRAV36DV7",     "TRA",
  "ENSG00000211763",    "TRAV38-1",      "TRA",
  "ENSG00000211762",    "TRAV38-2DV8",   "TRA",
  "ENSG00000211761",    "TRAV39",        "TRA",
  "ENSG00000211760",    "TRAV40",        "TRA",
  "ENSG00000211759",    "TRAV41",        "TRA",
  # Alpha joining
  "ENSG00000211658",    "TRAJ1",         "TRA",
  "ENSG00000211659",    "TRAJ2",         "TRA",
  "ENSG00000211660",    "TRAJ3",         "TRA",
  "ENSG00000211661",    "TRAJ4",         "TRA",
  "ENSG00000211662",    "TRAJ5",         "TRA",
  "ENSG00000211663",    "TRAJ6",         "TRA",
  "ENSG00000211664",    "TRAJ7",         "TRA",
  "ENSG00000211665",    "TRAJ8",         "TRA",
  # Beta variable
  "ENSG00000211813",    "TRBV2",         "TRB",
  "ENSG00000211814",    "TRBV3-1",       "TRB",
  "ENSG00000211815",    "TRBV4-1",       "TRB",
  "ENSG00000211816",    "TRBV4-2",       "TRB",
  "ENSG00000211817",    "TRBV4-3",       "TRB",
  "ENSG00000211818",    "TRBV5-1",       "TRB",
  "ENSG00000211819",    "TRBV5-4",       "TRB",
  "ENSG00000211820",    "TRBV5-5",       "TRB",
  "ENSG00000211821",    "TRBV5-6",       "TRB",
  "ENSG00000211822",    "TRBV5-8",       "TRB",
  "ENSG00000211823",    "TRBV6-1",       "TRB",
  "ENSG00000211824",    "TRBV6-2",       "TRB",
  "ENSG00000211825",    "TRBV6-5",       "TRB",
  "ENSG00000211826",    "TRBV6-6",       "TRB",
  "ENSG00000211827",    "TRBV6-8",       "TRB",
  "ENSG00000211828",    "TRBV6-9",       "TRB",
  "ENSG00000211829",    "TRBV7-2",       "TRB",
  "ENSG00000211830",    "TRBV9",         "TRB",
  "ENSG00000211831",    "TRBV10-1",      "TRB",
  "ENSG00000211832",    "TRBV10-2",      "TRB",
  "ENSG00000211833",    "TRBV10-3",      "TRB",
  "ENSG00000211834",    "TRBV11-1",      "TRB",
  "ENSG00000211835",    "TRBV11-2",      "TRB",
  "ENSG00000211836",    "TRBV12-3",      "TRB",
  "ENSG00000211837",    "TRBV12-4",      "TRB",
  "ENSG00000211838",    "TRBV13",        "TRB",
  "ENSG00000211839",    "TRBV14",        "TRB",
  "ENSG00000211840",    "TRBV15",        "TRB",
  "ENSG00000211841",    "TRBV16",        "TRB",
  "ENSG00000211842",    "TRBV18",        "TRB",
  "ENSG00000211843",    "TRBV19",        "TRB",
  "ENSG00000211844",    "TRBV20-1",      "TRB",
  "ENSG00000211845",    "TRBV25-1",      "TRB",
  "ENSG00000211846",    "TRBV27",        "TRB",
  "ENSG00000211847",    "TRBV28",        "TRB",
  "ENSG00000211848",    "TRBV29-1",      "TRB",
  "ENSG00000211849",    "TRBV30",        "TRB",
  # Beta joining
  "ENSG00000211750",    "TRBJ1-1",       "TRB",
  "ENSG00000211749",    "TRBJ1-2",       "TRB",
  "ENSG00000211748",    "TRBJ1-3",       "TRB",
  "ENSG00000211747",    "TRBJ1-4",       "TRB",
  "ENSG00000211746",    "TRBJ1-5",       "TRB",
  "ENSG00000211745",    "TRBJ1-6",       "TRB",
  "ENSG00000211744",    "TRBJ2-1",       "TRB",
  "ENSG00000211743",    "TRBJ2-2",       "TRB",
  "ENSG00000211742",    "TRBJ2-3",       "TRB",
  "ENSG00000211741",    "TRBJ2-4",       "TRB",
  "ENSG00000211740",    "TRBJ2-5",       "TRB",
  "ENSG00000211739",    "TRBJ2-6",       "TRB",
  "ENSG00000211738",    "TRBJ2-7",       "TRB"
)

# ── 1. Load and filter counts for TCR genes ───────────────────────────────────
counts_raw <- read_tsv("countsHTSeq.txt", show_col_types = FALSE)

tcr_counts <- counts_raw %>%
  filter(GENEID %in% tcr_ref$ENSG) %>%
  left_join(tcr_ref, by = c("GENEID" = "ENSG"))

message("TCR genes found in counts file: ", nrow(tcr_counts))

# ── 2. Compute library size–normalized values (CPM) ──────────────────────────
# Use total library size from all genes (not just TCR) for normalization
lib_sizes <- counts_raw %>%
  summarise(across(-GENEID, sum)) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "LibSize")

sample_meta <- tibble(
  Sample = c("F01counts", "C01counts", "E01counts",
             "G01counts", "H01counts", "D01counts"),
  Group  = c("4-1BB+", "4-1BB+", "4-1BB+",
             "4-1BB-", "4-1BB-", "4-1BB-")
)

# ── 3. Long format for plotting ───────────────────────────────────────────────
tcr_long <- tcr_counts %>%
  pivot_longer(cols = ends_with("counts"),
               names_to = "Sample", values_to = "RawCounts") %>%
  left_join(lib_sizes, by = "Sample") %>%
  left_join(sample_meta, by = "Sample") %>%
  mutate(
    CPM     = (RawCounts / LibSize) * 1e6,
    log2CPM = log2(CPM + 1),
    Group   = factor(Group, levels = c("4-1BB-", "4-1BB+"))
  )

# ── 4. PLOT 1 — Total TCR expression per sample (boxplot + jitter) ────────────
group_cols <- c("4-1BB+" = "#C0392B", "4-1BB-" = "#2980B9")

total_per_sample <- tcr_long %>%
  group_by(Sample, Group) %>%
  summarise(TotalCPM = sum(CPM), .groups = "drop")

p1 <- ggplot(total_per_sample, aes(x = Group, y = TotalCPM, fill = Group)) +
  geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.6,
               colour = "black", linewidth = 0.4) +
  geom_jitter(aes(colour = Group), width = 0.08, size = 2.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.5, label.y = max(total_per_sample$TotalCPM) * 1.08,
                     size = 3.2) +
  scale_fill_manual(values  = group_cols) +
  scale_colour_manual(values = group_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    title = "Total TCR gene expression",
    x = NULL,
    y = "Total CPM (TCR genes)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
  )

# ── 5. PLOT 2 — Expression per TCR chain (TRA / TRB) ─────────────────────────
# Focus on alpha and beta (most relevant for CD4 T cells)
chain_summary <- tcr_long %>%
  filter(Chain %in% c("TRA", "TRB")) %>%
  group_by(Sample, Group, Chain) %>%
  summarise(TotalCPM = sum(CPM), .groups = "drop")

p2 <- ggplot(chain_summary, aes(x = Group, y = TotalCPM, fill = Group)) +
  geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.6,
               colour = "black", linewidth = 0.4) +
  geom_jitter(aes(colour = Group), width = 0.08, size = 2.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.5, size = 3.0) +
  scale_fill_manual(values  = group_cols) +
  scale_colour_manual(values = group_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  facet_wrap(~ Chain, scales = "free_y") +
  labs(
    title = "TCR expression by chain",
    x = NULL,
    y = "Total CPM"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10)
  )

# ── 6. PLOT 3 — Heatmap of top expressed TCR genes ───────────────────────────
# Select genes expressed in at least 1 sample (log2CPM > 0)
top_genes <- tcr_long %>%
  group_by(Symbol) %>%
  summarise(MaxExpr = max(log2CPM)) %>%
  filter(MaxExpr > 0) %>%
  arrange(desc(MaxExpr)) %>%
  slice_head(n = 40) %>%
  pull(Symbol)

heatmap_mat <- tcr_long %>%
  filter(Symbol %in% top_genes) %>%
  select(Symbol, Sample, log2CPM) %>%
  pivot_wider(names_from = Sample, values_from = log2CPM) %>%
  column_to_rownames("Symbol") %>%
  as.matrix()

# Clean sample names for display
colnames(heatmap_mat) <- gsub("counts", "", colnames(heatmap_mat))

# Annotation for columns (groups)
col_annotation <- data.frame(
  Group = c("4-1BB+", "4-1BB+", "4-1BB+", "4-1BB-", "4-1BB-", "4-1BB-")
) %>% `rownames<-`(colnames(heatmap_mat))

# Reorder columns: 4-1BB- first, 4-1BB+ second (consistent with other plots)
col_order <- c("G01", "H01", "D01", "F01", "C01", "E01")
heatmap_mat <- heatmap_mat[, col_order]
col_annotation <- col_annotation[col_order, , drop = FALSE]

ann_colors <- list(Group = c("4-1BB+" = "#C0392B", "4-1BB-" = "#2980B9"))

# Row annotation for chain
row_annotation <- tcr_ref %>%
  filter(Symbol %in% rownames(heatmap_mat)) %>%
  select(Symbol, Chain) %>%
  distinct() %>%
  column_to_rownames("Symbol")

chain_cols <- c("TRA" = "#8E44AD", "TRB" = "#27AE60", "TRD" = "#E67E22", "TRG" = "#7F8C8D")
ann_colors$Chain <- chain_cols

pdf("Figure6C_heatmap_TCRgenes.pdf", width = 6, height = 8)
pheatmap(
  heatmap_mat,
  scale            = "row",           # z-score per gene
  cluster_cols     = FALSE,           # keep group order
  cluster_rows     = TRUE,
  annotation_col   = col_annotation,
  annotation_row   = row_annotation,
  annotation_colors = ann_colors,
  color            = colorRampPalette(c("#2980B9", "white", "#C0392B"))(100),
  border_color     = NA,
  fontsize_row     = 7,
  fontsize_col     = 9,
  main             = "TCR gene expression\n4-1BB+ vs 4-1BB− CD4 T cells",
  silent           = TRUE
)
dev.off()

# ── 7. Save boxplots ──────────────────────────────────────────────────────────
ggsave("Figure6C_boxplot_totalTCR.pdf",  p1, width = 3.0, height = 4.0)
ggsave("Figure6C_boxplot_byChain.pdf",   p2, width = 4.5, height = 4.0)

# Combined p1 + p2
combined <- p1 + p2 + plot_layout(widths = c(1, 2)) +
  plot_annotation(
    title = "TCR gene expression in 4-1BB+ vs 4-1BB\u2212 CD4 T cells",
    theme = theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))
  )
ggsave("Figure6C_combined_boxplots.pdf", combined, width = 7.5, height = 4.2)

# ── 8. Summary stats table ────────────────────────────────────────────────────
stats_table <- total_per_sample %>%
  group_by(Group) %>%
  summarise(
    Mean_TotalCPM  = round(mean(TotalCPM), 1),
    SD_TotalCPM    = round(sd(TotalCPM), 1),
    Median_TotalCPM = round(median(TotalCPM), 1),
    .groups = "drop"
  )

write_csv(stats_table, "Figure6C_stats_summary.csv")
print(stats_table)

message("\nDone. Output files:")
message("  Figure6C_boxplot_totalTCR.pdf")
message("  Figure6C_boxplot_byChain.pdf")
message("  Figure6C_combined_boxplots.pdf")
message("  Figure6C_heatmap_TCRgenes.pdf")
message("  Figure6C_stats_summary.csv")
