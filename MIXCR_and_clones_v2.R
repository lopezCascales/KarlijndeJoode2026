cd ~/LABMEMBERS/RENO

# Compara tamaños de archivo para el mismo sample en ambas carpetas
echo "=== Tamaños comparados para C01 ==="
ls -lh fastq/I25-1078-C01* mixcr_fastq/I25-1078-C01*

# Compara número de reads en R1 para ambas carpetas (rápido, solo primer archivo)
echo ""
echo "=== Reads en fastq/C01 ==="
zcat fastq/I25-1078-C01_S7_R1_001.fastq.gz | wc -l | awk '{print $1/4, "reads"}'

echo ""
echo "=== Reads en mixcr_fastq/C01 ==="
zcat mixcr_fastq/I25-1078-C01_S82_R1_001.fastq.gz | wc -l | awk '{print $1/4, "reads"}'
Dime los resultados y te digo exactamente qué hacer. Las opciones son:
SituaciónQué hacerMismos tamaños / mismo nº readsSon duplicados — usa solo una carpetaTamaños distintosSon runs complementarios — hay que concatenar R1+R1 y R2+R2 antes de correr TRUST4Un folder tiene más readsUsa el que tenga más (probablemente fastq/)

cd ~/LABMEMBERS/RENO

# Ver cuántos reads tiene cada carpeta para C01
echo "=== fastq/ (bulk?) ==="
zcat fastq/I25-1078-C01_S7_R1_001.fastq.gz | wc -l | awk '{print $1/4, "reads"}'

echo "=== mixcr_fastq/ (TCR enriched?) ==="
zcat mixcr_fastq/I25-1078-C01_S82_R1_001.fastq.gz | wc -l | awk '{print $1/4, "reads"}'

# Ver las primeras secuencias de ambos — si son TCR enriched deberían
# empezar con secuencias V(D)J, no con genes housekeeping
echo ""
echo "=== Primeras 8 secuencias fastq/ ==="
zcat fastq/I25-1078-C01_S7_R1_001.fastq.gz | head -32 | grep -v "^@" | grep -v "^+" | grep "^[ACGT]"

echo ""
echo "=== Primeras 8 secuencias mixcr_fastq/ ==="
zcat mixcr_fastq/I25-1078-C01_S82_R1_001.fastq.gz | head -32 | grep -v "^@" | grep -v "^+" | grep "^[ACGT]"

echo "=== mixcr_fastq/ (TCR enriched?) ==="
zcat mixcr_fastq/I25-1078-C01_S82_R1_001.fastq.gz | wc -l | awk '{print $1/4, "reads"}'
# Ver las primeras secuencias de ambos — si son TCR enriched deberían
# empezar con secuencias V(D)J, no con genes housekeeping
echo ""
echo "=== Primeras 8 secuencias fastq/ ==="
zcat fastq/I25-1078-C01_S7_R1_001.fastq.gz | head -32 | grep -v "^@" | grep -v "^+" | grep "^[ACGT]"
echo ""
echo "=== Primeras 8 secuencias mixcr_fastq/ ==="
zcat mixcr_fastq/I25-1078-C01_S82_R1_001.fastq.gz | head -32 | grep -v "^@" | grep -v "^+" | grep "^[ACGT]"
=== fastq/ (bulk?) ===
20161756 reads
=== mixcr_fastq/ (TCR enriched?) ===
2411089 reads
=== Primeras 8 secuencias fastq/ ===
GGCTTTGAGAGCATAAGGCTTCTGGGTCCCCTTCTGTTTGCATCTGTACAC
GGTCTTGGTGGTAGTAGCAAATATTCAAACGAGAACTTTGAAGGCCGAAGT
TTTTAGTAGAGAGGGGGTTTCACCATGTTGCCCAGGCTGGTCGTGAACTCC
TCGTATTACACCTGAGGAGGCCAAGTACAAGTTGTGCAAAGTGAGAAAGAT
GCCTTGGAAAGCGTCGCGGTTCCGGCGGCGTCCGGTGAGCTCTCGCTGGCC
GAGTATGAGTTCCTCACTCCAGTTGAGGAGGCTCCCAAGGGCATGCTGGCG
GGCTGGGAGTTTGGGACCAGCCTGACCAACATGGAGAAACCCCATCTCTAC
GTCTTTCGCCCCTATACCCAGGTCGGACGACCGATTTGCACGTCAGGACCG
=== Primeras 8 secuencias mixcr_fastq/ ===
TCGGCGATCAGAGGGCGATGAAGTTCTAGATCCATTGAGACAAGCTCTAGA
CAACAGACCGGAGCTTCTCAGCCGTCAGGGCCCCAAACCGTACGGGCACCT
TCTGGTAAAGGGCCCGAACAGCTTCACCACGGCGGTCATGTACTTCCTGTC
TTTATACCAAGAAGGTTGGGAAAGCACCAAAATCTGCATGTGGTGTGTGCC
ATAGAGATTAATACAACTACTTAAAAAATATAGTCAATAGGTTACTAAGAT
CCCATATCCGCAGCAGGTCTCCAAGGTGAACAGCCTCTGGCATGTTGGAAC
GTCCTGTGTGCCGGCAGCGCTTTCTGCTCAAGAATCTCCGGCCCAATCGAC
TCTTAAAGGTTTCCCTGTGAAATTACAAATGTACATGGATCAACATCTGTC
lopez_cascales@BM-581122AEE16D:~/LABMEMBERS/RENO$

###############################################################################################################

Muy claro ahora:

fastq/ → 20M reads → bulk RNA-seq (el que usaste para HTSeq/DESeq2)
mixcr_fastq/ → 2.4M reads → también bulk, solo que submuestreado o una lane menor

Las secuencias de ambas carpetas se ven igual — genes housekeeping mezclados, nada que indique enriquecimiento TCR. No hay librería TCR enriquecida.


Usar Figure6C_TCRexpression_v2.R con los counts y ser explícito en el título y Methods. Lo que tienes es científicamente válido y publicable — simplemente no son clonotypes sino expresión de genes TCR.
Título figura: "TCR gene expression in 4-1BB⁺ vs 4-1BB⁻ CD4 T cells"


Since bulk RNA-seq libraries were not enriched for V(D)J sequences, TCR repertoire activity was assessed by quantifying expression of TCR variable, joining and constant genes (TRAV, TRAJ, TRAC, TRBV, TRBJ, TRBC1/2) from HTSeq read counts, normalized as counts per million (CPM). Groups were compared using the Wilcoxon rank-sum test.

################################################################################
# Figure 6C — TCR gene expression in 4-1BB+ vs 4-1BB− CD4 T cells
# Fixed version: TCR reference as read_csv (avoids tribble() slowness)
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(ggpubr)
library(pheatmap)
library(patchwork)

# ── 0. TCR gene reference (fast: CSV string instead of tribble) ───────────────
tcr_ref <- read_csv(
"ENSG,Symbol,Chain
ENSG00000277734,TRAC,TRA
ENSG00000211751,TRBC1,TRB
ENSG00000211772,TRBC2,TRB
ENSG00000211666,TRDC,TRD
ENSG00000211673,TRGC1,TRG
ENSG00000211668,TRGC2,TRG
ENSG00000276308,TRAV1-1,TRA
ENSG00000229134,TRAV1-2,TRA
ENSG00000211801,TRAV2,TRA
ENSG00000211800,TRAV3,TRA
ENSG00000211799,TRAV4,TRA
ENSG00000211798,TRAV5,TRA
ENSG00000211797,TRAV6,TRA
ENSG00000211796,TRAV7,TRA
ENSG00000211795,TRAV8-1,TRA
ENSG00000211794,TRAV8-2,TRA
ENSG00000211793,TRAV8-3,TRA
ENSG00000211792,TRAV8-4,TRA
ENSG00000211791,TRAV9-1,TRA
ENSG00000211790,TRAV9-2,TRA
ENSG00000211789,TRAV10,TRA
ENSG00000211788,TRAV12-1,TRA
ENSG00000211787,TRAV12-2,TRA
ENSG00000211786,TRAV12-3,TRA
ENSG00000211785,TRAV13-1,TRA
ENSG00000211784,TRAV13-2,TRA
ENSG00000211783,TRAV14DV4,TRA
ENSG00000211782,TRAV16,TRA
ENSG00000211781,TRAV17,TRA
ENSG00000211780,TRAV18,TRA
ENSG00000211779,TRAV19,TRA
ENSG00000211778,TRAV20,TRA
ENSG00000211777,TRAV21,TRA
ENSG00000211776,TRAV22,TRA
ENSG00000211775,TRAV23DV6,TRA
ENSG00000211774,TRAV24,TRA
ENSG00000211773,TRAV25,TRA
ENSG00000211771,TRAV26-1,TRA
ENSG00000211770,TRAV26-2,TRA
ENSG00000211769,TRAV27,TRA
ENSG00000211768,TRAV29DV5,TRA
ENSG00000211767,TRAV30,TRA
ENSG00000211766,TRAV34,TRA
ENSG00000211765,TRAV35,TRA
ENSG00000211764,TRAV36DV7,TRA
ENSG00000211763,TRAV38-1,TRA
ENSG00000211762,TRAV38-2DV8,TRA
ENSG00000211761,TRAV39,TRA
ENSG00000211760,TRAV40,TRA
ENSG00000211759,TRAV41,TRA
ENSG00000211658,TRAJ1,TRA
ENSG00000211659,TRAJ2,TRA
ENSG00000211660,TRAJ3,TRA
ENSG00000211661,TRAJ4,TRA
ENSG00000211662,TRAJ5,TRA
ENSG00000211663,TRAJ6,TRA
ENSG00000211664,TRAJ7,TRA
ENSG00000211665,TRAJ8,TRA
ENSG00000211813,TRBV2,TRB
ENSG00000211814,TRBV3-1,TRB
ENSG00000211815,TRBV4-1,TRB
ENSG00000211816,TRBV4-2,TRB
ENSG00000211817,TRBV4-3,TRB
ENSG00000211818,TRBV5-1,TRB
ENSG00000211819,TRBV5-4,TRB
ENSG00000211820,TRBV5-5,TRB
ENSG00000211821,TRBV5-6,TRB
ENSG00000211822,TRBV5-8,TRB
ENSG00000211823,TRBV6-1,TRB
ENSG00000211824,TRBV6-2,TRB
ENSG00000211825,TRBV6-5,TRB
ENSG00000211826,TRBV6-6,TRB
ENSG00000211827,TRBV6-8,TRB
ENSG00000211828,TRBV6-9,TRB
ENSG00000211829,TRBV7-2,TRB
ENSG00000211830,TRBV9,TRB
ENSG00000211831,TRBV10-1,TRB
ENSG00000211832,TRBV10-2,TRB
ENSG00000211833,TRBV10-3,TRB
ENSG00000211834,TRBV11-1,TRB
ENSG00000211835,TRBV11-2,TRB
ENSG00000211836,TRBV12-3,TRB
ENSG00000211837,TRBV12-4,TRB
ENSG00000211838,TRBV13,TRB
ENSG00000211839,TRBV14,TRB
ENSG00000211840,TRBV15,TRB
ENSG00000211841,TRBV16,TRB
ENSG00000211842,TRBV18,TRB
ENSG00000211843,TRBV19,TRB
ENSG00000211844,TRBV20-1,TRB
ENSG00000211845,TRBV25-1,TRB
ENSG00000211846,TRBV27,TRB
ENSG00000211847,TRBV28,TRB
ENSG00000211848,TRBV29-1,TRB
ENSG00000211849,TRBV30,TRB
ENSG00000211750,TRBJ1-1,TRB
ENSG00000211749,TRBJ1-2,TRB
ENSG00000211748,TRBJ1-3,TRB
ENSG00000211747,TRBJ1-4,TRB
ENSG00000211746,TRBJ1-5,TRB
ENSG00000211745,TRBJ1-6,TRB
ENSG00000211744,TRBJ2-1,TRB
ENSG00000211743,TRBJ2-2,TRB
ENSG00000211742,TRBJ2-3,TRB
ENSG00000211741,TRBJ2-4,TRB
ENSG00000211740,TRBJ2-5,TRB
ENSG00000211739,TRBJ2-6,TRB
ENSG00000211738,TRBJ2-7,TRB",
show_col_types = FALSE)


head(tcr_ref)
# A tibble: 6 × 3
ENSG            Symbol Chain
<chr>           <chr>  <chr>
  1 ENSG00000277734 TRAC   TRA  
2 ENSG00000211751 TRBC1  TRB  
3 ENSG00000211772 TRBC2  TRB  
4 ENSG00000211666 TRDC   TRD  
5 ENSG00000211673 TRGC1  TRG  
6 ENSG00000211668 TRGC2  TRG  


message("TCR reference loaded: ", nrow(tcr_ref), " genes")

# TCR reference loaded: 108 genes


# ── 1. Load and filter counts ─────────────────────────────────────────────────
counts_raw <- read_tsv("countsHTSeq.txt", show_col_types = FALSE)

tcr_counts <- counts_raw %>%
  filter(GENEID %in% tcr_ref$ENSG) %>%
  left_join(tcr_ref, by = c("GENEID" = "ENSG"))

message("TCR genes found in counts file: ", nrow(tcr_counts))
# TCR genes found in counts file: 88 

# ── 2. Library sizes (from all genes, for CPM normalisation) ─────────────────
lib_sizes <- counts_raw %>%
  summarise(across(-GENEID, sum)) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "LibSize")

sample_meta <- tibble(
  Sample = c("F01counts", "C01counts", "E01counts",
             "G01counts", "H01counts", "D01counts"),
  Group  = c("4-1BB+", "4-1BB+", "4-1BB+",
             "4-1BB-", "4-1BB-", "4-1BB-")
) %>% mutate(Group = factor(Group, levels = c("4-1BB-", "4-1BB+")))

group_cols <- c("4-1BB+" = "#C0392B", "4-1BB-" = "#2980B9")

# ── 3. Long format + CPM ─────────────────────────────────────────────────────
tcr_long <- tcr_counts %>%
  pivot_longer(cols = ends_with("counts"),
               names_to = "Sample", values_to = "RawCounts") %>%
  left_join(lib_sizes,   by = "Sample") %>%
  left_join(sample_meta, by = "Sample") %>%
  mutate(
    CPM     = (RawCounts / LibSize) * 1e6,
    log2CPM = log2(CPM + 1)
  )

# ── 4. PLOT 1 — Total TCR expression per sample ───────────────────────────────
total_per_sample <- tcr_long %>%
  group_by(Sample, Group) %>%
  summarise(TotalCPM = sum(CPM), .groups = "drop")

p1 <- ggplot(total_per_sample, aes(x = Group, y = TotalCPM, fill = Group)) +
  geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.6,
               colour = "black", linewidth = 0.4) +
  geom_jitter(aes(colour = Group), width = 0.08, size = 2.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.5,
                     label.y = max(total_per_sample$TotalCPM) * 1.08,
                     size = 3.2) +
  scale_fill_manual(values  = group_cols) +
  scale_colour_manual(values = group_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)),
                     labels = scales::comma) +
  labs(title = "Total TCR gene expression",
       x = NULL, y = "Total CPM (TCR genes)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5))

# ── 5. PLOT 2 — By chain (TRA / TRB) ─────────────────────────────────────────
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
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)),
                     labels = scales::comma) +
  facet_wrap(~ Chain, scales = "free_y") +
  labs(title = "TCR expression by chain",
       x = NULL, y = "Total CPM") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10))

# ── 6. PLOT 3 — Heatmap top 40 TCR genes ─────────────────────────────────────
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

colnames(heatmap_mat) <- gsub("counts", "", colnames(heatmap_mat))

# Reorder: 4-1BB- first, 4-1BB+ second
col_order <- c("G01", "H01", "D01", "F01", "C01", "E01")
heatmap_mat <- heatmap_mat[, col_order]

col_ann <- data.frame(
  Group = c("4-1BB-", "4-1BB-", "4-1BB-", "4-1BB+", "4-1BB+", "4-1BB+"),
  row.names = col_order
)

row_ann <- tcr_ref %>%
  filter(Symbol %in% rownames(heatmap_mat)) %>%
  select(Symbol, Chain) %>%
  distinct() %>%
  column_to_rownames("Symbol")

ann_colors <- list(
  Group = c("4-1BB+" = "#C0392B", "4-1BB-" = "#2980B9"),
  Chain = c("TRA" = "#8E44AD", "TRB" = "#27AE60", "TRD" = "#E67E22", "TRG" = "#7F8C8D")
)

pdf("Figure6C_heatmap_TCRgenes.pdf", width = 6, height = 8)
pheatmap(
  heatmap_mat,
  scale             = "row",
  cluster_cols      = FALSE,
  cluster_rows      = TRUE,
  annotation_col    = col_ann,
  annotation_row    = row_ann,
  annotation_colors = ann_colors,
  color             = colorRampPalette(c("#2980B9", "white", "#C0392B"))(100),
  border_color      = NA,
  fontsize_row      = 7,
  fontsize_col      = 9,
  main              = "TCR gene expression\n4-1BB+ vs 4-1BB− CD4 T cells",
  silent            = TRUE
)
dev.off()

# ── 7. Save combined boxplots ─────────────────────────────────────────────────
ggsave("Figure6C_boxplot_totalTCR.pdf",  p1, width = 3.0, height = 4.0)
ggsave("Figure6C_boxplot_byChain.pdf",   p2, width = 4.5, height = 4.0)

combined <- p1 + p2 + plot_layout(widths = c(1, 2)) +
  plot_annotation(
    title = "TCR gene expression in 4-1BB\u207a vs 4-1BB\u207b CD4 T cells",
    theme = theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))
  )
ggsave("Figure6C_combined_boxplots.pdf", combined, width = 7.5, height = 4.2)

# ── 8. Summary stats ──────────────────────────────────────────────────────────
stats_table <- total_per_sample %>%
  group_by(Group) %>%
  summarise(Mean = round(mean(TotalCPM), 1),
            SD   = round(sd(TotalCPM), 1),
            .groups = "drop")

write_csv(stats_table, "Figure6C_stats_summary.csv")
print(stats_table)

message("\nDone. Output: Figure6C_boxplot_totalTCR.pdf, Figure6C_boxplot_byChain.pdf,")
message("              Figure6C_heatmap_TCRgenes.pdf, Figure6C_combined_boxplots.pdf")

########################################################################################################
./run_mixcr_rnaseq.sh
MiXCR version:
MiXCR v4.7.0 (built Wed Aug 07 21:19:48 CEST 2024; rev=976ba14139; branch=no_branch; host=fv-az1019-185)
RepSeq.IO v2.5.0 (rev=06fa1852ee)
MiLib v3.5.0 (rev=b6cfcdc2af)
Built-in V/D/J/C library: repseqio.v5.1

Library search path:
- built-in libraries
- /home/lopez_cascales/LABMEMBERS/RENO/.

========================================
Processing: C01 (1000)
========================================

>>>>>>>>>>>>>>>>>>>>>>> mixcr align <<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr align --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.align.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.align.report.json --preset rna-seq --save-output-file-names /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.align.list.tsv --species hsa /home/lopez_cascales/LABMEMBERS/RENO/fastq/I25-1078-C01_S7_R1_001.fastq.gz /home/lopez_cascales/LABMEMBERS/RENO/fastq/I25-1078-C01_S7_R2_001.fastq.gz /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.alignments.vdjca
Alignment: 0%
Alignment: 10.5%  ETA: 00:01:50
Alignment: 21.2%  ETA: 00:01:28
Alignment: 31.9%  ETA: 00:01:16
Alignment: 42.6%  ETA: 00:01:04
Alignment: 53.1%  ETA: 00:00:53
Alignment: 63.3%  ETA: 00:00:43
Alignment: 73.5%  ETA: 00:00:31
Alignment: 83.7%  ETA: 00:00:19
Alignment: 93.7%  ETA: 00:00:07
====================== report: align ======================
Analysis time: 1.95m
Total sequencing reads: 20161756
Successfully aligned reads: 4666 (0.02%)
Coverage (percent of successfully aligned):
  CDR3: 595 (12.75%)
  FR3_TO_FR4: 0 (0%)
  CDR2_TO_FR4: 0 (0%)
  FR2_TO_FR4: 0 (0%)
  CDR1_TO_FR4: 0 (0%)
  VDJRegion: 0 (0%)
Alignment failed: no hits (not TCR/IG?): 19595641 (97.19%)
Alignment failed: no CDR3 parts: 9476 (0.05%)
Alignment failed: low total score: 551973 (2.74%)
Overlapped: 9777063 (48.49%)
Overlapped and aligned: 1541 (0.01%)
Overlapped and not aligned: 9775522 (48.49%)
Alignment-aided overlaps, percent of overlapped and aligned: 225 (14.6%)
No CDR3 parts alignments, percent of successfully aligned: 113 (2.42%)
Partial aligned reads, percent of successfully aligned: 3958 (84.83%)
J gene chimeras: 11 (0%)
Paired-end alignment conflicts eliminated: 1 (0%)
Realigned with forced non-floating bound: 20769836 (103.02%)
Realigned with forced non-floating right bound in left read: 738 (0%)
Realigned with forced non-floating left bound in right read: 738 (0%)
TRA chains: 823 (17.64%)
TRA non-functional: 15 (1.82%)
TRB chains: 3583 (76.79%)
TRB non-functional: 7 (0.2%)
TRD chains: 2 (0.04%)
TRD non-functional: 0 (0%)
TRG chains: 37 (0.79%)
TRG non-functional: 6 (16.22%)
TRAD chains: 211 (4.52%)
TRAD non-functional: 0 (0%)
IGH chains: 4 (0.09%)
IGH non-functional: 0 (0%)
IGL chains: 6 (0.13%)
IGL non-functional: 0 (0%)
Trimming report:
  R1 reads trimmed left: 57216 (0.28%)
  R1 reads trimmed right: 2 (0%)
  Average R1 nucleotides trimmed left: 0.008501392438238019
  Average R1 nucleotides trimmed right: 9.919770877100189E-8
  R2 reads trimmed left: 1 (0%)
  R2 reads trimmed right: 13 (0%)
  Average R2 nucleotides trimmed left: 9.919770877100189E-8
  Average R2 nucleotides trimmed right: 2.9759312631300566E-6

>>>>>>>>>>>>>>>>>> mixcr assemblePartial <<<<<<<<<<<<<<<<<<
Running:
mixcr assemblePartial --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.assemblePartial.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.assemblePartial.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.alignments.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.passembled.1.vdjca
Building index: 0%
Searching for overlaps: 0%
================= report: assemblePartial =================
Analysis time: 266ms
Total alignments analysed: 4666
Number of output alignments: 4499 (96.42%)
Alignments already with CDR3 (no overlapping is performed): 898 (19.25%)
Successfully overlapped alignments: 167 (3.58%)
Left parts with too small N-region (failed to extract k-mer): 710 (15.22%)
Extracted k-mer diversity: 3837
Dropped due to wildcard in k-mer: 0 (0%)
Dropped due to too short NRegion parts in overlap: 126 (2.7%)
Dropped overlaps with empty N region due to no complete NDN coverage: 1 (0.02%)
Number of left-side alignments: 752 (16.12%)
Number of right-side alignments: 2490 (53.36%)
Complex overlaps: 0 (0%)
Over-overlaps: 0 (0%)
Partial alignments written to output: 3434 (73.6%)

>>>>>>>>>>>>>>>>>> mixcr assemblePartial <<<<<<<<<<<<<<<<<<
Running:
mixcr assemblePartial --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.assemblePartial.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.assemblePartial.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.passembled.1.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.passembled.2.vdjca
Building index: 0%
Searching for overlaps: 0%
================= report: assemblePartial =================
Analysis time: 71ms
Total alignments analysed: 4499
Number of output alignments: 4498 (99.98%)
Alignments already with CDR3 (no overlapping is performed): 1065 (23.67%)
Successfully overlapped alignments: 1 (0.02%)
Left parts with too small N-region (failed to extract k-mer): 702 (15.6%)
Extracted k-mer diversity: 2917
Dropped due to wildcard in k-mer: 0 (0%)
Dropped due to too short NRegion parts in overlap: 94 (2.09%)
Dropped overlaps with empty N region due to no complete NDN coverage: 1 (0.02%)
Number of left-side alignments: 583 (12.96%)
Number of right-side alignments: 2313 (51.41%)
Complex overlaps: 0 (0%)
Over-overlaps: 0 (0%)
Partial alignments written to output: 3432 (76.28%)

>>>>>>>>>>>>>>>>>>>>>>> mixcr extend <<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr extend --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.extend.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.extend.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.passembled.2.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.extended.vdjca
Extending alignments: 0%
====================== report: extend ======================
Analysis time: 55ms
Extended alignments count: 301 (6.69%)
V extensions total: 34 (0.76%)
V extensions with merged targets: 0 (0%)
J extensions total: 267 (5.94%)
J extensions with merged targets: 12 (0.27%)
V+J extensions: 0 (0%)
Mean V extension length: 2.764705882352941
Mean J extension length: 5.348314606741573

>>>>>>>>>>>>>>>>>>>>>> mixcr assemble <<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr assemble --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.assemble.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.assemble.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.extended.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.clna
Initialization: progress unknown
Initialization: progress unknown
===================== report: assemble =====================
Analysis time: 389ms
Final clonotype count: 348
Reads used in clonotypes, percent of total: 1223 (0.01%)
Average number of reads per clonotype: 3.51
Reads dropped due to the lack of a clone sequence, percent of total: 3436 (0.02%)
Reads dropped due to a too short clonal sequence, percent of total: 0 (0%)
Reads dropped due to low quality, percent of total: 0 (0%)
Reads dropped due to failed mapping, percent of total: 7 (0%)
Reads dropped with low quality clones, percent of total: 0 (0%)
Aligned reads processed: 1230
Reads used in clonotypes before clustering, percent of total: 1223 (0.01%)
Number of reads used as a core, percent of used: 1208 (98.77%)
Mapped low quality reads, percent of used: 15 (1.23%)
Reads clustered in PCR error correction, percent of used: 0 (0%)
Reads pre-clustered due to the similar VJC-lists, percent of used: 0 (0%)
Clonotypes dropped as low quality: 0
Clonotypes eliminated by PCR error correction: 0
Clonotypes pre-clustered due to the similar VJC-lists: 0
Clones dropped in post filtering: 0 (0%)
Reads dropped in post filtering: 0.0 (0%)
Reads filtered by tag prefix: 0 (0%)
TRA chains: 98 (28.16%)
TRA non-functional: 16 (16.33%)
TRB chains: 247 (70.98%)
TRB non-functional: 8 (3.24%)
TRG chains: 2 (0.57%)
TRG non-functional: 2 (100%)
IGL chains: 1 (0.29%)
IGL non-functional: 0 (0%)

>>>>>>>>>>>>>>>>>> mixcr assembleContigs <<<<<<<<<<<<<<<<<<
Running:
mixcr assembleContigs --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.assembleContigs.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.assembleContigs.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.clna /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.contigs.clns
Assembling contigs: 0%
================= report: assembleContigs =================
Analysis time: 280ms
Initial clonotype count: 348
Final clonotype count: 348 (100%)
Canceled assemblies: 0 (0%)
Clones filtered out after assemble: 0 (0%)
Number of premature termination assembly events, percent of number of initial clonotypes: 0.0 (0%)
Longest contig length: 254
Clustered variants: 0 (0%)
Reads in clustered variants: 0.0 (0%)
Reads in divided (newly created) clones: 0.0 (0%)
Clones with ambiguous letters in splitting region: 0 (0%)
Reads in clones with ambiguous letters in splitting region: 0.0 (0%)
Average number of ambiguous letters per clone with ambiguous letters in splitting region: NaN
Clones combined after clustering: 0 (0%)

>>>>>>>>>>>>>>>>>>>>>>>>> mixcr qc <<<<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr qc --print-to-stdout /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.contigs.clns /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.qc.txt /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.qc.json

  Alignments dropped due to low sequence quality: 0.0% [OK]

>>>>>>>>>>>>>>>>>>>> mixcr exportClones <<<<<<<<<<<<<<<<<<<<
Running:
mixcr exportClones /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.contigs.clns /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/C01_rnaseq.clones.tsv
Exporting TRB
Exporting clones: 0.4%
Filtered 247 of 348 clones (29.02%).
Filtered 857.0 of 1223.0 reads (29.93%).
Exporting TRA
Exporting clones: 1%
Filtered 98 of 348 clones (71.84%).
Filtered 355.0 of 1223.0 reads (70.97%).
Exporting TRG
Exporting clones: 100%
Filtered 2 of 348 clones (99.43%).
Filtered 9.0 of 1223.0 reads (99.26%).
Exporting IGL
Filtered 1 of 348 clones (99.71%).
Filtered 2.0 of 1223.0 reads (99.84%).
Analysis finished successfully.
  Done: C01

========================================
Processing: E01 (24)
========================================

>>>>>>>>>>>>>>>>>>>>>>> mixcr align <<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr align --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.align.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.align.report.json --preset rna-seq --save-output-file-names /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.align.list.tsv --species hsa /home/lopez_cascales/LABMEMBERS/RENO/fastq/I25-1078-E01_S9_R1_001.fastq.gz /home/lopez_cascales/LABMEMBERS/RENO/fastq/I25-1078-E01_S9_R2_001.fastq.gz /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.alignments.vdjca
Alignment: 0%
Alignment: 10.1%  ETA: 00:02:04
Alignment: 20.6%  ETA: 00:01:38
Alignment: 30.7%  ETA: 00:01:29
Alignment: 41%  ETA: 00:01:13
Alignment: 51.1%  ETA: 00:01:03
Alignment: 61.6%  ETA: 00:00:51
Alignment: 72.3%  ETA: 00:00:36
Alignment: 82.6%  ETA: 00:00:23
Alignment: 93%  ETA: 00:00:09
====================== report: align ======================
Analysis time: 2.21m
Total sequencing reads: 23966622
Successfully aligned reads: 2260 (0.01%)
Coverage (percent of successfully aligned):
  CDR3: 239 (10.58%)
  FR3_TO_FR4: 0 (0%)
  CDR2_TO_FR4: 0 (0%)
  FR2_TO_FR4: 0 (0%)
  CDR1_TO_FR4: 0 (0%)
  VDJRegion: 0 (0%)
Alignment failed: no hits (not TCR/IG?): 23404901 (97.66%)
Alignment failed: no CDR3 parts: 5308 (0.02%)
Alignment failed: low total score: 554153 (2.31%)
Overlapped: 10958832 (45.73%)
Overlapped and aligned: 542 (0%)
Overlapped and not aligned: 10958290 (45.72%)
Alignment-aided overlaps, percent of overlapped and aligned: 78 (14.39%)
No CDR3 parts alignments, percent of successfully aligned: 65 (2.88%)
Partial aligned reads, percent of successfully aligned: 1956 (86.55%)
Realigned with forced non-floating bound: 26015736 (108.55%)
Realigned with forced non-floating right bound in left read: 615 (0%)
Realigned with forced non-floating left bound in right read: 615 (0%)
TRA chains: 390 (17.26%)
TRA non-functional: 1 (0.26%)
TRB chains: 1545 (68.36%)
TRB non-functional: 0 (0%)
TRAD chains: 325 (14.38%)
TRAD non-functional: 0 (0%)
Trimming report:
  R1 reads trimmed left: 64993 (0.27%)
  R1 reads trimmed right: 4 (0%)
  Average R1 nucleotides trimmed left: 0.008123923346393997
  Average R1 nucleotides trimmed right: 2.5034817171981934E-7
  R2 reads trimmed left: 1 (0%)
  R2 reads trimmed right: 21 (0%)
  Average R2 nucleotides trimmed left: 8.344939057327311E-8
  Average R2 nucleotides trimmed right: 5.4659350825493886E-6

>>>>>>>>>>>>>>>>>> mixcr assemblePartial <<<<<<<<<<<<<<<<<<
Running:
mixcr assemblePartial --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.assemblePartial.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.assemblePartial.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.alignments.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.passembled.1.vdjca
Building index: 0%
Searching for overlaps: 0%
================= report: assemblePartial =================
Analysis time: 273ms
Total alignments analysed: 2260
Number of output alignments: 2115 (93.58%)
Alignments already with CDR3 (no overlapping is performed): 366 (16.19%)
Successfully overlapped alignments: 145 (6.42%)
Left parts with too small N-region (failed to extract k-mer): 465 (20.58%)
Extracted k-mer diversity: 300
Dropped due to wildcard in k-mer: 0 (0%)
Dropped due to too short NRegion parts in overlap: 53 (2.35%)
Dropped overlaps with empty N region due to no complete NDN coverage: 1 (0.04%)
Number of left-side alignments: 483 (21.37%)
Number of right-side alignments: 1007 (44.56%)
Complex overlaps: 0 (0%)
Over-overlaps: 3 (0.13%)
Partial alignments written to output: 1604 (70.97%)

>>>>>>>>>>>>>>>>>> mixcr assemblePartial <<<<<<<<<<<<<<<<<<
Running:
mixcr assemblePartial --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.assemblePartial.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.assemblePartial.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.passembled.1.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.passembled.2.vdjca
Building index: 0%
Searching for overlaps: 0%
================= report: assemblePartial =================
Analysis time: 40ms
Total alignments analysed: 2115
Number of output alignments: 2110 (99.76%)
Alignments already with CDR3 (no overlapping is performed): 511 (24.16%)
Successfully overlapped alignments: 5 (0.24%)
Left parts with too small N-region (failed to extract k-mer): 457 (21.61%)
Extracted k-mer diversity: 262
Dropped due to wildcard in k-mer: 0 (0%)
Dropped due to too short NRegion parts in overlap: 50 (2.36%)
Dropped overlaps with empty N region due to no complete NDN coverage: 0 (0%)
Number of left-side alignments: 338 (15.98%)
Number of right-side alignments: 857 (40.52%)
Complex overlaps: 0 (0%)
Over-overlaps: 0 (0%)
Partial alignments written to output: 1594 (75.37%)

>>>>>>>>>>>>>>>>>>>>>>> mixcr extend <<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr extend --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.extend.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.extend.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.passembled.2.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.extended.vdjca
Extending alignments: 0%
====================== report: extend ======================
Analysis time: 24ms
Extended alignments count: 131 (6.21%)
V extensions total: 3 (0.14%)
V extensions with merged targets: 0 (0%)
J extensions total: 128 (6.07%)
J extensions with merged targets: 0 (0%)
V+J extensions: 0 (0%)
Mean V extension length: 1.0
Mean J extension length: 3.2421875

>>>>>>>>>>>>>>>>>>>>>> mixcr assemble <<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr assemble --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.assemble.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.assemble.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.extended.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.clna
Initialization: progress unknown
Initialization: progress unknown
===================== report: assemble =====================
Analysis time: 161ms
Final clonotype count: 10
Reads used in clonotypes, percent of total: 665 (0%)
Average number of reads per clonotype: 66.5
Reads dropped due to the lack of a clone sequence, percent of total: 1594 (0.01%)
Reads dropped due to a too short clonal sequence, percent of total: 0 (0%)
Reads dropped due to low quality, percent of total: 0 (0%)
Reads dropped due to failed mapping, percent of total: 0 (0%)
Reads dropped with low quality clones, percent of total: 0 (0%)
Aligned reads processed: 666
Reads used in clonotypes before clustering, percent of total: 666 (0%)
Number of reads used as a core, percent of used: 659 (98.95%)
Mapped low quality reads, percent of used: 7 (1.05%)
Reads clustered in PCR error correction, percent of used: 1 (0.15%)
Reads pre-clustered due to the similar VJC-lists, percent of used: 0 (0%)
Clonotypes dropped as low quality: 0
Clonotypes eliminated by PCR error correction: 1
Clonotypes pre-clustered due to the similar VJC-lists: 0
Clones dropped in post filtering: 0 (0%)
Reads dropped in post filtering: 0.0 (0%)
Reads filtered by tag prefix: 0 (0%)
TRA chains: 6 (60%)
TRA non-functional: 1 (16.67%)
TRB chains: 4 (40%)
TRB non-functional: 0 (0%)

>>>>>>>>>>>>>>>>>> mixcr assembleContigs <<<<<<<<<<<<<<<<<<
Running:
mixcr assembleContigs --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.assembleContigs.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.assembleContigs.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.clna /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.contigs.clns
Assembling contigs: 0%
================= report: assembleContigs =================
Analysis time: 108ms
Initial clonotype count: 10
Final clonotype count: 10 (100%)
Canceled assemblies: 0 (0%)
Clones filtered out after assemble: 0 (0%)
Number of premature termination assembly events, percent of number of initial clonotypes: 0.0 (0%)
Longest contig length: 350
Clustered variants: 0 (0%)
Reads in clustered variants: 0.0 (0%)
Reads in divided (newly created) clones: 0.0 (0%)
Clones with ambiguous letters in splitting region: 0 (0%)
Reads in clones with ambiguous letters in splitting region: 0.0 (0%)
Average number of ambiguous letters per clone with ambiguous letters in splitting region: NaN
Clones combined after clustering: 0 (0%)

>>>>>>>>>>>>>>>>>>>>>>>>> mixcr qc <<<<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr qc --print-to-stdout /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.contigs.clns /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.qc.txt /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.qc.json

  Alignments dropped due to low sequence quality: 0.0% [OK]

>>>>>>>>>>>>>>>>>>>> mixcr exportClones <<<<<<<<<<<<<<<<<<<<
Running:
mixcr exportClones /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.contigs.clns /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/E01_rnaseq.clones.tsv
Exporting TRB
Exporting clones: 25%
Filtered 4 of 10 clones (60%).
Filtered 495.0 of 665.0 reads (25.56%).
Exporting TRA
Exporting clones: 16.7%
Filtered 6 of 10 clones (40%).
Filtered 170.0 of 665.0 reads (74.44%).
Analysis finished successfully.
  Done: E01

========================================
Processing: F01 (27)
========================================

>>>>>>>>>>>>>>>>>>>>>>> mixcr align <<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr align --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.align.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.align.report.json --preset rna-seq --save-output-file-names /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.align.list.tsv --species hsa /home/lopez_cascales/LABMEMBERS/RENO/fastq/I25-1078-F01_S10_R1_001.fastq.gz /home/lopez_cascales/LABMEMBERS/RENO/fastq/I25-1078-F01_S10_R2_001.fastq.gz /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.alignments.vdjca
Alignment: 0%
Alignment: 10%  ETA: 00:01:56
Alignment: 20.6%  ETA: 00:01:30
Alignment: 30.9%  ETA: 00:01:20
Alignment: 41%  ETA: 00:01:10
Alignment: 51.2%  ETA: 00:00:57
Alignment: 61.8%  ETA: 00:00:46
Alignment: 72.5%  ETA: 00:00:33
Alignment: 83%  ETA: 00:00:20
Alignment: 93.7%  ETA: 00:00:07
====================== report: align ======================
Analysis time: 2.03m
Total sequencing reads: 21201912
Successfully aligned reads: 4836 (0.02%)
Coverage (percent of successfully aligned):
  CDR3: 673 (13.92%)
  FR3_TO_FR4: 0 (0%)
  CDR2_TO_FR4: 0 (0%)
  FR2_TO_FR4: 0 (0%)
  CDR1_TO_FR4: 0 (0%)
  VDJRegion: 0 (0%)
Alignment failed: no hits (not TCR/IG?): 20632119 (97.31%)
Alignment failed: no CDR3 parts: 10152 (0.05%)
Alignment failed: low total score: 554805 (2.62%)
Overlapped: 11439433 (53.95%)
Overlapped and aligned: 1867 (0.01%)
Overlapped and not aligned: 11437566 (53.95%)
Alignment-aided overlaps, percent of overlapped and aligned: 230 (12.32%)
No CDR3 parts alignments, percent of successfully aligned: 134 (2.77%)
Partial aligned reads, percent of successfully aligned: 4029 (83.31%)
V gene chimeras: 1 (0%)
J gene chimeras: 6 (0%)
Paired-end alignment conflicts eliminated: 5 (0%)
Realigned with forced non-floating bound: 19525418 (92.09%)
Realigned with forced non-floating right bound in left read: 679 (0%)
Realigned with forced non-floating left bound in right read: 679 (0%)
TRA chains: 1118 (23.12%)
TRA non-functional: 47 (4.2%)
TRB chains: 3223 (66.65%)
TRB non-functional: 4 (0.12%)
TRD chains: 9 (0.19%)
TRD non-functional: 0 (0%)
TRG chains: 69 (1.43%)
TRG non-functional: 1 (1.45%)
TRAD chains: 388 (8.02%)
TRAD non-functional: 0 (0%)
IGH chains: 25 (0.52%)
IGH non-functional: 0 (0%)
IGL chains: 4 (0.08%)
IGL non-functional: 0 (0%)
Trimming report:
  R1 reads trimmed left: 59247 (0.28%)
  R1 reads trimmed right: 8 (0%)
  Average R1 nucleotides trimmed left: 0.00837160346670621
  Average R1 nucleotides trimmed right: 3.7732446017132795E-7
  R2 reads trimmed left: 1 (0%)
  R2 reads trimmed right: 16 (0%)
  Average R2 nucleotides trimmed left: 9.433111504283199E-8
  Average R2 nucleotides trimmed right: 4.7165557521416E-6

>>>>>>>>>>>>>>>>>> mixcr assemblePartial <<<<<<<<<<<<<<<<<<
Running:
mixcr assemblePartial --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.assemblePartial.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.assemblePartial.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.alignments.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.passembled.1.vdjca
Building index: 0%
Searching for overlaps: 0%
================= report: assemblePartial =================
Analysis time: 235ms
Total alignments analysed: 4836
Number of output alignments: 4651 (96.17%)
Alignments already with CDR3 (no overlapping is performed): 1050 (21.71%)
Successfully overlapped alignments: 185 (3.83%)
Left parts with too small N-region (failed to extract k-mer): 731 (15.12%)
Extracted k-mer diversity: 2356
Dropped due to wildcard in k-mer: 0 (0%)
Dropped due to too short NRegion parts in overlap: 151 (3.12%)
Dropped overlaps with empty N region due to no complete NDN coverage: 1 (0.02%)
Number of left-side alignments: 708 (14.64%)
Number of right-side alignments: 2535 (52.42%)
Complex overlaps: 0 (0%)
Over-overlaps: 0 (0%)
Partial alignments written to output: 3416 (70.64%)

>>>>>>>>>>>>>>>>>> mixcr assemblePartial <<<<<<<<<<<<<<<<<<
Running:
mixcr assemblePartial --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.assemblePartial.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.assemblePartial.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.passembled.1.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.passembled.2.vdjca
Building index: 0%
Searching for overlaps: 0%
================= report: assemblePartial =================
Analysis time: 67ms
Total alignments analysed: 4651
Number of output alignments: 4650 (99.98%)
Alignments already with CDR3 (no overlapping is performed): 1235 (26.55%)
Successfully overlapped alignments: 1 (0.02%)
Left parts with too small N-region (failed to extract k-mer): 722 (15.52%)
Extracted k-mer diversity: 1716
Dropped due to wildcard in k-mer: 0 (0%)
Dropped due to too short NRegion parts in overlap: 113 (2.43%)
Dropped overlaps with empty N region due to no complete NDN coverage: 0 (0%)
Number of left-side alignments: 520 (11.18%)
Number of right-side alignments: 2336 (50.23%)
Complex overlaps: 0 (0%)
Over-overlaps: 0 (0%)
Partial alignments written to output: 3414 (73.4%)

>>>>>>>>>>>>>>>>>>>>>>> mixcr extend <<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr extend --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.extend.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.extend.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.passembled.2.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.extended.vdjca
Extending alignments: 0%
====================== report: extend ======================
Analysis time: 51ms
Extended alignments count: 371 (7.98%)
V extensions total: 58 (1.25%)
V extensions with merged targets: 0 (0%)
J extensions total: 316 (6.8%)
J extensions with merged targets: 16 (0.34%)
V+J extensions: 3 (0.06%)
Mean V extension length: 3.586206896551724
Mean J extension length: 6.006329113924051

>>>>>>>>>>>>>>>>>>>>>> mixcr assemble <<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr assemble --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.assemble.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.assemble.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.extended.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.clna
Initialization: progress unknown
Initialization: progress unknown
===================== report: assemble =====================
Analysis time: 290ms
Final clonotype count: 209
Reads used in clonotypes, percent of total: 1409 (0.01%)
Average number of reads per clonotype: 6.74
Reads dropped due to the lack of a clone sequence, percent of total: 3424 (0.02%)
Reads dropped due to a too short clonal sequence, percent of total: 0 (0%)
Reads dropped due to low quality, percent of total: 0 (0%)
Reads dropped due to failed mapping, percent of total: 3 (0%)
Reads dropped with low quality clones, percent of total: 0 (0%)
Aligned reads processed: 1412
Reads used in clonotypes before clustering, percent of total: 1409 (0.01%)
Number of reads used as a core, percent of used: 1402 (99.5%)
Mapped low quality reads, percent of used: 7 (0.5%)
Reads clustered in PCR error correction, percent of used: 0 (0%)
Reads pre-clustered due to the similar VJC-lists, percent of used: 0 (0%)
Clonotypes dropped as low quality: 0
Clonotypes eliminated by PCR error correction: 0
Clonotypes pre-clustered due to the similar VJC-lists: 0
Clones dropped in post filtering: 0 (0%)
Reads dropped in post filtering: 0.0 (0%)
Reads filtered by tag prefix: 0 (0%)
TRA chains: 83 (39.71%)
TRA non-functional: 23 (27.71%)
TRB chains: 124 (59.33%)
TRB non-functional: 4 (3.23%)
TRG chains: 2 (0.96%)
TRG non-functional: 1 (50%)

>>>>>>>>>>>>>>>>>> mixcr assembleContigs <<<<<<<<<<<<<<<<<<
Running:
mixcr assembleContigs --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.assembleContigs.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.assembleContigs.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.clna /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.contigs.clns
Assembling contigs: 0%
================= report: assembleContigs =================
Analysis time: 274ms
Initial clonotype count: 209
Final clonotype count: 209 (100%)
Canceled assemblies: 0 (0%)
Clones filtered out after assemble: 0 (0%)
Number of premature termination assembly events, percent of number of initial clonotypes: 0.0 (0%)
Longest contig length: 318
Clustered variants: 0 (0%)
Reads in clustered variants: 0.0 (0%)
Reads in divided (newly created) clones: 0.0 (0%)
Clones with ambiguous letters in splitting region: 0 (0%)
Reads in clones with ambiguous letters in splitting region: 0.0 (0%)
Average number of ambiguous letters per clone with ambiguous letters in splitting region: NaN
Clones combined after clustering: 0 (0%)

>>>>>>>>>>>>>>>>>>>>>>>>> mixcr qc <<<<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr qc --print-to-stdout /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.contigs.clns /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.qc.txt /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.qc.json

  Alignments dropped due to low sequence quality: 0.0% [OK]

>>>>>>>>>>>>>>>>>>>> mixcr exportClones <<<<<<<<<<<<<<<<<<<<
Running:
mixcr exportClones /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.contigs.clns /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/F01_rnaseq.clones.tsv
Exporting TRB
Exporting clones: 0.8%
Filtered 124 of 209 clones (40.67%).
Filtered 859.0 of 1409.0 reads (39.03%).
Exporting TRA
Exporting clones: 1.2%
Filtered 83 of 209 clones (60.29%).
Filtered 525.0 of 1409.0 reads (62.74%).
Exporting TRG
Exporting clones: 100%
Filtered 2 of 209 clones (99.04%).
Filtered 25.0 of 1409.0 reads (98.23%).
Analysis finished successfully.
  Done: F01

========================================
Processing: D01 (30)
========================================

>>>>>>>>>>>>>>>>>>>>>>> mixcr align <<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr align --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.align.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.align.report.json --preset rna-seq --save-output-file-names /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.align.list.tsv --species hsa /home/lopez_cascales/LABMEMBERS/RENO/fastq/I25-1078-D01_S8_R1_001.fastq.gz /home/lopez_cascales/LABMEMBERS/RENO/fastq/I25-1078-D01_S8_R2_001.fastq.gz /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.alignments.vdjca
Alignment: 0%
Alignment: 10.5%  ETA: 00:02:08
Alignment: 20.6%  ETA: 00:01:42
Alignment: 31.2%  ETA: 00:01:30
Alignment: 41.8%  ETA: 00:01:16
Alignment: 52.4%  ETA: 00:01:02
Alignment: 62.6%  ETA: 00:00:51
Alignment: 72.7%  ETA: 00:00:37
Alignment: 82.8%  ETA: 00:00:23
Alignment: 92.8%  ETA: 00:00:09
====================== report: align ======================
Analysis time: 2.28m
Total sequencing reads: 24042301
Successfully aligned reads: 4667 (0.02%)
Coverage (percent of successfully aligned):
  CDR3: 482 (10.33%)
  FR3_TO_FR4: 0 (0%)
  CDR2_TO_FR4: 0 (0%)
  FR2_TO_FR4: 0 (0%)
  CDR1_TO_FR4: 0 (0%)
  VDJRegion: 0 (0%)
Alignment failed: no hits (not TCR/IG?): 23459397 (97.58%)
Alignment failed: no CDR3 parts: 8660 (0.04%)
Alignment failed: low total score: 569577 (2.37%)
Overlapped: 11007159 (45.78%)
Overlapped and aligned: 1200 (0%)
Overlapped and not aligned: 11005959 (45.78%)
Alignment-aided overlaps, percent of overlapped and aligned: 169 (14.08%)
No CDR3 parts alignments, percent of successfully aligned: 82 (1.76%)
Partial aligned reads, percent of successfully aligned: 4103 (87.92%)
J gene chimeras: 16 (0%)
Paired-end alignment conflicts eliminated: 9 (0%)
Realigned with forced non-floating bound: 26070622 (108.44%)
Realigned with forced non-floating right bound in left read: 539 (0%)
Realigned with forced non-floating left bound in right read: 539 (0%)
TRA chains: 786 (16.84%)
TRA non-functional: 0 (0%)
TRB chains: 3780 (80.99%)
TRB non-functional: 35 (0.93%)
TRAD chains: 101 (2.16%)
TRAD non-functional: 0 (0%)
Trimming report:
  R1 reads trimmed left: 64353 (0.27%)
  R1 reads trimmed right: 3 (0%)
  Average R1 nucleotides trimmed left: 0.008020322181308687
  Average R1 nucleotides trimmed right: 1.2478006992758305E-7
  R2 reads trimmed left: 0 (0%)
  R2 reads trimmed right: 18 (0%)
  Average R2 nucleotides trimmed left: 0.0
  Average R2 nucleotides trimmed right: 4.242522377537824E-6

>>>>>>>>>>>>>>>>>> mixcr assemblePartial <<<<<<<<<<<<<<<<<<
Running:
mixcr assemblePartial --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.assemblePartial.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.assemblePartial.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.alignments.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.passembled.1.vdjca
Building index: 0%
Searching for overlaps: 0%
================= report: assemblePartial =================
Analysis time: 411ms
Total alignments analysed: 4667
Number of output alignments: 4334 (92.86%)
Alignments already with CDR3 (no overlapping is performed): 804 (17.23%)
Successfully overlapped alignments: 333 (7.14%)
Left parts with too small N-region (failed to extract k-mer): 506 (10.84%)
Extracted k-mer diversity: 429
Dropped due to wildcard in k-mer: 0 (0%)
Dropped due to too short NRegion parts in overlap: 304 (6.51%)
Dropped overlaps with empty N region due to no complete NDN coverage: 0 (0%)
Number of left-side alignments: 918 (19.67%)
Number of right-side alignments: 2597 (55.65%)
Complex overlaps: 0 (0%)
Over-overlaps: 0 (0%)
Partial alignments written to output: 3197 (68.5%)

>>>>>>>>>>>>>>>>>> mixcr assemblePartial <<<<<<<<<<<<<<<<<<
Running:
mixcr assemblePartial --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.assemblePartial.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.assemblePartial.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.passembled.1.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.passembled.2.vdjca
Building index: 0%
Searching for overlaps: 0%
================= report: assemblePartial =================
Analysis time: 108ms
Total alignments analysed: 4334
Number of output alignments: 4334 (100%)
Alignments already with CDR3 (no overlapping is performed): 1137 (26.23%)
Successfully overlapped alignments: 0 (0%)
Left parts with too small N-region (failed to extract k-mer): 494 (11.4%)
Extracted k-mer diversity: 335
Dropped due to wildcard in k-mer: 0 (0%)
Dropped due to too short NRegion parts in overlap: 221 (5.1%)
Dropped overlaps with empty N region due to no complete NDN coverage: 0 (0%)
Number of left-side alignments: 578 (13.34%)
Number of right-side alignments: 2250 (51.92%)
Complex overlaps: 0 (0%)
Over-overlaps: 0 (0%)
Partial alignments written to output: 3197 (73.77%)

>>>>>>>>>>>>>>>>>>>>>>> mixcr extend <<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr extend --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.extend.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.extend.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.passembled.2.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.extended.vdjca
Extending alignments: 0%
====================== report: extend ======================
Analysis time: 57ms
Extended alignments count: 324 (7.48%)
V extensions total: 43 (0.99%)
V extensions with merged targets: 1 (0.02%)
J extensions total: 281 (6.48%)
J extensions with merged targets: 5 (0.12%)
V+J extensions: 0 (0%)
Mean V extension length: 3.2093023255813953
Mean J extension length: 4.480427046263345

>>>>>>>>>>>>>>>>>>>>>> mixcr assemble <<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr assemble --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.assemble.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.assemble.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.extended.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.clna
Initialization: progress unknown
Initialization: progress unknown
===================== report: assemble =====================
Analysis time: 207ms
Final clonotype count: 20
Reads used in clonotypes, percent of total: 1463 (0.01%)
Average number of reads per clonotype: 73.15
Reads dropped due to the lack of a clone sequence, percent of total: 3196 (0.01%)
Reads dropped due to a too short clonal sequence, percent of total: 0 (0%)
Reads dropped due to low quality, percent of total: 0 (0%)
Reads dropped due to failed mapping, percent of total: 1 (0%)
Reads dropped with low quality clones, percent of total: 0 (0%)
Aligned reads processed: 1471
Reads used in clonotypes before clustering, percent of total: 1470 (0.01%)
Number of reads used as a core, percent of used: 1454 (98.91%)
Mapped low quality reads, percent of used: 16 (1.09%)
Reads clustered in PCR error correction, percent of used: 7 (0.48%)
Reads pre-clustered due to the similar VJC-lists, percent of used: 0 (0%)
Clonotypes dropped as low quality: 0
Clonotypes eliminated by PCR error correction: 5
Clonotypes pre-clustered due to the similar VJC-lists: 0
Clones dropped in post filtering: 0 (0%)
Reads dropped in post filtering: 0.0 (0%)
Reads filtered by tag prefix: 0 (0%)
TRA chains: 5 (25%)
TRA non-functional: 0 (0%)
TRB chains: 15 (75%)
TRB non-functional: 2 (13.33%)

>>>>>>>>>>>>>>>>>> mixcr assembleContigs <<<<<<<<<<<<<<<<<<
Running:
mixcr assembleContigs --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.assembleContigs.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.assembleContigs.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.clna /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.contigs.clns
Assembling contigs: 0%
================= report: assembleContigs =================
Analysis time: 189ms
Initial clonotype count: 20
Final clonotype count: 20 (100%)
Canceled assemblies: 0 (0%)
Clones filtered out after assemble: 0 (0%)
Number of premature termination assembly events, percent of number of initial clonotypes: 0.0 (0%)
Longest contig length: 390
Clustered variants: 0 (0%)
Reads in clustered variants: 0.0 (0%)
Reads in divided (newly created) clones: 0.0 (0%)
Clones with ambiguous letters in splitting region: 0 (0%)
Reads in clones with ambiguous letters in splitting region: 0.0 (0%)
Average number of ambiguous letters per clone with ambiguous letters in splitting region: NaN
Clones combined after clustering: 0 (0%)

>>>>>>>>>>>>>>>>>>>>>>>>> mixcr qc <<<<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr qc --print-to-stdout /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.contigs.clns /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.qc.txt /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.qc.json

  Alignments dropped due to low sequence quality: 0.0% [OK]

>>>>>>>>>>>>>>>>>>>> mixcr exportClones <<<<<<<<<<<<<<<<<<<<
Running:
mixcr exportClones /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.contigs.clns /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/D01_rnaseq.clones.tsv
Exporting TRB
Exporting clones: 6.7%
Filtered 15 of 20 clones (25%).
Filtered 1157.0 of 1463.0 reads (20.92%).
Exporting TRA
Exporting clones: 20%
Filtered 5 of 20 clones (75%).
Filtered 306.0 of 1463.0 reads (79.08%).
Analysis finished successfully.
  Done: D01

========================================
Processing: G01 (46)
========================================

>>>>>>>>>>>>>>>>>>>>>>> mixcr align <<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr align --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.align.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.align.report.json --preset rna-seq --save-output-file-names /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.align.list.tsv --species hsa /home/lopez_cascales/LABMEMBERS/RENO/fastq/I25-1078-G01_S11_R1_001.fastq.gz /home/lopez_cascales/LABMEMBERS/RENO/fastq/I25-1078-G01_S11_R2_001.fastq.gz /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.alignments.vdjca
Alignment: 0%
Alignment: 10.8%  ETA: 00:02:12
Alignment: 21.5%  ETA: 00:01:42
Alignment: 32.2%  ETA: 00:01:29
Alignment: 42.9%  ETA: 00:01:14
Alignment: 53.4%  ETA: 00:01:02
Alignment: 63.5%  ETA: 00:00:50
Alignment: 73.7%  ETA: 00:00:36
Alignment: 84.5%  ETA: 00:00:21
Alignment: 95.1%  ETA: 00:00:06
====================== report: align ======================
Analysis time: 2.29m
Total sequencing reads: 21726191
Successfully aligned reads: 4844 (0.02%)
Coverage (percent of successfully aligned):
  CDR3: 573 (11.83%)
  FR3_TO_FR4: 0 (0%)
  CDR2_TO_FR4: 0 (0%)
  FR2_TO_FR4: 0 (0%)
  CDR1_TO_FR4: 0 (0%)
  VDJRegion: 0 (0%)
Alignment failed: no hits (not TCR/IG?): 21074733 (97%)
Alignment failed after alignment-aided overlap: 3 (0%)
Alignment failed: no CDR3 parts: 10370 (0.05%)
Alignment failed: low total score: 636241 (2.93%)
Overlapped: 9665007 (44.49%)
Overlapped and aligned: 1442 (0.01%)
Overlapped and not aligned: 9663565 (44.48%)
Alignment-aided overlaps, percent of overlapped and aligned: 212 (14.7%)
No CDR3 parts alignments, percent of successfully aligned: 130 (2.68%)
Partial aligned reads, percent of successfully aligned: 4141 (85.49%)
V gene chimeras: 4 (0%)
J gene chimeras: 6 (0%)
Paired-end alignment conflicts eliminated: 9 (0%)
Realigned with forced non-floating bound: 24122792 (111.03%)
Realigned with forced non-floating right bound in left read: 878 (0%)
Realigned with forced non-floating left bound in right read: 878 (0%)
TRA chains: 1311 (27.06%)
TRA non-functional: 25 (1.91%)
TRB chains: 3067 (63.32%)
TRB non-functional: 15 (0.49%)
TRG chains: 6 (0.12%)
TRG non-functional: 0 (0%)
TRAD chains: 449 (9.27%)
TRAD non-functional: 0 (0%)
IGH chains: 9 (0.19%)
IGH non-functional: 0 (0%)
IGL chains: 2 (0.04%)
IGL non-functional: 0 (0%)
Trimming report:
  R1 reads trimmed left: 59487 (0.27%)
  R1 reads trimmed right: 6 (0%)
  Average R1 nucleotides trimmed left: 0.00820166774746664
  Average R1 nucleotides trimmed right: 3.221917730540066E-7
  R2 reads trimmed left: 0 (0%)
  R2 reads trimmed right: 23 (0%)
  Average R2 nucleotides trimmed left: 0.0
  Average R2 nucleotides trimmed right: 6.673972441832993E-6

>>>>>>>>>>>>>>>>>> mixcr assemblePartial <<<<<<<<<<<<<<<<<<
Running:
mixcr assemblePartial --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.assemblePartial.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.assemblePartial.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.alignments.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.passembled.1.vdjca
Building index: 0%
Searching for overlaps: 0%
================= report: assemblePartial =================
Analysis time: 225ms
Total alignments analysed: 4844
Number of output alignments: 4692 (96.86%)
Alignments already with CDR3 (no overlapping is performed): 981 (20.25%)
Successfully overlapped alignments: 152 (3.14%)
Left parts with too small N-region (failed to extract k-mer): 822 (16.97%)
Extracted k-mer diversity: 4020
Dropped due to wildcard in k-mer: 0 (0%)
Dropped due to too short NRegion parts in overlap: 57 (1.18%)
Dropped overlaps with empty N region due to no complete NDN coverage: 7 (0.14%)
Number of left-side alignments: 790 (16.31%)
Number of right-side alignments: 2451 (50.6%)
Complex overlaps: 0 (0%)
Over-overlaps: 1 (0.02%)
Partial alignments written to output: 3559 (73.47%)

>>>>>>>>>>>>>>>>>> mixcr assemblePartial <<<<<<<<<<<<<<<<<<
Running:
mixcr assemblePartial --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.assemblePartial.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.assemblePartial.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.passembled.1.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.passembled.2.vdjca
Building index: 0%
Searching for overlaps: 0.1%
================= report: assemblePartial =================
Analysis time: 59ms
Total alignments analysed: 4692
Number of output alignments: 4691 (99.98%)
Alignments already with CDR3 (no overlapping is performed): 1132 (24.13%)
Successfully overlapped alignments: 1 (0.02%)
Left parts with too small N-region (failed to extract k-mer): 813 (17.33%)
Extracted k-mer diversity: 3249
Dropped due to wildcard in k-mer: 0 (0%)
Dropped due to too short NRegion parts in overlap: 47 (1%)
Dropped overlaps with empty N region due to no complete NDN coverage: 7 (0.15%)
Number of left-side alignments: 636 (13.55%)
Number of right-side alignments: 2291 (48.83%)
Complex overlaps: 0 (0%)
Over-overlaps: 1 (0.02%)
Partial alignments written to output: 3558 (75.83%)

>>>>>>>>>>>>>>>>>>>>>>> mixcr extend <<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr extend --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.extend.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.extend.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.passembled.2.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.extended.vdjca
Extending alignments: 0%
====================== report: extend ======================
Analysis time: 59ms
Extended alignments count: 418 (8.91%)
V extensions total: 69 (1.47%)
V extensions with merged targets: 0 (0%)
J extensions total: 349 (7.44%)
J extensions with merged targets: 18 (0.38%)
V+J extensions: 0 (0%)
Mean V extension length: 3.652173913043478
Mean J extension length: 5.6189111747851

>>>>>>>>>>>>>>>>>>>>>> mixcr assemble <<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr assemble --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.assemble.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.assemble.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.extended.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.clna
Initialization: progress unknown
Initialization: progress unknown
===================== report: assemble =====================
Analysis time: 417ms
Final clonotype count: 429
Reads used in clonotypes, percent of total: 1266 (0.01%)
Average number of reads per clonotype: 2.95
Reads dropped due to the lack of a clone sequence, percent of total: 3556 (0.02%)
Reads dropped due to a too short clonal sequence, percent of total: 0 (0%)
Reads dropped due to low quality, percent of total: 0 (0%)
Reads dropped due to failed mapping, percent of total: 22 (0%)
Reads dropped with low quality clones, percent of total: 0 (0%)
Aligned reads processed: 1288
Reads used in clonotypes before clustering, percent of total: 1266 (0.01%)
Number of reads used as a core, percent of used: 1236 (97.63%)
Mapped low quality reads, percent of used: 30 (2.37%)
Reads clustered in PCR error correction, percent of used: 0 (0%)
Reads pre-clustered due to the similar VJC-lists, percent of used: 0 (0%)
Clonotypes dropped as low quality: 0
Clonotypes eliminated by PCR error correction: 0
Clonotypes pre-clustered due to the similar VJC-lists: 0
Clones dropped in post filtering: 0 (0%)
Reads dropped in post filtering: 0.0 (0%)
Reads filtered by tag prefix: 0 (0%)
TRA chains: 189 (44.06%)
TRA non-functional: 21 (11.11%)
TRB chains: 240 (55.94%)
TRB non-functional: 9 (3.75%)

>>>>>>>>>>>>>>>>>> mixcr assembleContigs <<<<<<<<<<<<<<<<<<
Running:
mixcr assembleContigs --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.assembleContigs.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.assembleContigs.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.clna /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.contigs.clns
Assembling contigs: 0%
================= report: assembleContigs =================
Analysis time: 285ms
Initial clonotype count: 429
Final clonotype count: 429 (100%)
Canceled assemblies: 0 (0%)
Clones filtered out after assemble: 0 (0%)
Number of premature termination assembly events, percent of number of initial clonotypes: 0.0 (0%)
Longest contig length: 277
Clustered variants: 0 (0%)
Reads in clustered variants: 0.0 (0%)
Reads in divided (newly created) clones: 0.0 (0%)
Clones with ambiguous letters in splitting region: 0 (0%)
Reads in clones with ambiguous letters in splitting region: 0.0 (0%)
Average number of ambiguous letters per clone with ambiguous letters in splitting region: NaN
Clones combined after clustering: 0 (0%)

>>>>>>>>>>>>>>>>>>>>>>>>> mixcr qc <<<<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr qc --print-to-stdout /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.contigs.clns /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.qc.txt /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.qc.json

  Alignments dropped due to low sequence quality: 0.0% [OK]

>>>>>>>>>>>>>>>>>>>> mixcr exportClones <<<<<<<<<<<<<<<<<<<<
Running:
mixcr exportClones /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.contigs.clns /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/G01_rnaseq.clones.tsv
Exporting TRB
Exporting clones: 0.4%
Filtered 240 of 429 clones (44.06%).
Filtered 731.0 of 1266.0 reads (42.26%).
Exporting TRA
Exporting clones: 1.1%
Filtered 189 of 429 clones (55.94%).
Filtered 535.0 of 1266.0 reads (57.74%).
Analysis finished successfully.
  Done: G01

========================================
Processing: H01 (122)
========================================

>>>>>>>>>>>>>>>>>>>>>>> mixcr align <<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr align --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.align.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.align.report.json --preset rna-seq --save-output-file-names /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.align.list.tsv --species hsa /home/lopez_cascales/LABMEMBERS/RENO/fastq/I25-1078-H01_S12_R1_001.fastq.gz /home/lopez_cascales/LABMEMBERS/RENO/fastq/I25-1078-H01_S12_R2_001.fastq.gz /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.alignments.vdjca
Alignment: 0%
Alignment: 10.7%  ETA: 00:02:14
Alignment: 21.5%  ETA: 00:01:49
Alignment: 32%  ETA: 00:01:37
Alignment: 42.1%  ETA: 00:01:20
Alignment: 52.6%  ETA: 00:01:07
Alignment: 62.8%  ETA: 00:00:54
Alignment: 73.2%  ETA: 00:00:38
Alignment: 83.4%  ETA: 00:00:24
Alignment: 93.6%  ETA: 00:00:09
====================== report: align ======================
Analysis time: 2.42m
Total sequencing reads: 23679418
Successfully aligned reads: 5961 (0.03%)
Coverage (percent of successfully aligned):
  CDR3: 747 (12.53%)
  FR3_TO_FR4: 0 (0%)
  CDR2_TO_FR4: 0 (0%)
  FR2_TO_FR4: 0 (0%)
  CDR1_TO_FR4: 0 (0%)
  VDJRegion: 0 (0%)
Alignment failed: no hits (not TCR/IG?): 23016626 (97.2%)
Alignment failed after alignment-aided overlap: 2 (0%)
Alignment failed: no CDR3 parts: 12226 (0.05%)
Alignment failed: low total score: 644603 (2.72%)
Overlapped: 12733095 (53.77%)
Overlapped and aligned: 2190 (0.01%)
Overlapped and not aligned: 12730905 (53.76%)
Alignment-aided overlaps, percent of overlapped and aligned: 294 (13.42%)
No CDR3 parts alignments, percent of successfully aligned: 116 (1.95%)
Partial aligned reads, percent of successfully aligned: 5098 (85.52%)
J gene chimeras: 2 (0%)
Paired-end alignment conflicts eliminated: 9 (0%)
Realigned with forced non-floating bound: 21893234 (92.46%)
Realigned with forced non-floating right bound in left read: 1083 (0%)
Realigned with forced non-floating left bound in right read: 1083 (0%)
TRA chains: 1724 (28.92%)
TRA non-functional: 34 (1.97%)
TRB chains: 3791 (63.6%)
TRB non-functional: 30 (0.79%)
TRG chains: 6 (0.1%)
TRG non-functional: 0 (0%)
TRAD chains: 439 (7.36%)
TRAD non-functional: 0 (0%)
IGH chains: 1 (0.02%)
IGH non-functional: 0 (0%)
Trimming report:
  R1 reads trimmed left: 65361 (0.28%)
  R1 reads trimmed right: 3 (0%)
  Average R1 nucleotides trimmed left: 0.00826933330878318
  Average R1 nucleotides trimmed right: 1.2669230299494691E-7
  R2 reads trimmed left: 1 (0%)
  R2 reads trimmed right: 8 (0%)
  Average R2 nucleotides trimmed left: 8.446153532996461E-8
  Average R2 nucleotides trimmed right: 2.9983845042137438E-6

>>>>>>>>>>>>>>>>>> mixcr assemblePartial <<<<<<<<<<<<<<<<<<
Running:
mixcr assemblePartial --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.assemblePartial.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.assemblePartial.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.alignments.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.passembled.1.vdjca
Building index: 0%
Searching for overlaps: 0%
================= report: assemblePartial =================
Analysis time: 347ms
Total alignments analysed: 5961
Number of output alignments: 5619 (94.26%)
Alignments already with CDR3 (no overlapping is performed): 1232 (20.67%)
Successfully overlapped alignments: 342 (5.74%)
Left parts with too small N-region (failed to extract k-mer): 942 (15.8%)
Extracted k-mer diversity: 1415
Dropped due to wildcard in k-mer: 0 (0%)
Dropped due to too short NRegion parts in overlap: 237 (3.98%)
Dropped overlaps with empty N region due to no complete NDN coverage: 9 (0.15%)
Number of left-side alignments: 947 (15.89%)
Number of right-side alignments: 3086 (51.77%)
Complex overlaps: 0 (0%)
Over-overlaps: 0 (0%)
Partial alignments written to output: 4045 (67.86%)

>>>>>>>>>>>>>>>>>> mixcr assemblePartial <<<<<<<<<<<<<<<<<<
Running:
mixcr assemblePartial --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.assemblePartial.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.assemblePartial.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.passembled.1.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.passembled.2.vdjca
Building index: 0%
Searching for overlaps: 0.1%
================= report: assemblePartial =================
Analysis time: 104ms
Total alignments analysed: 5619
Number of output alignments: 5619 (100%)
Alignments already with CDR3 (no overlapping is performed): 1574 (28.01%)
Successfully overlapped alignments: 0 (0%)
Left parts with too small N-region (failed to extract k-mer): 926 (16.48%)
Extracted k-mer diversity: 969
Dropped due to wildcard in k-mer: 0 (0%)
Dropped due to too short NRegion parts in overlap: 173 (3.08%)
Dropped overlaps with empty N region due to no complete NDN coverage: 4 (0.07%)
Number of left-side alignments: 602 (10.71%)
Number of right-side alignments: 2726 (48.51%)
Complex overlaps: 0 (0%)
Over-overlaps: 0 (0%)
Partial alignments written to output: 4045 (71.99%)

>>>>>>>>>>>>>>>>>>>>>>> mixcr extend <<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr extend --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.extend.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.extend.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.passembled.2.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.extended.vdjca
Extending alignments: 0%
====================== report: extend ======================
Analysis time: 76ms
Extended alignments count: 484 (8.61%)
V extensions total: 82 (1.46%)
V extensions with merged targets: 0 (0%)
J extensions total: 402 (7.15%)
J extensions with merged targets: 12 (0.21%)
V+J extensions: 0 (0%)
Mean V extension length: 3.0365853658536586
Mean J extension length: 5.977611940298507

>>>>>>>>>>>>>>>>>>>>>> mixcr assemble <<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr assemble --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.assemble.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.assemble.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.extended.vdjca /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.clna
Initialization: progress unknown
Initialization: progress unknown
===================== report: assemble =====================
Analysis time: 279ms
Final clonotype count: 93
Reads used in clonotypes, percent of total: 1910 (0.01%)
Average number of reads per clonotype: 20.54
Reads dropped due to the lack of a clone sequence, percent of total: 4046 (0.02%)
Reads dropped due to a too short clonal sequence, percent of total: 0 (0%)
Reads dropped due to low quality, percent of total: 0 (0%)
Reads dropped due to failed mapping, percent of total: 2 (0%)
Reads dropped with low quality clones, percent of total: 1 (0%)
Aligned reads processed: 1915
Reads used in clonotypes before clustering, percent of total: 1912 (0.01%)
Number of reads used as a core, percent of used: 1879 (98.27%)
Mapped low quality reads, percent of used: 33 (1.73%)
Reads clustered in PCR error correction, percent of used: 2 (0.1%)
Reads pre-clustered due to the similar VJC-lists, percent of used: 0 (0%)
Clonotypes dropped as low quality: 1
Clonotypes eliminated by PCR error correction: 2
Clonotypes pre-clustered due to the similar VJC-lists: 0
Clones dropped in post filtering: 0 (0%)
Reads dropped in post filtering: 0.0 (0%)
Reads filtered by tag prefix: 0 (0%)
TRA chains: 35 (37.63%)
TRA non-functional: 5 (14.29%)
TRB chains: 58 (62.37%)
TRB non-functional: 5 (8.62%)

>>>>>>>>>>>>>>>>>> mixcr assembleContigs <<<<<<<<<<<<<<<<<<
Running:
mixcr assembleContigs --report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.assembleContigs.report.txt --json-report /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.assembleContigs.report.json /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.clna /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.contigs.clns
Assembling contigs: 0%
================= report: assembleContigs =================
Analysis time: 244ms
Initial clonotype count: 93
Final clonotype count: 93 (100%)
Canceled assemblies: 0 (0%)
Clones filtered out after assemble: 0 (0%)
Number of premature termination assembly events, percent of number of initial clonotypes: 0.0 (0%)
Longest contig length: 388
Clustered variants: 0 (0%)
Reads in clustered variants: 0.0 (0%)
Reads in divided (newly created) clones: 0.0 (0%)
Clones with ambiguous letters in splitting region: 0 (0%)
Reads in clones with ambiguous letters in splitting region: 0.0 (0%)
Average number of ambiguous letters per clone with ambiguous letters in splitting region: NaN
Clones combined after clustering: 0 (0%)

>>>>>>>>>>>>>>>>>>>>>>>>> mixcr qc <<<<<<<<<<<<<<<<<<<<<<<<<
Running:
mixcr qc --print-to-stdout /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.contigs.clns /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.qc.txt /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.qc.json

  Alignments dropped due to low sequence quality: 0.0% [OK]

>>>>>>>>>>>>>>>>>>>> mixcr exportClones <<<<<<<<<<<<<<<<<<<<
Running:
mixcr exportClones /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.contigs.clns /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output/H01_rnaseq.clones.tsv
Exporting TRB
Exporting clones: 1.7%
Filtered 58 of 93 clones (37.63%).
Filtered 1150.0 of 1910.0 reads (39.79%).
Exporting TRA
Exporting clones: 2.9%
Filtered 35 of 93 clones (62.37%).
Filtered 760.0 of 1910.0 reads (60.21%).
Analysis finished successfully.
  Done: H01

========================================
All samples processed.
Output in: /home/lopez_cascales/LABMEMBERS/RENO/mixcr_rnaseq_output
========================================

=== Alignment summary ===

Next step: check alignment rates above.
  If > 0% aligned → run Figure6C_mixcr_rnaseq_analysis.R
  If still 0%     → bulk RNA-seq truly has no recoverable TCR reads
