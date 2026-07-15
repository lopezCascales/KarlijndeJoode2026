################################################################################
# Figure 6C — TCR clonotype analysis from MiXCR rna-seq output
# Real clonotypes detected with preset rna-seq (not vergani)
#
# 4-1BB+ : C01, E01, F01
# 4-1BB- : D01, G01, H01
#
# Input files: mixcr_rnaseq_output/*_rnaseq.clones.tsv
# Run from: ~/LABMEMBERS/RENO/
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(patchwork)
library(vegan)

OUTDIR  <- "mixcr_rnaseq_output"

sample_meta <- tibble(
  Sample = c("C01", "E01", "F01", "D01", "G01", "H01"),
  Group  = c("4-1BB+", "4-1BB+", "4-1BB+", "4-1BB-", "4-1BB-", "4-1BB-")
) %>% mutate(Group = factor(Group, levels = c("4-1BB-", "4-1BB+")))

group_cols <- c("4-1BB+" = "#C0392B", "4-1BB-" = "#2980B9")

# ── 1. Load clonotype files ───────────────────────────────────────────────────
# MiXCR exports separate files per chain: clones_TRA.tsv and clones_TRB.tsv
# We load both chains and combine

load_clones <- function(sample_id, chain, outdir) {
  f <- file.path(outdir, paste0(sample_id, "_rnaseq.clones_", chain, ".tsv"))
  if (!file.exists(f)) {
    message("  Not found: ", f)
    return(NULL)
  }
  df <- read_tsv(f, show_col_types = FALSE)
  if (nrow(df) == 0) return(NULL)
  df %>%
    mutate(across(everything(), as.character)) %>%  # fuerza todo a character
    mutate(Sample = sample_id, Chain = chain)
}

clones_list <- list()
for (s in sample_meta$Sample) {
  for (ch in c("TRA", "TRB")) {
    cl <- load_clones(s, ch, OUTDIR)
    if (!is.null(cl)) clones_list[[paste(s, ch)]] <- cl
  }
}

clones_df <- bind_rows(clones_list) %>%
  left_join(sample_meta, by = "Sample")

# Busca las columnas correctas
count_col <- colnames(clones_df)[grep("readCount|Read count|reads", 
                                      colnames(clones_df), ignore.case = TRUE)[1]]
freq_col  <- colnames(clones_df)[grep("readFraction|Read fraction|frequency", 
                                      colnames(clones_df), ignore.case = TRUE)[1]]

message("count col: ", count_col)
message("freq col:  ", freq_col)

clones_df <- clones_df %>%
  rename(ReadCount = all_of(count_col),
         Frequency = all_of(freq_col)) %>%
  mutate(ReadCount = as.numeric(ReadCount),
         Frequency = as.numeric(Frequency))

message("Total clonotypes loaded: ", nrow(clones_df))
message("Per sample:")
print(clones_df %>% group_by(Sample, Group, Chain) %>% summarise(N = n(), .groups = "drop"))

# Total clonotypes loaded: 1104
Per sample:
  # A tibble: 12 × 4
  Sample Group  Chain     N
<chr>  <fct>  <chr> <int>
  1 C01    4-1BB+ TRA      98
2 C01    4-1BB+ TRB     247
3 D01    4-1BB- TRA       5
4 D01    4-1BB- TRB      15
5 E01    4-1BB+ TRA       6
6 E01    4-1BB+ TRB       4
7 F01    4-1BB+ TRA      83
8 F01    4-1BB+ TRB     124
9 G01    4-1BB- TRA     189
10 G01    4-1BB- TRB     240
11 H01    4-1BB- TRA      35
12 H01    4-1BB- TRB      58

# ── 2. Get read count column (MiXCR uses 'readCount' or 'Read count') ─────────
# Check column names
message("\nColumn names: ", paste(colnames(clones_df)[1:10], collapse = ", "))

# Column names: cloneId, ReadCount, Frequency, targetSequences, targetQualities, allVHitsWithScore, allDHitsWithScore, allJHitsWithScore, allCHitsWithScore, allVAlignments
 
  
  
# Standardize: find the read count column
count_col <- colnames(clones_df)[grep("readCount|Read count|reads", 
                                       colnames(clones_df), ignore.case = TRUE)[1]]
freq_col  <- colnames(clones_df)[grep("readFraction|Read fraction|frequency", 
                                       colnames(clones_df), ignore.case = TRUE)[1]]

message("Using count column: ", count_col)
message("Using freq column:  ", freq_col)

clones_df <- clones_df %>%
  rename(ReadCount = all_of(count_col),
         Frequency = all_of(freq_col))

# ── 3. Frequency categories ───────────────────────────────────────────────────
cat_levels <- c("rare", "low", "medium", "high", "hyperexpanded")
cat_cols   <- c(
  "rare"          = "#F4D03F",
  "low"           = "#82E0AA",
  "medium"        = "#5DADE2",
  "high"          = "#E59866",
  "hyperexpanded" = "#C0392B"
)

clones_df <- clones_df %>%
  mutate(Category = case_when(
    Frequency > 0      & Frequency <= 1e-05 ~ "rare",
    Frequency > 1e-05  & Frequency <= 1e-04 ~ "low",
    Frequency > 1e-04  & Frequency <= 1e-03 ~ "medium",
    Frequency > 1e-03  & Frequency <= 1e-02 ~ "high",
    Frequency > 1e-02                        ~ "hyperexpanded",
    TRUE                                     ~ NA_character_
  )) %>%
  filter(!is.na(Category)) %>%
  mutate(Category = factor(Category, levels = cat_levels))

# ── 4. PLOT 1 — Number of clonotypes per sample ───────────────────────────────
clone_counts <- clones_df %>%
  group_by(Sample, Group) %>%
  summarise(NumClones = n(), .groups = "drop")

# Summary for error bars
clone_summary <- clone_counts %>%
  group_by(Group) %>%
  summarise(Mean = mean(NumClones), SD = sd(NumClones), .groups = "drop")

p1 <- ggplot(clone_counts, aes(x = Group, y = NumClones, colour = Group)) +
  geom_errorbar(data = clone_summary,
                aes(x = Group, ymin = Mean - SD, ymax = Mean + SD),
                width = 0.15, linewidth = 0.6, colour = "grey40",
                inherit.aes = FALSE) +
  geom_crossbar(data = clone_summary,
                aes(x = Group, y = Mean, ymin = Mean, ymax = Mean),
                width = 0.3, linewidth = 0.6, colour = "grey30",
                inherit.aes = FALSE) +
  geom_jitter(aes(fill = Group), width = 0.08, size = 3,
              shape = 21, colour = "white", stroke = 0.4) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.5,
                     label.y = max(clone_counts$NumClones) * 1.12,
                     size = 3.5) +
  scale_fill_manual(values   = group_cols) +
  scale_colour_manual(values = group_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18)),
                     labels = scales::comma) +
  labs(title = "Number of TCR clonotypes",
       x = NULL, y = "Number of clonotypes") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10, face = "bold",
                                   colour = c("#2980B9", "#C0392B")))

# ── 5. PLOT 2 — Clonotype distribution (%) per group ─────────────────────────
cat_dist <- clones_df %>%
  group_by(Group, Category) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(Pct = N / sum(N) * 100)

p2 <- ggplot(cat_dist, aes(x = Group, y = Pct, fill = Category)) +
  geom_col(width = 0.55, colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = cat_cols, breaks = rev(cat_levels)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                     labels = scales::percent_format(scale = 1)) +
  labs(title = "Clonotype distribution per group",
       x = NULL, y = "% of clonotypes", fill = "Category") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10, face = "bold",
                                   colour = c("#2980B9", "#C0392B")),
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8))

# ── 6. Shannon diversity ──────────────────────────────────────────────────────
shannon_df <- clones_df %>%
  group_by(Sample, Group) %>%
  summarise(Shannon = diversity(ReadCount, index = "shannon"), .groups = "drop")

message("\nShannon diversity:")
print(shannon_df)

shannon_summary <- shannon_df %>%
  group_by(Group) %>%
  summarise(Mean = mean(Shannon), SD = sd(Shannon), .groups = "drop")

p3 <- ggplot(shannon_df, aes(x = Group, y = Shannon, colour = Group)) +
  geom_errorbar(data = shannon_summary,
                aes(x = Group, ymin = Mean - SD, ymax = Mean + SD),
                width = 0.15, linewidth = 0.6, colour = "grey40",
                inherit.aes = FALSE) +
  geom_crossbar(data = shannon_summary,
                aes(x = Group, y = Mean, ymin = Mean, ymax = Mean),
                width = 0.3, linewidth = 0.6, colour = "grey30",
                inherit.aes = FALSE) +
  geom_jitter(aes(fill = Group), width = 0.08, size = 3,
              shape = 21, colour = "white", stroke = 0.4) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.5, size = 3.5) +
  scale_fill_manual(values   = group_cols) +
  scale_colour_manual(values = group_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(title = "TCR clonotype diversity",
       x = NULL, y = "Shannon diversity index") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10, face = "bold",
                                   colour = c("#2980B9", "#C0392B")))
Shannon diversity:
  # A tibble: 6 × 3
  Sample Group  Shannon
<chr>  <fct>    <dbl>
  1 C01    4-1BB+    5.61
2 D01    4-1BB-    2.61
3 E01    4-1BB+    1.44
4 F01    4-1BB+    5.03
5 G01    4-1BB-    5.84
6 H01    4-1BB-    4.20

# ── 7. Save figures ───────────────────────────────────────────────────────────
ggsave("Figure6C_MiXCR_nclones.pdf",      p1, width = 3.0, height = 4.0)
ggsave("Figure6C_MiXCR_distribution.pdf", p2, width = 3.5, height = 4.0)
ggsave("Figure6C_MiXCR_diversity.pdf",    p3, width = 3.0, height = 4.0)

# Main combined figure: clonotype number + distribution (as Reno requested)
combined <- p1 + p2 +
  plot_annotation(
    title = "TCR clonotypes in 4-1BB+ vs 4-1BB- CD4 T cells",
    theme = theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))
  )
ggsave("Figure6C_MiXCR_combined.pdf", combined, width = 6.5, height = 4.2)
ggsave("Figure6C_MiXCR_combined.svg", combined, width = 6.5, height = 4.2)

# ── 8. Summary table ──────────────────────────────────────────────────────────
summary_table <- clone_counts %>%
  left_join(shannon_df, by = c("Sample", "Group")) %>%
  arrange(Group, Sample)

write_csv(summary_table, "Figure6C_MiXCR_summary.csv")
print(summary_table)

message("\nDone!")
message("Main figure: Figure6C_MiXCR_combined.pdf")
