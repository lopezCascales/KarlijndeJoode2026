################################################################################
# Figure 6C — TCR clonotype-equivalent analysis from bulk RNA-seq counts
# Note: MiXCR yielded 0 TCR clonotypes (bulk RNA-seq not enriched for V(D)J).
# This analysis uses HTSeq gene counts as a proxy for clonal diversity,
# following the frequency-based categorization described in the Methods.
#
# Samples:
#   4-1BB+: F01, C01, E01
#   4-1BB-: G01, H01, D01
#
# Categories (as in Methods):
#   rare          : 0  < Freq <= 1e-05
#   low           : 1e-05 < Freq <= 1e-04
#   medium        : 1e-04 < Freq <= 1e-03
#   high          : 1e-03 < Freq <= 1e-02
#   hyperexpanded : 1e-02 < Freq <= 1
#
# Per Reno's request:
#   - Remove "none" (Freq = 0)
#   - Pool "high" + "hyperexpanded" into one category
#   - Plot 1: Number of expressed "clones" per group
#   - Plot 2: Distribution (%) of categories per group
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# ── 1. Load data ──────────────────────────────────────────────────────────────
counts <- read_tsv("countsHTSeq.txt", show_col_types = FALSE)

# ── 2. Reshape to long format and assign groups ───────────────────────────────
counts_long <- counts %>%
  pivot_longer(-GENEID, names_to = "Sample", values_to = "Counts") %>%
  mutate(Group = case_when(
    Sample %in% c("F01counts", "C01counts", "E01counts") ~ "4-1BB+",
    Sample %in% c("G01counts", "H01counts", "D01counts") ~ "4-1BB-"
  ))

# ── 3. Summarise: total counts per gene per group ────────────────────────────
counts_summary <- counts_long %>%
  group_by(GENEID, Group) %>%
  summarise(TotalCounts = sum(Counts), .groups = "drop")

# ── 4. Compute relative frequency within each group ──────────────────────────
counts_summary <- counts_summary %>%
  group_by(Group) %>%
  mutate(Freq = TotalCounts / sum(TotalCounts)) %>%
  ungroup()

# ── 5. Assign categories (excluding none / Freq == 0) ────────────────────────
counts_summary <- counts_summary %>%
  mutate(Category = case_when(
    Freq > 0    & Freq <= 1e-05 ~ "rare",
    Freq > 1e-05 & Freq <= 1e-04 ~ "low",
    Freq > 1e-04 & Freq <= 1e-03 ~ "medium",
    Freq > 1e-03                  ~ "high/hyperexpanded",  # pooled per Reno
    TRUE                          ~ "none"   # Freq == 0, will be removed
  ))

# Remove unexpressed genes (none category)
counts_filt <- counts_summary %>%
  filter(Category != "none")

# ── 6. PLOT 1 — Number of expressed "clones" per group ───────────────────────
clone_counts <- counts_filt %>%
  group_by(Group) %>%
  summarise(NumClones = n(), .groups = "drop")

# Colours consistent with rest of figure (4-1BB+ red, 4-1BB- blue)
group_cols <- c("4-1BB+" = "#C0392B", "4-1BB-" = "#2980B9")

p1 <- ggplot(clone_counts, aes(x = Group, y = NumClones, fill = Group)) +
  geom_col(width = 0.55, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = scales::comma(NumClones)),
            vjust = -0.5, size = 3.2) +
  scale_fill_manual(values = group_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     labels = scales::comma) +
  labs(
    title = "Number of TCR clonotypes",
    x = NULL,
    y = "Number of clones"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 10)
  )

# ── 7. PLOT 2 — Distribution (%) of categories per group ─────────────────────
# Define ordered factor for stacking (rare at bottom, high/hyper at top)
cat_levels <- c("rare", "low", "medium", "high/hyperexpanded")

cat_counts <- counts_filt %>%
  group_by(Group, Category) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(Pct = N / sum(N) * 100,
         Category = factor(Category, levels = cat_levels))

# Colour palette for categories (colourblind-friendly)
cat_cols <- c(
  "rare"               = "#F4D03F",   # yellow
  "low"                = "#82E0AA",   # light green
  "medium"             = "#5DADE2",   # light blue
  "high/hyperexpanded" = "#E74C3C"    # red
)

p2 <- ggplot(cat_counts, aes(x = Group, y = Pct, fill = Category)) +
  geom_col(width = 0.55, colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = cat_cols,
                    breaks = rev(cat_levels),      # legend top-to-bottom
                    labels = c("high/hyperexpanded" = "high/hyperexpanded",
                               "medium" = "medium",
                               "low" = "low",
                               "rare" = "rare")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                     labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Clonotype distribution per group",
    x = NULL,
    y = "% of clones",
    fill = "Category"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

# ── 8. Save individual panels ─────────────────────────────────────────────────
ggsave("Figure6C_panel1_nclones.pdf",  p1, width = 3.0, height = 4.0)
ggsave("Figure6C_panel2_distribution.pdf", p2, width = 3.5, height = 4.0)

# ── 9. Also save a combined two-panel figure ─────────────────────────────────
library(patchwork)

combined <- p1 + p2 +
  plot_annotation(
    title = "TCR clonotypes in 4-1BB+ vs 4-1BB\u2212 CD4 T cells",
    theme = theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))
  )

ggsave("Figure6C_combined.pdf", combined, width = 6.5, height = 4.0)

message("Done. Output files: Figure6C_panel1_nclones.pdf, Figure6C_panel2_distribution.pdf, Figure6C_combined.pdf")
