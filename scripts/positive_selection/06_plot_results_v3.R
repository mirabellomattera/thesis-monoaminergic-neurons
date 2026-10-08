.libPaths(c("~/R/library", .libPaths()))
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(ggrepel)

df <- read_tsv("~/tesi_bioinfo/positive_selection_orthodb/positive_selection_results_v3.tsv",
               show_col_types = FALSE)

df <- df %>%
  mutate(
    condition = ifelse(grepl("FGF4", direction), "FGF4 (serotonergic)", "FGF8 (dopaminergic)"),
    is_significant = n_branches_significant > 0
  )

col_fgf4 <- "#1C7293"
col_fgf8 <- "#E07B39"

# ---- FIGURA 1: Scatter con ggrepel ----
genes_to_label <- df %>%
  filter(n_branches_significant >= 2) %>%
  pull(gene)

p1 <- ggplot(df, aes(x = abs(log2FC), y = n_branches_significant,
                      color = condition, shape = condition)) +
  geom_jitter(width = 0.05, height = 0.08, size = 2.5, alpha = 0.7) +
  geom_label_repel(
    data = df %>% filter(gene %in% genes_to_label),
    aes(label = gene),
    size = 2.8, fontface = "italic",
    box.padding = 0.4, point.padding = 0.3,
    max.overlaps = 30, segment.size = 0.3,
    segment.color = "gray60", fill = "white", alpha = 0.9,
    label.size = 0.2
  ) +
  scale_color_manual(values = c("FGF4 (serotonergic)" = col_fgf4,
                                "FGF8 (dopaminergic)" = col_fgf8)) +
  scale_shape_manual(values = c("FGF4 (serotonergic)" = 16,
                                "FGF8 (dopaminergic)" = 17)) +
  scale_y_continuous(breaks = 0:4) +
  labs(
    title = "Positive selection vs. differential expression",
    subtitle = "233 DEGs analyzed with aBSREL (HyPhy, primate tree)",
    x = "Absolute log2 fold change (Week4, FGF4 vs FGF8)",
    y = "Branches with positive selection (p ≤ 0.05)",
    color = "Condition", shape = "Condition"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40"))

ggsave("~/tesi_bioinfo/positive_selection_orthodb/v3_fig1_scatter.png",
       p1, width = 9, height = 6.5, dpi = 300)
cat("Fig 1 salvata\n")

# ---- FIGURA 2: Barplot ----
sig_df <- df %>% filter(n_branches_significant > 0)

p2 <- ggplot(sig_df, aes(x = reorder(gene, n_branches_significant),
                          y = n_branches_significant, fill = condition)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("FGF4 (serotonergic)" = col_fgf4,
                               "FGF8 (dopaminergic)" = col_fgf8)) +
  labs(
    title = "Genes with positive selection in primates",
    subtitle = "67 significant genes (aBSREL, p ≤ 0.05)",
    x = NULL,
    y = "Number of branches with positive selection",
    fill = "Condition"
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 7, face = "italic"))

ggsave("~/tesi_bioinfo/positive_selection_orthodb/v3_fig2_barplot.png",
       p2, width = 7, height = 12, dpi = 300)
cat("Fig 2 salvata\n")

# ---- FIGURA 3: Heatmap ----
top_genes <- df %>%
  filter(n_branches_significant > 0) %>%
  arrange(desc(n_branches_significant), gene) %>%
  slice_head(n = 30) %>%
  pull(gene)

all_species <- c("Homo_sapiens","Pan_troglodytes","Gorilla_gorilla_gorilla",
                 "Pongo_abelii","Nomascus_leucogenys","Macaca_mulatta",
                 "Macaca_fascicularis","Rhinopithecus_roxellana",
                 "Callithrix_jacchus","Cebus_imitator",
                 "Microcebus_murinus","Otolemur_garnettii",
                 "Node2","Node3","Node4","Node5","Node6",
                 "Node12","Node17","Node20")

heat_df <- df %>%
  filter(gene %in% top_genes) %>%
  rowwise() %>%
  mutate(branch_list = list(
    gsub("\\(p=[0-9.]+\\)", "", unlist(strsplit(significant_branches, ";")))
  )) %>%
  unnest(branch_list) %>%
  mutate(branch_list = trimws(branch_list)) %>%
  filter(branch_list != "none", branch_list != "") %>%
  mutate(present = 1) %>%
  select(gene, branch_list, present) %>%
  complete(gene = top_genes, branch_list = all_species, fill = list(present = 0))

heat_df <- heat_df %>%
  left_join(df %>% select(gene, condition), by = "gene")

heat_df$branch_list <- factor(heat_df$branch_list, levels = rev(all_species))
heat_df$gene <- factor(heat_df$gene,
                        levels = df %>%
                          filter(gene %in% top_genes) %>%
                          arrange(condition, desc(n_branches_significant)) %>%
                          pull(gene))

p3 <- ggplot(heat_df, aes(x = gene, y = branch_list, fill = factor(present))) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("0" = "gray92", "1" = "#D85A30"),
                    labels = c("not significant", "significant (p ≤ 0.05)")) +
  facet_grid(. ~ condition, scales = "free_x", space = "free") +
  labs(title = "Positive selection by branch and gene (top 30)",
       subtitle = "Each cell = a branch tested by aBSREL",
       x = NULL, y = NULL, fill = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7, face = "italic"),
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave("~/tesi_bioinfo/positive_selection_orthodb/v3_fig3_heatmap.png",
       p3, width = 14, height = 8, dpi = 300)
cat("Fig 3 salvata\n")

# ---- FIGURA 4: Pie chart riepilogo ----
pie_df <- df %>%
  group_by(condition) %>%
  summarise(
    sig = sum(n_branches_significant > 0),
    not_sig = sum(n_branches_significant == 0),
    total = n()
  ) %>%
  tidyr::pivot_longer(cols = c(sig, not_sig), names_to = "type", values_to = "count") %>%
  mutate(
    label = ifelse(type == "sig", "With positive\nselection", "No positive\nselection"),
    pct = round(100 * count / total, 1),
    label_full = paste0(label, "\n", count, " (", pct, "%)")
  )

p4 <- ggplot(pie_df, aes(x = "", y = count, fill = interaction(type, condition))) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
  coord_polar("y") +
  facet_wrap(~ condition) +
  geom_text(aes(label = label_full), 
            position = position_stack(vjust = 0.5),
            size = 3.2, color = "white", fontface = "bold") +
  scale_fill_manual(values = c(
    "sig.FGF4 (serotonergic)"     = col_fgf4,
    "not_sig.FGF4 (serotonergic)" = "#9FC8D8",
    "sig.FGF8 (dopaminergic)"     = col_fgf8,
    "not_sig.FGF8 (dopaminergic)" = "#F2C4A0"
  )) +
  labs(
    title = "Proportion of DEGs with positive selection",
    subtitle = "233 genes analyzed with aBSREL (primate tree, p ≤ 0.05)"
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(color = "gray40", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11)
  )

ggsave("~/tesi_bioinfo/positive_selection_orthodb/v3_fig4_pie.png",
       p4, width = 8, height = 5, dpi = 300)
cat("Fig 4 salvata\n")

cat("\nTutte le figure salvate!\n")
