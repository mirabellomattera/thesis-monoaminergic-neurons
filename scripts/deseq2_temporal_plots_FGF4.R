# ============================================================
# Plots temporal analysis FGF4
# Heatmap top 100 + Line plot split by trend (up/down)
# ============================================================

.libPaths(path.expand("~/R/library"))
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(gridExtra)

CATEGORIES <- c("protein_coding", "lncRNA", "other_ncRNA")
RES_DIR <- path.expand("~/tesi_bioinfo/results_deseq2/temporal_FGF4")
MAT_DIR <- path.expand("~/tesi_bioinfo/expression_matrices")
FIG_DIR <- path.expand("~/tesi_bioinfo/figures/DESeq2/temporal_FGF4")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

# Annotation table
annot <- read.table(path.expand("~/tesi_bioinfo/reference/annotation_table.tsv"),
                    header=TRUE, sep="\t")
annot <- unique(annot[, c("gene_id", "gene_name")])

# FGF4 samples in chronological order (without A1)
samples_FGF4 <- c("FEU8234A2",
                   "FEU8234A4","FEU8234A5","FEU8234A6",
                   "FEU8234A7","FEU8234A8","FEU8234A9",
                   "FEU8234A10","FEU8234A11","FEU8234A12",
                   "FEU8234A16","FEU8234A17","FEU8234A18",
                   "FEU8234A22","FEU8234A23","FEU8234A24",
                   "FEU8234A28","FEU8234A29","FEU8234A30")

timepoint_labels <- c("Week1",
                       "Week2","Week2","Week2",
                       "Week2","Week2","Week2",
                       "Week3","Week3","Week3",
                       "Week3","Week3","Week3",
                       "Week4","Week4","Week4",
                       "Week4","Week4","Week4")

sample_info <- data.frame(
  Timepoint = timepoint_labels,
  row.names = samples_FGF4
)

ann_colors <- list(
  Timepoint = c(Week1="#ffffe0", Week2="#ffa500", Week3="#762a83", Week4="#00441b")
)

for (category in CATEGORIES) {
  cat("================================================\n")
  cat("Temporal plots FGF4:", category, "\n")
  cat("================================================\n")

  # Load temporal results
  res_path <- file.path(RES_DIR, paste0(category, "_temporal_FGF4.tsv"))
  res <- read.table(res_path, header=TRUE, sep="\t")

  # Top 100 genes by padj
  top100 <- head(res[!is.na(res$padj) & res$padj < 0.05, ], 100)

  if (nrow(top100) == 0) {
    cat("No significant genes - skip\n")
    next
  }

  # Bridge gene_id -> gene_name
  top100 <- merge(top100, annot, by="gene_id", all.x=TRUE)
  top100$label <- ifelse(is.na(top100$gene_name) | top100$gene_name == "",
                         top100$gene_id, top100$gene_name)

  # Load TPM matrix
  tpm_path <- file.path(MAT_DIR, paste0(category, "_tpm_matrix.tsv"))
  tpm <- read.table(tpm_path, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

  tpm_sel <- tpm[rownames(tpm) %in% top100$label, samples_FGF4]

  cat("  Genes found in TPM matrix:", nrow(tpm_sel), "\n")

  if (nrow(tpm_sel) < 2) {
    cat("  Too few genes - skip\n\n")
    next
  }

  # ------------------------------------------------------------
  # 1. HEATMAP top 100
  # ------------------------------------------------------------
  mat_z <- t(scale(t(tpm_sel)))
  mat_z <- mat_z[complete.cases(mat_z), ]

  heat_path <- file.path(FIG_DIR, paste0("heatmap_temporal_FGF4_", category, "_top100.png"))
  png(heat_path, width=1600, height=2400, res=150)
  pheatmap(mat_z,
           annotation_col    = sample_info,
           annotation_colors = ann_colors,
           color             = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
           show_rownames     = TRUE,
           show_colnames     = TRUE,
           cluster_cols      = FALSE,
           cluster_rows      = TRUE,
           main              = paste("Top 100 genes temporal trend - FGF4 -", category),
           fontsize_row      = 6,
           fontsize_col      = 8)
  dev.off()
  cat("  Saved:", heat_path, "\n")

  # ------------------------------------------------------------
  # 2. LINE PLOT — split by trend up/down
  # ------------------------------------------------------------
  tpm_log <- log2(tpm_sel + 1)

  tp_order <- c("Week1","Week2","Week3","Week4")

  mean_by_tp <- data.frame(
    gene  = rownames(tpm_log),
    Week1 = rowMeans(tpm_log[, timepoint_labels == "Week1", drop=FALSE]),
    Week2 = rowMeans(tpm_log[, timepoint_labels == "Week2", drop=FALSE]),
    Week3 = rowMeans(tpm_log[, timepoint_labels == "Week3", drop=FALSE]),
    Week4 = rowMeans(tpm_log[, timepoint_labels == "Week4", drop=FALSE])
  )

  # Classify genes as up or down based on Week4 vs Week1
  mean_by_tp$trend <- ifelse(mean_by_tp$Week4 >= mean_by_tp$Week1, "Up", "Down")

  up_genes   <- mean_by_tp$gene[mean_by_tp$trend == "Up"]
  down_genes <- mean_by_tp$gene[mean_by_tp$trend == "Down"]

  cat("  Genes up in Week4:", length(up_genes), "\n")
  cat("  Genes down in Week4:", length(down_genes), "\n")

  mean_long <- pivot_longer(mean_by_tp, cols=c("Week1","Week2","Week3","Week4"),
                             names_to="Timepoint", values_to="log2TPM")
  mean_long$Timepoint <- factor(mean_long$Timepoint, levels=tp_order)

  # Up genes plot
  p_up <- ggplot(mean_long[mean_long$gene %in% up_genes, ],
                 aes(x=Timepoint, y=log2TPM, group=gene, color=gene)) +
    geom_line(linewidth=0.6) +
    geom_point(size=1.5) +
    labs(title=paste("Upregulated genes - FGF4 -", category),
         x="Timepoint", y="log2(TPM + 1)", color="Gene") +
    theme_bw() +
    theme(legend.position="right",
          legend.text=element_text(size=5),
          legend.key.size=unit(0.3, "cm"))

  # Down genes plot
  p_down <- ggplot(mean_long[mean_long$gene %in% down_genes, ],
                   aes(x=Timepoint, y=log2TPM, group=gene, color=gene)) +
    geom_line(linewidth=0.6) +
    geom_point(size=1.5) +
    labs(title=paste("Downregulated genes - FGF4 -", category),
         x="Timepoint", y="log2(TPM + 1)", color="Gene") +
    theme_bw() +
    theme(legend.position="right",
          legend.text=element_text(size=5),
          legend.key.size=unit(0.3, "cm"))

  # Save combined plot
  line_path <- file.path(FIG_DIR, paste0("lineplot_temporal_FGF4_", category, "_updown.png"))
  png(line_path, width=2400, height=800, res=120)
  grid.arrange(p_up, p_down, ncol=2)
  dev.off()
  cat("  Saved:", line_path, "\n\n")
}

cat("=== DONE ===\n")
