# ============================================================
# DESeq2 VST PCA — FGF4 and FGF8 temporal analysis
# v2: A1 included (exploratory stage, per Methods 3.1.4)
# ============================================================
.libPaths(path.expand("~/R/library"))
library(DESeq2)
library(ggplot2)

MAT_DIR <- path.expand("~/tesi_bioinfo/expression_matrices")
FIG_DIR <- path.expand("~/tesi_bioinfo/figures/DESeq2/PCA")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

# ------------------------------------------------------------
# METADATA — all 29 samples (A1 included for exploratory PCA)
# ------------------------------------------------------------
all_samples <- data.frame(
  sample    = c("FEU8234A1","FEU8234A2",
                "FEU8234A4","FEU8234A5","FEU8234A6",
                "FEU8234A7","FEU8234A8","FEU8234A9",
                "FEU8234A10","FEU8234A11","FEU8234A12",
                "FEU8234A13","FEU8234A14","FEU8234A15",
                "FEU8234A16","FEU8234A17","FEU8234A18",
                "FEU8234A19","FEU8234A20","FEU8234A21",
                "FEU8234A22","FEU8234A23","FEU8234A24",
                "FEU8234A25","FEU8234A26","FEU8234A27",
                "FEU8234A28","FEU8234A29","FEU8234A30"),
  condition = c("Undiff","Undiff",
                "Undiff","Undiff","Undiff",
                "Undiff","Undiff","Undiff",
                "FGF4","FGF4","FGF4",
                "FGF8","FGF8","FGF8",
                "FGF4","FGF4","FGF4",
                "FGF8","FGF8","FGF8",
                "FGF4","FGF4","FGF4",
                "FGF8","FGF8","FGF8",
                "FGF4","FGF4","FGF4"),
  timepoint = c("Week1","Week1",
                "Week2","Week2","Week2",
                "Week2","Week2","Week2",
                "Week3","Week3","Week3",
                "Week3","Week3","Week3",
                "Week3","Week3","Week3",
                "Week4","Week4","Week4",
                "Week4","Week4","Week4",
                "Week4","Week4","Week4",
                "Week4","Week4","Week4"),
  stringsAsFactors = FALSE
)
rownames(all_samples) <- all_samples$sample

condition_colors <- c(Undiff="grey60", FGF4="steelblue", FGF8="tomato")
timepoint_colors <- c(Week1="#ffffe0", Week2="#ffa500", Week3="#9970ab", Week4="#006837")
timepoint_shapes <- c(Week1=15, Week2=16, Week3=17, Week4=18)

CATEGORIES <- c("protein_coding", "lncRNA", "other_ncRNA")

for (category in CATEGORIES) {
  cat("================================================\n")
  cat("DESeq2 VST PCA (v2, A1 included):", category, "\n")
  cat("================================================\n")

  mat_path <- file.path(MAT_DIR, paste0(category, "_counts_matrix.tsv"))
  counts <- read.table(mat_path, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

  counts <- counts[, all_samples$sample]
  counts <- counts[rowSums(counts) >= 10, ]
  cat("  Genes after filter:", nrow(counts), "\n")

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = all_samples[colnames(counts), ],
    design    = ~ condition
  )

  vst <- vst(dds, blind=TRUE)
  cat("  VST done\n")

  pca <- prcomp(t(assay(vst)), scale=FALSE, center=TRUE)
  var_exp <- round(summary(pca)$importance[2, 1:3] * 100, 1)
  cat("  Variance explained PC1/PC2/PC3:", var_exp, "\n")

  pca_df <- data.frame(
    Sample    = rownames(pca$x),
    PC1       = pca$x[, 1],
    PC2       = pca$x[, 2],
    PC3       = pca$x[, 3],
    Condition = all_samples[rownames(pca$x), "condition"],
    Timepoint = all_samples[rownames(pca$x), "timepoint"]
  )
  pca_df$Timepoint <- factor(pca_df$Timepoint, levels=c("Week1","Week2","Week3","Week4"))

  p1 <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Condition, shape=Timepoint, label=Sample)) +
    geom_point(size=4) +
    geom_text(size=2.5, vjust=-0.8, hjust=0.5) +
    scale_color_manual(values=condition_colors) +
    scale_shape_manual(values=timepoint_shapes) +
    labs(title=paste("PCA VST (v2, A1 included) -", category, "- PC1 vs PC2"),
         x=paste0("PC1 (", var_exp[1], "%)"),
         y=paste0("PC2 (", var_exp[2], "%)")) +
    theme_bw() +
    theme(legend.position="right")
  p1_path <- file.path(FIG_DIR, paste0("PCA_VST_v2_", category, "_PC1_PC2_condition.png"))
  ggsave(p1_path, p1, width=10, height=7, dpi=150)
  cat("  Saved:", p1_path, "\n")

  p2 <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Timepoint, shape=Condition, label=Sample)) +
    geom_point(size=4) +
    geom_text(size=2.5, vjust=-0.8, hjust=0.5) +
    scale_color_manual(values=timepoint_colors) +
    labs(title=paste("PCA VST (v2, A1 included) -", category, "- PC1 vs PC2 by Timepoint"),
         x=paste0("PC1 (", var_exp[1], "%)"),
         y=paste0("PC2 (", var_exp[2], "%)")) +
    theme_bw() +
    theme(legend.position="right")
  p2_path <- file.path(FIG_DIR, paste0("PCA_VST_v2_", category, "_PC1_PC2_timepoint.png"))
  ggsave(p2_path, p2, width=10, height=7, dpi=150)
  cat("  Saved:", p2_path, "\n")

  p3 <- ggplot(pca_df, aes(x=PC1, y=PC3, color=Condition, shape=Timepoint, label=Sample)) +
    geom_point(size=4) +
    geom_text(size=2.5, vjust=-0.8, hjust=0.5) +
    scale_color_manual(values=condition_colors) +
    scale_shape_manual(values=timepoint_shapes) +
    labs(title=paste("PCA VST (v2, A1 included) -", category, "- PC1 vs PC3"),
         x=paste0("PC1 (", var_exp[1], "%)"),
         y=paste0("PC3 (", var_exp[3], "%)")) +
    theme_bw() +
    theme(legend.position="right")
  p3_path <- file.path(FIG_DIR, paste0("PCA_VST_v2_", category, "_PC1_PC3_condition.png"))
  ggsave(p3_path, p3, width=10, height=7, dpi=150)
  cat("  Saved:", p3_path, "\n\n")
}

cat("=== DONE ===\n")
