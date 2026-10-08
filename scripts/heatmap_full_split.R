library(pheatmap)
library(RColorBrewer)

categories <- c("protein_coding", "lncRNA", "other_ncRNA")

undiff_samples <- c("FEU8234A1","FEU8234A2","FEU8234A4","FEU8234A5","FEU8234A6",
                     "FEU8234A7","FEU8234A8","FEU8234A9")

fgf4_only <- c("FEU8234A10","FEU8234A11","FEU8234A12","FEU8234A16","FEU8234A17","FEU8234A18",
               "FEU8234A22","FEU8234A23","FEU8234A24","FEU8234A28","FEU8234A29","FEU8234A30")

fgf8_only <- c("FEU8234A13","FEU8234A14","FEU8234A15","FEU8234A19","FEU8234A20","FEU8234A21",
               "FEU8234A25","FEU8234A26","FEU8234A27")

all_samples <- c(undiff_samples, fgf4_only, fgf8_only)
fgf4_samples <- c(undiff_samples, fgf4_only)
fgf8_samples <- c(undiff_samples, fgf8_only)

fgf4_timepoints <- c("Week1","Week1","Week2_01","Week2_01","Week2_01",
                      "Week2_03","Week2_03","Week2_03",
                      "Week3_07_FGF4","Week3_07_FGF4","Week3_07_FGF4",
                      "Week3_10_FGF4","Week3_10_FGF4","Week3_10_FGF4",
                      "Week4_14_FGF4","Week4_14_FGF4","Week4_14_FGF4",
                      "Week4_17_FGF4","Week4_17_FGF4","Week4_17_FGF4")

fgf8_timepoints <- c("Week1","Week1","Week2_01","Week2_01","Week2_01",
                      "Week2_03","Week2_03","Week2_03",
                      "Week3_10_FGF8","Week3_10_FGF8","Week3_10_FGF8",
                      "Week4_14_FGF8","Week4_14_FGF8","Week4_14_FGF8",
                      "Week4_17_FGF8","Week4_17_FGF8","Week4_17_FGF8")

ann_colors <- list(
  Condition = c(Undifferentiated="grey80", FGF4="steelblue", FGF8="tomato"),
  Timepoint = c(
    Week1="#ffffe0", Week2_01="#ffd700", Week2_03="#ffa500",
    Week3_07_FGF4="#c6dbef", Week3_10_FGF4="#6baed6",
    Week4_14_FGF4="#2171b5", Week4_17_FGF4="#08306b",
    Week3_10_FGF8="#fdbe85", Week4_14_FGF8="#fd8d3c", Week4_17_FGF8="#a63603"
  )
)

for (cat in categories) {
  mat <- read.table(
    sprintf("/home/mirab/tesi_bioinfo/expression_matrices/%s_tpm_matrix.tsv", cat),
    header=TRUE, sep="\t", row.names=1, check.names=FALSE
  )

  mat_all <- mat[, all_samples]
  gene_means_all <- rowMeans(mat_all)
  keep_genes <- rownames(mat_all)[gene_means_all >= 1]
  mat_filtered <- mat[keep_genes, ]

  cat(sprintf("%s: master gene set = %d genes\n", cat, length(keep_genes)))

  for (cond in c("FGF4", "FGF8")) {
    samples <- if (cond == "FGF4") fgf4_samples else fgf8_samples
    timepoints <- if (cond == "FGF4") fgf4_timepoints else fgf8_timepoints
    condition_labels <- ifelse(grepl("^Week1$|^Week2", timepoints), "Undifferentiated", cond)

    mat_sub <- mat_filtered[, samples]

    row_sd <- apply(mat_sub, 1, sd)
    n_dropped <- sum(row_sd == 0)
    mat_sub <- mat_sub[row_sd > 0, ]

    mat_zscore <- t(scale(t(mat_sub)))

    sample_info <- data.frame(
      Timepoint = timepoints,
      Condition = condition_labels,
      row.names = samples
    )

    cat(sprintf("Generating %s - %s (%d genes, %d dropped as zero-variance)...\n",
                cat, cond, nrow(mat_zscore), n_dropped))

    png(sprintf("/home/mirab/tesi_bioinfo/figures/heatmaps/heatmap_full_%s_%s.png", cond, cat),
        width=1400, height=1800, res=150)
    pheatmap(
      mat_zscore,
      annotation_col = sample_info,
      annotation_colors = ann_colors,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      show_rownames = FALSE,
      show_colnames = TRUE,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      main = sprintf("%s Genes - %s - Z-score TPM (n = %d)", cat, cond, nrow(mat_zscore)),
      fontsize_col = 8
    )
    dev.off()
  }
}
cat("Done!\n")
