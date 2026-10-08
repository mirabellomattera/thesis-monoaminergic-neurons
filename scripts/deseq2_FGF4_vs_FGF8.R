# ============================================================
# DESeq2 — Analisi 1: FGF4 vs FGF8
# Solo campioni differenziati (Week3 + Week4)
# ============================================================

.libPaths(path.expand("~/R/library"))
library(DESeq2)

# ------------------------------------------------------------
# PARAMETRI
# ------------------------------------------------------------
OUTLIERS <- c("FEU8234A1")
CATEGORIES <- c("protein_coding", "lncRNA", "other_ncRNA")
MAT_DIR  <- path.expand("~/tesi_bioinfo/expression_matrices")
OUT_DIR  <- path.expand("~/tesi_bioinfo/results_deseq2/FGF4_vs_FGF8")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# ------------------------------------------------------------
# METADATI
# ------------------------------------------------------------
sample_info <- data.frame(
  sample    = c("FEU8234A10","FEU8234A11","FEU8234A12",
                "FEU8234A13","FEU8234A14","FEU8234A15",
                "FEU8234A16","FEU8234A17","FEU8234A18",
                "FEU8234A19","FEU8234A20","FEU8234A21",
                "FEU8234A22","FEU8234A23","FEU8234A24",
                "FEU8234A25","FEU8234A26","FEU8234A27",
                "FEU8234A28","FEU8234A29","FEU8234A30"),
  condition = c("FGF4","FGF4","FGF4",
                "FGF8","FGF8","FGF8",
                "FGF4","FGF4","FGF4",
                "FGF8","FGF8","FGF8",
                "FGF4","FGF4","FGF4",
                "FGF8","FGF8","FGF8",
                "FGF4","FGF4","FGF4"),
  timepoint = c("Week3","Week3","Week3",
                "Week3","Week3","Week3",
                "Week3","Week3","Week3",
                "Week4","Week4","Week4",
                "Week4","Week4","Week4",
                "Week4","Week4","Week4",
                "Week4","Week4","Week4"),
  stringsAsFactors = FALSE
)

sample_info <- sample_info[!sample_info$sample %in% OUTLIERS, ]
sample_info$condition <- factor(sample_info$condition, levels=c("FGF4","FGF8"))
rownames(sample_info) <- sample_info$sample

cat("Campioni inclusi:", nrow(sample_info), "\n")
cat("FGF4:", sum(sample_info$condition=="FGF4"), "\n")
cat("FGF8:", sum(sample_info$condition=="FGF8"), "\n\n")

# ------------------------------------------------------------
# LOOP SU CATEGORIE
# ------------------------------------------------------------
for (category in CATEGORIES) {
  cat("================================================\n")
  cat("Analisi:", category, "\n")
  cat("================================================\n")

  mat_path <- file.path(MAT_DIR, paste0(category, "_counts_matrix.tsv"))
  counts <- read.table(mat_path, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

  samples_use <- sample_info$sample
  counts <- counts[, colnames(counts) %in% samples_use]
  counts <- counts[, samples_use[samples_use %in% colnames(counts)]]

  cat("Dimensioni matrice:", nrow(counts), "geni x", ncol(counts), "campioni\n")

  counts <- counts[rowSums(counts) >= 10, ]
  cat("Geni dopo filtro (rowSums >= 10):", nrow(counts), "\n")

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = sample_info[colnames(counts), ],
    design    = ~ condition
  )

  cat("Eseguo DESeq2...\n")
  dds <- DESeq(dds)

  res <- results(dds,
                 contrast = c("condition", "FGF8", "FGF4"),
                 alpha    = 0.05)

  cat("Summary risultati:\n")
  summary(res)

  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[order(res_df$padj, na.last=TRUE), ]

  deg_up   <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1,  na.rm=TRUE)
  deg_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm=TRUE)
  cat("DEG upregolati in FGF8:", deg_up, "\n")
  cat("DEG downregolati in FGF8:", deg_down, "\n\n")

  out_path <- file.path(OUT_DIR, paste0(category, "_DESeq2_FGF4vsFGF8.tsv"))
  write.table(res_df, out_path, sep="\t", quote=FALSE, row.names=FALSE)
  cat("Salvato:", out_path, "\n\n")
}

cat("=== DONE ===\n")
