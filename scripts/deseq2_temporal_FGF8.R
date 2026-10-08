# ============================================================
# DESeq2 — Analisi temporale FGF8
# Campioni: Week1 + Week2 + Week3_FGF8 + Week4_FGF8
# Modello: ~ timepoint (LRT)
# ============================================================

.libPaths(path.expand("~/R/library"))
library(DESeq2)

CATEGORIES <- c("protein_coding", "lncRNA", "other_ncRNA")
MAT_DIR  <- path.expand("~/tesi_bioinfo/expression_matrices")
OUT_DIR  <- path.expand("~/tesi_bioinfo/results_deseq2/temporal_FGF8")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# ------------------------------------------------------------
# METADATI — FGF8 + baseline
# ------------------------------------------------------------
sample_info <- data.frame(
  sample    = c("FEU8234A1","FEU8234A2",              # Week1
                "FEU8234A4","FEU8234A5","FEU8234A6",  # Week2 day1
                "FEU8234A7","FEU8234A8","FEU8234A9",  # Week2 day3
                "FEU8234A13","FEU8234A14","FEU8234A15", # Week3 FGF8
                "FEU8234A19","FEU8234A20","FEU8234A21", # Week4 FGF8
                "FEU8234A25","FEU8234A26","FEU8234A27"), # Week4 FGF8
  timepoint = c(1, 1,
                2, 2, 2,
                2, 2, 2,
                3, 3, 3,
                4, 4, 4,
                4, 4, 4),
  stringsAsFactors = FALSE
)

# Escludi A1 (outlier estremo)
sample_info <- sample_info[sample_info$sample != "FEU8234A1", ]
rownames(sample_info) <- sample_info$sample

cat("Campioni inclusi:", nrow(sample_info), "\n")
cat("Timepoint:", table(sample_info$timepoint), "\n\n")

# ------------------------------------------------------------
# LOOP SU CATEGORIE
# ------------------------------------------------------------
for (category in CATEGORIES) {
  cat("================================================\n")
  cat("Analisi temporale FGF8:", category, "\n")
  cat("================================================\n")

  mat_path <- file.path(MAT_DIR, paste0(category, "_counts_matrix.tsv"))
  counts <- read.table(mat_path, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

  samples_use <- sample_info$sample
  counts <- counts[, samples_use]

  cat("Dimensioni matrice:", nrow(counts), "geni x", ncol(counts), "campioni\n")

  counts <- counts[rowSums(counts) >= 10, ]
  cat("Geni dopo filtro:", nrow(counts), "\n")

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = sample_info[colnames(counts), ],
    design    = ~ timepoint
  )

  cat("Eseguo DESeq2 LRT...\n")
  dds <- DESeq(dds, test="LRT", reduced= ~ 1)

  res <- results(dds, alpha=0.05)

  cat("Summary risultati:\n")
  summary(res)

  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[order(res_df$padj, na.last=TRUE), ]

  deg_n <- sum(res_df$padj < 0.05, na.rm=TRUE)
  cat("Geni con trend temporale significativo (padj<0.05):", deg_n, "\n\n")

  out_path <- file.path(OUT_DIR, paste0(category, "_temporal_FGF8.tsv"))
  write.table(res_df, out_path, sep="\t", quote=FALSE, row.names=FALSE)
  cat("Salvato:", out_path, "\n\n")
}

cat("=== DONE ===\n")
