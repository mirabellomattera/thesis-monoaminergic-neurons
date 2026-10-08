# ============================================================
# Export TF temporal analysis tables - FGF4 and FGF8
# ============================================================

.libPaths(path.expand("~/R/library"))

# Load TF list
tf <- read.table(path.expand("~/tesi_bioinfo/reference/Homo_sapiens_TF.txt"),
                 header=TRUE, sep="\t")
tf_clean <- tf[tf$Symbol != "", c("Symbol", "Family")]

# Load annotation table
annot <- read.table(path.expand("~/tesi_bioinfo/reference/annotation_table.tsv"),
                    header=TRUE, sep="\t")
annot <- unique(annot[, c("gene_id", "gene_name")])

OUT_DIR <- path.expand("~/tesi_bioinfo/results_deseq2/TF_tables")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

for (condition in c("FGF4", "FGF8")) {
  cat("Processing", condition, "\n")

  res <- read.table(
    path.expand(paste0("~/tesi_bioinfo/results_deseq2/temporal_", condition,
                       "/protein_coding_temporal_", condition, ".tsv")),
    header=TRUE, sep="\t")

  # Bridge gene_id -> gene_name
  res <- merge(res, annot, by="gene_id", all.x=TRUE)
  res$gene_name <- ifelse(is.na(res$gene_name) | res$gene_name == "",
                          res$gene_id, res$gene_name)

  # Filter significant
  res_sig <- res[!is.na(res$padj) & res$padj < 0.05, ]

  # Keep only TF
  res_tf <- res_sig[res_sig$gene_name %in% tf_clean$Symbol, ]

  # Add TF family
  res_tf <- merge(res_tf, tf_clean, by.x="gene_name", by.y="Symbol", all.x=TRUE)

  # Order by padj
  res_tf <- res_tf[order(res_tf$padj), ]

  # Select and rename columns
  res_tf <- res_tf[, c("gene_name", "gene_id", "Family", "baseMean",
                        "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  colnames(res_tf)[1:3] <- c("TF_name", "Ensembl_ID", "TF_family")

  cat("  TF significativi:", nrow(res_tf), "\n")

  # Save
  out_path <- file.path(OUT_DIR, paste0("TF_temporal_", condition, "_significant.tsv"))
  write.table(res_tf, out_path, sep="\t", quote=FALSE, row.names=FALSE)
  cat("  Salvato:", out_path, "\n")
}

cat("=== DONE ===\n")
