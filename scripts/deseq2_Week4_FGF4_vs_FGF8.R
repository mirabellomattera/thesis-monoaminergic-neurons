# Week4-only DESeq2 analysis: FGF4 vs FGF8 (protein_coding)
# Design: ~condition, solo campioni Week4 (fase differenziata/stabilizzata)
# Output: lista DEG per analisi di selezione positiva (Selectome, primati)

.libPaths(path.expand("~/R/library"))
library(DESeq2)

counts <- read.table("~/tesi_bioinfo/expression_matrices/protein_coding_counts_matrix.tsv",
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

week4_fgf4 <- c("FEU8234A22", "FEU8234A23", "FEU8234A24",
                "FEU8234A28", "FEU8234A29", "FEU8234A30")
week4_fgf8 <- c("FEU8234A19", "FEU8234A20", "FEU8234A21",
                "FEU8234A25", "FEU8234A26", "FEU8234A27")

week4_samples <- c(week4_fgf4, week4_fgf8)
counts_week4 <- counts[, week4_samples]

condition <- factor(c(rep("FGF4", length(week4_fgf4)), rep("FGF8", length(week4_fgf8))),
                     levels = c("FGF4", "FGF8"))
coldata <- data.frame(condition = condition, row.names = week4_samples)

dds <- DESeqDataSetFromMatrix(countData = round(counts_week4),
                               colData = coldata,
                               design = ~ condition)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "FGF8", "FGF4"), alpha = 0.05)
res <- res[order(res$padj), ]

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
write.table(res_df, "~/tesi_bioinfo/results_deseq2/Week4_FGF4_vs_FGF8/Week4_FGF4_vs_FGF8_protein_coding_full.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
write.table(sig, "~/tesi_bioinfo/results_deseq2/Week4_FGF4_vs_FGF8/Week4_FGF4_vs_FGF8_protein_coding_DEG.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Geni testati:", nrow(res_df), "\n")
cat("DEG significativi (padj<0.05, |log2FC|>1):", nrow(sig), "\n")
cat("Up in FGF8 (dopaminergico):", sum(sig$log2FoldChange > 1), "\n")
cat("Up in FGF4 (serotoninergico):", sum(sig$log2FoldChange < -1), "\n")
