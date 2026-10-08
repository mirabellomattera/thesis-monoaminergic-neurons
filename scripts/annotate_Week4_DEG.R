# Annotazione DEG Week4 FGF4 vs FGF8 con gene_name e Ensembl ID senza versione
# Output: tabella pronta per ricerca su Selectome (primati) e per discussione col prof

deg <- read.table("~/tesi_bioinfo/results_deseq2/Week4_FGF4_vs_FGF8/Week4_FGF4_vs_FGF8_protein_coding_DEG.tsv",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

ann <- read.table("~/tesi_bioinfo/reference/annotation_table.tsv",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "", fill = TRUE)

# Una riga per gene_id (gene_name e gene_id non cambiano tra trascritti dello stesso gene)
ann_gene <- unique(ann[, c("gene_id", "gene_name")])

# Merge
deg_ann <- merge(deg, ann_gene, by = "gene_id", all.x = TRUE)

# ID Ensembl senza versione (serve per Selectome/Ensembl Compara)
deg_ann$gene_id_clean <- sub("\\..*", "", deg_ann$gene_id)

# Direzione biologica
deg_ann$direction <- ifelse(deg_ann$log2FoldChange > 1, "up_FGF8_dopaminergic",
                      ifelse(deg_ann$log2FoldChange < -1, "up_FGF4_serotonergic", "NA"))

# Riordino colonne e ordino per padj
deg_ann <- deg_ann[order(deg_ann$padj),
                    c("gene_name", "gene_id_clean", "gene_id", "direction",
                      "log2FoldChange", "padj", "baseMean")]

write.table(deg_ann, "~/tesi_bioinfo/results_deseq2/Week4_FGF4_vs_FGF8/Week4_FGF4_vs_FGF8_protein_coding_DEG_annotated.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Geni totali:", nrow(deg_ann), "\n")
cat("Geni con gene_name trovato:", sum(!is.na(deg_ann$gene_name) & deg_ann$gene_name != ""), "\n")
cat("Geni senza match (NA):", sum(is.na(deg_ann$gene_name) | deg_ann$gene_name == ""), "\n")
cat("\nTop 10 geni per significatività:\n")
print(head(deg_ann, 10))
