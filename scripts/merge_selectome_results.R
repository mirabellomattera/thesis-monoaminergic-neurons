# Unione dei 5 batch Selectome + filtro per taxon Primates
# Aggiunge info di direzione/significatività dal nostro DEG Week4

col_names <- c("tree_id", "selection_in_tree", "tree_name", "tree_description",
               "tree_taxon_root", "selection_in_gene", "gene_id", "gene_name",
               "gene_description", "species", "taxid", "match_score")

read_selectome <- function(path) {
  read.table(path, header = FALSE, skip = 1, sep = "\t",
             quote = "", stringsAsFactors = FALSE, fill = TRUE,
             col.names = col_names)
}

base_path <- "~/tesi_bioinfo/results_deseq2/Week4_FGF4_vs_FGF8/"
files <- file.path(base_path, paste0("selectome_batch_0", 0:4, ".tsv"))

all_sel <- do.call(rbind, lapply(files, read_selectome))

cat("Righe totali combinate:", nrow(all_sel), "\n")
cat("Geni unici (qualsiasi taxon):", length(unique(all_sel$gene_id)), "\n")
cat("Distribuzione per taxon:\n")
print(table(all_sel$tree_taxon_root))

# Carica il nostro DEG annotato per aggiungere direzione/log2FC/padj
deg <- read.table(paste0(base_path, "Week4_FGF4_vs_FGF8_protein_coding_DEG_annotated.tsv"),
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Risultati confermati specificamente sui Primati
primates_hits <- subset(all_sel, tree_taxon_root == "Primates" & selection_in_gene == "YES")
primates_final <- merge(primates_hits[, c("gene_id", "gene_name", "tree_id")],
                         deg[, c("gene_id_clean", "direction", "log2FoldChange", "padj")],
                         by.x = "gene_id", by.y = "gene_id_clean")
primates_final <- primates_final[order(primates_final$padj), ]

# Risultati solo a livello Euteleostomi (da controllare a mano se interessanti)
euteleostomi_hits <- subset(all_sel, tree_taxon_root == "Euteleostomi" & selection_in_gene == "YES")
euteleostomi_final <- merge(euteleostomi_hits[, c("gene_id", "gene_name", "tree_id")],
                             deg[, c("gene_id_clean", "direction", "log2FoldChange", "padj")],
                             by.x = "gene_id", by.y = "gene_id_clean")
euteleostomi_final <- euteleostomi_final[order(euteleostomi_final$padj), ]

write.table(primates_final, paste0(base_path, "POSITIVE_SELECTION_PRIMATES_CONFIRMED.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(euteleostomi_final, paste0(base_path, "positive_selection_euteleostomi_to_review.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== RISULTATO FINALE ===\n")
cat("Geni con selezione positiva CONFERMATA su ramo Primati:", nrow(primates_final), "\n")
cat("Geni con selezione positiva solo a livello Euteleostomi (da verificare):", nrow(euteleostomi_final), "\n")
cat("\nTop hit Primates:\n")
print(head(primates_final, 15))
