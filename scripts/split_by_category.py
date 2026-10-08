#!/usr/bin/env python3

import os
import csv

# Percorsi
ANNOTATION = "/home/mirab/tesi_bioinfo/reference/annotation_table.tsv"
RESULTS = "/home/mirab/tesi_bioinfo/results"

print("Caricamento tabella di annotazione...")

# Carica la tabella di annotazione in un dizionario
annotation = {}
with open(ANNOTATION, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        # Prendiamo solo l'ID base del trascritto (senza versione)
        t_id = row['transcript_id']
        annotation[t_id] = {
            'gene_id': row['gene_id'],
            'gene_name': row['gene_name'],
            'gene_type': row['gene_type'],
            'transcript_type': row['transcript_type']
        }

print(f"Annotazione caricata: {len(annotation)} trascritti")

# Categorie
def get_category(gene_type):
    if gene_type == 'protein_coding':
        return 'protein_coding'
    elif gene_type == 'lncRNA':
        return 'lncRNA'
    else:
        return 'other_ncRNA'

# Processa ogni campione
samples = sorted(os.listdir(RESULTS))

for sample in samples:
    sample_dir = os.path.join(RESULTS, sample)
    abundance_file = os.path.join(sample_dir, 'abundance.tsv')
    
    if not os.path.exists(abundance_file):
        continue
    
    print(f"Processando {sample}...")
    
    # Dizionario per raccogliere i dati per categoria
    categories = {
        'protein_coding': [],
        'lncRNA': [],
        'other_ncRNA': []
    }
    
    with open(abundance_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # L'ID nel file abundance ha formato lungo, prendiamo la parte prima del |
            target_id = row['target_id']
            # Estrai il transcript_id (prima parte prima del |)
            t_id = target_id.split('|')[0]
            
            if t_id in annotation:
                gene_type = annotation[t_id]['gene_type']
                category = get_category(gene_type)
                
                # Aggiungi info annotazione alla riga
                row['gene_name'] = annotation[t_id]['gene_name']
                row['gene_id'] = annotation[t_id]['gene_id']
                row['gene_type'] = gene_type
                categories[category].append(row)
            else:
                categories['other_ncRNA'].append(row)
    
    # Scrivi i file per categoria
    for cat, rows in categories.items():
        out_file = os.path.join(sample_dir, f'{cat}.tsv')
        if rows:
            with open(out_file, 'w', newline='') as f:
                fieldnames = ['target_id', 'gene_id', 'gene_name', 'gene_type', 'length', 'eff_length', 'est_counts', 'tpm']
                writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
                writer.writeheader()
                writer.writerows(rows)
            print(f"  {cat}: {len(rows)} trascritti")

print("Spacchettamento completato!")
