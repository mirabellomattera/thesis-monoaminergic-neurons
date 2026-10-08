#!/usr/bin/env python3

import os
import csv
from collections import defaultdict

# Percorsi
RESULTS = "/home/mirab/tesi_bioinfo/results"

# Categorie da collassare
CATEGORIES = ['protein_coding', 'lncRNA', 'other_ncRNA']

print("Inizio collasso per gene...")

samples = sorted(os.listdir(RESULTS))

for sample in samples:
    sample_dir = os.path.join(RESULTS, sample)
    print(f"Processando {sample}...")

    for cat in CATEGORIES:
        input_file = os.path.join(sample_dir, f'{cat}.tsv')
        output_file = os.path.join(sample_dir, f'{cat}_by_gene.tsv')

        if not os.path.exists(input_file):
            continue

        # Dizionario per collassare per gene
        genes = defaultdict(lambda: {'est_counts': 0.0, 'tpm': 0.0, 'gene_id': '', 'gene_type': ''})

        with open(input_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_name = row['gene_name']
                genes[gene_name]['gene_id'] = row['gene_id']
                genes[gene_name]['gene_type'] = row['gene_type']
                genes[gene_name]['est_counts'] += float(row['est_counts'])
                genes[gene_name]['tpm'] += float(row['tpm'])

        # Scrivi il file collassato
        with open(output_file, 'w', newline='') as f:
            fieldnames = ['gene_name', 'gene_id', 'gene_type', 'est_counts', 'tpm']
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for gene_name, info in sorted(genes.items()):
                writer.writerow({
                    'gene_name': gene_name,
                    'gene_id': info['gene_id'],
                    'gene_type': info['gene_type'],
                    'est_counts': round(info['est_counts'], 4),
                    'tpm': round(info['tpm'], 4)
                })

        print(f"  {cat}: {len(genes)} geni unici")

print("Collasso completato!")
