#!/usr/bin/env python3

import os
import csv

# Percorsi
RESULTS = "/home/mirab/tesi_bioinfo/results"
OUT_DIR = "/home/mirab/tesi_bioinfo/expression_matrices"

os.makedirs(OUT_DIR, exist_ok=True)

SAMPLES = [
    'FEU8234A1', 'FEU8234A2',
    'FEU8234A4', 'FEU8234A5', 'FEU8234A6',
    'FEU8234A7', 'FEU8234A8', 'FEU8234A9',
    'FEU8234A10', 'FEU8234A11', 'FEU8234A12',
    'FEU8234A13', 'FEU8234A14', 'FEU8234A15',
    'FEU8234A16', 'FEU8234A17', 'FEU8234A18',
    'FEU8234A19', 'FEU8234A20', 'FEU8234A21',
    'FEU8234A22', 'FEU8234A23', 'FEU8234A24',
    'FEU8234A25', 'FEU8234A26', 'FEU8234A27',
    'FEU8234A28', 'FEU8234A29', 'FEU8234A30'
]

CATEGORIES = ['protein_coding', 'lncRNA', 'other_ncRNA']

for cat in CATEGORIES:
    print(f"Costruendo matrice per {cat}...")

    # Step 1: raccogli tutti i geni unici
    all_genes = set()
    for sample in SAMPLES:
        f_path = os.path.join(RESULTS, sample, f'{cat}_by_gene.tsv')
        if not os.path.exists(f_path):
            print(f"  ATTENZIONE: file mancante per {sample}")
            continue
        with open(f_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                all_genes.add(row['gene_name'])

    all_genes = sorted(all_genes)
    print(f"  Geni totali: {len(all_genes)}")

    # Step 2: costruisci dizionario sample -> {gene -> tpm}
    data = {}
    for sample in SAMPLES:
        f_path = os.path.join(RESULTS, sample, f'{cat}_by_gene.tsv')
        if not os.path.exists(f_path):
            continue
        data[sample] = {}
        with open(f_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                data[sample][row['gene_name']] = row['tpm']
        print(f"  Letti {len(data[sample])} geni da {sample}")

    # Step 3: scrivi la matrice
    out_file = os.path.join(OUT_DIR, f'{cat}_tpm_matrix.tsv')
    print(f"  Scrittura in {out_file}...")

    with open(out_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['gene_name'] + SAMPLES)
        for gene in all_genes:
            row = [gene]
            for sample in SAMPLES:
                row.append(data.get(sample, {}).get(gene, '0'))
            writer.writerow(row)

    # Verifica
    lines = sum(1 for _ in open(out_file))
    print(f"  Righe scritte: {lines}")
    print(f"  Matrice completata: {out_file}")

print("Tutto completato!")
