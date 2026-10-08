#!/usr/bin/env python3

import re
import csv

# Percorsi
GTF = "/home/mirab/tesi_bioinfo/reference/gencode.v46.annotation.gtf"
OUT = "/home/mirab/tesi_bioinfo/reference/annotation_table.tsv"

print("Lettura del file GTF...")

annotation = {}

with open(GTF, 'r') as f:
    for line in f:
        # Salta le righe di commento
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        
        # Prendiamo solo le righe di tipo 'transcript'
        if fields[2] != 'transcript':
            continue
        
        info = fields[8]
        
        # Estrai transcript_id
        t_id = re.search(r'transcript_id "([^"]+)"', info)
        # Estrai gene_id
        g_id = re.search(r'gene_id "([^"]+)"', info)
        # Estrai gene_name
        g_name = re.search(r'gene_name "([^"]+)"', info)
        # Estrai gene_type
        g_type = re.search(r'gene_type "([^"]+)"', info)
        # Estrai transcript_type
        t_type = re.search(r'transcript_type "([^"]+)"', info)

        if t_id and g_id and g_name and g_type:
            annotation[t_id.group(1)] = {
                'gene_id': g_id.group(1),
                'gene_name': g_name.group(1),
                'gene_type': g_type.group(1),
                'transcript_type': t_type.group(1) if t_type else 'NA'
            }

print(f"Trovati {len(annotation)} trascritti annotati")

# Scrivi la tabella di output
print("Scrittura tabella di annotazione...")

with open(OUT, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['transcript_id', 'gene_id', 'gene_name', 'gene_type', 'transcript_type'])
    for t_id, info in annotation.items():
        writer.writerow([t_id, info['gene_id'], info['gene_name'], info['gene_type'], info['transcript_type']])

print(f"Tabella salvata in {OUT}")
print("Done!")
