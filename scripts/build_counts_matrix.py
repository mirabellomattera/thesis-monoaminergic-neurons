import os
import pandas as pd

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

BASE_DIR = os.path.expanduser('~/tesi_bioinfo/results')
OUT_DIR  = os.path.expanduser('~/tesi_bioinfo/expression_matrices')

for category in CATEGORIES:
    print(f"\n=== Costruzione matrice est_counts: {category} ===")

    matrix = {}

    for sample in SAMPLES:
        filepath = os.path.join(BASE_DIR, sample, f'{category}_by_gene.tsv')

        if not os.path.exists(filepath):
            print(f"  ATTENZIONE: file mancante per {sample} — {filepath}")
            continue

        df = pd.read_csv(filepath, sep='\t')

        # Arrotonda est_counts all'intero (DESeq2 richiede interi)
        matrix[sample] = df.set_index('gene_id')['est_counts'].round().astype(int)

        print(f"  {sample}: {len(matrix[sample])} geni letti")

    # Costruisci dataframe finale
    mat_df = pd.DataFrame(matrix)
    mat_df.index.name = 'gene_id'

    # Ordine cronologico (segue lista SAMPLES)
    mat_df = mat_df[SAMPLES]

    # Salva
    out_path = os.path.join(OUT_DIR, f'{category}_counts_matrix.tsv')
    mat_df.to_csv(out_path, sep='\t')

    # Verifica
    print(f"  -> Salvato: {out_path}")
    print(f"  -> Dimensioni: {mat_df.shape[0]} geni x {mat_df.shape[1]} campioni")

    if (mat_df < 0).any().any():
        print(f"  ERRORE: valori negativi trovati!")
    else:
        print(f"  Check OK: nessun valore negativo")

    # Verifica valori nulli
    nulls = mat_df.isnull().sum().sum()
    if nulls > 0:
        print(f"  ATTENZIONE: {nulls} valori NaN trovati!")
    else:
        print(f"  Check OK: nessun valore NaN")

print("\n=== DONE — tutte le matrici generate ===")
