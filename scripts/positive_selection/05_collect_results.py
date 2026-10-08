#!/usr/bin/env python3
import os, json

RES_DIR = os.path.expanduser("~/tesi_bioinfo/positive_selection_orthodb/results")
OUT = os.path.expanduser("~/tesi_bioinfo/positive_selection_orthodb/absrel_summary.tsv")

results = []
for fname in sorted(os.listdir(RES_DIR)):
    if not fname.endswith('_absrel.json'):
        continue
    gene = fname.replace('_absrel.json','')
    fpath = os.path.join(RES_DIR, fname)
    try:
        with open(fpath) as f:
            data = json.load(f)
        ba = data['branch attributes']['0']
        sig_branches = []
        for branch, attrs in ba.items():
            pc = attrs.get('Corrected P-value')
            if pc is not None and pc <= 0.05:
                sig_branches.append(f"{branch}(p={pc:.4f})")
        n_tested = data['test results']['tested']
        n_sig = data['test results']['positive test results']
        results.append({
            'gene': gene,
            'n_tested': n_tested,
            'n_sig': n_sig,
            'sig_branches': ';'.join(sig_branches) if sig_branches else 'none'
        })
    except Exception as e:
        print(f"Errore {gene}: {e}")

# Scrivi tabella
with open(OUT, 'w') as f:
    f.write('gene\tn_branches_tested\tn_branches_significant\tsignificant_branches\n')
    for r in results:
        f.write(f"{r['gene']}\t{r['n_tested']}\t{r['n_sig']}\t{r['sig_branches']}\n")

sig = [r for r in results if r['n_sig'] > 0]
print(f"Geni analizzati: {len(results)}")
print(f"Geni con almeno 1 ramo significativo (p<=0.05): {len(sig)}")
print(f"\nTop geni sotto selezione positiva:")
for r in sorted(sig, key=lambda x: -x['n_sig'])[:20]:
    print(f"  {r['gene']}: {r['n_sig']} rami - {r['sig_branches']}")
