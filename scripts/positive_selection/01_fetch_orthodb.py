#!/usr/bin/env python3
"""
Script 01: Scarica sequenze CDS da OrthoDB v12 per tutti i DEG Week4
Per ciascun gene, cerca l'OG a livello Primates (taxid 9443) e scarica le CDS.
Output: una cartella per gene in positive_selection_orthodb/sequences/
"""

import requests
import os
import time
import json

# Percorsi
BASE_DIR = os.path.expanduser("~/tesi_bioinfo/positive_selection_orthodb")
SEQ_DIR = os.path.join(BASE_DIR, "sequences")
LOG_DIR = os.path.join(BASE_DIR, "logs")
DEG_FILE = os.path.expanduser("~/tesi_bioinfo/results_deseq2/Week4_FGF4_vs_FGF8/Week4_FGF4_vs_FGF8_protein_coding_DEG_annotated.tsv")

API_BASE = "https://data.orthodb.org/v12"
PRIMATE_TAXID = 9443
MIN_SPECIES = 5  # numero minimo di specie per considerare il gene

def search_og(gene_name, ensembl_id):
    """Cerca l'OG di un gene a livello Primates."""
    # Prima prova con Ensembl ID
    url = f"{API_BASE}/search?query={ensembl_id}&level={PRIMATE_TAXID}&limit=3"
    r = requests.get(url, timeout=30)
    if r.status_code == 200:
        data = r.json()
        if data.get("count", "0") != "0" and data.get("data"):
            return data["data"][0], data["bigdata"][0]
    
    # Poi prova con gene name
    url = f"{API_BASE}/search?query={gene_name}&level={PRIMATE_TAXID}&limit=3"
    r = requests.get(url, timeout=30)
    if r.status_code == 200:
        data = r.json()
        if data.get("count", "0") != "0" and data.get("data"):
            return data["data"][0], data["bigdata"][0]
    
    return None, None

def download_cds(og_id, output_file):
    """Scarica le sequenze CDS per un OG."""
    url = f"{API_BASE}/fasta?id={og_id}&seqtype=cds"
    r = requests.get(url, timeout=60)
    if r.status_code == 200 and r.text.startswith(">"):
        with open(output_file, "w") as f:
            f.write(r.text)
        return True
    return False

def main():
    # Leggi lista DEG
    genes = []
    with open(DEG_FILE) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                gene_name = parts[0]
                gene_id = parts[1]  # gene_id_clean (ENSG...)
                genes.append((gene_name, gene_id))
    
    print(f"Geni da processare: {len(genes)}")
    
    results = []
    not_found = []
    too_few_species = []
    
    for i, (gene_name, gene_id) in enumerate(genes):
        print(f"[{i+1}/{len(genes)}] {gene_name} ({gene_id})", end=" ... ")
        
        gene_dir = os.path.join(SEQ_DIR, gene_name)
        os.makedirs(gene_dir, exist_ok=True)
        output_fasta = os.path.join(gene_dir, f"{gene_name}_primates_cds.fasta")
        
        # Salta se già scaricato
        if os.path.exists(output_fasta) and os.path.getsize(output_fasta) > 0:
            print("già presente, salto")
            continue
        
        og_id, og_info = search_og(gene_name, gene_id)
        
        if og_id is None:
            print("NON TROVATO")
            not_found.append(gene_name)
            time.sleep(0.5)
            continue
        
        n_species = int(og_info.get("present_in", 0))
        if n_species < MIN_SPECIES:
            print(f"poche specie ({n_species}), salto")
            too_few_species.append(gene_name)
            time.sleep(0.5)
            continue
        
        ok = download_cds(og_id, output_fasta)
        if ok:
            print(f"OK ({n_species} specie, OG: {og_id})")
            results.append({"gene": gene_name, "og_id": og_id, "n_species": n_species})
        else:
            print("ERRORE download")
            not_found.append(gene_name)
        
        time.sleep(0.3)  # pausa per non sovraccaricare l'API
    
    # Salva log
    with open(os.path.join(LOG_DIR, "01_fetch_results.json"), "w") as f:
        json.dump({
            "found": results,
            "not_found": not_found,
            "too_few_species": too_few_species
        }, f, indent=2)
    
    print(f"\nCompletato: {len(results)} trovati, {len(not_found)} non trovati, {len(too_few_species)} con poche specie")

if __name__ == "__main__":
    main()
