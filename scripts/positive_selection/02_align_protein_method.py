#!/usr/bin/env python3
import os, json, subprocess, requests, time

BASE_DIR = os.path.expanduser("~/tesi_bioinfo/positive_selection_orthodb")
SEQ_DIR  = os.path.join(BASE_DIR, "sequences")
ALN_DIR  = os.path.join(BASE_DIR, "alignments_v3")
LOG_DIR  = os.path.join(BASE_DIR, "logs")
DEG_FILE = os.path.expanduser("~/tesi_bioinfo/results_deseq2/Week4_FGF4_vs_FGF8/Week4_FGF4_vs_FGF8_protein_coding_DEG_annotated.tsv")
PAL2NAL  = os.path.expanduser("~/miniconda3/bin/pal2nal.pl")
API_BASE = "https://data.orthodb.org/v12"

os.makedirs(ALN_DIR, exist_ok=True)

TREE_SPECIES = {
    "Homo sapiens":"Homo_sapiens","Pan troglodytes":"Pan_troglodytes",
    "Gorilla gorilla gorilla":"Gorilla_gorilla_gorilla","Pongo abelii":"Pongo_abelii",
    "Nomascus leucogenys":"Nomascus_leucogenys","Macaca mulatta":"Macaca_mulatta",
    "Macaca fascicularis":"Macaca_fascicularis","Rhinopithecus roxellana":"Rhinopithecus_roxellana",
    "Callithrix jacchus":"Callithrix_jacchus","Cebus imitator":"Cebus_imitator",
    "Microcebus murinus":"Microcebus_murinus","Otolemur garnettii":"Otolemur_garnettii",
}
STOP = {'TAA','TAG','TGA'}
MIN_SP = 6

# Carica log esistente per OG IDs già trovati
log_file = os.path.join(LOG_DIR, '01_fetch_results.json')
og_cache = {}
with open(log_file) as f:
    data = json.load(f)
for entry in data.get('found', []):
    og_cache[entry['gene']] = entry['og_id']

# Carica mappa gene_name -> ensembl_id dal file DEG
gene_to_ensembl = {}
with open(DEG_FILE) as f:
    f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            gene_to_ensembl[parts[0]] = parts[1]

def get_og_id(gene, ensembl_id):
    """Recupera OG ID da cache o da API."""
    if gene in og_cache:
        return og_cache[gene]
    # Cerca via API
    for query in [ensembl_id, gene]:
        try:
            url = f"{API_BASE}/search?query={query}&level=9443&limit=3"
            r = requests.get(url, timeout=15)
            if r.status_code == 200:
                data = r.json()
                if data.get('count','0') != '0' and data.get('data'):
                    og_id = data['data'][0]
                    og_cache[gene] = og_id
                    time.sleep(0.3)
                    return og_id
        except:
            pass
        time.sleep(0.3)
    return None

def parse_fasta(text):
    seqs = {}
    current_sp = None
    current_seq = []
    seen = set()
    for line in text.strip().split('\n'):
        if line.startswith('>'):
            if current_sp and current_sp in TREE_SPECIES and TREE_SPECIES[current_sp] not in seen:
                sp = TREE_SPECIES[current_sp]
                seqs[sp] = ''.join(current_seq)
                seen.add(sp)
            try:
                info = json.loads(line[line.index('{'):])
                current_sp = info.get('organism_name','')
            except:
                current_sp = None
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_sp and current_sp in TREE_SPECIES and TREE_SPECIES[current_sp] not in seen:
        sp = TREE_SPECIES[current_sp]
        seqs[sp] = ''.join(current_seq)
    return seqs

def rm_stop(seq):
    if len(seq) >= 3 and seq[-3:].upper() in STOP:
        seq = seq[:-3]
    remainder = len(seq) % 3
    if remainder:
        seq = seq[:len(seq)-remainder]
    return seq

done, skipped_few, skipped_pal2nal, skipped_no_og = [], [], [], []
genes = sorted(os.listdir(SEQ_DIR))
total = len(genes)

for i, gene in enumerate(genes):
    aln_out = os.path.join(ALN_DIR, f"{gene}_aligned.fasta")
    if os.path.exists(aln_out) and os.path.getsize(aln_out) > 0:
        print(f"[{i+1}/{total}] {gene}: già presente, salto")
        done.append(gene); continue

    ensembl_id = gene_to_ensembl.get(gene, gene)
    og_id = get_og_id(gene, ensembl_id)
    if not og_id:
        print(f"[{i+1}/{total}] {gene}: OG ID non trovato, salto")
        skipped_no_og.append(gene); continue

    try:
        r_cds  = requests.get(f"{API_BASE}/fasta?id={og_id}&seqtype=cds", timeout=30)
        r_prot = requests.get(f"{API_BASE}/fasta?id={og_id}&seqtype=protein", timeout=30)
        time.sleep(0.3)
    except:
        print(f"[{i+1}/{total}] {gene}: errore download, salto")
        skipped_no_og.append(gene); continue

    cds_seqs  = parse_fasta(r_cds.text)
    prot_seqs = parse_fasta(r_prot.text)
    common = set(cds_seqs.keys()) & set(prot_seqs.keys())

    if len(common) < MIN_SP:
        print(f"[{i+1}/{total}] {gene}: solo {len(common)} specie, salto")
        skipped_few.append(gene); continue

    cds_clean = {sp: rm_stop(cds_seqs[sp]) for sp in common}

    tmp_prot = os.path.join(ALN_DIR, f"{gene}_tmp_prot.fasta")
    tmp_cds  = os.path.join(ALN_DIR, f"{gene}_tmp_cds.fasta")
    with open(tmp_prot,'w') as f:
        for sp in common: f.write(f'>{sp}\n{prot_seqs[sp]}\n')
    with open(tmp_cds,'w') as f:
        for sp in common: f.write(f'>{sp}\n{cds_clean[sp]}\n')

    r = subprocess.run(['mafft','--auto','--quiet', tmp_prot], capture_output=True, text=True)
    if r.returncode != 0 or not r.stdout.startswith('>'):
        print(f"[{i+1}/{total}] {gene}: MAFFT fallito")
        for fp in [tmp_prot, tmp_cds]:
            if os.path.exists(fp): os.remove(fp)
        skipped_pal2nal.append(gene); continue

    tmp_prot_aln = os.path.join(ALN_DIR, f"{gene}_tmp_prot_aln.fasta")
    with open(tmp_prot_aln,'w') as f: f.write(r.stdout)

    r2 = subprocess.run(
        ['perl', PAL2NAL, tmp_prot_aln, tmp_cds, '-output', 'fasta', '-nogap'],
        capture_output=True, text=True
    )
    for fp in [tmp_prot, tmp_cds, tmp_prot_aln]:
        if os.path.exists(fp): os.remove(fp)

    if r2.returncode != 0 or not r2.stdout.startswith('>'):
        print(f"[{i+1}/{total}] {gene}: PAL2NAL fallito")
        skipped_pal2nal.append(gene); continue

    with open(aln_out,'w') as f: f.write(r2.stdout)
    print(f"[{i+1}/{total}] {gene}: OK ({len(common)} specie)")
    done.append(gene)

with open(os.path.join(LOG_DIR,'02_align_v4.json'),'w') as f:
    json.dump({'done':done,'few':skipped_few,'pal2nal':skipped_pal2nal,'no_og':skipped_no_og},f,indent=2)
print(f"\nFine: {len(done)} OK, {len(skipped_few)} poche specie, {len(skipped_pal2nal)} PAL2NAL, {len(skipped_no_og)} no OG")
