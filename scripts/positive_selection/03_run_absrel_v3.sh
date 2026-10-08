#!/bin/bash
BASE_DIR="$HOME/tesi_bioinfo/positive_selection_orthodb"
ALN_DIR="$BASE_DIR/alignments_v3"
RES_DIR="$BASE_DIR/results_v3"
LOG_DIR="$BASE_DIR/logs"
TREE="$BASE_DIR/primates_tree.nwk"
HYPHY="$HOME/miniconda3/bin/hyphy"

mkdir -p "$RES_DIR"

total=$(ls "$ALN_DIR"/*_aligned.fasta 2>/dev/null | wc -l)
count=0
done=0
failed=0

echo "Geni da processare: $total"
echo "Inizio: $(date)"

for aln in "$ALN_DIR"/*_aligned.fasta; do
    gene=$(basename "$aln" _aligned.fasta)
    result_json="${aln}.ABSREL.json"
    final_json="$RES_DIR/${gene}_absrel.json"
    count=$((count + 1))

    if [ -f "$final_json" ] && [ -s "$final_json" ]; then
        echo "[$count/$total] $gene: già presente, salto"
        done=$((done + 1))
        continue
    fi

    echo -n "[$count/$total] $gene ... "

    printf "${aln}\n${TREE}\nAll\nNone\nNo\nYes\n" | \
        "$HYPHY" absrel > "$LOG_DIR/${gene}_hyphy.log" 2>&1

    if [ -f "$result_json" ] && [ -s "$result_json" ]; then
        mv "$result_json" "$final_json"
        echo "OK"
        done=$((done + 1))
    else
        echo "ERRORE"
        failed=$((failed + 1))
    fi
done

echo ""
echo "Fine: $(date)"
echo "Completati: $done / Errori: $failed"
