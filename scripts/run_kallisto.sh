#!/bin/bash

# Percorsi
DATA="/mnt/c/Users/mirab/Desktop/RNA-seq_data"
IDX="/home/mirab/tesi_bioinfo/index/gencode.v46.idx"
OUT="/home/mirab/tesi_bioinfo/results"
THREADS=4

echo "=== Inizio mappatura Kallisto ==="

# ─── WEEK1 ───────────────────────────────────────────────────────────────────
for REP in Rep1 Rep2; do
    DIR="$DATA/Week1/$REP"
    SAMPLE=$(ls $DIR/*.fastq.gz | head -1 | xargs basename | cut -d'_' -f1)

    R1_L006=$(ls $DIR/${SAMPLE}_*_L006_R1_001.fastq.gz)
    R2_L006=$(ls $DIR/${SAMPLE}_*_L006_R2_001.fastq.gz)
    R1_L007=$(ls $DIR/${SAMPLE}_*_L007_R1_001.fastq.gz)
    R2_L007=$(ls $DIR/${SAMPLE}_*_L007_R2_001.fastq.gz)

    mkdir -p "$OUT/$SAMPLE"
    echo "Mappando $SAMPLE (Week1/$REP)..."

    kallisto quant \
        -i "$IDX" \
        -o "$OUT/$SAMPLE" \
        -t $THREADS \
        $R1_L006 $R2_L006 $R1_L007 $R2_L007

    echo "✓ $SAMPLE completato"
done

# ─── WEEK2+ ──────────────────────────────────────────────────────────────────
for WEEK_DIR in "$DATA"/Week2_* "$DATA"/Week3_* "$DATA"/Week4_*; do
    WEEK=$(basename $WEEK_DIR)
    echo "--- Processando $WEEK ---"

    SAMPLES=$(ls $WEEK_DIR/*.fastq.gz | xargs -I{} basename {} | cut -d'_' -f1 | sort -u)

    for SAMPLE in $SAMPLES; do
        R1_L006=$(ls $WEEK_DIR/${SAMPLE}_*_L006_R1_001.fastq.gz)
        R2_L006=$(ls $WEEK_DIR/${SAMPLE}_*_L006_R2_001.fastq.gz)
        R1_L007=$(ls $WEEK_DIR/${SAMPLE}_*_L007_R1_001.fastq.gz)
        R2_L007=$(ls $WEEK_DIR/${SAMPLE}_*_L007_R2_001.fastq.gz)

        mkdir -p "$OUT/$SAMPLE"
        echo "Mappando $SAMPLE ($WEEK)..."

        kallisto quant \
            -i "$IDX" \
            -o "$OUT/$SAMPLE" \
            -t $THREADS \
            $R1_L006 $R2_L006 $R1_L007 $R2_L007

        echo "✓ $SAMPLE completato"
    done
done

echo "=== Mappatura completata! ==="
