#!/bin/bash
################################################################################
# MiXCR re-run with correct preset for bulk RNA-seq (no UMI enrichment)
#
# Previous run used: vergani-et-al-2017-full-length  → 0% aligned (wrong preset)
# This run uses:     rnaseq                           → designed for bulk RNA-seq
#
# Samples:
#   4-1BB+ : C01 (S7), E01 (S9), F01 (S10)
#   4-1BB- : D01 (S8), G01 (S11), H01 (S12)
#
# Run from: ~/LABMEMBERS/RENO/
# Usage:
#   chmod +x run_mixcr_rnaseq.sh
#   ./run_mixcr_rnaseq.sh
################################################################################

set -euo pipefail

WORKDIR="$HOME/LABMEMBERS/RENO"
FASTQDIR="$WORKDIR/fastq"          # use the 20M reads folder (bulk)
OUTDIR="$WORKDIR/mixcr_rnaseq_output"
MIXCR="mixcr"                      # assumes mixcr is in PATH; adjust if needed

mkdir -p "$OUTDIR"

# Check mixcr is available
if ! command -v mixcr &> /dev/null; then
    echo "ERROR: mixcr not found in PATH."
    echo "If installed elsewhere, set MIXCR variable at top of this script."
    exit 1
fi

echo "MiXCR version:"
mixcr --version

# ── Sample definitions ────────────────────────────────────────────────────────
SAMPLES=("C01" "E01" "F01" "D01" "G01" "H01")
SNUMS=("S7" "S9" "S10" "S8" "S11" "S12")
GROUPS=("4-1BBpos" "4-1BBpos" "4-1BBpos" "4-1BBneg" "4-1BBneg" "4-1BBneg")

# ── Run MiXCR per sample ──────────────────────────────────────────────────────
for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"
    SNUM="${SNUMS[$i]}"
    GROUP="${GROUPS[$i]}"

    R1="$FASTQDIR/I25-1078-${SAMPLE}_${SNUM}_R1_001.fastq.gz"
    R2="$FASTQDIR/I25-1078-${SAMPLE}_${SNUM}_R2_001.fastq.gz"
    PREFIX="$OUTDIR/${SAMPLE}_rnaseq"

    echo ""
    echo "========================================"
    echo "Processing: $SAMPLE ($GROUP)"
    echo "========================================"

    if [ ! -f "$R1" ]; then
        echo "  [WARNING] R1 not found: $R1 — skipping"
        continue
    fi

    if [ -f "${PREFIX}.clones_TRB.tsv" ]; then
        echo "  [SKIP] Already processed: $SAMPLE"
        continue
    fi

    # Run MiXCR with rnaseq preset
    # --species hsa = human
    # rnaseq preset: handles fragmented reads, no UMI assumption
    mixcr analyze rnaseq \
        --species hsa \
        --report "${PREFIX}.report.txt" \
        --json-report "${PREFIX}.report.json" \
        "$R1" "$R2" \
        "$PREFIX" \
        2>&1 | tee "$OUTDIR/${SAMPLE}_mixcr.log"

    echo "  Done: $SAMPLE"
done

echo ""
echo "========================================"
echo "All samples processed."
echo "Output in: $OUTDIR"
echo "========================================"

# ── Quick summary of alignment rates ─────────────────────────────────────────
echo ""
echo "=== Alignment summary ==="
for SAMPLE in "${SAMPLES[@]}"; do
    REPORT="$OUTDIR/${SAMPLE}_rnaseq.report.txt"
    if [ -f "$REPORT" ]; then
        echo ""
        echo "--- $SAMPLE ---"
        grep -E "Total sequencing reads|Successfully aligned|Alignment failed: no hits" "$REPORT" || true
    fi
done

echo ""
echo "Next step: check alignment rates above."
echo "  If > 0% aligned → run Figure6C_mixcr_rnaseq_analysis.R"
echo "  If still 0%     → bulk RNA-seq truly has no recoverable TCR reads"
