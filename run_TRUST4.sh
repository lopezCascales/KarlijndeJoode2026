#!/bin/bash
################################################################################
# TRUST4 installation + TCR clonotype analysis
# 4-1BB+ : C01, F01, E01
# 4-1BB- : D01, G01, H01
#
# Run from: ~/LABMEMBERS/RENO/
# FASTQs in: ~/LABMEMBERS/RENO/fastq/
#
# Usage:
#   chmod +x run_TRUST4.sh
#   ./run_TRUST4.sh
################################################################################

set -euo pipefail

WORKDIR="$HOME/LABMEMBERS/RENO"
FASTQDIR="$WORKDIR/fastq"
OUTDIR="$WORKDIR/trust4_output"
TRUST4DIR="$WORKDIR/TRUST4"

mkdir -p "$OUTDIR"

# ── 1. Install TRUST4 ─────────────────────────────────────────────────────────
echo "=== Installing TRUST4 ==="

if [ ! -f "$TRUST4DIR/run-trust4" ]; then
    cd "$WORKDIR"
    git clone https://github.com/liulab-dfci/TRUST4.git
    cd TRUST4
    make
    echo "TRUST4 compiled successfully."
else
    echo "TRUST4 already installed, skipping."
fi

TRUST4="$TRUST4DIR/run-trust4"

# ── 2. Download human TCR/BCR reference (hg38) ───────────────────────────────
echo ""
echo "=== Downloading TRUST4 human reference ==="

REFDIR="$TRUST4DIR/hg38_ref"
mkdir -p "$REFDIR"

if [ ! -f "$REFDIR/hg38_bcrtcr.fa" ]; then
    # TRUST4 provides a ready-made human reference
    cd "$REFDIR"
    wget -q https://raw.githubusercontent.com/liulab-dfci/TRUST4/master/hg38/hg38_bcrtcr.fa
    wget -q https://raw.githubusercontent.com/liulab-dfci/TRUST4/master/hg38/IMGT+C.fa
    echo "Reference files downloaded."
else
    echo "Reference files already present, skipping."
fi

BCRTCR_REF="$REFDIR/hg38_bcrtcr.fa"
IMGT_REF="$REFDIR/IMGT+C.fa"

# ── 3. Define samples ─────────────────────────────────────────────────────────
declare -A SAMPLES=(
    # Sample_ID   Group
    ["C01"]="4-1BBpos"
    ["F01"]="4-1BBpos"
    ["E01"]="4-1BBpos"
    ["D01"]="4-1BBneg"
    ["G01"]="4-1BBneg"
    ["H01"]="4-1BBneg"
)

# Map sample ID to FASTQ suffix (S number)
declare -A SNUM=(
    ["A01"]="S5"  ["B01"]="S6"  ["C01"]="S7"  ["D01"]="S8"
    ["E01"]="S9"  ["F01"]="S10" ["G01"]="S11" ["H01"]="S12"
    ["J01"]="S13" ["K01"]="S14"
)

# ── 4. Run TRUST4 per sample ──────────────────────────────────────────────────
echo ""
echo "=== Running TRUST4 ==="

for SAMPLE in "${!SAMPLES[@]}"; do
    SNUM_VAL="${SNUM[$SAMPLE]}"
    R1="$FASTQDIR/I25-1078-${SAMPLE}_${SNUM_VAL}_R1_001.fastq.gz"
    R2="$FASTQDIR/I25-1078-${SAMPLE}_${SNUM_VAL}_R2_001.fastq.gz"
    OUTPREFIX="$OUTDIR/${SAMPLE}_trust4"

    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "  [WARNING] FASTQs not found for $SAMPLE: $R1"
        continue
    fi

    if [ -f "${OUTPREFIX}_report.tsv" ]; then
        echo "  [SKIP] $SAMPLE already processed."
        continue
    fi

    echo "  Processing $SAMPLE (${SAMPLES[$SAMPLE]})..."
    "$TRUST4" \
        -u "$R1" \
        --secondary-read "$R2" \
        -f "$BCRTCR_REF" \
        --ref "$IMGT_REF" \
        -o "$OUTPREFIX" \
        --od "$OUTDIR" \
        -t 4 \
        2> "$OUTDIR/${SAMPLE}_trust4.log"

    echo "  Done: $SAMPLE"
done

echo ""
echo "=== TRUST4 runs complete ==="
echo "Output files in: $OUTDIR"
echo ""

# ── 5. Check output ───────────────────────────────────────────────────────────
echo "=== Report files generated ==="
ls -lh "$OUTDIR"/*_report.tsv 2>/dev/null || echo "No report files found — check logs in $OUTDIR"

echo ""
echo "Next step: run Figure6C_TRUST4_analysis.R"
