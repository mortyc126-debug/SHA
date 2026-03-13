#!/bin/bash
# run_experiments.sh
# Launcher for SHA-256 GPU differential experiments.
# Run from the directory containing birthday_search_17 and p_adic_tower.

set -euo pipefail

RESULTS_DIR="results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS_DIR"

echo "==================================================================="
echo "  SHA-256 GPU Differential Cryptanalysis — Experiment Runner"
echo "  Results dir: $RESULTS_DIR"
echo "==================================================================="

# Detect GPU
nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv,noheader 2>/dev/null \
    | head -1 | awk -F, '{printf "GPU: %s | VRAM: %s | Compute: %s\n", $1, $2, $3}'
echo ""

# ===================================================================
# EXPERIMENT 1: A1 — GPU Birthday Search (collect 50 pairs)
# ===================================================================
echo "--- Experiment A1: Birthday Search (De3..De17=0) ---"
echo "Scanning 50 random W0 values (1 full 2^32 scan each)"
echo "Expected: ~0.04s per pair on A100, ~0.3s on RTX 3090"
echo ""

time ./birthday_search_17 0 1 50 "$RESULTS_DIR/pairs_a1.csv"

NPAIRS=$(grep -c "^0x" "$RESULTS_DIR/pairs_a1.csv" 2>/dev/null || echo 0)
echo ""
echo "A1 result: $NPAIRS pairs found"
echo "Output: $RESULTS_DIR/pairs_a1.csv"
echo ""

# ===================================================================
# EXPERIMENT 2: H3 — DW0 scan (all 32 bit positions)
# ===================================================================
echo "--- Experiment H3: DW0 = 2^j scan (j=0..31) ---"
echo "Testing all 32 values of DW0 with fixed W0"
echo "Hypothesis: P(De17=0) should be ~2^{-32} for all j"
echo "If any j gives significantly more hits -> structure found!"
echo ""

W0_FIXED="e82222c7"
for j in $(seq 0 31); do
    DW0=$(python3 -c "print(hex(1 << $j)[2:])")
    OUTF="$RESULTS_DIR/dw0_j${j}.csv"
    # Scan 1 W0 for each j (fast, ~0.04s per j on A100)
    ./birthday_search_17 $W0_FIXED $DW0 1 "$OUTF" 2>/dev/null
    N=$(grep -c "^0x" "$OUTF" 2>/dev/null || echo 0)
    printf "  j=%2d  DW0=2^%2d=0x%08x  pairs=%d\n" $j $j $((1 << j)) $N
done
echo ""

# ===================================================================
# EXPERIMENT 3: H1 — p-adic Tower Height Test
# ===================================================================
echo "--- Experiment H1: p-adic Tower Height ---"
echo "Testing height_2(SHA-256) >= 11 (from П-53)"
echo "Using 1M seeds, testing up to Sol_25"
echo ""

time ./p_adic_tower 1000000 25 "$RESULTS_DIR/tower_h1.txt"

echo ""
echo "H1 result:"
cat "$RESULTS_DIR/tower_h1.txt" | grep -v "^#" | head -30
echo ""

# ===================================================================
# EXPERIMENT 4: Statistical analysis of found pairs (A1 results)
# ===================================================================
if [ -s "$RESULTS_DIR/pairs_a1.csv" ]; then
    echo "--- Statistical analysis of A1 pairs ---"
    python3 - <<'PYEOF' "$RESULTS_DIR/pairs_a1.csv"
import sys, struct

fname = sys.argv[1]
pairs = []
with open(fname) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'): continue
        parts = line.split(',')
        if len(parts) < 5: continue
        W0,w1,dw0,Da13,DW16 = [int(x,16) for x in parts]
        pairs.append({'W0':W0,'w1':w1,'dw0':dw0,'Da13':Da13,'DW16':DW16})

print(f"  Total pairs: {len(pairs)}")
if pairs:
    # HW distribution of Da13
    hws = [bin(p['Da13']).count('1') for p in pairs]
    print(f"  HW(Da13): min={min(hws)} max={max(hws)} mean={sum(hws)/len(hws):.2f}")
    print(f"  (Expected ~16 for uniform distribution)")

    # Check Da13 + DW16 == 0 for all
    valid = sum(1 for p in pairs if (p['Da13'] + p['DW16']) & 0xFFFFFFFF == 0)
    print(f"  Verification Da13+DW16=0: {valid}/{len(pairs)}")

    # Distribution of low bits of Da13
    low4 = [p['Da13'] & 0xF for p in pairs]
    print(f"  Low nibble Da13 distribution: {sorted(set(low4))} (expect uniform)")
PYEOF
fi

echo ""
echo "==================================================================="
echo "  All experiments complete. Results in: $RESULTS_DIR/"
echo "==================================================================="
ls -la "$RESULTS_DIR/"
