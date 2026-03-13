#!/usr/bin/env python3
"""
h5_analysis.py — анализ результатов H5 State Differential.

Читает CSV из h5_state_diff (или birthday_search_17),
выводит статистику v2(Da17), v2(De18), v2(DW17) и ищет структуру.

Использование:
  python3 h5_analysis.py h5_results.csv
  python3 h5_analysis.py pairs.txt         # обычный birthday output
"""
import sys, collections, math

def v2(x):
    x = x & 0xFFFFFFFF
    if x == 0: return 32
    v = 0
    while x % 2 == 0: v += 1; x //= 2
    return v

def hw(x):
    return bin(x & 0xFFFFFFFF).count('1')

def parse_csv(fname):
    rows = []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split(',')
            if not parts[0].startswith('0x'): continue
            try:
                row = [int(x, 16) for x in parts if x.startswith('0x')]
                rows.append(row)
            except ValueError:
                pass
    return rows

def analyze_birthday(fname):
    """Analyze standard birthday_search output: W0,w1,DW0,Da13,DW16"""
    rows = parse_csv(fname)
    if not rows:
        print("No data found.")
        return
    print(f"=== Birthday pairs analysis: {len(rows)} pairs ===\n")

    da13_list = [r[3] for r in rows if len(r) >= 5]
    dw16_list = [r[4] for r in rows if len(r) >= 5]

    # v2(Da13)
    v2_da13 = [v2(x) for x in da13_list]
    cnt = collections.Counter(v2_da13)
    print("v2(Da13) distribution (expected Geometric: P(k)=1/2^{k+1}):")
    for k in sorted(cnt):
        exp = len(da13_list) / 2**(k+1)
        print(f"  v2={k:2d}: {cnt[k]:5d}  obs={cnt[k]/len(da13_list):.4f}  exp={exp/len(da13_list):.4f}")

    # HW(Da13)
    hws = [hw(x) for x in da13_list]
    print(f"\nHW(Da13): min={min(hws)} max={max(hws)} mean={sum(hws)/len(hws):.3f} (exp=16)")

    # Da13 + DW16 = 0 check
    ok = sum(1 for a,d in zip(da13_list, dw16_list) if (a+d) & 0xFFFFFFFF == 0)
    print(f"Da13+DW16=0: {ok}/{len(da13_list)}")

def analyze_h5(fname):
    """Analyze h5_state_diff output: full state diff at rounds 17-20"""
    rows = parse_csv(fname)
    if not rows:
        print("No data.")
        return

    # Columns: W0,w1,DW0,Da13,DW16,Da17,Db17,Dc17,Dd17,Df17,Dg17,Dh17,De18,De19,De20,DW17,DW18
    # Plus 3 v2 values at end (integers, not hex — parse separately)
    hex_rows = []
    int_rows = []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split(',')
            if not parts[0].startswith('0x'): continue
            try:
                hex_part = [int(x, 16) for x in parts if x.startswith('0x')]
                int_part = [int(x) for x in parts if not x.startswith('0x')]
                hex_rows.append(hex_part)
                int_rows.append(int_part)
            except ValueError:
                pass

    N = len(hex_rows)
    print(f"=== H5 State Differential Analysis: {N} pairs ===\n")
    if N == 0: return

    fields = ['W0','w1','DW0','Da13','DW16','Da17','Db17','Dc17','Dd17','Df17','Dg17','Dh17','De18','De19','De20','DW17','DW18']
    data = {f: [r[i] for r in hex_rows if i < len(r)] for i, f in enumerate(fields)}

    # --- Verify De17=0 indirectly (Da13+DW16=0) ---
    ok = sum(1 for a,d in zip(data['Da13'], data['DW16']) if (a+d)&0xFFFFFFFF == 0)
    print(f"De17=0 check (Da13+DW16=0): {ok}/{N}")

    # --- Round 17 state differentials ---
    print("\nRound 17 state differentials (De17=0 by construction):")
    for name in ['Da17','Db17','Dc17','Dd17','Df17','Dg17','Dh17']:
        vals = data.get(name, [])
        if not vals: continue
        hws = [hw(x) for x in vals]
        v2s = [v2(x) for x in vals]
        print(f"  {name}: HW={sum(hws)/len(hws):.2f}  v2_avg={sum(v2s)/len(v2s):.2f}  "
              f"v2=0: {sum(1 for v in v2s if v==0)/len(v2s)*100:.1f}%")

    # --- Key H5 question: De18 distribution ---
    de18 = data.get('De18', [])
    if de18:
        v2_de18 = [v2(x) for x in de18]
        cnt = collections.Counter(v2_de18)
        print(f"\n--- De18 distribution (H5 key question) ---")
        print(f"De18=0 exactly: {sum(1 for x in de18 if x==0)}/{N}")
        print(f"v2(De18) distribution (exp Geometric if random):")
        for k in sorted(cnt)[:8]:
            exp_p = 1/2**(k+1)
            obs_p = cnt[k]/N
            diff = obs_p - exp_p
            flag = "*** BIAS ***" if abs(diff) > 0.05 else ""
            print(f"  v2={k}: {cnt[k]:5d}  obs={obs_p:.4f}  exp={exp_p:.4f}  {flag}")
        print(f"HW(De18): mean={sum(hw(x) for x in de18)/N:.2f} (exp=16)")

    # --- De19, De20 ---
    for name in ['De19','De20']:
        vals = data.get(name, [])
        if not vals: continue
        v2s = [v2(x) for x in vals]
        print(f"\n{name}: De=0: {sum(1 for x in vals if x==0)}/{N}  "
              f"v2_avg={sum(v2s)/len(v2s):.2f}  HW={sum(hw(x) for x in vals)/N:.2f}")

    # --- DW17: schedule difference ---
    dw17 = data.get('DW17', [])
    if dw17:
        v2_dw17 = [v2(x) for x in dw17]
        print(f"\nDW17 (schedule diff at r=18):")
        print(f"  DW17=0: {sum(1 for x in dw17 if x==0)}/{N}")
        cnt = collections.Counter(v2_dw17)
        for k in sorted(cnt)[:6]:
            print(f"  v2={k}: {cnt[k]:5d}  ({cnt[k]/N*100:.1f}%)")

    # --- Correlation: v2(Da13) vs v2(De18) ---
    if de18 and data.get('Da13'):
        print("\nCorrelation v2(Da13) vs v2(De18):")
        pairs_v2 = list(zip([v2(x) for x in data['Da13']], [v2(x) for x in de18]))
        for v_a in range(5):
            subset = [p[1] for p in pairs_v2 if p[0] == v_a]
            if subset:
                print(f"  v2(Da13)={v_a}: n={len(subset)}  avg_v2(De18)={sum(subset)/len(subset):.2f}")

    # --- Summary ---
    print("\n=== Summary ===")
    de18_bias = False
    if de18:
        v2s = [v2(x) for x in de18]
        obs_v0 = sum(1 for v in v2s if v==0)/N
        if abs(obs_v0 - 0.5) > 0.05:
            de18_bias = True
            print(f"⚠ De18 v2=0: {obs_v0:.3f} (exp=0.500) → BIAS DETECTED")
        else:
            print(f"✓ De18: no significant bias (v2=0: {obs_v0:.3f})")
    if not de18_bias:
        print("  De18 ≈ Uniform(Z/2^32) → no exploitable structure after round 17")
        print("  T_B1_DIFFUSION confirmed: differential fully randomized after 1 round")
    else:
        print("  *** POTENTIAL STRUCTURE FOUND — investigate! ***")

if __name__ == '__main__':
    fname = sys.argv[1] if len(sys.argv) > 1 else 'h5_results.csv'
    rows = parse_csv(fname)
    if not rows:
        print("Empty file or no data.")
        sys.exit(1)
    # Detect file type by column count
    with open(fname) as f:
        for line in f:
            if line.startswith('#'): continue
            cols = line.strip().split(',')
            n_hex = sum(1 for c in cols if c.startswith('0x'))
            break

    if n_hex >= 12:
        analyze_h5(fname)
    else:
        analyze_birthday(fname)
        print("\n(For full H5 analysis, use h5_state_diff binary)")
