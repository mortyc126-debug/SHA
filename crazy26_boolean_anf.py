#!/usr/bin/env python3
"""
crazy26_boolean_anf.py — Algebraic Normal Form analysis of De17(W[14])

Determines the algebraic degree of De17 (XOR of SHA-256 register e after
round 17 between original and Wang-twin) as a Boolean function of W[14] bits.

Approach: restrict De17 to 16-bit subsets of W[14], compute ANF via Möbius
transform, measure degree and term distribution.
"""

import random
import time
import numpy as np

MASK = 0xFFFFFFFF
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,
     0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,
     0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,
     0x0fc19dc6,0x240ca1cc]
H0 = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Ch(e,f,g): return ((e&f)^(~e&g))&MASK
def Maj(a,b,c): return ((a&b)^(a&c)^(b&c))&MASK

def add32(*a):
    s=0
    for x in a: s=(s+x)&MASK
    return s

def sha_round(st,w,k):
    a,b,c,d,e,f,g,h=st
    T1=add32(h,Sig1(e),Ch(e,f,g),k,w)
    T2=add32(Sig0(a),Maj(a,b,c))
    return [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

def expand_W(W16):
    W=list(W16)
    for i in range(16,20):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def get_de17(msg, iv=list(H0)):
    W=list(msg); Wp=list(W); Wp[0]^=0x80000000
    s=list(iv); sp=list(iv)
    s=sha_round(s,W[0],K[0]); sp=sha_round(sp,Wp[0],K[0])
    for t in range(1,16):
        a,b,c,d,e,f,g,h=s
        a2,b2,c2,d2,e2,f2,g2,h2=sp
        tp=add32(h,Sig1(e),Ch(e,f,g),K[t])
        tp2=add32(h2,Sig1(e2),Ch(e2,f2,g2),K[t])
        target=add32(d,tp,W[t])
        Wp[t]=(target-d2-tp2)&MASK
        s=sha_round(s,W[t],K[t]); sp=sha_round(sp,Wp[t],K[t])
    We=expand_W(msg); Wpe=expand_W(Wp)
    s1=list(iv); s2=list(iv)
    for t in range(17):
        s1=sha_round(s1,We[t],K[t]); s2=sha_round(s2,Wpe[t],K[t])
    return s1[4]^s2[4]

def popcount(x):
    c = 0
    while x:
        c += 1
        x &= x - 1
    return c

def mobius_transform(tt, n):
    """In-place Möbius transform on truth table array of size 2^n."""
    for i in range(n):
        step = 1 << i
        for j in range(1 << n):
            if j & step:
                tt[j] ^= tt[j ^ step]

def analyze_anf(anf, n):
    """Analyze ANF: degree, term counts per degree, density."""
    degree_counts = [0] * (n + 1)
    max_degree = 0
    total_terms = 0
    for s in range(1 << n):
        if anf[s]:
            d = popcount(s)
            degree_counts[d] += 1
            total_terms += 1
            if d > max_degree:
                max_degree = d
    return max_degree, degree_counts, total_terms

def main():
    print("=" * 70)
    print("CRAZY-26: Algebraic Normal Form Analysis of De17(W[14])")
    print("=" * 70)

    # Base message from seed
    rng = random.Random(0xC026)
    base_msg = [rng.randint(0, MASK) for _ in range(16)]
    print(f"\nBase message W[14] = 0x{base_msg[14]:08x}")

    # Verify De17 at base
    de17_base = get_de17(base_msg)
    print(f"De17 at base = 0x{de17_base:08x}  (hw={bin(de17_base).count('1')})")

    N_SUBSETS = 5
    OUTPUT_BITS = [0, 7, 15]  # 3 different output bits to study
    n = 16  # 16-bit subsets

    all_results = []

    for subset_idx in range(N_SUBSETS):
        # Pick random 16 bit positions from the 32 bits of W[14]
        rng_sub = random.Random(0xC026 + subset_idx * 137)
        bit_positions = sorted(rng_sub.sample(range(32), n))

        print(f"\n{'=' * 70}")
        print(f"SUBSET {subset_idx}: bit positions {bit_positions}")
        print(f"{'=' * 70}")

        # Build truth table: enumerate all 2^16 values for selected bits
        # The other 16 bits stay at their base_msg[14] values
        base_w14 = base_msg[14]
        # Mask for the fixed bits
        var_mask = 0
        for bp in bit_positions:
            var_mask |= (1 << bp)
        fixed_val = base_w14 & (~var_mask & MASK)

        t0 = time.time()

        # Compute all 2^16 De17 values
        # For each input index i (0..2^16-1), set the selected bits accordingly
        de17_values = np.zeros(1 << n, dtype=np.uint32)

        for i in range(1 << n):
            # Map the 16-bit index to actual W[14] value
            w14 = fixed_val
            for k, bp in enumerate(bit_positions):
                if i & (1 << k):
                    w14 |= (1 << bp)
            msg_copy = list(base_msg)
            msg_copy[14] = w14
            de17_values[i] = get_de17(msg_copy)

        t1 = time.time()
        print(f"  Computed {1<<n} De17 values in {t1-t0:.1f}s")

        for out_bit in OUTPUT_BITS:
            # Extract single output bit as truth table
            tt = np.zeros(1 << n, dtype=np.uint8)
            for i in range(1 << n):
                tt[i] = (int(de17_values[i]) >> out_bit) & 1

            # Weight of the Boolean function
            weight = int(np.sum(tt))
            balance = weight / (1 << n)

            # Möbius transform (in-place)
            anf = tt.copy()
            for bit_i in range(n):
                step = 1 << bit_i
                for j in range(1 << n):
                    if j & step:
                        anf[j] ^= anf[j ^ step]

            # Analyze
            max_deg, deg_counts, total_terms = analyze_anf(anf, n)

            print(f"\n  Output bit {out_bit}:")
            print(f"    Weight: {weight}/{1<<n} = {balance:.4f} (ideal 0.5)")
            print(f"    ANF degree: {max_deg}")
            print(f"    Total nonzero ANF terms: {total_terms} / {1<<n}")
            print(f"    Term density: {total_terms/(1<<n):.4f}")
            print(f"    Terms by degree:")

            from math import comb
            for d in range(n + 1):
                possible = comb(n, d)
                if deg_counts[d] > 0 or d <= 3 or d == n:
                    frac = deg_counts[d] / possible if possible > 0 else 0
                    print(f"      deg {d:2d}: {deg_counts[d]:6d} / {possible:6d} ({frac:.4f})")

            all_results.append({
                'subset': subset_idx,
                'out_bit': out_bit,
                'degree': max_deg,
                'total_terms': total_terms,
                'density': total_terms / (1 << n),
                'weight': weight,
                'balance': balance,
                'deg_counts': deg_counts
            })

    # Summary
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}")

    degrees = [r['degree'] for r in all_results]
    densities = [r['density'] for r in all_results]
    balances = [r['balance'] for r in all_results]

    print(f"\nAcross {len(all_results)} trials ({N_SUBSETS} subsets x {len(OUTPUT_BITS)} output bits):")
    print(f"  Algebraic degree: min={min(degrees)}, max={max(degrees)}, "
          f"mean={sum(degrees)/len(degrees):.1f}")
    print(f"  ANF density:      min={min(densities):.4f}, max={max(densities):.4f}, "
          f"mean={sum(densities)/len(densities):.4f}")
    print(f"  Balance (weight):  min={min(balances):.4f}, max={max(balances):.4f}, "
          f"mean={sum(balances)/len(balances):.4f}")

    # Degree distribution summary
    print(f"\n  Average terms by degree (across all trials):")
    from math import comb
    for d in range(n + 1):
        avg_count = sum(r['deg_counts'][d] for r in all_results) / len(all_results)
        possible = comb(n, d)
        if avg_count > 0.1 or d <= 3 or d >= n - 1:
            avg_frac = avg_count / possible if possible > 0 else 0
            print(f"    deg {d:2d}: avg {avg_count:8.1f} / {possible:6d} "
                  f"(fraction {avg_frac:.4f})")

    print(f"\n  Interpretation:")
    mean_deg = sum(degrees) / len(degrees)
    if mean_deg >= n - 1:
        print(f"    De17 restricted to 16-bit subsets has FULL algebraic degree ({mean_deg:.1f}/16).")
        print(f"    This means De17 is algebraically complex — no low-degree shortcut.")
        print(f"    The function behaves like a random Boolean function (expected degree ~n).")
    elif mean_deg >= n * 0.7:
        print(f"    De17 has HIGH algebraic degree ({mean_deg:.1f}/16).")
        print(f"    Algebraic attacks are unlikely to be efficient.")
    elif mean_deg <= 3:
        print(f"    De17 has LOW algebraic degree ({mean_deg:.1f}/16)!")
        print(f"    This could enable Gröbner basis / linearization attacks.")
    else:
        print(f"    De17 has MODERATE algebraic degree ({mean_deg:.1f}/16).")

    mean_density = sum(densities) / len(densities)
    print(f"\n    ANF density {mean_density:.4f} ", end="")
    if abs(mean_density - 0.5) < 0.05:
        print("≈ 0.5 (consistent with random Boolean function).")
    elif mean_density < 0.3:
        print("is LOW — function has sparse ANF, potentially exploitable.")
    else:
        print(f"(reference: random function ≈ 0.5).")

    print(f"\nDone.")

if __name__ == "__main__":
    main()
