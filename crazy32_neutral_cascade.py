#!/usr/bin/env python3
"""
crazy32_neutral_cascade.py — Neutral bit cascading experiment for SHA-256 collision research.

Tests whether neutral bits (for De17) form a linear subspace by checking
pairs, triples, and larger k-subsets of simultaneously flipped bits.
"""

import random
import time

MASK = 0xFFFFFFFF

K_ext = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
]

H0 = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

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
def hw(x): return bin(x&MASK).count('1')

def sha_round(st,w,k):
    a,b,c,d,e,f,g,h=st
    T1=add32(h,Sig1(e),Ch(e,f,g),k,w)
    T2=add32(Sig0(a),Maj(a,b,c))
    return [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

def expand_W(W16, need):
    W=list(W16)
    for i in range(16, need):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def get_de(msg, target_round, iv=None):
    """Compute De_{target_round} = e_reg(msg) XOR e_reg(msg') after Wang twin construction."""
    if iv is None:
        iv = list(H0)
    W = list(msg)
    Wp = list(W)
    Wp[0] ^= 0x80000000  # flip MSB of W[0]

    # Run round 0 with original/flipped W[0]
    s = list(iv)
    sp = list(iv)
    s = sha_round(s, W[0], K_ext[0])
    sp = sha_round(sp, Wp[0], K_ext[0])

    # Rounds 1..15: correct Wp[t] so that twin states match after each round
    for t in range(1, 16):
        a,b,c,d,e,f,g,h = s
        a2,b2,c2,d2,e2,f2,g2,h2 = sp
        tp = add32(h, Sig1(e), Ch(e,f,g), K_ext[t])
        tp2 = add32(h2, Sig1(e2), Ch(e2,f2,g2), K_ext[t])
        target = add32(d, tp, W[t])
        Wp[t] = (target - d2 - tp2) & MASK
        s = sha_round(s, W[t], K_ext[t])
        sp = sha_round(sp, Wp[t], K_ext[t])

    # Expand both schedules
    need = max(target_round + 1, 17)
    We = expand_W(W, need)
    Wpe = expand_W(Wp, need)

    # Re-run from IV through target_round
    s1 = list(iv)
    s2 = list(iv)
    for t in range(target_round):
        s1 = sha_round(s1, We[t], K_ext[t])
        s2 = sha_round(s2, Wpe[t], K_ext[t])

    return s1[4] ^ s2[4]  # e register XOR


def flip_bits(msg, positions):
    """Return a copy of msg with specified bit positions flipped.
    Position = word_idx * 32 + bit_idx (bit 0 = LSB)."""
    m = list(msg)
    for pos in positions:
        word = pos // 32
        bit = pos % 32
        m[word] ^= (1 << bit)
    return m


def main():
    random.seed(0xC032)
    print("=" * 72)
    print("NEUTRAL BIT CASCADING EXPERIMENT")
    print("Seed: 0xC032 | Target: De17, De18, De19")
    print("=" * 72)

    # Step 1: Generate base message
    msg = [random.getrandbits(32) for _ in range(16)]
    print(f"\nBase message W[0..3]: {[hex(w) for w in msg[:4]]}")

    de17_base = get_de(msg, 17)
    de18_base = get_de(msg, 18)
    de19_base = get_de(msg, 19)
    print(f"Base De17 = {hex(de17_base)}, hw = {hw(de17_base)}")
    print(f"Base De18 = {hex(de18_base)}, hw = {hw(de18_base)}")
    print(f"Base De19 = {hex(de19_base)}, hw = {hw(de19_base)}")

    # Step 2: Find all neutral bits for De17
    print("\n" + "-" * 72)
    print("STEP 2: Finding neutral bits (De17 unchanged when flipped individually)")
    print("-" * 72)

    t0 = time.time()
    neutral_bits = []
    for pos in range(512):
        m2 = flip_bits(msg, [pos])
        de17_flip = get_de(m2, 17)
        if de17_flip == de17_base:
            neutral_bits.append(pos)

    elapsed = time.time() - t0
    print(f"Found {len(neutral_bits)} neutral bits out of 512 (in {elapsed:.2f}s)")
    # Show by word
    by_word = {}
    for pos in neutral_bits:
        w = pos // 32
        by_word.setdefault(w, []).append(pos % 32)
    for w in sorted(by_word.keys()):
        bits = sorted(by_word[w])
        print(f"  W[{w:2d}]: {len(bits):2d} neutral bits -> bit positions {bits}")

    if len(neutral_bits) < 2:
        print("Too few neutral bits found. Aborting.")
        return

    N = len(neutral_bits)

    # Step 3: Test ALL PAIRS
    print("\n" + "-" * 72)
    print(f"STEP 3: Testing ALL pairs of {N} neutral bits ({N*(N-1)//2} pairs)")
    print("-" * 72)

    t0 = time.time()
    total_pairs = 0
    compatible_pairs = 0
    interfering_pairs = []

    for i in range(N):
        for j in range(i+1, N):
            m2 = flip_bits(msg, [neutral_bits[i], neutral_bits[j]])
            de17_flip = get_de(m2, 17)
            total_pairs += 1
            if de17_flip == de17_base:
                compatible_pairs += 1
            else:
                interfering_pairs.append((neutral_bits[i], neutral_bits[j], hw(de17_flip ^ de17_base)))

    elapsed = time.time() - t0
    frac_pairs = compatible_pairs / total_pairs if total_pairs > 0 else 0
    print(f"Compatible pairs: {compatible_pairs} / {total_pairs} = {frac_pairs*100:.2f}%")
    print(f"Interfering pairs: {len(interfering_pairs)} (elapsed {elapsed:.2f}s)")
    if interfering_pairs:
        print(f"First 10 interfering pairs (pos_a, pos_b, hamming_dist_of_De17_change):")
        for a, b, h in interfering_pairs[:10]:
            print(f"    ({a:3d}, {b:3d})  hw_change={h}")

    # Step 4: Test RANDOM TRIPLES
    print("\n" + "-" * 72)
    print("STEP 4: Testing random TRIPLES of neutral bits (sample 10000)")
    print("-" * 72)

    t0 = time.time()
    n_triple_samples = min(10000, N*(N-1)*(N-2)//6)
    triple_ok = 0
    for _ in range(n_triple_samples):
        subset = random.sample(neutral_bits, 3)
        m2 = flip_bits(msg, subset)
        if get_de(m2, 17) == de17_base:
            triple_ok += 1
    elapsed = time.time() - t0
    frac_triple = triple_ok / n_triple_samples
    print(f"Triples preserving De17: {triple_ok} / {n_triple_samples} = {frac_triple*100:.2f}% ({elapsed:.2f}s)")

    # Step 5: Test RANDOM k-SUBSETS for k = 5, 10, 20, 40
    print("\n" + "-" * 72)
    print("STEP 5: Random k-subsets — fraction preserving De17")
    print("-" * 72)

    k_values = [5, 10, 20, 40]
    n_samples_per_k = 1000
    k_results = {}

    for k in k_values:
        if k > N:
            print(f"  k={k:3d}: SKIP (only {N} neutral bits)")
            continue
        t0 = time.time()
        ok = 0
        for _ in range(n_samples_per_k):
            subset = random.sample(neutral_bits, k)
            m2 = flip_bits(msg, subset)
            if get_de(m2, 17) == de17_base:
                ok += 1
        elapsed = time.time() - t0
        frac = ok / n_samples_per_k
        k_results[k] = frac
        print(f"  k={k:3d}: {ok}/{n_samples_per_k} = {frac*100:.1f}% preserving De17  ({elapsed:.2f}s)")

    # Determine effective subspace dimension
    print("\n  Effective neutral subspace dimension (>50% threshold):")
    eff_dim = 0
    # Also check k=2, k=3 from earlier
    all_k_frac = {2: frac_pairs, 3: frac_triple}
    all_k_frac.update(k_results)
    for k in sorted(all_k_frac.keys()):
        f = all_k_frac[k]
        marker = " <-- threshold" if f < 0.5 and (k == min(kk for kk in all_k_frac if all_k_frac[kk] < 0.5)) else ""
        if f >= 0.5:
            eff_dim = k
        print(f"    k={k:3d}: {f*100:.1f}%{marker}")

    # Step 6: Test finer granularity around threshold
    print("\n" + "-" * 72)
    print("STEP 5b: Finer granularity around threshold")
    print("-" * 72)
    fine_ks = [2, 3, 5, 8, 10, 15, 20, 25, 30, 40, 50, 60]
    fine_ks = [k for k in fine_ks if k <= N and k not in k_results and k not in [2, 3]]
    for k in fine_ks:
        t0 = time.time()
        ok = 0
        for _ in range(n_samples_per_k):
            subset = random.sample(neutral_bits, k)
            m2 = flip_bits(msg, subset)
            if get_de(m2, 17) == de17_base:
                ok += 1
        elapsed = time.time() - t0
        frac = ok / n_samples_per_k
        all_k_frac[k] = frac
        print(f"  k={k:3d}: {ok}/{n_samples_per_k} = {frac*100:.1f}%  ({elapsed:.2f}s)")

    # Step 7: CASCADE effect on De18, De19
    print("\n" + "-" * 72)
    print("STEP 7: Cascade effect on De18 and De19")
    print("Testing: do neutral bits for De17 preserve De18? De19?")
    print("-" * 72)

    for target_r, label in [(18, "De18"), (19, "De19")]:
        print(f"\n  --- {label} ---")
        # Single bits
        single_ok = 0
        de_base = get_de(msg, target_r)
        for pos in neutral_bits:
            m2 = flip_bits(msg, [pos])
            if get_de(m2, target_r) == de_base:
                single_ok += 1
        print(f"  Singles preserving {label}: {single_ok}/{N} = {single_ok*100/N:.1f}%")

        # k-subsets
        for k in [2, 5, 10, 20, 40]:
            if k > N:
                continue
            ok = 0
            ns = min(1000, n_samples_per_k)
            for _ in range(ns):
                subset = random.sample(neutral_bits, k)
                m2 = flip_bits(msg, subset)
                if get_de(m2, target_r) == de_base:
                    ok += 1
            frac = ok / ns
            print(f"  k={k:3d}: {ok}/{ns} = {frac*100:.1f}% preserving {label}")

    # Step 8: Hamming weight analysis for non-preserved cases
    print("\n" + "-" * 72)
    print("STEP 8: When De17 is NOT preserved, how bad is it?")
    print("Hamming weight of De17_flipped XOR De17_base for k-subset flips")
    print("-" * 72)

    for k in [5, 10, 20, 40]:
        if k > N:
            continue
        hws = []
        for _ in range(1000):
            subset = random.sample(neutral_bits, k)
            m2 = flip_bits(msg, subset)
            de = get_de(m2, 17)
            diff = de ^ de17_base
            if diff != 0:
                hws.append(hw(diff))
        if hws:
            avg_hw = sum(hws) / len(hws)
            print(f"  k={k:3d}: {len(hws)} non-zero cases, avg hw(De17 change) = {avg_hw:.2f}, max = {max(hws)}")
        else:
            print(f"  k={k:3d}: ALL preserved De17 (0 non-zero cases)")

    # Final verdict
    print("\n" + "=" * 72)
    print("VERDICT")
    print("=" * 72)

    # Recompute effective dimension
    eff_dim = 0
    for k in sorted(all_k_frac.keys()):
        if all_k_frac[k] >= 0.5:
            eff_dim = k

    print(f"\nTotal neutral bits found: {N}")
    print(f"Pair compatibility: {frac_pairs*100:.1f}%")
    print(f"Triple compatibility: {frac_triple*100:.1f}%")
    print(f"\nFraction preserving De17 by k:")
    for k in sorted(all_k_frac.keys()):
        bar = "#" * int(all_k_frac[k] * 50)
        print(f"  k={k:3d}: {all_k_frac[k]*100:5.1f}% |{bar}")

    print(f"\nEffective neutral subspace dimension (>50% preserved): k = {eff_dim}")

    if frac_pairs > 0.95:
        print("\nCONCLUSION: Neutral bits are HIGHLY COMPATIBLE in pairs.")
        if frac_triple > 0.90:
            print("Triples also highly compatible — strong evidence for near-linear subspace.")
        else:
            print("But triples show degradation — higher-order interactions exist.")
    elif frac_pairs > 0.50:
        print("\nCONCLUSION: Neutral bits are PARTIALLY COMPATIBLE. Some interactions reduce space.")
    else:
        print("\nCONCLUSION: Neutral bits are STRONGLY INTERACTING. Not a linear subspace.")

    if eff_dim >= 20:
        print(f"Effective dimension {eff_dim} is LARGE — cascaded flips remain powerful.")
    elif eff_dim >= 5:
        print(f"Effective dimension {eff_dim} is MODERATE — some cascading works.")
    else:
        print(f"Effective dimension {eff_dim} is SMALL — cascading breaks down quickly.")

    print()


if __name__ == "__main__":
    main()
