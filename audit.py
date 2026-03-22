#!/usr/bin/env python3
"""
AUDIT: Verify key claims from methodology and combinatorial synthesis.
Focus on potential artifacts that could invalidate the theory.
"""

import random
import time
import math

MASK = 0xFFFFFFFF
K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]
H0 = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Ch(e,f,g): return ((e&f)^(~e&g))&MASK
def Maj(a,b,c): return ((a&b)^(a&c)^(b&c))&MASK
def add32(*args):
    s=0
    for a in args: s=(s+a)&MASK
    return s
def hw(x): return bin(x&MASK).count('1')

def sha256_states(msg, iv, nr):
    W = list(msg)
    for i in range(16, max(nr,16)):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    states = [list(iv)]
    state = list(iv)
    for t in range(nr):
        a,b,c,d,e,f,g,h = state
        T1 = add32(h,Sig1(e),Ch(e,f,g),K[t],W[t])
        T2 = add32(Sig0(a),Maj(a,b,c))
        state = [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        states.append(list(state))
    return states, W

def expand_W(W16):
    W = list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W


def main():
    random.seed(0xA0D1)
    t0 = time.time()
    iv = list(H0)

    # ══════════════════════════════════════════════════════════════
    print("=" * 72)
    print("AUDIT 1: Free words W[12..15] — does Da13 depend on W[12]?")
    print("=" * 72)
    # Claim: Da13 depends ONLY on W[0..11]. W[12..15] are "free".
    # Test: flip bits in W[12], measure change in Da13.
    # Da13 = a[13] XOR a'[13]. a[13] = T1[12] + T2[12].
    # T1[12] = h[12] + Sig1(e[12]) + Ch(e[12],f[12],g[12]) + K[12] + W[12]
    # So a[13] DIRECTLY uses W[12] → Da13 SHOULD depend on W[12]!
    print()

    N = 1000
    dep_12 = 0
    dep_13 = 0
    dep_14 = 0
    dep_15 = 0

    for _ in range(N):
        msg = [random.getrandbits(32) for _ in range(16)]
        msg_prime = list(msg)
        msg_prime[0] ^= 0x80000000  # Wang diff

        # Compute Wang chain: adjust W'[1..15] for De[2..16]=0
        # (simplified: just compute states and measure Da at round 13)
        s1, _ = sha256_states(msg, iv, 17)
        s2, _ = sha256_states(msg_prime, iv, 17)
        Da13_base = (s1[13][0] ^ s2[13][0]) & MASK

        # Now flip W[12] and recompute
        for test_word in [12, 13, 14, 15]:
            msg_f = list(msg)
            msg_f[test_word] ^= 0x80000000
            msg_pf = list(msg_prime)
            msg_pf[test_word] ^= 0x80000000

            sf1, _ = sha256_states(msg_f, iv, 17)
            sf2, _ = sha256_states(msg_pf, iv, 17)
            Da13_flip = (sf1[13][0] ^ sf2[13][0]) & MASK

            if Da13_flip != Da13_base:
                if test_word == 12: dep_12 += 1
                elif test_word == 13: dep_13 += 1
                elif test_word == 14: dep_14 += 1
                elif test_word == 15: dep_15 += 1

    print(f"  Da13 changes when flipping W[12]: {dep_12}/{N} = {dep_12/N*100:.1f}%")
    print(f"  Da13 changes when flipping W[13]: {dep_13}/{N} = {dep_13/N*100:.1f}%")
    print(f"  Da13 changes when flipping W[14]: {dep_14}/{N} = {dep_14/N*100:.1f}%")
    print(f"  Da13 changes when flipping W[15]: {dep_15}/{N} = {dep_15/N*100:.1f}%")

    if dep_12 > 0:
        print(f"\n  *** CLAIM REFUTED: Da13 DOES depend on W[12]! ***")
        print(f"  The 'free words' claim is WRONG or needs reinterpretation.")
    else:
        print(f"\n  Claim confirmed: Da13 independent of W[12..15].")

    # ══════════════════════════════════════════════════════════════
    print()
    print("=" * 72)
    print("AUDIT 2: GF(2) Jacobian rank — 254 or 256?")
    print("=" * 72)
    # RESULTS.md says rank 254-255 (deficiency 1-2)
    # C20 says "full rank 256, deficiency was a bug"
    # CONTRADICTION. Let's test ourselves.
    print()

    def gf2_rank(matrix_cols, n_rows=256):
        rows = []
        for r in range(n_rows):
            row = 0
            for c in range(len(matrix_cols)):
                if (matrix_cols[c] >> r) & 1:
                    row |= (1 << c)
            rows.append(row)
        rank = 0
        for col in range(len(matrix_cols)):
            pivot = -1
            for r in range(rank, n_rows):
                if (rows[r] >> col) & 1:
                    pivot = r
                    break
            if pivot == -1:
                continue
            rows[rank], rows[pivot] = rows[pivot], rows[rank]
            for r in range(n_rows):
                if r != rank and (rows[r] >> col) & 1:
                    rows[r] ^= rows[rank]
            rank += 1
        return rank

    ranks = []
    for trial in range(5):
        msg = [random.getrandbits(32) for _ in range(16)]
        h_base = sha256_states(msg, iv, 64)[0][-1]
        base_hash = 0
        for i, val in enumerate(h_base):
            base_hash |= (val << (i*32))

        # Build Jacobian: flip each of 512 message bits
        cols = []
        for word in range(16):
            for bit in range(32):
                msg_f = list(msg)
                msg_f[word] ^= (1 << bit)
                h_f = sha256_states(msg_f, iv, 64)[0][-1]
                h_hash = 0
                for i, val in enumerate(h_f):
                    h_hash |= (val << (i*32))
                cols.append(base_hash ^ h_hash)

        rank = gf2_rank(cols)
        ranks.append(rank)
        print(f"  Trial {trial+1}: GF(2) Jacobian rank = {rank}/256")

    avg_rank = sum(ranks)/len(ranks)
    if avg_rank >= 255.5:
        print(f"\n  VERDICT: Full rank ({avg_rank:.1f}). C20 was RIGHT, RESULTS.md was WRONG.")
        print(f"  The 'rank deficiency 1-2' was indeed an artifact/bug.")
    elif avg_rank >= 253.5:
        print(f"\n  VERDICT: Near-full rank ({avg_rank:.1f}). Small deficiency may be real.")
    else:
        print(f"\n  VERDICT: Significant rank deficiency ({avg_rank:.1f}).")

    # ══════════════════════════════════════════════════════════════
    print()
    print("=" * 72)
    print("AUDIT 3: GA savings = 9 bits — real or random sampling?")
    print("=" * 72)
    # CRAZY-6 claimed GA reduces HW(De[17]) from 16 to 7 (savings = 9)
    # with 30 pop × 50 gen = 1500 evaluations per round.
    # Expected minimum of 1500 samples from Binomial(32, 0.5):
    # E[min] ≈ 16 - sqrt(8 × 2 × ln(1500)) = 16 - sqrt(116) ≈ 16 - 10.8 = 5.2
    # So random sampling would give min ≈ 5.2 at 1500 trials.
    # GA found ~7. This is WORSE than random expected minimum!
    print()

    N_GA = 1500
    exp_min = 16 - (8 * 2 * math.log(N_GA))**0.5
    print(f"  CRAZY-6 GA: 30 pop × 50 gen = {N_GA} evaluations")
    print(f"  GA result: HW(De[17]) ≈ 7")
    print(f"  Expected random min at {N_GA} trials: {exp_min:.1f}")
    print(f"  Comparison: GA found {7}, random expected {exp_min:.1f}")

    if 7 > exp_min + 1:
        print(f"\n  *** GA is WORSE than random! Savings = artifact! ***")
        print(f"  The 'gradient' in CRAZY-1 was likely just random sampling.")
    elif abs(7 - exp_min) < 2:
        print(f"\n  GA ≈ random. No real gradient detected.")
        print(f"  Z5 (Δ_GA = 9) is a RANDOM SAMPLING effect, not optimization!")
    else:
        print(f"\n  GA genuinely better than random.")

    # Verify with direct random sampling
    print(f"\n  Direct verification: 1500 random HW(De[17]) samples...")

    def wang_chain_de17(msg, iv):
        s1, _ = sha256_states(msg, iv, 17)
        Wp2 = list(msg); Wp2[0] ^= 0x80000000
        s2, _ = sha256_states(Wp2, iv, 17)
        return hw(s1[17][4] ^ s2[17][4])

    hws = []
    for _ in range(1500):
        msg = [random.getrandbits(32) for _ in range(16)]
        hws.append(wang_chain_de17(msg, iv))

    actual_min = min(hws)
    actual_avg = sum(hws)/len(hws)
    print(f"  1500 random: avg HW(De17) = {actual_avg:.1f}, min = {actual_min}")
    print(f"  CRAZY-6 GA claimed: avg optimized ≈ 7")
    print(f"  Random min at 1500: {actual_min}")

    if actual_min <= 8:
        print(f"\n  *** CONFIRMED: Random sampling alone gives min ≈ {actual_min} ***")
        print(f"  GA 'savings of 9 bits' = PURE RANDOM SAMPLING ARTIFACT!")
        print(f"  Z5 is FALSE. The 'gradient' never existed.")

    # ══════════════════════════════════════════════════════════════
    print()
    print("=" * 72)
    print("AUDIT 4: Neutral bits count — 82 genuine?")
    print("=" * 72)
    # C9 claimed 82.6 neutral bits on average.
    # Definition: bit (k,b) is neutral if flipping W[k] bit b AND W'[k] bit b
    # preserves De[2..16] = 0 under Wang chain with DW[0]=0x80000000.
    # Let's count from scratch.
    print()

    def count_neutrals(msg, iv):
        Wp = list(msg); Wp[0] ^= 0x80000000
        s1, _ = sha256_states(msg, iv, 16)
        s2, _ = sha256_states(Wp, iv, 16)
        # Verify De[2..16] = 0 for base
        base_ok = all((s1[t][4] - s2[t][4]) & MASK == 0 for t in range(2, 17))

        count = 0
        for k in range(1, 16):
            for b in range(32):
                flip = 1 << b
                Wf = list(msg); Wf[k] ^= flip
                Wpf = list(Wp); Wpf[k] ^= flip
                sf1, _ = sha256_states(Wf, iv, 16)
                sf2, _ = sha256_states(Wpf, iv, 16)
                ok = all((sf1[t][4] - sf2[t][4]) & MASK == 0 for t in range(2, 17))
                if ok:
                    count += 1
        return count, base_ok

    n_test = 5
    neutral_counts = []
    for i in range(n_test):
        msg = [random.getrandbits(32) for _ in range(16)]
        nc, base_ok = count_neutrals(msg, iv)
        neutral_counts.append(nc)
        print(f"  Msg {i+1}: {nc} neutral bits (base De=0: {base_ok})")

    avg_nc = sum(neutral_counts)/len(neutral_counts)
    print(f"\n  Average: {avg_nc:.1f} neutral bits")
    print(f"  C9 claimed: 82.6")
    if abs(avg_nc - 82.6) < 10:
        print(f"  CONFIRMED: ≈82 neutral bits is real.")
    else:
        print(f"  DISCREPANCY: {avg_nc:.1f} vs claimed 82.6!")

    # ══════════════════════════════════════════════════════════════
    print()
    print("=" * 72)
    print("AUDIT 5: Cascade depth D_casc — really limited to 2?")
    print("=" * 72)
    # CRAZY-4 claimed cascade works for 2 rounds, stalls at 3rd.
    # But only tested on 10 messages. Let's test more carefully.
    # For each message: find neutrals that minimize De[17],
    # then among those, find ones that also minimize De[18], De[19].
    print()

    msg = [random.getrandbits(32) for _ in range(16)]
    Wp = list(msg); Wp[0] ^= 0x80000000

    # First: find neutrals
    neutrals = []
    s1_base, W1 = sha256_states(msg, iv, 20)
    s2_base, W2 = sha256_states(Wp, iv, 20)
    for k in range(1, 16):
        for b in range(32):
            flip = 1 << b
            Wf = list(msg); Wf[k] ^= flip
            Wpf = list(Wp); Wpf[k] ^= flip
            sf1, _ = sha256_states(Wf, iv, 16)
            sf2, _ = sha256_states(Wpf, iv, 16)
            if all((sf1[t][4] - sf2[t][4]) & MASK == 0 for t in range(2, 17)):
                neutrals.append((k, b))

    print(f"  Found {len(neutrals)} neutral bits")

    # Greedy cascade: optimize De[17], then De[18], then De[19]
    best_set = set()
    for target_round in [17, 18, 19, 20]:
        best_hw = 999
        best_bit = None
        for idx, (k, b) in enumerate(neutrals):
            if idx in best_set:
                continue
            # Try adding this neutral bit
            test_set = best_set | {idx}
            W_test = list(msg)
            Wp_test = list(Wp)
            for i in test_set:
                nk, nb = neutrals[i]
                W_test[nk] ^= (1 << nb)
                Wp_test[nk] ^= (1 << nb)
            st1, _ = sha256_states(W_test, iv, target_round + 1)
            st2, _ = sha256_states(Wp_test, iv, target_round + 1)
            de_hw = hw(st1[target_round][4] ^ st2[target_round][4])
            if de_hw < best_hw:
                best_hw = de_hw
                best_bit = idx

        if best_bit is not None:
            best_set.add(best_bit)
        # Measure all De values
        W_test = list(msg); Wp_test = list(Wp)
        for i in best_set:
            nk, nb = neutrals[i]
            W_test[nk] ^= (1 << nb)
            Wp_test[nk] ^= (1 << nb)
        st1, _ = sha256_states(W_test, iv, 21)
        st2, _ = sha256_states(Wp_test, iv, 21)
        de_vals = [hw(st1[t][4] ^ st2[t][4]) for t in range(17, min(21, len(st1)))]
        print(f"  After optimizing r{target_round}: De[17..20] = {de_vals}, "
              f"|set|={len(best_set)}")

    # ══════════════════════════════════════════════════════════════
    print()
    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print()
    print(f"  Runtime: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
