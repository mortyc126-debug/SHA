#!/usr/bin/env python3
"""
CRAZY-17: Coupled Barrier System — the REAL question after audit

AUDIT PROVED: free words W[12..15] affect ALL barriers De17-De20.
Old claim "independently controllable at O(2^34)" was WRONG.

NEW QUESTION: 128 bits freedom (4 words) for 128 bits constraints (4 barriers).
The system is EXACTLY DETERMINED. What is its actual cost?

Options:
A) 2^128 — if system acts like random function (no structure)
B) 2^64  — birthday on 128-bit space
C) < 2^64 — if algebraic structure of SHA-256 helps
D) << 2^64 — if carry-AND identity enables algebraic solving

This experiment measures the ACTUAL structure of the 4-barrier system.
"""

import random
import time

MASK = 0xFFFFFFFF
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc]
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

def sha_round(state, W_t, K_t):
    a,b,c,d,e,f,g,h = state
    T1 = add32(h,Sig1(e),Ch(e,f,g),K_t,W_t)
    T2 = add32(Sig0(a),Maj(a,b,c))
    return [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

def expand_W(W16):
    W = list(W16)
    for i in range(16, 21):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def wang_chain(msg, iv):
    W = list(msg); Wp = list(W); Wp[0] ^= 0x80000000
    s = list(iv); sp = list(iv)
    s = sha_round(s,W[0],K[0]); sp = sha_round(sp,Wp[0],K[0])
    for t in range(1,16):
        a,b,c,d,e,f,g,h = s; a2,b2,c2,d2,e2,f2,g2,h2 = sp
        tp = add32(h,Sig1(e),Ch(e,f,g),K[t])
        tp2 = add32(h2,Sig1(e2),Ch(e2,f2,g2),K[t])
        target = add32(d,tp,W[t])
        Wp[t] = (target - d2 - tp2) & MASK
        s = sha_round(s,W[t],K[t]); sp = sha_round(sp,Wp[t],K[t])
    return Wp

def get_de(msg, iv, nr):
    """Get De at rounds 17..nr given Wang chain."""
    Wp = wang_chain(msg, iv)
    W_exp = expand_W(msg); Wp_exp = expand_W(Wp)
    s1 = list(iv); s2 = list(iv)
    for t in range(nr):
        s1 = sha_round(s1, W_exp[t], K[t])
        s2 = sha_round(s2, Wp_exp[t], K[t])
    return s1[4] ^ s2[4]  # De at round nr


def main():
    print("=" * 72)
    print("CRAZY-17: Coupled Barrier System Analysis")
    print("=" * 72)

    random.seed(0xC017)
    t0 = time.time()
    iv = list(H0)

    # ── Phase 1: How coupled are the barriers? ────────────────────
    print("\nPhase 1: Coupling matrix — how each free word affects each barrier")
    print("-" * 72)

    # Fix W[0..11], vary W[12..15] one at a time, measure De17..De20
    base_msg = [random.getrandbits(32) for _ in range(16)]
    N = 200

    for free_word in [12, 13, 14, 15]:
        de_changes = {17: 0, 18: 0, 19: 0, 20: 0}
        for _ in range(N):
            msg1 = list(base_msg)
            msg2 = list(base_msg)
            msg2[free_word] = random.getrandbits(32)  # randomize ONE free word

            for r in [17, 18, 19, 20]:
                d1 = get_de(msg1, iv, r)
                d2 = get_de(msg2, iv, r)
                if d1 != d2:
                    de_changes[r] += 1

        pcts = {r: de_changes[r]/N*100 for r in [17,18,19,20]}
        print(f"  W[{free_word}]: De17={pcts[17]:.0f}% De18={pcts[18]:.0f}% "
              f"De19={pcts[19]:.0f}% De20={pcts[20]:.0f}%")

    # ── Phase 2: Joint search — cost of De17=De18=0 ──────────────
    print(f"\nPhase 2: Joint birthday search De17=De18=0 via W[14]+W[15]")
    print("-" * 72)

    # Fix W[0..13], vary W[14] and W[15] randomly
    # Measure: how often De17=0 AND De18=0?
    base_msg = [random.getrandbits(32) for _ in range(16)]
    N_joint = 500_000
    de17_zeros = 0
    de18_zeros = 0
    both_zeros = 0
    best_sum = 999

    for trial in range(N_joint):
        msg = list(base_msg)
        msg[14] = random.getrandbits(32)
        msg[15] = random.getrandbits(32)

        d17 = get_de(msg, iv, 17)
        d18 = get_de(msg, iv, 18)

        hw17 = hw(d17)
        hw18 = hw(d18)
        s = hw17 + hw18

        if d17 == 0: de17_zeros += 1
        if d18 == 0: de18_zeros += 1
        if d17 == 0 and d18 == 0: both_zeros += 1
        if s < best_sum: best_sum = s

        if trial % 100000 == 0 and trial > 0:
            print(f"  {trial}: De17=0: {de17_zeros}, De18=0: {de18_zeros}, "
                  f"both=0: {both_zeros}, best_sum={best_sum}")

    print(f"\n  Results ({N_joint} trials, varying W[14]+W[15]):")
    print(f"  De17=0: {de17_zeros} (P={de17_zeros/N_joint:.6f}, expected 2^-32)")
    print(f"  De18=0: {de18_zeros} (P={de18_zeros/N_joint:.6f}, expected 2^-32)")
    print(f"  Both=0: {both_zeros}")
    print(f"  Best HW(De17)+HW(De18): {best_sum}")

    if de17_zeros > 0:
        print(f"\n  De17=0 found {de17_zeros} times! P = {de17_zeros/N_joint:.2e}")
        print(f"  Expected at 2^-32: {N_joint/2**32:.2e}")
        ratio = (de17_zeros/N_joint) / (1/2**32)
        print(f"  Ratio observed/expected: {ratio:.1f}x")

    # ── Phase 3: How many free word bits affect De17? ─────────────
    print(f"\nPhase 3: Sensitivity — how many bits of W[14] affect De17?")
    print("-" * 72)

    base_msg = [random.getrandbits(32) for _ in range(16)]
    base_de17 = get_de(base_msg, iv, 17)

    sensitive_bits = 0
    for bit in range(32):
        msg_f = list(base_msg)
        msg_f[14] ^= (1 << bit)
        de17_f = get_de(msg_f, iv, 17)
        if de17_f != base_de17:
            sensitive_bits += 1

    print(f"  W[14]: {sensitive_bits}/32 bits affect De17")
    print(f"  Expected if random function: 32/32")

    # Same for W[12], W[13], W[15]
    for w in [12, 13, 15]:
        sens = 0
        for bit in range(32):
            msg_f = list(base_msg)
            msg_f[w] ^= (1 << bit)
            de17_f = get_de(msg_f, iv, 17)
            if de17_f != base_de17: sens += 1
        print(f"  W[{w}]: {sens}/32 bits affect De17")

    # ── Phase 4: GF(2) Jacobian of the 4-barrier system ──────────
    print(f"\nPhase 4: GF(2) Jacobian of barrier system (128 × 128)")
    print("-" * 72)

    # Build Jacobian: 128 output bits (De17..De20 × 32) vs 128 input bits (W[12..15] × 32)
    base_msg = [random.getrandbits(32) for _ in range(16)]
    base_de = []
    for r in [17, 18, 19, 20]:
        base_de.append(get_de(base_msg, iv, r))

    # Base output as 128-bit value
    base_out = 0
    for i, d in enumerate(base_de):
        base_out |= (d << (i * 32))

    # Build 128 columns (one per input bit)
    cols = []
    for w_idx, w in enumerate([12, 13, 14, 15]):
        for bit in range(32):
            msg_f = list(base_msg)
            msg_f[w] ^= (1 << bit)
            de_f = []
            for r in [17, 18, 19, 20]:
                de_f.append(get_de(msg_f, iv, r))
            out_f = 0
            for i, d in enumerate(de_f):
                out_f |= (d << (i * 32))
            cols.append(base_out ^ out_f)

    # Compute rank
    rows = []
    for r in range(128):
        row = 0
        for c in range(128):
            if (cols[c] >> r) & 1:
                row |= (1 << c)
        rows.append(row)

    rank = 0
    for col in range(128):
        pivot = -1
        for r in range(rank, 128):
            if (rows[r] >> col) & 1:
                pivot = r; break
        if pivot == -1: continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for r in range(128):
            if r != rank and (rows[r] >> col) & 1:
                rows[r] ^= rows[rank]
        rank += 1

    print(f"  GF(2) Jacobian rank: {rank}/128")
    print(f"  Kernel dimension: {128 - rank}")

    if rank < 128:
        print(f"  RANK DEFICIENT! {128-rank} free directions in barrier space.")
        print(f"  These directions change free words WITHOUT changing De17..De20 (mod 2).")
    else:
        print(f"  Full rank. System is non-degenerate over GF(2).")
        print(f"  No 'free lunch' — every bit flip changes some barrier.")

    # ── Verdict ───────────────────────────────────────────────────
    print(f"\n{'='*72}")
    print("VERDICT")
    print(f"{'='*72}\n")

    if de17_zeros > 0:
        p_de17 = de17_zeros / N_joint
        cost_single = int(1 / p_de17)
        print(f"  De17=0 probability: {p_de17:.2e}")
        print(f"  Single barrier cost: ~{cost_single} ≈ 2^{cost_single.bit_length()}")
    else:
        print(f"  De17=0: 0/{N_joint} → P < {1/N_joint:.2e}")
        print(f"  Consistent with P = 2^-32 (need 2^32 trials)")

    print(f"  GF(2) Jacobian rank: {rank}/128")

    if rank == 128:
        print(f"\n  System is EXACTLY DETERMINED (128 eq, 128 var, full rank).")
        print(f"  Theoretical cost: 2^64 (birthday on 128-bit space).")
        print(f"  Can we do better? Only if carry structure helps.")
        print(f"  Carry-AND identity: De = linear + 2*(quadratic)")
        print(f"  If quadratic correction is small → linearize → solve → 2^0!")
        print(f"  If correction is large → stuck at 2^64.")
    else:
        deficiency = 128 - rank
        print(f"\n  RANK DEFICIENT by {deficiency}!")
        print(f"  Cost could be 2^{64 - deficiency//2} (birthday on reduced space).")

    print(f"\n  Runtime: {time.time()-t0:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
