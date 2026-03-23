#!/usr/bin/env python3
"""CRAZY-29: Round-split MITM — split computation, not bits."""
import random, time

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
def sha_round(st,w,k):
    a,b,c,d,e,f,g,h=st
    T1=add32(h,Sig1(e),Ch(e,f,g),k,w); T2=add32(Sig0(a),Maj(a,b,c))
    return [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
def expand_W(W16):
    W=list(W16)
    for i in range(16,20): W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def wang_twin(msg):
    W=list(msg); Wp=list(W); Wp[0]^=0x80000000
    s=list(H0); sp=list(H0)
    s=sha_round(s,W[0],K[0]); sp=sha_round(sp,Wp[0],K[0])
    for t in range(1,16):
        a,b,c,d,e,f,g,h=s; a2,b2,c2,d2,e2,f2,g2,h2=sp
        tp=add32(h,Sig1(e),Ch(e,f,g),K[t]); tp2=add32(h2,Sig1(e2),Ch(e2,f2,g2),K[t])
        target=add32(d,tp,W[t]); Wp[t]=(target-d2-tp2)&MASK
        s=sha_round(s,W[t],K[t]); sp=sha_round(sp,Wp[t],K[t])
    return Wp

def get_de17(msg):
    Wp=wang_twin(msg); We=expand_W(msg); Wpe=expand_W(Wp)
    s1=list(H0); s2=list(H0)
    for t in range(17): s1=sha_round(s1,We[t],K[t]); s2=sha_round(s2,Wpe[t],K[t])
    return s1[4]^s2[4]

def main():
    random.seed(0xC029)
    t0 = time.time()
    base_msg = [random.getrandbits(32) for _ in range(16)]

    # ── 1. W[14] enters at round 14 AND via W[16] = sig1(W[14]) + W[9] + sig0(W[1]) + W[0]
    # These are TWO paths. Measure how much each contributes.
    print("="*72)
    print("CRAZY-29: Round-split MITM")
    print("="*72)

    # ── 2. Decomposition test: f(direct_effect) vs g(schedule_effect)
    # Direct: W[14] enters T1 at round 14
    # Schedule: W[14] → sig1(W[14]) → W[16] enters T1 at round 16
    # Can we decompose De17 ≈ f(T1_14) XOR g(T1_16)?

    print("\nPART 1: Two-path decomposition")
    print("-"*72)

    # For each W[14], compute:
    # path_A = state change due to W[14] at round 14 (direct)
    # path_B = state change due to sig1(W[14]) at round 16 (schedule)
    # XOR decomposition residual

    N = 50000
    base_de17 = get_de17(base_msg)
    We_base = expand_W(base_msg)

    # Measure: correlation between De17 and sig1(W[14])
    # If De17 ≈ h(sig1(W14)), then sig1 is the bottleneck
    sig1_to_de17 = {}
    collisions = 0
    for _ in range(N):
        w14 = random.getrandbits(32)
        msg = list(base_msg); msg[14] = w14
        de17 = get_de17(msg)
        s1v = sig1(w14)
        if s1v in sig1_to_de17:
            if sig1_to_de17[s1v] != de17:
                collisions += 1
        else:
            sig1_to_de17[s1v] = de17

    print(f"  sig1(W14) unique values: {len(sig1_to_de17)}/{N}")
    print(f"  sig1 collisions with different De17: {collisions}")
    print(f"  → De17 depends on MORE than just sig1(W14)" if collisions > 0 else "  → De17 = f(sig1(W14))!")

    # ── 3. Round-15 intermediate state analysis
    print(f"\nPART 2: Intermediate state at round 15")
    print("-"*72)

    # Compute (e[15], a[15]) for many W[14] values
    # These are the two registers affected by W[14] at round 14
    # e[15] = d[14] + T1[14], a[15] = T1[14] + T2[14]
    # T1[14] depends on W[14] linearly (one addition)

    Wp_base = wang_twin(base_msg)
    We_b = expand_W(base_msg); Wpe_b = expand_W(Wp_base)

    # Run rounds 0-13 (fixed, independent of W[14])
    s1_r13 = list(H0); s2_r13 = list(H0)
    for t in range(14):
        s1_r13 = sha_round(s1_r13, We_b[t], K[t])
        s2_r13 = sha_round(s2_r13, Wpe_b[t], K[t])

    print(f"  State after round 13 (fixed): computed")
    print(f"  Registers: a={s1_r13[0]:08x} e={s1_r13[4]:08x}")

    # For round 14: T1 = h + Sig1(e) + Ch(e,f,g) + K[14] + W[14]
    # Only W[14] varies! So T1 = CONST + W[14] (mod 2^32)
    a,b,c,d,e,f,g,h = s1_r13
    const_T1 = add32(h, Sig1(e), Ch(e,f,g), K[14])
    const_T2 = add32(Sig0(a), Maj(a,b,c))

    print(f"  T1[14] = {const_T1:08x} + W[14]  (single addition!)")
    print(f"  T2[14] = {const_T2:08x}  (constant, no W[14] dependence)")
    print(f"  e[15] = d[14] + T1 = {d:08x} + {const_T1:08x} + W[14]")
    print(f"  a[15] = T1 + T2 = ({const_T1:08x} + W[14]) + {const_T2:08x}")

    # ── 4. How does W[14] enter round 16 via W[16]?
    print(f"\nPART 3: W[16] schedule path")
    print("-"*72)

    # W[16] = sig1(W[14]) + W[9] + sig0(W[1]) + W[0]
    const_W16 = add32(We_b[9], sig0(We_b[1]), We_b[0])
    print(f"  W[16] = sig1(W[14]) + {const_W16:08x}")
    print(f"  sig1 is NONLINEAR (rotations + shift + XOR)")

    # Key insight: W[14] enters TWICE:
    # Path A: round 14 as W[14] (LINEAR addition to T1)
    # Path B: round 16 as sig1(W[14]) (NONLINEAR transform + addition)
    # These two paths COUPLE at the intermediate state

    # ── 5. Decomposition residual
    print(f"\nPART 4: Decomposition test")
    print("-"*72)

    # If we could separate: De17 = f(W14_direct) XOR g(sig1(W14))
    # Then MITM: precompute f table, invert g
    # Test: fix sig1(W14) [by fixing W14], vary only the direct path somehow
    # Problem: W14 determines BOTH paths. No separation.

    # But: sig1 is 2-to-1? No, it's a bijection on 32 bits (invertible via linear algebra over GF(2)... wait, no: sig1 has the >>10 shift which loses bits)
    # sig1(x) = rotr(x,17) XOR rotr(x,19) XOR (x >> 10)
    # This is LINEAR over GF(2)! It's a 32×32 GF(2) matrix.

    # So sig1 is GF(2)-linear. What's its rank?
    sig1_cols = []
    for j in range(32):
        sig1_cols.append(sig1(1 << j))
    # Build matrix and compute rank
    rows = list(sig1_cols)  # each is a column of the matrix
    mat = list(rows)
    rank = 0
    for col in range(32):
        pivot = -1
        for r in range(rank, 32):
            if (mat[r] >> col) & 1: pivot = r; break
        if pivot == -1: continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        for r in range(32):
            if r != rank and (mat[r] >> col) & 1: mat[r] ^= mat[rank]
        rank += 1

    print(f"  sig1 GF(2) rank: {rank}/32")
    if rank == 32:
        print(f"  sig1 is a GF(2) bijection! Every sig1(W14) maps to unique W14.")
        print(f"  → No information loss through schedule path")
        print(f"  → W14 and sig1(W14) carry SAME 32 bits of info")
        print(f"  → CANNOT separate into independent paths for MITM")
    else:
        print(f"  sig1 has kernel dimension {32-rank}")
        print(f"  → {32-rank} bits lost → potential MITM opportunity!")

    # ── 6. The REAL question: is there ANY decomposition?
    print(f"\nPART 5: Additive decomposition test")
    print("-"*72)

    # De17(W14) = function of (W14, sig1(W14))
    # Since sig1 is bijective, this is just De17(W14) — one variable.
    # No MITM possible on a single 32-bit variable.

    # HOWEVER: what if we split at the ROUND level?
    # Define: F(W14) = state after round 15 (depends on W14 at round 14)
    #         G(W14) = W[16] = sig1(W14) + const (enters at round 16)
    # De17 = H(F(W14), G(W14)) where H combines round 15-16-17 effects

    # If F and G were independent, we could MITM.
    # They're NOT independent (both functions of W14).
    # But maybe the COUPLING through H is weak?

    # Test: fix F (by fixing e[15], a[15]) and vary G independently
    # This is impossible with a single W14, but we can MEASURE the coupling:
    # For pairs (W14_a, W14_b) where F(W14_a) ≈ F(W14_b) but G differs,
    # how different is De17?

    from collections import defaultdict

    # e[15] = d + const_T1 + W14 = (d + const_T1) + W14
    e15_const = add32(d, const_T1)
    # So e[15] = e15_const + W14. Two W14 values give same e[15] only if they're equal.
    # a[15] = const_T1 + W14 + const_T2. Same story.
    # Both are bijective in W14! So F is a bijection → F(W14_a) = F(W14_b) iff W14_a = W14_b.
    print(f"  F(W14) = (e15_const + W14, const_T1 + W14 + const_T2)")
    print(f"  F is BIJECTIVE in W14 (mod-add is bijective)")
    print(f"  G(W14) = sig1(W14) + const (sig1 is GF(2)-bijective)")
    print(f"  Both F and G are bijections of W14")
    print(f"  → F and G are in 1-1 correspondence")
    print(f"  → NO independent variation possible")
    print(f"  → Round-split MITM is IMPOSSIBLE for single-word barrier")

    # ── 7. What about MULTI-WORD?
    print(f"\nPART 6: Multi-word round-split possibility")
    print("-"*72)
    print(f"  W[14] and W[15] both enter rounds 14-17.")
    print(f"  W[14] enters round 14 (direct) and round 16 (via W[16]=sig1(W14)+...)")
    print(f"  W[15] enters round 15 (direct) and round 17 (via W[17]=sig1(W15)+...)")
    print(f"  Total: 64 bits of freedom → 4 entry points into rounds 14-17")
    print(f"  Potential MITM: split at round 15.5")
    print(f"    Forward: W[14] → state after round 15 (32-bit subspace in 256-dim)")
    print(f"    Backward: W[15] + W[17]=sig1(W15)+... → what round-15 state needed")

    # Test: vary W[14] and W[15] independently
    # Forward: for 2^16 W14 values, compute state after round 15
    # The state is determined by e[15] and a[15] (other regs shift from previous)
    # e[15] = d + const + W14, a[15] = const + W14 + const2
    # So the "state fingerprint" from W14 is just W14 itself (bijection)

    # Backward: rounds 16-17 use W[16]=sig1(W14)+const and W[17]=sig1(W15)+const
    # These depend on BOTH W14 and W15.

    # MITM split: choose W14 (forward, 2^32), choose W15 (backward, 2^32)
    # Match at round-15 boundary: needs e[15] from forward to feed into round 16
    # But round 16 uses W[16] which depends on W14 (forward variable!)
    # So backward computation ALSO needs W14 → not independent.

    print(f"\n  Round 16 uses W[16] = sig1(W[14]) + const")
    print(f"  → Backward (round 16-17) depends on W[14] (forward variable!)")
    print(f"  → Split at round 15 doesn't decouple W14 from backward path")
    print(f"  → Multi-word MITM at round boundary: IMPOSSIBLE")
    print(f"     (message schedule creates cross-round coupling)")

    # ── VERDICT
    print(f"\n{'='*72}")
    print(f"VERDICT: Round-split MITM")
    print(f"{'='*72}")
    print(f"  sig1 rank: {rank}/32 (bijection)")
    print(f"  F (round-14 path): bijective in W14")
    print(f"  G (round-16 path via schedule): bijective in W14")
    print(f"  F and G are 1-1 locked → no independent variation")
    print(f"  Multi-word split: blocked by message schedule coupling")
    print(f"")
    print(f"  **DEAD**: Round-split MITM impossible because:")
    print(f"  1. Single-word: both paths are bijections of same variable")
    print(f"  2. Multi-word: message schedule couples forward/backward")
    print(f"  3. This is BY DESIGN in SHA-256 — schedule prevents MITM")
    print(f"")
    print(f"  Runtime: {time.time()-t0:.1f}s")
    print(f"{'='*72}")

if __name__ == "__main__":
    main()
