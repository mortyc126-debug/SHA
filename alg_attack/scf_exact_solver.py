#!/usr/bin/env python3
"""
SCF: EXACT NONLINEAR SOLVER — наш собственный инструмент.

Newton не работает потому что аппроксимирует f(x+δ) ≈ f(x) + J·δ.
Но мы ЗНАЕМ f ТОЧНО. Мы можем вычислить ОБРАТНУЮ функцию раунда.

SHA round:
  T1 = h + Σ1(e) + Ch(e,f,g) + K[r] + W[r]
  T2 = Σ0(a) + Maj(a,b,c)
  a_new = T1 + T2
  e_new = d + T1

ОБРАТНЫЙ РАУНД (при известном state):
  T1 = e_new - d        (вычислимо!)
  W[r] = T1 - h - Σ1(e) - Ch(e,f,g) - K[r]  (вычислимо!)

Для ПАРЫ (W1, W2) и целевого δe[r+1]=0:
  e_new_1 = d_1 + T1_1
  e_new_2 = d_2 + T1_2
  δe_new = 0 → d_1 + T1_1 = d_2 + T1_2 (mod 2^32)
  → T1_2 = T1_1 + d_1 - d_2
  → W2[r] = T1_2 - h_2 - Σ1(e_2) - Ch(e_2,f_2,g_2) - K[r]

Это ТОЧНАЯ формула! Не приближение!

Проблема: W2[r] определяется schedule от W2[0..15].
Но мы ЗНАЕМ W2[0..15] → W2[r] фиксировано → формула даёт ТРЕБУЕМОЕ W2[r].
Если требуемое ≠ фактическое → нужно менять W2[0..15].

НОВЫЙ ПОДХОД: Wang chain через gap, используя ОБРАТНЫЙ РАУНД.
"""
import os, sys

MASK32 = 0xFFFFFFFF
K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK32
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def add32(*args):
    s=0
    for x in args: s=(s+x)&MASK32
    return s
def sub32(a,b): return (a-b)&MASK32
def hw(x): return bin(x).count('1')

def expand_real(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def sha_all_states(W16):
    W=expand_real(W16); s=list(IV); states=[list(s)]
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        states.append(list(s))
    return states, W

def sha_compress_full(W16):
    st, _ = sha_all_states(W16)
    return [add32(IV[i], st[64][i]) for i in range(8)]


def required_W_for_de_zero(state1_r, state2_r, state1_r1, r):
    """Given states at round r for both messages, and target δe[r+1]=0,
    compute the REQUIRED W2[r].

    state1_r1 = state of message 1 at round r+1 (already computed).
    We want state2 at r+1 to have e2[r+1] = e1[r+1].

    e1[r+1] = d1[r] + T1_1[r]
    We want e2[r+1] = e1[r+1]
    → d2[r] + T1_2[r] = d1[r] + T1_1[r]
    → T1_2 = T1_1 + d1 - d2 = T1_1 + δd (additive diff of d register)

    T1_2 = h2 + Σ1(e2) + Ch(e2,f2,g2) + K[r] + W2[r]
    → W2[r] = T1_2 - h2 - Σ1(e2) - Ch(e2,f2,g2) - K[r]
    """
    a1,b1,c1,d1,e1,f1,g1,h1 = state1_r
    a2,b2,c2,d2,e2,f2,g2,h2 = state2_r

    # T1_1 from message 1
    T1_1 = sub32(state1_r1[4], d1)  # e1_new = d1 + T1_1 → T1_1 = e1_new - d1

    # Required T1_2
    T1_2_required = add32(T1_1, sub32(d1, d2))

    # Required W2[r]
    W2_r = sub32(sub32(sub32(sub32(T1_2_required, h2), Sig1(e2)), Ch(e2,f2,g2)), K[r])

    return W2_r


# ============================================================
# EXP 1: For each gap round, compute required W2[r] for δe=0
# ============================================================
def exp1_required_W(N):
    print("="*70)
    print("EXP 1: REQUIRED W2[r] FOR δe[r+1]=0 AT EACH GAP ROUND")
    print("  Compute exact W2[r] needed, compare with actual schedule W2[r]")
    print("="*70)

    for trial in range(min(N, 5)):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        states1, Wexp1 = sha_all_states(W1)
        states2, Wexp2 = sha_all_states(W2)

        print(f"\n  Trial {trial}:")
        print(f"  {'r':>3} | {'W2_required':>12} {'W2_actual':>12} {'diff':>12} {'HW(diff)':>9}")
        print("  " + "-"*55)

        for r in range(17, 52):
            W2_req = required_W_for_de_zero(states1[r], states2[r], states1[r+1], r)
            W2_act = Wexp2[r]
            diff = W2_req ^ W2_act

            de_actual = states1[r+1][4] ^ states2[r+1][4]
            marker = " ★" if diff == 0 else ""

            if r < 22 or r > 47 or diff == 0:
                print(f"  {r:3d} | 0x{W2_req:08x} 0x{W2_act:08x} 0x{diff:08x} {hw(diff):9d}{marker}")

        # Summary: if we COULD set W2[r] = W2_required for all r,
        # we'd have δe=0 through the entire gap.
        # The COST is the diff between required and actual schedule.
        total_diff = 0
        for r in range(17, 52):
            W2_req = required_W_for_de_zero(states1[r], states2[r], states1[r+1], r)
            total_diff += hw(W2_req ^ Wexp2[r])

        print(f"  Total HW(required - actual): {total_diff} bits")
        print(f"  Average per round: {total_diff/35:.1f} bits")


# ============================================================
# EXP 2: The MISMATCH — required W vs schedule W
# ============================================================
def exp2_mismatch_structure(N):
    print("\n" + "="*70)
    print("EXP 2: MISMATCH STRUCTURE")
    print("  mismatch[r] = W2_required[r] ⊕ W2_actual[r]")
    print("  If mismatch is STRUCTURED → can be absorbed by kernel")
    print("="*70)

    mismatches = {r: [] for r in range(17, 52)}

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= (trial + 1) | 1

        states1, Wexp1 = sha_all_states(W1)
        states2, Wexp2 = sha_all_states(W2)

        for r in range(17, 52):
            W2_req = required_W_for_de_zero(states1[r], states2[r], states1[r+1], r)
            diff = W2_req ^ Wexp2[r]
            mismatches[r].append(diff)

    # HW distribution of mismatch
    print(f"\n  HW(mismatch) per round (N={N}):")
    print(f"  {'r':>3} | avg_HW | min_HW | max_HW | P(HW≤5)")
    print("  " + "-"*50)
    for r in range(17, 52, 2):
        hws = [hw(m) for m in mismatches[r]]
        avg = sum(hws)/len(hws)
        mn = min(hws)
        mx = max(hws)
        p_low = sum(1 for h in hws if h <= 5) / len(hws)
        marker = " ★" if mn <= 3 else ""
        print(f"  {r:3d} | {avg:6.1f} | {mn:6d} | {mx:6d} | {p_low:.4f}{marker}")


# ============================================================
# EXP 3: Can the kernel ABSORB the mismatch?
# ============================================================
def exp3_kernel_absorb(N):
    print("\n" + "="*70)
    print("EXP 3: CAN SCHEDULE KERNEL ABSORB THE MISMATCH?")
    print("  Required: δW2_schedule = mismatch pattern")
    print("  Available: kernel produces specific δW patterns")
    print("  Match = mismatch lies in image of kernel")
    print("="*70)

    # The mismatch at round r is a specific 32-bit value.
    # The schedule kernel at R=52 gives δW_xor patterns.
    # Question: does the kernel's IMAGE at round r contain the mismatch?

    # For each basis vector, compute what δW_real[r] it produces
    from scf_newton_scale import extract_kernel_basis, bits_to_words
    basis = extract_kernel_basis(52)
    n_basis = len(basis)
    basis_words = [bits_to_words(v) for v in basis]

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    W2 = list(W1); W2[0] ^= 0x80000000

    states1, Wexp1 = sha_all_states(W1)
    states2, Wexp2 = sha_all_states(W2)

    # For one gap round, check if mismatch is achievable
    for target_r in [20, 30, 40]:
        W2_req = required_W_for_de_zero(states1[target_r], states2[target_r],
                                         states1[target_r+1], target_r)
        mismatch = W2_req ^ Wexp2[target_r]

        # What δW_real[target_r] does each basis vector produce?
        basis_effects = []
        for i in range(n_basis):
            W2_mod = [(W2[j] ^ basis_words[i][j]) for j in range(16)]
            Wexp_mod = expand_real(W2_mod)
            effect = Wexp_mod[target_r] ^ Wexp2[target_r]
            basis_effects.append(effect)

        # Check: is mismatch in the span of basis_effects (over Z)?
        # Over GF(2): check if mismatch XOR is in span
        # Build matrix: each column = basis_effect (32 bits)
        # Check if mismatch is in column space

        from copy import deepcopy
        mat = [[((basis_effects[i]>>bit)&1) for i in range(n_basis)] for bit in range(32)]
        target_bits = [(mismatch>>bit)&1 for bit in range(32)]

        # Augmented system
        aug = [mat[b] + [target_bits[b]] for b in range(32)]
        # RREF
        m = deepcopy(aug); pivots=[]; ri=0
        for col in range(n_basis):
            pv=-1
            for r in range(ri,32):
                if m[r][col]: pv=r; break
            if pv==-1: continue
            m[ri],m[pv]=m[pv],m[ri]
            for r in range(32):
                if r!=ri and m[r][col]:
                    for c in range(n_basis+1): m[r][c]^=m[ri][c]
            pivots.append(col); ri+=1

        rank = len(pivots)
        consistent = all(m[r][n_basis]==0 for r in range(rank, 32))

        if consistent:
            print(f"  r={target_r}: mismatch IN kernel image (GF2)! rank={rank} ★")
            # Solve
            sol = [0]*n_basis
            for i,pc in enumerate(pivots): sol[pc] = m[i][n_basis]
            # Build correction
            bits = [0]*512
            for i in range(n_basis):
                if sol[i]:
                    for j in range(512): bits[j] ^= basis[i][j]
            dW = bits_to_words(bits)
            if all(w==0 for w in dW):
                print(f"    BUT: correction is trivial (δW=0)")
            else:
                # Apply and check
                W2_corrected = [(W2[j]^dW[j]) for j in range(16)]
                sc, Wexpc = sha_all_states(W2_corrected)
                de = states1[target_r+1][4] ^ sc[target_r+1][4]
                print(f"    Correction applied: HW(δe)={hw(de)}")
                if de == 0:
                    print(f"    ★★★ δe[{target_r+1}] = 0 EXACTLY!")
                    # Check hash
                    H1 = [add32(IV[i],states1[64][i]) for i in range(8)]
                    H2 = sha_compress_full(W2_corrected)
                    hhw = sum(hw(H1[i]^H2[i]) for i in range(8))
                    print(f"    Hash diff: HW = {hhw}")
        else:
            print(f"  r={target_r}: mismatch NOT in kernel image. rank={rank}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100

    exp1_required_W(3)
    exp2_mismatch_structure(min(N, 200))
    exp3_kernel_absorb(1)
