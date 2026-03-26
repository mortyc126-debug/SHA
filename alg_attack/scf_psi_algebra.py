#!/usr/bin/env python3
"""
SCF: Ψ-ALGEBRA — пятый фреймворк. Carry как ОСНОВА, не враг.

Все предыдущие подходы: "carry мешает, как обойти?"
Новый подход: "carry = структура, как ИСПОЛЬЗОВАТЬ?"

Ψ(a,b) = (a+b) ⊕ (a⊕b) = carry << 1.

Ψ — НЕ случайный. Ψ ДЕТЕРМИНИРОВАН при известных a,b.
Мы знаем: Ψ[bit_i] = OR_{j<i}(a[j] AND b[j] AND AND_{k=j+1..i-1}(a[k] XOR b[k]))

Определение Ψ-алгебры:
  Элемент: тройка (x, y, ψ) ∈ Z^3_{2^32} где ψ = Ψ(x, y)
  Сложение: (x1,y1,ψ1) ⊕ (x2,y2,ψ2) = (x1⊕x2, y1⊕y2, Ψ(x1⊕x2, y1⊕y2))
  Carry-op: (x,y,ψ) → (x+y, 0, 0) = "resolve carry"

Ключ: в Ψ-алгебре, CARRY ВИДИМ. Мы можем его отслеживать,
предсказывать и компенсировать.

SHA в Ψ-алгебре: каждое сложение порождает Ψ.
Collision в Ψ-алгебре: Σ всех Ψ = 0 (carry-corrections cancel).
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

def psi(a, b):
    """Ψ(a,b) = (a+b) ⊕ (a⊕b) = carry correction."""
    return ((a+b) & MASK32) ^ (a ^ b)

def expand_real(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W


def sha_with_psi_trace(W16):
    """Run SHA and record EVERY Ψ produced at each addition."""
    W = expand_real(W16); s = list(IV)
    psi_trace = []  # List of (round, operation, psi_value)

    for r in range(64):
        a,b,c,d,e,f,g,h = s

        # T1 = h + Σ1(e) + Ch(e,f,g) + K[r] + W[r]
        # 4 additions → 4 Ψ values
        p1 = psi(h, Sig1(e)); t = (h + Sig1(e)) & MASK32
        p2 = psi(t, Ch(e,f,g)); t = (t + Ch(e,f,g)) & MASK32
        p3 = psi(t, K[r]); t = (t + K[r]) & MASK32
        p4 = psi(t, W[r]); T1 = (t + W[r]) & MASK32

        # T2 = Σ0(a) + Maj(a,b,c)
        p5 = psi(Sig0(a), Maj(a,b,c)); T2 = add32(Sig0(a), Maj(a,b,c))

        # a_new = T1 + T2
        p6 = psi(T1, T2); a_new = add32(T1, T2)

        # e_new = d + T1
        p7 = psi(d, T1); e_new = add32(d, T1)

        psi_trace.append({
            'round': r,
            'psi_T1': [p1, p2, p3, p4],
            'psi_T2': p5,
            'psi_a': p6,
            'psi_e': p7,
            'total_hw': hw(p1)+hw(p2)+hw(p3)+hw(p4)+hw(p5)+hw(p6)+hw(p7),
        })

        s = [a_new, a, b, c, e_new, e, f, g]

    H = [add32(IV[i], s[i]) for i in range(8)]
    # MD addition Ψ
    psi_md = [psi(IV[i], s[i]) for i in range(8)]

    return H, s, psi_trace, psi_md


# ============================================================
# EXP 1: Ψ-DIFFERENTIAL — δΨ between two messages
# ============================================================
def exp1_psi_differential(N):
    print("="*70)
    print("EXP 1: Ψ-DIFFERENTIAL — δΨ between two messages")
    print("  δΨ_total = XOR of all Ψ differences across 64 rounds")
    print("  If δΨ_total has STRUCTURE → collision condition on Ψ")
    print("="*70)

    # For collision: δH = δ(IV + state) = δstate (since IV fixed)
    # δstate = δstate_xor ⊕ δΨ_total
    # For collision: δstate = 0 → δstate_xor = δΨ_total
    # So: collision ⟺ XOR-differential = Ψ-differential

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        H1, s1, trace1, psi_md1 = sha_with_psi_trace(W1)
        H2, s2, trace2, psi_md2 = sha_with_psi_trace(W2)

        # Total δΨ across all rounds
        total_dpsi_hw = 0
        per_round_dpsi = []
        for r in range(64):
            round_dpsi = 0
            for i in range(4):
                round_dpsi += hw(trace1[r]['psi_T1'][i] ^ trace2[r]['psi_T1'][i])
            round_dpsi += hw(trace1[r]['psi_T2'] ^ trace2[r]['psi_T2'])
            round_dpsi += hw(trace1[r]['psi_a'] ^ trace2[r]['psi_a'])
            round_dpsi += hw(trace1[r]['psi_e'] ^ trace2[r]['psi_e'])
            total_dpsi_hw += round_dpsi
            per_round_dpsi.append(round_dpsi)

        # MD addition δΨ
        md_dpsi = sum(hw(psi_md1[i] ^ psi_md2[i]) for i in range(8))

        hash_hw = sum(hw(H1[i]^H2[i]) for i in range(8))

        if trial < 3:
            print(f"\n  Trial {trial}: HW(δH)={hash_hw}")
            print(f"    Total δΨ across rounds: {total_dpsi_hw} bits")
            print(f"    MD addition δΨ: {md_dpsi} bits")
            print(f"    δΨ per phase:")
            for phase, start, end in [("early(0-15)", 0, 16),
                                       ("gap(16-51)", 16, 52),
                                       ("tail(52-63)", 52, 64)]:
                phase_dpsi = sum(per_round_dpsi[start:end])
                print(f"      {phase}: {phase_dpsi} ({phase_dpsi/(end-start):.1f}/round)")


# ============================================================
# EXP 2: Ψ-CANCELLATION — can Ψ from different rounds cancel?
# ============================================================
def exp2_psi_cancellation(N):
    print("\n" + "="*70)
    print("EXP 2: Ψ-CANCELLATION")
    print("  Do Ψ corrections from different rounds CANCEL each other?")
    print("  Collision requires: total carry correction = 0 (in XOR)")
    print("="*70)

    # Track cumulative Ψ (XOR accumulation)
    for trial in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        _, _, trace, psi_md = sha_with_psi_trace(W)

        # XOR all Ψ values together (cumulative)
        cumulative_psi = [0]*8  # One per "channel"

        # Accumulate per-round Ψ
        for r in range(64):
            for p in trace[r]['psi_T1']:
                cumulative_psi[0] ^= p
            cumulative_psi[1] ^= trace[r]['psi_T2']
            cumulative_psi[2] ^= trace[r]['psi_a']
            cumulative_psi[3] ^= trace[r]['psi_e']

        cum_hw = sum(hw(c) for c in cumulative_psi[:4])

        if trial < 3:
            print(f"\n  Trial {trial}: cumulative Ψ XOR:")
            labels = ['ΣΨ_T1', 'ΣΨ_T2', 'ΣΨ_a', 'ΣΨ_e']
            for i in range(4):
                print(f"    {labels[i]}: HW={hw(cumulative_psi[i])}")
            print(f"    Total cumulative HW: {cum_hw}")
            print(f"    Expected random: ~{4*16}={4*16}")

            if cum_hw < 4*16 - 10:
                print(f"    ★ Below random — partial cancellation!")


# ============================================================
# EXP 3: Ψ-COLLISION CONDITION — what Ψ pattern = collision?
# ============================================================
def exp3_psi_collision_condition(N):
    print("\n" + "="*70)
    print("EXP 3: Ψ-COLLISION CONDITION")
    print("  SHA(W) = SHA_xor(W) ⊕ Ψ_total(W)")
    print("  Collision: SHA(W1) = SHA(W2)")
    print("  ⟺ SHA_xor(W1) ⊕ Ψ_total(W1) = SHA_xor(W2) ⊕ Ψ_total(W2)")
    print("  ⟺ δSHA_xor = δΨ_total")
    print("  This is the MATCHING CONDITION in Ψ-algebra!")
    print("="*70)

    # Compute SHA_xor and Ψ_total separately
    def sha_xor_only(W16):
        W = list(W16)
        for i in range(16,64):
            W.append(sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16])
        s = list(IV)
        for r in range(64):
            a,b,c,d,e,f,g,h=s
            T1=h^Sig1(e)^Ch(e,f,g)^K[r]^W[r]
            T2=Sig0(a)^Maj(a,b,c)
            s=[T1^T2,a,b,c,d^T1,e,f,g]
        return [IV[i]^s[i] for i in range(8)]

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        H1_real = sha_with_psi_trace(W1)[0]
        H2_real = sha_with_psi_trace(W2)[0]

        H1_xor = sha_xor_only(W1)
        H2_xor = sha_xor_only(W2)

        # Ψ_total = SHA_real ⊕ SHA_xor
        psi1 = [H1_real[i] ^ H1_xor[i] for i in range(8)]
        psi2 = [H2_real[i] ^ H2_xor[i] for i in range(8)]

        # δSHA_xor
        d_xor = [H1_xor[i] ^ H2_xor[i] for i in range(8)]
        # δΨ_total
        d_psi = [psi1[i] ^ psi2[i] for i in range(8)]
        # δH_real
        d_real = [H1_real[i] ^ H2_real[i] for i in range(8)]

        # Verify: d_real = d_xor ⊕ d_psi
        verify = all(d_real[i] == (d_xor[i] ^ d_psi[i]) for i in range(8))

        hw_real = sum(hw(d) for d in d_real)
        hw_xor = sum(hw(d) for d in d_xor)
        hw_psi = sum(hw(d) for d in d_psi)

        # Matching: HW(d_xor ⊕ d_psi) = HW(d_real). For collision: d_xor = d_psi.
        matching_hw = sum(hw(d_xor[i] ^ d_psi[i]) for i in range(8))

        if trial < 5:
            print(f"\n  Trial {trial}: verify={verify}")
            print(f"    HW(δH_real)={hw_real}, HW(δH_xor)={hw_xor}, HW(δΨ)={hw_psi}")
            print(f"    HW(δH_xor ⊕ δΨ) = {matching_hw} (should = HW(δH_real) = {hw_real})")
            print(f"    For collision: need HW(δH_xor ⊕ δΨ) = 0")
            print(f"    i.e., δH_xor EXACTLY EQUALS δΨ (256-bit match)")


# ============================================================
# EXP 4: Can we PREDICT Ψ from message difference?
# ============================================================
def exp4_psi_prediction(N):
    print("\n" + "="*70)
    print("EXP 4: Ψ-PREDICTION — can we predict δΨ from δW?")
    print("  If yes → set δW to make δΨ = δH_xor → collision!")
    print("="*70)

    # δΨ depends on BOTH messages (not just δW).
    # But: Ψ(a,b) depends on a AND b (not just a⊕b).
    # So δΨ is state-dependent — same barrier as carry.

    # HOWEVER: in Ψ-algebra, we TRACK Ψ explicitly.
    # The Ψ-Jacobian: ∂(Ψ_total)/∂(W[word])
    # Can we compute it?

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    def psi_total(W):
        H_real = sha_with_psi_trace(W)[0]
        H_xor = [0]*8
        Wx = list(W)
        for i in range(16,64):
            Wx.append(sig1(Wx[i-2])^Wx[i-7]^sig0(Wx[i-15])^Wx[i-16])
        s=list(IV)
        for r in range(64):
            a,b,c,d,e,f,g,h=s
            T1=h^Sig1(e)^Ch(e,f,g)^K[r]^Wx[r]
            T2=Sig0(a)^Maj(a,b,c)
            s=[T1^T2,a,b,c,d^T1,e,f,g]
        H_xor = [IV[i]^s[i] for i in range(8)]
        return [H_real[i]^H_xor[i] for i in range(8)]

    psi_base = psi_total(W1)

    # Jacobian: for each word, how does Ψ_total change?
    print(f"  Ψ-Jacobian: ∂Ψ_total/∂W[word]")
    for word in range(16):
        W_test = list(W1); W_test[word] = add32(W1[word], 1)
        psi_test = psi_total(W_test)
        dpsi = sum(hw(psi_base[i] ^ psi_test[i]) for i in range(8))
        print(f"    W[{word:2d}]: HW(δΨ_total) = {dpsi}")

    avg_dpsi = sum(
        sum(hw(psi_total(list(W1[:w])+[add32(W1[w],1)]+W1[w+1:])[i]^psi_base[i]) for i in range(8))
        for w in range(16)) / 16
    print(f"\n  Average: {avg_dpsi:.1f}")
    print(f"  Random Ψ change: ~128 bits")
    print(f"  If avg << 128 → Ψ is predictable → exploitable!")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    exp1_psi_differential(3)
    exp2_psi_cancellation(3)
    exp3_psi_collision_condition(5)
    exp4_psi_prediction(1)
