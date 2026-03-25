#!/usr/bin/env python3
"""
AI-MATH: SHA как ПОТОК с точкой компенсации.

Видение: T1 — единственная точка входа message.
T1 идёт в ОБА выхода (a_new, e_new).
Если T1_pair1 = T1_pair2 → состояния СЛИВАЮТСЯ.

T1 = h + Σ1(e) + Ch(e,f,g) + K + W
ΔT1 = Δh + ΔΣ1 + ΔCh + ΔW = 0

→ ΔW = -(Δh + ΔΣ1 + ΔCh)

W КОМПЕНСИРУЕТ state-разницу. Не обнуляем δ — W ОТВЕЧАЕТ на δ.

Новое видение: collision = W₂ на каждом раунде ТОЧНО компенсирует
разницу state, порождённую предыдущими раундами.
Это ЦЕПОЧКА КОМПЕНСАЦИЙ, не цепочка нулей.
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


def flow_trace(W16):
    """Trace the FLOW: T1 at each round, and what feeds into it."""
    W = expand_real(W16)
    s = list(IV)
    trace = []

    for r in range(64):
        a,b,c,d,e,f,g,h = s

        # T1 components
        part_h = h
        part_sig1 = Sig1(e)
        part_ch = Ch(e,f,g)
        part_k = K[r]
        part_w = W[r]

        T1 = add32(part_h, part_sig1, part_ch, part_k, part_w)
        T2 = add32(Sig0(a), Maj(a,b,c))

        trace.append({
            'r': r,
            'T1': T1,
            'T2': T2,
            'h': part_h, 'sig1': part_sig1, 'ch': part_ch,
            'k': part_k, 'w': part_w,
            'state': list(s),
        })

        s = [add32(T1,T2), a, b, c, add32(d,T1), e, f, g]

    return trace, s


# ============================================================
# VISION 1: T1 matching — if ΔT1=0 at round r, what happens?
# ============================================================
def vision1_t1_matching(N):
    print("="*70)
    print("VISION 1: T1 MATCHING")
    print("  If T1₁[r] = T1₂[r] AND state₁[r] = state₂[r]")
    print("  → state₁[r+1] = state₂[r+1] → collision propagates!")
    print()
    print("  But state₁[r] ≠ state₂[r] in general.")
    print("  So: what if only T1 matches but state differs?")
    print("="*70)

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    trace1, final1 = flow_trace(W1)

    # For each round: compute REQUIRED W2[r] to make T1₁=T1₂
    # Given state₂ at round r, W2[r] must satisfy:
    # T1₂ = h₂ + Σ1(e₂) + Ch(e₂,f₂,g₂) + K[r] + W₂[r] = T1₁[r]
    # → W₂[r] = T1₁[r] - h₂ - Σ1(e₂) - Ch(e₂,f₂,g₂) - K[r]

    # Start: W2 = W1 + small δ
    W2 = list(W1); W2[0] ^= 0x80000000
    trace2, final2 = flow_trace(W2)

    print(f"\n  T1 comparison (W1 vs W2=W1⊕MSB):")
    print(f"  {'r':>3} | {'ΔT1':>10} | {'Δstate_hw':>10} | notes")
    print("  " + "-"*45)

    for r in range(64):
        dT1 = hw(trace1[r]['T1'] ^ trace2[r]['T1'])
        ds = sum(hw(trace1[r]['state'][i] ^ trace2[r]['state'][i]) for i in range(8))
        notes = ""
        if dT1 == 0: notes = " ★ T1 MATCH!"
        if ds == 0: notes = " ★★ STATE MATCH!"
        if r < 5 or r > 58 or dT1 == 0:
            print(f"  {r:3d} | {dT1:10d} | {ds:10d} |{notes}")


# ============================================================
# VISION 2: COMPENSATION CHAIN — W answers state difference
# ============================================================
def vision2_compensation():
    print("\n" + "="*70)
    print("VISION 2: COMPENSATION CHAIN")
    print("  At each round: W₂[r] = T1₁[r] - (h₂ + Σ1(e₂) + Ch₂ + K)")
    print("  This is the EXACT W₂ that makes T1 match.")
    print("  Then: state₂[r+1] = (T1+T2₂, a₂, b₂, c₂, d₂+T1, e₂, f₂, g₂)")
    print("  Even with T1 matching, T2 DIFFERS (because a₂≠a₁)")
    print("  → a_new₂ = T1 + T2₂ ≠ T1 + T2₁ = a_new₁")
    print("  → BUT e_new₂ = d₂ + T1 and e_new₁ = d₁ + T1")
    print("  → δe_new = δd. e inherits d's difference!")
    print("="*70)

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    W2 = list(W1); W2[0] ^= 0x80000000

    trace1, _ = flow_trace(W1)
    Wexp2 = expand_real(W2)

    # Run W2 with FORCED T1 matching
    s2 = list(IV)
    Wexp1 = expand_real(W1)

    print(f"\n  Running W2 with T1-forced matching:")
    print(f"  {'r':>3} | {'δa':>4} {'δe':>4} | {'δT2':>5} | {'forced_W2':>12} {'sched_W2':>12} {'gap':>5}")
    print("  " + "-"*60)

    for r in range(64):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2

        # Force T1₂ = T1₁
        T1_target = trace1[r]['T1']
        W2_forced = sub32(sub32(sub32(sub32(T1_target, h2), Sig1(e2)), Ch(e2,f2,g2)), K[r])
        W2_actual = Wexp2[r]
        gap = hw(W2_forced ^ W2_actual)

        T2_2 = add32(Sig0(a2), Maj(a2,b2,c2))
        a_new2 = add32(T1_target, T2_2)
        e_new2 = add32(d2, T1_target)

        # Compare with message 1
        a_new1 = trace1[r]['state'][0] if r == 0 else add32(trace1[r]['T1'], add32(Sig0(trace1[r]['state'][0]), Maj(trace1[r]['state'][0],trace1[r]['state'][1],trace1[r]['state'][2])))
        # Actually get from next round
        if r < 63:
            a1_next = trace1[r+1]['state'][0]
            e1_next = trace1[r+1]['state'][4]
        else:
            a1_next = 0; e1_next = 0

        da = hw(a_new2 ^ a1_next) if r < 63 else 0
        de = hw(e_new2 ^ e1_next) if r < 63 else 0

        if r < 10 or r > 55 or gap < 5:
            print(f"  {r:3d} | {da:4d} {de:4d} | {hw(T2_2 ^ add32(Sig0(trace1[r]['state'][0]),Maj(trace1[r]['state'][0],trace1[r]['state'][1],trace1[r]['state'][2]))):5d} | 0x{W2_forced:08x} 0x{W2_actual:08x} {gap:5d}")

        s2 = [a_new2, a2, b2, c2, e_new2, e2, f2, g2]


# ============================================================
# VISION 3: THE SPLIT — T1 goes to both a and e
# ============================================================
def vision3_split():
    print("\n" + "="*70)
    print("VISION 3: THE SPLIT — T1 feeds BOTH a_new and e_new")
    print("  a_new = T1 + T2  (T2 depends on a,b,c — 'left branch')")
    print("  e_new = d + T1   (d = a[-3] — 'right branch')")
    print()
    print("  T1 is the SHARED TRUNK. a and e are TWO BRANCHES.")
    print("  δa_new = δT1 + δT2")
    print("  δe_new = δd + δT1")
    print("  Subtract: δa_new - δe_new = δT2 - δd")
    print("  → δ(a-e) = δT2 - δd (INDEPENDENT of T1!)")
    print()
    print("  THIS IS THE D-COUPLING THEOREM seen as FLOW:")
    print("  The DIFFERENCE between two branches doesn't depend")
    print("  on the shared trunk (T1/message).")
    print("  Only the LEFT-BRANCH-ONLY part (T2) and the")
    print("  RIGHT-BRANCH-ONLY part (d) matter.")
    print("="*70)

    # This means: to make δ(a-e)=0, we need δT2 = δd.
    # T2 = Σ0(a) + Maj(a,b,c)
    # d = a[-3]
    # So δT2 = δd means:
    # δ(Σ0(a) + Maj(a,b,c)) = δa[-3]
    # This connects CURRENT a-state to 3-rounds-ago a-state.

    # For the LAST 4 rounds (collision target):
    # δ(a-e)[64] = δT2[63] - δd[63] = δT2[63] - δa[60]
    # For this to be 0: δT2[63] = δa[60]

    print(f"\n  Collision condition in FLOW terms:")
    print(f"    δ(a-e)[r+1] = δT2[r] - δa[r-3]")
    print(f"    = δ(Σ0(a[r]) + Maj(a[r],a[r-1],a[r-2])) - δa[r-3]")
    print(f"")
    print(f"  For collision at r=64: need δa[61..64]=0 AND δe[61..64]=0")
    print(f"    δa[61..64]=0: from T1+T2 matching")
    print(f"    δe[61..64]=0: from δd + δT1 = 0, i.e. δa[58..61] + δT1 = 0")
    print(f"")
    print(f"  The FLOW INSIGHT: message (W) only controls T1.")
    print(f"  T2 is message-FREE (depends only on state a,b,c).")
    print(f"  So the LEFT BRANCH (a-chain via T2) flows AUTONOMOUSLY.")
    print(f"  Only the RIGHT BRANCH (e-chain via d+T1) is steerable by W.")
    print(f"")
    print(f"  COLLISION = left branch NATURALLY converges +")
    print(f"             right branch STEERED by W to match.")


if __name__ == '__main__':
    vision1_t1_matching(1)
    vision2_compensation()
    vision3_split()
