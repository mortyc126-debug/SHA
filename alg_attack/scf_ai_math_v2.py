#!/usr/bin/env python3
"""
AI-MATH v2: W ЛИНЕЙНО входит в a[r+1]. Schedule — единственное ограничение.

a[r+1] = STUFF(state[r]) + W[r]  ← W просто прибавляется!
δa[r+1] = δSTUFF + δW[r]
Collision: δW[r] = -δSTUFF для ВСЕХ 64 раундов.

STUFF(state) = h + Σ1(e) + Ch(e,f,g) + K + Σ0(a) + Maj(a,b,c)
δSTUFF = вычислимо из state₁ и state₂.

Задача: 64 ТРЕБУЕМЫХ δW[r] → уместить в schedule из 16 слов.
Schedule: W[r] = σ1(W[r-2]) + W[r-7] + σ0(W[r-15]) + W[r-16]

Это задача: найти W₂[0..15] такие что schedule(W₂)[r] ≈ required[r].
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

def run_sha_with_W(W_expanded):
    """Run SHA with pre-expanded W (64 words)."""
    s = list(IV)
    for r in range(64):
        a,b,c,d,e,f,g,h = s
        T1 = add32(h,Sig1(e),Ch(e,f,g),K[r],W_expanded[r])
        T2 = add32(Sig0(a),Maj(a,b,c))
        s = [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return [add32(IV[i],s[i]) for i in range(8)]


def compute_required_W(W1):
    """For message W1: if W2 matches T1 exactly, what W2[r] is needed?
    This uses T1-forced matching where state₂ = state₁ at all rounds.
    Required W2[r] = W1[r] (trivially, since states are identical).

    BUT: for W2 ≠ W1, states WILL differ. So we need iterative approach.
    Start with W2 = W1 + δ, compute states, compute required W."""

    Wexp1 = expand_real(W1)
    # The "required W" for collision is simply Wexp1 (since T1-matching
    # with identical states = identical messages). Not useful.

    # Instead: start W2 with a DIFFERENT W[0..15] but compute
    # what W_expanded SHOULD be to maintain T1-matching.

    # The key: if we run BOTH messages in PARALLEL, at each step
    # computing required W₂[r] from current state₂ and target T1₁:
    return Wexp1  # Return for now


def parallel_flow(W1, W2_init):
    """Run both messages in parallel.
    For message 1: normal SHA with W1.
    For message 2: at each round, compute what W₂[r] SHOULD BE
    to make T1₂ = T1₁, then check if actual schedule W₂[r] matches."""

    Wexp1 = expand_real(W1)
    Wexp2 = expand_real(W2_init)

    s1 = list(IV); s2 = list(IV)
    required_W2 = []
    actual_W2 = Wexp2

    flow = []

    for r in range(64):
        a1,b1,c1,d1,e1,f1,g1,h1 = s1
        a2,b2,c2,d2,e2,f2,g2,h2 = s2

        # T1 for message 1
        T1_1 = add32(h1, Sig1(e1), Ch(e1,f1,g1), K[r], Wexp1[r])

        # Required W2[r] to make T1₂ = T1₁
        W2_req = sub32(sub32(sub32(sub32(T1_1, h2), Sig1(e2)), Ch(e2,f2,g2)), K[r])
        required_W2.append(W2_req)

        # Gap between required and actual
        gap = hw(W2_req ^ actual_W2[r])

        # Advance both states (message 2 uses ACTUAL schedule, not forced)
        T1_2 = add32(h2, Sig1(e2), Ch(e2,f2,g2), K[r], actual_W2[r])
        T2_1 = add32(Sig0(a1), Maj(a1,b1,c1))
        T2_2 = add32(Sig0(a2), Maj(a2,b2,c2))

        s1 = [add32(T1_1,T2_1),a1,b1,c1,add32(d1,T1_1),e1,f1,g1]
        s2 = [add32(T1_2,T2_2),a2,b2,c2,add32(d2,T1_2),e2,f2,g2]

        da = hw(s1[0] ^ s2[0])
        de = hw(s1[4] ^ s2[4])

        flow.append({
            'r': r, 'gap': gap, 'da': da, 'de': de,
            'W2_req': W2_req, 'W2_act': actual_W2[r],
        })

    H1 = [add32(IV[i],s1[i]) for i in range(8)]
    H2 = [add32(IV[i],s2[i]) for i in range(8)]
    dH = sum(hw(H1[i]^H2[i]) for i in range(8))

    return flow, required_W2, dH


def iterative_W_fit(W1, W2_init, max_iter=20):
    """Iteratively adjust W2[0..15] to minimize gap between
    required_W and schedule_W across all 64 rounds.

    At each iteration:
    1. Run parallel flow → get required_W2[0..63] and gaps
    2. For r=0..15: just SET W2[r] = required_W2[r] (free!)
    3. For r=16..63: gap = required - schedule. Adjust W2[0..15]
       to reduce the gap at the WORST round."""

    W2 = list(W2_init)
    best_dH = 256
    best_W2 = list(W2)

    for iteration in range(max_iter):
        flow, required, dH = parallel_flow(W1, W2)

        if dH < best_dH:
            best_dH = dH
            best_W2 = list(W2)

        # For r=0..15: set W2[r] = required[r] directly!
        for r in range(16):
            W2[r] = required[r]

        # Recompute after fixing r=0..15
        flow2, required2, dH2 = parallel_flow(W1, W2)

        if dH2 < best_dH:
            best_dH = dH2
            best_W2 = list(W2)

        # Check: did fixing r=0..15 help?
        total_gap = sum(flow2[r]['gap'] for r in range(16, 64))

        if iteration < 5 or dH2 < 100:
            print(f"  Iter {iteration}: dH={dH2}, gap[16..63]={total_gap} "
                  f"({total_gap/48:.1f}/round)")

        if dH2 == 0:
            print(f"  ★★★ COLLISION!")
            break

        # For r≥16: find the round with SMALLEST gap and try to fix it
        # by adjusting the W2[0..15] word that has most influence on that round
        worst_r = max(range(16, 64), key=lambda r: flow2[r]['gap'])

        # Try each word perturbation
        for word in range(16):
            W2_test = list(W2)
            W2_test[word] = add32(W2[word], 1)
            Wexp_test = expand_real(W2_test)
            # Did gap at worst_r decrease?
            old_gap = flow2[worst_r]['gap']
            new_val = Wexp_test[worst_r]
            new_gap = hw(required2[worst_r] ^ new_val)
            if new_gap < old_gap:
                # Check: does full flow improve?
                _, _, dH_test = parallel_flow(W1, W2_test)
                if dH_test < best_dH:
                    W2 = W2_test
                    best_dH = dH_test
                    best_W2 = list(W2)
                    break

    return best_W2, best_dH


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    print("="*70)
    print("AI-MATH v2: W IS LINEAR CORRECTOR + SCHEDULE FIT")
    print("="*70)

    results = []

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2_init = list(W1); W2_init[0] ^= 0x80000000

        print(f"\n--- Trial {trial} ---")

        # Phase 1: show the parallel flow
        flow, required, dH_raw = parallel_flow(W1, W2_init)
        print(f"  Raw: dH={dH_raw}")
        print(f"  Gap profile [required vs schedule]:")
        for r in [0, 5, 10, 15, 16, 20, 30, 40, 50, 60, 63]:
            print(f"    r={r:2d}: gap={flow[r]['gap']:2d}, δa={flow[r]['da']:3d}, δe={flow[r]['de']:3d}")

        # Phase 2: iterative W-fit
        W2_opt, dH_opt = iterative_W_fit(W1, W2_init, max_iter=10)
        results.append(dH_opt)

        # Check nontrivial
        n_diff = sum(1 for i in range(16) if W1[i] != W2_opt[i])
        print(f"  Optimized: dH={dH_opt}, words_diff={n_diff}")

    print(f"\n{'='*70}")
    print(f"SUMMARY ({N} trials)")
    print(f"  dH results: {results}")
    print(f"  avg={sum(results)/len(results):.1f}, min={min(results)}")
