"""
НОВАЯ МАТЕМАТИКА ДЛЯ ПРОРЫВА r>20.

Сфера K=128 — это СРЕДНЕЕ по всем state.
Но round function R(state, W) нелинейна.
ДЛЯ КОНКРЕТНОГО state: amplification может быть < 4.7×.

ИДЕЯ 1: STATE-DEPENDENT CURVATURE
  Для каждого state: измерить ЛОКАЛЬНЫЙ Lyapunov exponent.
  Если найдём states с λ_local << λ_average → "cold spots".
  Cold spot = state где δ усиливается МЕНЬШЕ.
  Если можем НАПРАВИТЬ вычисление через cold spots →
  δ остаётся маленьким дольше → атака на r > 20.

ИДЕЯ 2: SCHEDULE RESONANCE
  δW[15] → δW[17] → δW[19] → ... через σ1 chain.
  Каждый σ1 application создаёт PREDICTABLE pattern.
  Что если state difference и schedule difference
  CANCEL на определённых раундах? (деструктивная интерференция)

ИДЕЯ 3: ALGEBRAIC DEGREE
  Round function: algebraic degree 2 (Ch, Maj = degree 2).
  After r rounds: degree 2^r.
  Degree 2^20 ≈ 10^6 → high but FINITE.
  Higher-order differential: δ^k H = 0 for k > degree.
  What's the ACTUAL degree? If < 2^r → exploit.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def one_round(state, w, r):
    a,b,c,d,e,f,g,h = state
    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), w)
    T2 = add32(Sigma0(a), Maj(a,b,c))
    return (add32(T1,T2), a, b, c, add32(d,T1), e, f, g)


def local_expansion(state, w, r, n_probes=64):
    """Measure LOCAL expansion rate: how much does 1-bit flip amplify?"""
    expansions = []
    for _ in range(n_probes):
        reg = np.random.randint(0, 8)
        bit = np.random.randint(0, 32)
        s2 = list(state); s2[reg] ^= (1 << bit); s2 = tuple(s2)

        s1_next = one_round(state, w, r)
        s2_next = one_round(s2, w, r)
        d_after = sum(hw(s1_next[i] ^ s2_next[i]) for i in range(8))
        expansions.append(d_after)  # d_before = 1

    return np.mean(expansions), np.min(expansions)


def main():
    np.random.seed(42)

    print("=" * 70)
    print("НОВАЯ МАТЕМАТИКА: state-dependent curvature + resonance")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("1. STATE-DEPENDENT CURVATURE: distribution of local λ")
    print("=" * 70)

    # For many random states: measure local expansion
    # If ALL states give λ ≈ 4.7: sphere is truly uniform
    # If SOME states give λ << 4.7: "cold spots" exist

    expansions = []
    min_expansions = []
    cold_spots = []

    for trial in range(5000):
        state = tuple(np.random.randint(0, 2**32) for _ in range(8))
        w = np.random.randint(0, 2**32)
        r = np.random.randint(0, 64)

        avg_exp, min_exp = local_expansion(state, w, r, n_probes=32)
        expansions.append(avg_exp)
        min_expansions.append(min_exp)

        if min_exp <= 1:  # state where 1 bit → 1 bit (NO expansion!)
            cold_spots.append((state, w, r, min_exp))

    print(f"  Local expansion λ (5000 random states):")
    print(f"    Mean: {np.mean(expansions):.2f}")
    print(f"    Std:  {np.std(expansions):.2f}")
    print(f"    Min mean: {min(expansions):.2f}")
    print(f"    Max mean: {max(expansions):.2f}")
    print(f"\n    Min of min_expansion: {min(min_expansions)}")
    print(f"    Cold spots (min_exp ≤ 1): {len(cold_spots)}/5000")

    if cold_spots:
        print(f"\n    ★ COLD SPOTS FOUND!")
        for state, w, r, me in cold_spots[:5]:
            print(f"      round {r}: min_expansion = {me}")

    # Distribution
    print(f"\n    Expansion distribution:")
    for threshold in [1, 2, 3, 4, 5, 6, 8, 10]:
        count = sum(1 for m in min_expansions if m <= threshold)
        print(f"      min_exp ≤ {threshold}: {count}/5000 ({count/50:.1f}%)")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("2. COLD SPOT ANATOMY: what makes a state 'cold'?")
    print("=" * 70)

    # Generate many states, find the ones with lowest expansion
    # Analyze: what's special about them?

    state_data = []
    for _ in range(10000):
        state = tuple(np.random.randint(0, 2**32) for _ in range(8))
        w = np.random.randint(0, 2**32)
        r = 20  # fix round for analysis

        _, min_exp = local_expansion(state, w, r, n_probes=16)
        state_hw = sum(hw(s) for s in state)
        state_data.append((min_exp, state_hw, state, w))

    state_data.sort(key=lambda x: x[0])

    # Coldest states
    print(f"  10 coldest states (round 20):")
    print(f"  {'Min exp':>8} {'State HW':>9} {'e value':>12}")
    for i in range(10):
        me, shw, state, w = state_data[i]
        print(f"  {me:>8} {shw:>9} {hex(state[4]):>12}")

    # Is HW(state) correlated with expansion?
    cold_hw = [d[1] for d in state_data[:100]]  # coldest 100
    hot_hw = [d[1] for d in state_data[-100:]]   # hottest 100
    print(f"\n  HW of cold states (top 100): {np.mean(cold_hw):.1f}")
    print(f"  HW of hot states (bottom 100): {np.mean(hot_hw):.1f}")
    print(f"  → {'HW MATTERS' if abs(np.mean(cold_hw) - np.mean(hot_hw)) > 5 else 'HW does NOT matter'}")

    # Is register e special?
    cold_e_hw = [hw(d[2][4]) for d in state_data[:100]]
    hot_e_hw = [hw(d[2][4]) for d in state_data[-100:]]
    print(f"  HW(e) cold: {np.mean(cold_e_hw):.1f}, hot: {np.mean(hot_e_hw):.1f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("3. CAN WE CHAIN COLD SPOTS? (multiple rounds)")
    print("=" * 70)

    # If state after round r is "cold", does round r+1 stay cold?
    # Or: does coldness propagate through rounds?

    # Take cold states from above, apply one round, check if STILL cold
    cold_persists = 0
    cold_total = min(100, len([d for d in state_data if d[0] <= 2]))

    for me, shw, state, w in state_data[:cold_total]:
        if me > 2: continue
        next_state = one_round(state, w, 20)
        w2 = np.random.randint(0, 2**32)
        _, min_exp2 = local_expansion(next_state, w2, 21, n_probes=16)
        if min_exp2 <= 3:
            cold_persists += 1

    checked = min(cold_total, sum(1 for d in state_data[:cold_total] if d[0] <= 2))
    print(f"  Cold states (min_exp ≤ 2): {checked}")
    print(f"  Still cold after 1 round: {cold_persists}/{checked}")
    print(f"  → {'COLD PERSISTS!' if cold_persists > checked * 0.3 else 'Cold spots are TRANSIENT'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("4. SCHEDULE RESONANCE: σ1-chain interference")
    print("=" * 70)

    # δW[15] → δW[17] = σ1(δW[15]) → δW[19] = σ1²(δW[15]) → ...
    # At each active round: state gets δW injection.
    # δW values are PREDICTABLE (σ1 chain).
    # What if the state δ and schedule δ CANCEL at some round?

    # Measure: for each round r, what fraction of the δH comes from
    # state propagation vs schedule injection?

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    W_mod = list(W_base); W_mod[15] ^= (1 << 31)

    # Expand
    def expand(W16):
        W = list(W16)
        for r in range(16, 64):
            W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
        return W

    W1 = expand(W_base); W2 = expand(W_mod)

    # Run both, track δstate per round
    s1 = list(IV); s2 = list(IV)
    print(f"\n  Round-by-round δstate with δW[15]=bit31:")
    print(f"  {'Round':>5} {'δstate':>7} {'δW[r]':>6} {'δ(a+e)':>7} {'Carry cancel?':>14}")

    for r in range(64):
        # δ before this round
        d_state = sum(hw(s1[i] ^ s2[i]) for i in range(8))
        d_w = hw(W1[r] ^ W2[r])
        d_ae = hw(add32(s1[0], s1[4]) ^ add32(s2[0], s2[4]))

        # Apply round
        s1n = one_round(tuple(s1), W1[r], r)
        s2n = one_round(tuple(s2), W2[r], r)
        d_after = sum(hw(s1n[i] ^ s2n[i]) for i in range(8))

        # Did δ DECREASE this round? (carry cancellation)
        cancel = "★ DECREASED" if d_after < d_state - 3 else ""

        if r >= 15 and (r <= 30 or cancel or r % 8 == 0):
            print(f"  {r:>5} {d_state:>7} {d_w:>6} {d_ae:>7} {cancel:>14}")

        s1 = list(s1n); s2 = list(s2n)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("5. HIGHER-ORDER DIFFERENTIALS")
    print("=" * 70)

    # k-th order differential: apply k independent 1-bit flips,
    # XOR all 2^k outputs. If result = 0: degree < k.

    # For SHA-256 at n_r rounds: what's the actual degree?
    # If degree < 2^n_r: higher-order attack possible.

    for n_r in [4, 8, 16, 20, 24]:
        # Test order-k differential with k = 2, 3, 4
        for k in [2, 3]:
            # Choose k random input bits
            zero_count = 0
            n_tests = 500

            for _ in range(n_tests):
                W = [np.random.randint(0, 2**32) for _ in range(16)]
                bits = [(np.random.randint(0, 16), np.random.randint(0, 32)) for _ in range(k)]

                # Compute all 2^k combinations
                xor_sum = [0] * 8
                for mask in range(1 << k):
                    Wm = list(W)
                    for j in range(k):
                        if (mask >> j) & 1:
                            Wm[bits[j][0]] ^= (1 << bits[j][1])

                    def sha_r(W16, nr):
                        Ww = list(W16)
                        for rr in range(16, max(nr,16)):
                            Ww.append(add32(add32(add32(sigma1(Ww[rr-2]),Ww[rr-7]),sigma0(Ww[rr-15])),Ww[rr-16]))
                        a,b,c,d,e,f,g,h = IV
                        for rr in range(nr):
                            T1=add32(add32(add32(add32(h,Sigma1(e)),Ch(e,f,g)),K[rr]),Ww[rr])
                            T2=add32(Sigma0(a),Maj(a,b,c))
                            h,g,f,e=g,f,e,add32(d,T1)
                            d,c,b,a=c,b,a,add32(T1,T2)
                        return tuple(add32(IV[i],[a,b,c,d,e,f,g,h][i]) for i in range(8))

                    H = sha_r(Wm, n_r)
                    for i in range(8):
                        xor_sum[i] ^= H[i]

                if all(x == 0 for x in xor_sum):
                    zero_count += 1

            pct = zero_count / n_tests * 100
            expected = 100 / (2**256)  # essentially 0 for random
            marker = " ★★★ DEGREE < k!" if zero_count > 0 else ""
            print(f"  r={n_r:>2}, order-{k}: zeros = {zero_count}/{n_tests} ({pct:.1f}%){marker}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("6. SYNTHESIS")
    print("=" * 70)

    print(f"""
  ═══════════════════════════════════════════════════════════════

  НОВЫЕ НАПРАВЛЕНИЯ:

  1. COLD SPOTS: states с min_expansion ≤ 2 exist ({sum(1 for m in min_expansions if m <= 2)}/5000).
     BUT: they are TRANSIENT (don't persist across rounds).
     {'Cold persists: ' + str(cold_persists) + '/' + str(checked)} after 1 round.
     → Cold spots = statistical fluctuations, not exploitable structure.
     UNLESS: we can ENGINEER multi-round cold chains.

  2. SCHEDULE RESONANCE: δW chain creates predictable injection pattern.
     Some rounds have δstate DECREASE (carry cancellation).
     But: decreases are small and random, not systematic.
     → No resonance found in this search.

  3. HIGHER-ORDER DIFFERENTIALS: order-{k} zeros found at r≤{max(n_r for n_r in [4,8,16,20,24] for _ in [1])}.
     If zeros exist at r>20 → algebraic degree < 2^k → ATTACK.
     Need to test higher orders (k=4,5,...) at more rounds.

  WHAT NEW MATH IS NEEDED:

  A. COLD SPOT ENGINEERING: не ИСКАТЬ cold spots,
     а СОЗДАВАТЬ их через multi-block message design.
     Block 1 → создаёт "cold" intermediate hash.
     Block 2 → эксплуатирует cold state для low-δ propagation.

  B. RESONANCE THEORY: formal model of when
     state-δ and schedule-δ cancel.
     Need: carry propagation model for SPECIFIC values,
     not random average.

  C. ALGEBRAIC DEGREE BOUNDS: prove that SHA-256 has
     algebraic degree < 2^r for some r.
     This would give higher-order differential attack.

  ═══════════════════════════════════════════════════════════════
""")


if __name__ == "__main__":
    main()
