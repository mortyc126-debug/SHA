#!/usr/bin/env python3
"""
SCF: SEQUENCE SOLVER — инструмент на основе Sequence Mathematics.

Теорема: SHA = одна последовательность {a[r]} порядка 8.
Collision = δa[57..64] = 0.

Новый инструмент: BACKWARD SEQUENCE INVERSION.
1. Зафиксировать цель: a2[57..64] = a1[57..64]
2. Обратить рекурренцию: из a[r+1] вычислить НУЖНЫЙ W[r]
3. Сравнить нужный W[r] с фактическим schedule W[r]
4. Mismatch = разница → оптимизировать W[0..15]

Это ПРЯМОЕ использование D[r] формулы:
  e[r] = a[r] - D[r] где D зависит только от a-sequence.
  Обратный раунд: W[r] = a[r+1] - T2[r] - h[r] - Σ1(e[r]) - Ch(...) - K[r]
  Всё выражается через a-sequence!
"""
import os, sys, math

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

def sha_a_sequence(W16):
    """Extract just the a-sequence. ALL state info is here."""
    W=expand_real(W16); s=list(IV)
    a_seq = [s[0]]
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        a_seq.append(s[0])
    return a_seq, W

def D_from_a(a_seq, r):
    """D[r] = Σ0(a[r-1]) + Maj(a[r-1],a[r-2],a[r-3]) - a[r-4]."""
    return sub32(add32(Sig0(a_seq[r-1]), Maj(a_seq[r-1],a_seq[r-2],a_seq[r-3])),
                 a_seq[r-4])

def e_from_a(a_seq, r):
    """e[r] = a[r] - D[r]."""
    return sub32(a_seq[r], D_from_a(a_seq, r))

def required_W_from_a_seq(a_seq, r):
    """Given a-sequence, compute what W[r] must be.
    From the round function: a[r+1] = T1 + T2
    T2 = Σ0(a[r]) + Maj(a[r],a[r-1],a[r-2])
    T1 = a[r+1] - T2
    T1 = h + Σ1(e) + Ch(e,f,g) + K[r] + W[r]
    where h=e[r-3], e=e[r]=a[r]-D[r], f=e[r-1], g=e[r-2]
    W[r] = T1 - h - Σ1(e) - Ch(e,f,g) - K[r]
    """
    T2 = add32(Sig0(a_seq[r]), Maj(a_seq[r], a_seq[r-1], a_seq[r-2]))
    T1 = sub32(a_seq[r+1], T2)

    e_r = e_from_a(a_seq, r)
    e_r1 = e_from_a(a_seq, r-1) if r >= 5 else sub32(a_seq[r-1], D_from_a(a_seq, r-1)) if r >= 4 else IV[4+3-(r-1)] # Approximate for early rounds
    e_r2 = e_from_a(a_seq, r-2) if r >= 6 else IV[4+3-(r-2)] if r-2 >= 0 else IV[4]
    h_r = e_from_a(a_seq, r-3) if r >= 7 else IV[4+3-(r-3)] if r-3 >= 0 else IV[7]

    W_r = sub32(sub32(sub32(sub32(T1, h_r), Sig1(e_r)), Ch(e_r, e_r1, e_r2)), K[r])
    return W_r


# ============================================================
# SEQUENCE SOLVER: target δa[R..R+7]=0
# ============================================================
def sequence_mismatch(W1, W2, target_start=57):
    """Compute mismatch between required and actual schedule
    for δa[target_start..target_start+7]=0 in the a-sequence."""
    a1, Wexp1 = sha_a_sequence(W1)
    a2, Wexp2 = sha_a_sequence(W2)

    total_mm = 0
    per_round = []

    for r in range(target_start, min(target_start + 8, 64)):
        # For δa[r+1]=0: we need a2[r+1] = a1[r+1]
        # Build the a-sequence where a2[target..64] = a1[target..64]
        # Required W2[r] is computed from the DESIRED a-sequence

        # The desired a2-sequence: matches a1 from target_start onwards
        # But to compute required_W, we need the ACTUAL a2 up to round r
        # and then FORCE a2[r+1] = a1[r+1]

        # Required W2[r] from actual a2[r..r-7] + target a2[r+1]=a1[r+1]
        # This is complex because we need hybrid a-sequence.
        # Simplify: just compute mismatch at the hash level for now.

        W_req = required_W_from_a_seq(a1, r)
        mm = hw(W_req ^ Wexp2[r])
        total_mm += mm
        per_round.append(mm)

    return total_mm, per_round


def sequence_score(W1, W2):
    """How close is W2's a-sequence to W1's in last 8 rounds?"""
    a1, _ = sha_a_sequence(W1)
    a2, _ = sha_a_sequence(W2)
    return sum(hw(a1[r] ^ a2[r]) for r in range(57, 65))


def solve_sequence(W1, budget=3000):
    """Find W2 minimizing δa in last 8 rounds."""

    # Screen deltas
    a1, _ = sha_a_sequence(W1)

    screen = []
    for word in range(8):
        for bit in range(32):
            W2 = list(W1); W2[word] ^= (1<<bit)
            s = sequence_score(W1, W2)
            screen.append((s, word, bit))
    screen.sort()

    global_best = 256; global_W2 = None
    per_start = (budget - 256) // 3

    for s, dw, db in screen[:3]:
        W2 = list(W1); W2[dw] ^= (1<<db)

        # SA on sequence score
        cur = list(W2); best = sequence_score(W1, cur)
        for it in range(per_start):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            if w==dw and b==db: continue
            t[w] ^= (1<<b)
            if t == W1: continue
            sc = sequence_score(W1, t)
            T = max(0.01, 1-it/per_start)
            if sc < best or math.exp(-(sc-best)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur = t
                if sc < best: best = sc

        if best < global_best:
            global_best = best
            global_W2 = list(cur)

    return global_W2, global_best


def solve_hash(W1, budget=3000):
    """Standard: minimize HW(δH) with same budget. Baseline."""
    def sha_compress(W16):
        W=expand_real(W16); s=list(IV)
        for r in range(64):
            a,b,c,d,e,f,g,h=s
            T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
            T2=add32(Sig0(a),Maj(a,b,c))
            s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        return [add32(IV[i],s[i]) for i in range(8)]

    H1 = sha_compress(W1)

    screen = []
    for word in range(8):
        for bit in range(32):
            W2 = list(W1); W2[word] ^= (1<<bit)
            H2 = sha_compress(W2)
            s = sum(hw(H1[i]^H2[i]) for i in range(8))
            screen.append((s, word, bit))
    screen.sort()

    global_best = 256; global_W2 = None
    per_start = (budget - 256) // 3

    for s, dw, db in screen[:3]:
        W2 = list(W1); W2[dw] ^= (1<<db)
        cur = list(W2); best = s
        for it in range(per_start):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            if w==dw and b==db: continue
            t[w] ^= (1<<b)
            if t == W1: continue
            H2 = sha_compress(t)
            sc = sum(hw(H1[i]^H2[i]) for i in range(8))
            T = max(0.01, 1-it/per_start)
            if sc < best or math.exp(-(sc-best)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur = t
                if sc < best: best = sc
        if best < global_best:
            global_best = best; global_W2 = list(cur)

    return global_W2, global_best


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    print("="*70)
    print("SCF: SEQUENCE SOLVER vs HASH SOLVER")
    print("  Sequence: minimize δa[57..64]")
    print("  Hash:     minimize HW(δH)")
    print("="*70)

    seq_results = []; hash_results = []
    seq_hash_hw = []  # HW(δH) when optimizing sequence metric

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        # Sequence solver
        W2_seq, seq_score = solve_sequence(W1, budget=3000)

        # Hash solver
        W2_hash, hash_score = solve_hash(W1, budget=3000)

        # Cross-check: what's the hash HW for sequence-optimized result?
        def sha_compress_fn(W16):
            W=expand_real(W16); s=list(IV)
            for r in range(64):
                a,b,c,d,e,f,g,h=s
                T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
                T2=add32(Sig0(a),Maj(a,b,c))
                s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
            return [add32(IV[i],s[i]) for i in range(8)]

        H1 = sha_compress_fn(W1)
        H2_seq = sha_compress_fn(W2_seq)
        seq_hw = sum(hw(H1[i]^H2_seq[i]) for i in range(8))

        seq_results.append(seq_score)
        hash_results.append(hash_score)
        seq_hash_hw.append(seq_hw)

        print(f"  Trial {trial:2d}: seq_metric={seq_score:3d} (→HW={seq_hw:3d})  hash_metric={hash_score:3d}")

    print(f"\n{'='*70}")
    print(f"SUMMARY ({N} trials)")
    print(f"  Sequence solver: avg_metric={sum(seq_results)/N:.1f}")
    print(f"    → avg HW(δH)={sum(seq_hash_hw)/N:.1f}, min={min(seq_hash_hw)}")
    print(f"  Hash solver:     avg HW={sum(hash_results)/N:.1f}, min={min(hash_results)}")
    print(f"  Correlation: seq_metric vs hash_HW")
