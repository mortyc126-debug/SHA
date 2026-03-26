#!/usr/bin/env python3
"""
SCF: SEQUENCE MATHEMATICS — новая математика SHA-256.

ОТКРЫТИЕ: SHA-256 = ОДНА последовательность {a[r]}, r=0..64.
  - e[r] определяется через a[r]: e[r] = a[r] - D[r]
  - D[r] = Σ0(a[r-1]) + Maj(a[r-1],a[r-2],a[r-3]) - a[r-4]
  - Collision = δa[57..64] = 0 (8 consecutive zeros)

Это РАДИКАЛЬНО меняет задачу:
  Старое: обнулить 8 независимых 32-бит регистров (256 бит)
  Новое:  обнулить 8 последовательных значений ОДНОЙ рекурренции

Для рекурренции порядка k: если k+1 последовательных значений
совпадают, ВСЕ последующие тоже (для ЛИНЕЙНОЙ рекурренции).
SHA — нелинейная, но может быть "почти линейная" в хвосте.

Вопросы:
1. Верифицировать: D[r] зависит только от a-последовательности?
2. Какой ПОРЯДОК рекурренции для a[r]?
3. Если δa[r..r+k]=0 для k значений, сколько ещё нулей бесплатно?
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
def Ch(e,f,g): return (e&f)^(~e&g)&MASK32
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def add32(*args):
    s=0
    for x in args: s=(s+x)&MASK32
    return s
def sub32(a,b): return (a-b)&MASK32
def hw(x): return bin(x).count('1')

def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)

def expand_real(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def sha_sequences(W16):
    """Extract the a-sequence and e-sequence."""
    W=expand_real(W16); s=list(IV)
    a_seq = [s[0]]  # a[0] = IV[0]
    e_seq = [s[4]]  # e[0] = IV[4]

    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        a_seq.append(s[0])
        e_seq.append(s[4])

    return a_seq, e_seq, W


# ============================================================
# EXP 1: Verify D[r] = a[r] - e[r] depends only on a-sequence
# ============================================================
def exp1_verify_D():
    print("="*70)
    print("EXP 1: VERIFY D[r] = a[r] - e[r] FORMULA")
    print("  D[r+1] should equal Σ0(a[r]) + Maj(a[r],a[r-1],a[r-2]) - a[r-3]")
    print("="*70)

    W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    a_seq, e_seq, _ = sha_sequences(W)

    D_actual = [sub32(a_seq[r], e_seq[r]) for r in range(65)]

    errors = 0
    for r in range(4, 65):
        # D[r] should be Σ0(a[r-1]) + Maj(a[r-1],a[r-2],a[r-3]) - a[r-4]
        # = T2[r-1] - a[r-4]
        D_predicted = sub32(
            add32(Sig0(a_seq[r-1]), Maj(a_seq[r-1], a_seq[r-2], a_seq[r-3])),
            a_seq[r-4]
        )

        if D_actual[r] != D_predicted:
            errors += 1
            if errors <= 3:
                print(f"  r={r}: MISMATCH actual=0x{D_actual[r]:08x} pred=0x{D_predicted:08x}")

    print(f"\n  Errors: {errors}/61")
    if errors == 0:
        print(f"  ✓ D[r] = Σ0(a[r-1]) + Maj(a[r-1],a[r-2],a[r-3]) - a[r-4]")
        print(f"  ✓ e[r] is FULLY DETERMINED by a-sequence!")
        print(f"  ✓ SHA-256 = ONE sequence {{a[r]}} of order 8")


# ============================================================
# EXP 2: Collision in sequence space
# ============================================================
def exp2_sequence_collision(N):
    print("\n" + "="*70)
    print("EXP 2: COLLISION AS SEQUENCE MATCHING")
    print("  Collision = δa[57..64] = 0 (8 consecutive zeros)")
    print("  How many consecutive δa=0 can we get near the end?")
    print("="*70)

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        a1, _, _ = sha_sequences(W1)
        a2, _, _ = sha_sequences(W2)

        # Count consecutive δa=0 from the END
        consec_end = 0
        for r in range(64, -1, -1):
            if a1[r] == a2[r]:
                consec_end += 1
            else:
                break

        # Count ALL δa=0 in last 20 rounds
        zeros_last20 = sum(1 for r in range(45, 65) if a1[r] == a2[r])

        # SA to maximize consecutive δa=0 from end
        best_consec = consec_end
        cur = list(W2)
        for it in range(500):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            if w==0 and b==31: continue
            t[w] ^= (1<<b)
            if t == W1: continue

            at, _, _ = sha_sequences(t)
            c = 0
            for r in range(64, -1, -1):
                if a1[r] == at[r]: c += 1
                else: break

            import math
            T = max(0.01, 1-it/500)
            if c > best_consec or (c == best_consec and
               int.from_bytes(os.urandom(1),'big') < 50):
                cur = t; best_consec = max(best_consec, c)

        if trial < 5:
            print(f"  Trial {trial}: raw consec_end={consec_end}, "
                  f"optimized={best_consec}, zeros_last20={zeros_last20}")
            if best_consec > 0:
                print(f"    ★ {best_consec} consecutive δa=0 from end!")


# ============================================================
# EXP 3: Recurrence order — how many δa=0 give free propagation?
# ============================================================
def exp3_propagation(N):
    print("\n" + "="*70)
    print("EXP 3: FREE PROPAGATION — δa[r..r+k]=0 → δa[r+k+1]=0?")
    print("  If k consecutive zeros in a-seq → does next one follow?")
    print("  For LINEAR recurrence of order m: k≥m → guaranteed")
    print("="*70)

    # For each k=1..8: if δa[r..r+k-1]=0, what's P(δa[r+k]=0)?
    for k in range(1, 9):
        total = 0; propagated = 0

        for trial in range(N * 50):
            W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            W2 = list(W1)
            W2[trial % 16] ^= ((trial // 16) + 1) & MASK32

            a1, _, _ = sha_sequences(W1)
            a2, _, _ = sha_sequences(W2)

            # Find stretches of k consecutive δa=0
            for r in range(20, 60 - k):
                if all(a1[r+i] == a2[r+i] for i in range(k)):
                    total += 1
                    if a1[r+k] == a2[r+k]:
                        propagated += 1

            if total >= 100:
                break

        if total > 0:
            p = propagated / total
            random_p = 2**(-32)
            ratio = p / random_p if random_p > 0 else 0
            marker = " ★★★" if p > 0.01 else (" ★" if ratio > 10 else "")
            print(f"  k={k}: P(δa[r+k]=0 | δa[r..r+k-1]=0) = "
                  f"{p:.6f} ({propagated}/{total}){marker}")
            if p > 0.01:
                print(f"    → {k} zeros PROPAGATE with P={p:.3f}!")
        else:
            print(f"  k={k}: no stretches of {k} consecutive zeros found")


# ============================================================
# EXP 4: New collision metric — consecutive δa zeros
# ============================================================
def exp4_new_metric(N):
    print("\n" + "="*70)
    print("EXP 4: NEW COLLISION METRIC — max consecutive δa=0")
    print("  Old metric: HW(δH) (Hamming weight of hash diff)")
    print("  New metric: max_k such that δa[r..r+k]=0 for some r≥50")
    print("="*70)

    import math

    old_results = []
    new_results = []

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        a1, _, _ = sha_sequences(W1)
        a2, _, _ = sha_sequences(W2)
        H1 = [add32(IV[i], [a1[64],a1[63],a1[62],a1[61],
               sub32(a1[64],sub32(add32(Sig0(a1[63]),Maj(a1[63],a1[62],a1[61])),a1[60])),0,0,0][i])
               for i in range(1)]  # Skip complex e computation

        # Just use standard hash for old metric
        from scf_ultimate import sha_hash, hash_diff
        H1_full = sha_hash(W1)

        # SA optimizing NEW metric: consecutive δa zeros near end
        def new_score(w2):
            at, _, _ = sha_sequences(w2)
            best_k = 0
            for r in range(50, 65):
                k = 0
                for rr in range(r, 65):
                    if a1[rr] == at[rr]: k += 1
                    else: break
                if k > best_k: best_k = k
            return -best_k  # Negative because we minimize

        def old_score(w2):
            H2 = sha_hash(w2)
            return sum(hw(H1_full[i]^H2[i]) for i in range(8))

        # Optimize new metric
        cur = list(W2); best_new = new_score(cur)
        for it in range(1000):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            if w==0 and b==31: continue
            t[w] ^= (1<<b)
            if t == W1: continue
            s = new_score(t)
            T = max(0.01, 1-it/1000)
            if s < best_new or math.exp(-(s-best_new)/(T*0.5)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur = t
                if s < best_new: best_new = s

        # What's the old metric for the new-metric-optimized result?
        final_old = old_score(cur)
        new_results.append(-best_new)  # k consecutive zeros
        old_results.append(final_old)

        if trial < 5:
            print(f"  Trial {trial}: max_consec_δa_zeros={-best_new}, HW(δH)={final_old}")

    avg_k = sum(new_results)/N
    avg_hw = sum(old_results)/N
    print(f"\n  Summary:")
    print(f"    Avg max consecutive δa=0: {avg_k:.1f}")
    print(f"    Avg HW(δH) at best δa-point: {avg_hw:.1f}")
    print(f"    Need: 8 consecutive for collision")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    exp1_verify_D()
    exp3_propagation(N)
    exp2_sequence_collision(5)
    # exp4 needs scf_ultimate import — skip if unavailable
    try:
        exp4_new_metric(5)
    except ImportError:
        print("\n  (EXP4 skipped — needs scf_ultimate)")
