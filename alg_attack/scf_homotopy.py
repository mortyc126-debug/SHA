#!/usr/bin/env python3
"""
SCF: ГОМОТОПИЯ — деформация от решаемой системы к SHA-256.

EXP2 показал: прямой путь t·δW не работает (дискретно, oscillates).
Новая идея: homotopy НЕ по t, а по СТРУКТУРЕ SHA.

H(λ): система уравнений, параметризованная λ ∈ [0,1].
  λ=0: SHA_xor (без carry) — ЛИНЕЙНАЯ, решаем точно
  λ=1: SHA_real (полная) — цель

Промежуточные λ: SHA с ЧАСТИЧНЫМ carry.
Carry включается постепенно, от бита 0 к биту 31.

λ=k/32: carry работает на битах 0..k-1, остальные = XOR.

На каждом шаге: берём решение от предыдущего λ,
корректируем для нового λ. Если шаг мал → коррекция мала.

Это CONTINUATION METHOD — стандарт в численном анализе,
но применённый к ДИСКРЕТНОЙ задаче SHA.
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


def add_truncated(a, b, n_carry_bits):
    """Addition with carry only on lower n_carry_bits.
    Bits 0..n_carry_bits-1: real addition (with carry).
    Bits n_carry_bits..31: XOR (no carry)."""
    if n_carry_bits >= 32:
        return add32(a, b)
    if n_carry_bits <= 0:
        return a ^ b

    mask_low = (1 << n_carry_bits) - 1
    mask_high = MASK32 ^ mask_low

    low = (a & mask_low) + (b & mask_low)  # Real add on low bits
    carry_out = (low >> n_carry_bits) & 1
    low &= mask_low

    high = (a & mask_high) ^ (b & mask_high)  # XOR on high bits
    # Carry from low part propagates ONE bit into high? NO — we cut it.

    return (high | low) & MASK32


def add_truncated_multi(n_carry_bits, *args):
    s = 0
    for x in args:
        s = add_truncated(s, x, n_carry_bits)
    return s


def expand_lambda(W16, n_carry_bits):
    W = list(W16)
    for i in range(16, 64):
        W.append(add_truncated_multi(n_carry_bits,
                 sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W


def sha_lambda(W16, n_carry_bits, iv=None):
    """SHA-256 with carry truncated to n_carry_bits.
    n_carry_bits=0: pure XOR-SHA (linear)
    n_carry_bits=32: real SHA-256"""
    if iv is None: iv = IV
    W = expand_lambda(W16, n_carry_bits)
    s = list(iv)
    for r in range(64):
        a,b,c,d,e,f,g,h = s
        T1 = add_truncated_multi(n_carry_bits, h, Sig1(e), Ch(e,f,g), K[r], W[r])
        T2 = add_truncated_multi(n_carry_bits, Sig0(a), Maj(a,b,c))
        a_new = add_truncated(T1, T2, n_carry_bits)
        e_new = add_truncated(d, T1, n_carry_bits)
        s = [a_new, a, b, c, e_new, e, f, g]
    return [add_truncated(iv[i], s[i], n_carry_bits) for i in range(8)]


def hash_diff_lambda(W1, W2, n_carry_bits):
    H1 = sha_lambda(W1, n_carry_bits)
    H2 = sha_lambda(W2, n_carry_bits)
    return sum(hw(H1[i]^H2[i]) for i in range(8))


# ============================================================
# EXP 1: Collision at λ=0 (XOR-SHA) → lift to λ=k
# ============================================================
def exp1_xor_collision():
    """Find collision in XOR-SHA, then track it as carry increases."""
    print("="*70)
    print("EXP 1: XOR-SHA COLLISION → CARRY CONTINUATION")
    print("  Step 1: Find collision in XOR-SHA (linear → exact)")
    print("  Step 2: Track HW as carry bits increase 0→32")
    print("="*70)

    # XOR-SHA collision: since it's GF(2)-linear (after Ch/Maj),
    # we can find collisions by solving linear system.
    # Simpler: brute-force for reduced rounds, or use kernel.

    # For full 64-round XOR-SHA: we know rank=256 (from Carleman).
    # So XOR-collision requires 2^{256/2} = 2^128. Can't brute force.

    # BUT: for REDUCED round XOR-SHA (R rounds), rank < 256 for small R.
    # Find collision for R rounds, then see what happens with carry.

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    # For λ=0 (XOR-SHA), find pair with δH=0 via birthday on truncated hash
    print(f"\n  Step 1: Find XOR-SHA near-collision via birthday")

    # Use small truncation for feasibility
    seen = {}
    best_xor = 256
    best_W2 = None

    for trial in range(50000):
        W2 = list(W1)
        W2[trial % 16] ^= (trial + 1)
        H_xor = sha_lambda(W2, 0)
        H1_xor = sha_lambda(W1, 0)

        d = sum(hw(H_xor[i]^H1_xor[i]) for i in range(8))
        if d < best_xor:
            best_xor = d
            best_W2 = list(W2)
            if d == 0:
                print(f"  ★ XOR-SHA COLLISION at trial {trial}!")
                break

        # Truncated birthday (match first 32 bits)
        trunc = H_xor[0]
        if trunc in seen:
            other_W = seen[trunc]
            H_other = sha_lambda(other_W, 0)
            full_d = sum(hw(H_xor[i]^H_other[i]) for i in range(8))
            if full_d < best_xor:
                best_xor = full_d
                W1_use = other_W
                best_W2 = list(W2)
                if full_d == 0:
                    W1 = W1_use
                    print(f"  ★ XOR-SHA COLLISION via birthday at trial {trial}!")
                    break
        else:
            seen[trunc] = list(W2)

    print(f"  Best XOR-SHA diff: HW = {best_xor}")

    if best_xor > 0:
        print(f"  No exact XOR-collision found. Using best near-collision.")

    if best_W2 is None:
        return

    # Step 2: Track HW as carry increases
    print(f"\n  Step 2: Continuation — increase carry bits 0→32")
    print(f"  {'carry_bits':>10} | HW(δH) | notes")
    print("  " + "-"*40)

    for n_carry in range(33):
        d = hash_diff_lambda(W1, best_W2, n_carry)
        marker = ""
        if d == 0: marker = " ★ COLLISION!"
        elif n_carry == 0: marker = " (XOR-SHA)"
        elif n_carry == 32: marker = " (real SHA)"
        print(f"  {n_carry:10d} | {d:6d} |{marker}")


# ============================================================
# EXP 2: Continuation with correction at each step
# ============================================================
def exp2_continuation_corrected(N):
    """At each λ step: adjust W2 to maintain low δH.
    This is the CORE of homotopy continuation."""
    print("\n" + "="*70)
    print("EXP 2: CONTINUATION WITH CORRECTION")
    print("  At each carry_bit increase, SA-correct W2")
    print("="*70)

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        # Start: find good W2 for XOR-SHA (carry=0)
        def score_xor(w2): return hash_diff_lambda(W1, w2, 0)

        best_xor = score_xor(W2); cur = list(W2)
        for it in range(500):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            if w==0 and b==31: continue
            t[w] ^= (1<<b)
            if t == W1: continue
            s = score_xor(t)
            T = max(0.01, 1-it/500)
            if s < best_xor or math.exp(-(s-best_xor)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur = t
                if s < best_xor: best_xor = s
        W2 = list(cur)

        print(f"\n  Trial {trial}: XOR-SHA optimized: HW={best_xor}")

        # Now: increase carry bits one at a time, correcting W2 each step
        hw_path = [best_xor]

        for n_carry in range(1, 33):
            def score_n(w2): return hash_diff_lambda(W1, w2, n_carry)

            # Correct W2 for new carry level (small SA)
            cur = list(W2); best_n = score_n(cur)
            for it in range(100):
                t = list(cur)
                w = int.from_bytes(os.urandom(1),'big') % 16
                b = int.from_bytes(os.urandom(1),'big') % 32
                if w==0 and b==31: continue
                t[w] ^= (1<<b)
                if t == W1: continue
                s = score_n(t)
                T = max(0.01, 1-it/100)
                if s < best_n or math.exp(-(s-best_n)/(T*1.5)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                    cur = t
                    if s < best_n: best_n = s

            W2 = list(cur)
            hw_path.append(best_n)

        # Show path
        print(f"  Continuation path (carry_bits → HW):")
        for n_carry in [0, 1, 2, 4, 8, 16, 24, 32]:
            marker = " ★" if hw_path[n_carry] < 100 else ""
            print(f"    carry={n_carry:2d}: HW={hw_path[n_carry]:3d}{marker}")

        # Compare with direct SA at carry=32
        print(f"  Direct SA at carry=32: HW≈100")
        print(f"  Continuation at carry=32: HW={hw_path[32]}")
        print(f"  Advantage: {100 - hw_path[32]:+d}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3

    exp1_xor_collision()
    exp2_continuation_corrected(N)

    print(f"\n{'='*70}")
    print(f"ИТОГ ГОМОТОПИИ")
    print(f"{'='*70}")
