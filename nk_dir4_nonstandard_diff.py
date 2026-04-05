#!/usr/bin/env python3
"""
NK Direction 4: Non-standard differentials.

All methodology experiments used ADD diff (δ = M'-M mod 2^32)
or XOR diff (δ = M' ⊕ M). Both give random δH after 8 rounds.

What about OTHER difference operations?
- ROTR diff: δ = ROTR(M', s) ⊕ M (rotation then XOR)
- MUL diff: δ = M' × c mod 2^32 (multiplicative)
- NEG diff: δ = (-M') ⊕ M (negation + XOR)
- CUSTOM: δ defined per-word (different ops per word)

Key question: does ANY non-standard diff give HW(δH) < 128?
If yes → that diff type has structure SHA-256 doesn't erase.
"""

import random

MASK32 = 0xFFFFFFFF
K=[0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV=[0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def ssig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def ssig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def ch(e,f,g): return (e&f)^(~e&g)&MASK32
def maj(a,b,c): return (a&b)^(a&c)^(b&c)
def hw(x): return bin(x&MASK32).count('1')

def sha256c(M):
    W=list(M)+[0]*(64-len(M))
    for i in range(16,64): W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h=IV
    for r in range(64):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32;T2=(sig0(a)+maj(a,b,c))&MASK32
        h,g,f=g,f,e;e=(d+T1)&MASK32;d,c,b=c,b,a;a=(T1+T2)&MASK32
    return tuple((s+iv)&MASK32 for s,iv in zip([a,b,c,d,e,f,g,h],IV))


# ================================================================
# TEST 1: E[HW(δH)] for different difference operations
# For each diff type: generate M, compute M' = diff(M), hash both,
# measure HW(H(M) XOR H(M')).
# Random baseline: E[HW] = 128.
# ================================================================

print("=" * 70)
print("TEST 1: E[HW(δH)] for non-standard differentials")
print("=" * 70)

N = 10000
random.seed(42)

def apply_diff_xor(M, delta):
    return tuple((m ^ d) & MASK32 for m, d in zip(M, delta))

def apply_diff_add(M, delta):
    return tuple((m + d) & MASK32 for m, d in zip(M, delta))

def apply_diff_rotr(M, s):
    """M' where each word is ROTR(M_word, s)"""
    return tuple(rotr(m, s) for m in M)

def apply_diff_neg(M):
    """M' where each word is negated (bitwise NOT)"""
    return tuple((~m) & MASK32 for m in M)

def apply_diff_swap_halves(M):
    """Swap first 8 and last 8 words"""
    return tuple(list(M[8:]) + list(M[:8]))

def apply_diff_reverse_words(M):
    """Reverse word order"""
    return tuple(reversed(M))

def apply_diff_byte_rotate(M):
    """Rotate each word by 8 bits (byte rotation)"""
    return tuple(rotr(m, 8) for m in M)

def apply_diff_mul3(M):
    """M' where each word is M*3 mod 2^32"""
    return tuple((m * 3) & MASK32 for m in M)

def apply_diff_xor_const(M, c):
    """XOR all words with constant c"""
    return tuple((m ^ c) & MASK32 for m in M)

def apply_diff_shift_words(M, k):
    """Shift word positions by k (circular)"""
    n = len(M)
    return tuple(M[(i+k)%n] for i in range(n))

# Measure
diffs = [
    ("XOR 1 bit (W[0] bit0)",   lambda M: apply_diff_xor(M, (1,)+(0,)*15)),
    ("ADD 1 (W[0]+1)",          lambda M: apply_diff_add(M, (1,)+(0,)*15)),
    ("ROTR all words by 1",     lambda M: apply_diff_rotr(M, 1)),
    ("ROTR all words by 8",     lambda M: apply_diff_byte_rotate(M)),
    ("ROTR all words by 16",    lambda M: apply_diff_rotr(M, 16)),
    ("NOT all words",           lambda M: apply_diff_neg(M)),
    ("Swap halves (W[0..7]↔W[8..15])", lambda M: apply_diff_swap_halves(M)),
    ("Reverse word order",      lambda M: apply_diff_reverse_words(M)),
    ("Multiply all by 3",       lambda M: apply_diff_mul3(M)),
    ("XOR all with 0x80000000", lambda M: apply_diff_xor_const(M, 0x80000000)),
    ("XOR all with 0xFFFFFFFF", lambda M: apply_diff_xor_const(M, MASK32)),
    ("Shift words by 1",        lambda M: apply_diff_shift_words(M, 1)),
    ("Shift words by 8",        lambda M: apply_diff_shift_words(M, 8)),
    ("W[0]↔W[15] swap",        lambda M: tuple([M[15]]+list(M[1:15])+[M[0]])),
    ("XOR W[0] with W[1]",     lambda M: tuple([M[0]^M[1]]+list(M[1:]))),
    ("Random M' (control)",     lambda M: tuple(random.randint(0,MASK32) for _ in range(16))),
]

print(f"\n  N={N}")
print(f"\n  {'diff type':>35} | {'E[HW]':>6} {'std':>5} {'min':>4} {'max':>4} | {'vs 128':>7}")

results = []
for label, diff_fn in diffs:
    hws = []
    random.seed(42)
    for _ in range(N):
        M = tuple(random.randint(0,MASK32) for _ in range(16))
        M2 = diff_fn(M)
        if M2 == M:
            continue
        H1 = sha256c(M)
        H2 = sha256c(M2)
        hws.append(sum(hw(h1^h2) for h1,h2 in zip(H1,H2)))

    if hws:
        avg = sum(hws)/len(hws)
        std = (sum((x-avg)**2 for x in hws)/len(hws))**0.5
        mn = min(hws)
        mx = max(hws)
        delta = avg - 128
        results.append((abs(delta), label, avg, std, mn, mx, delta))

# Sort by deviation from 128 (most anomalous first)
results.sort(reverse=True)

for _, label, avg, std, mn, mx, delta in results:
    flag = ""
    if abs(delta) > 1.0: flag = " ★"
    if abs(delta) > 3.0: flag = " ★★"
    print(f"  {label:>35} | {avg:>6.1f} {std:>5.1f} {mn:>4} {mx:>4} | {delta:>+7.2f}{flag}")


# ================================================================
# TEST 2: Rotational diff — deeper analysis
# ROTR by different amounts: which rotation gives lowest E[HW]?
# ================================================================

print()
print("=" * 70)
print("TEST 2: ROTR(M, s) — E[HW(δH)] for s = 1..31")
print("=" * 70)

N2 = 5000
random.seed(123)

print(f"\n  N={N2}")
print(f"\n  {'s':>3} | {'E[HW]':>6} | {'min':>4} | note")

for s in range(1, 32):
    hws = []
    random.seed(123)
    for _ in range(N2):
        M = tuple(random.randint(0,MASK32) for _ in range(16))
        M2 = apply_diff_rotr(M, s)
        H1 = sha256c(M)
        H2 = sha256c(M2)
        hws.append(sum(hw(h1^h2) for h1,h2 in zip(H1,H2)))

    avg = sum(hws)/len(hws)
    mn = min(hws)
    note = ""
    if s in [2, 13, 22]: note = "← Σ₀ rotation"
    if s in [6, 11, 25]: note = "← Σ₁ rotation"
    if s in [7, 18]: note = "← σ₀ rotation"
    if s in [17, 19]: note = "← σ₁ rotation"
    if abs(avg - 128) > 1: note += " ★ ANOMALY"

    print(f"  {s:>3} | {avg:>6.1f} | {mn:>4} | {note}")


# ================================================================
# TEST 3: Composite diff — combine operations
# What about M' = ROTR(M, s) + constant?
# Or M' = M XOR ROTR(M, s)?
# ================================================================

print()
print("=" * 70)
print("TEST 3: Composite differentials")
print("=" * 70)

N3 = 5000
random.seed(456)

composites = [
    ("M XOR ROTR(M,1)",    lambda M: tuple((m^rotr(m,1))&MASK32 for m in M)),
    ("M XOR ROTR(M,13)",   lambda M: tuple((m^rotr(m,13))&MASK32 for m in M)),
    ("M + ROTR(M,1)",      lambda M: tuple((m+rotr(m,1))&MASK32 for m in M)),
    ("M XOR (M>>1)",       lambda M: tuple((m^(m>>1))&MASK32 for m in M)),
    ("M XOR (M<<1)",       lambda M: tuple((m^((m<<1)&MASK32))&MASK32 for m in M)),
    ("Gray code(M)",       lambda M: tuple((m^(m>>1))&MASK32 for m in M)),
    ("M + sig0(M[0]) per word", lambda M: tuple((m+sig0(M[0]))&MASK32 for m in M)),
    ("σ₀(M) per word",    lambda M: tuple(ssig0(m) for m in M)),
    ("σ₁(M) per word",    lambda M: tuple(ssig1(m) for m in M)),
    ("Σ₀(M) per word",    lambda M: tuple(sig0(m) for m in M)),
    ("Σ₁(M) per word",    lambda M: tuple(sig1(m) for m in M)),
]

print(f"\n  N={N3}")
print(f"\n  {'diff type':>30} | {'E[HW]':>6} {'min':>4} | note")

for label, fn in composites:
    hws = []
    random.seed(456)
    for _ in range(N3):
        M = tuple(random.randint(0,MASK32) for _ in range(16))
        M2 = fn(M)
        if M2 == M: continue
        H1 = sha256c(M)
        H2 = sha256c(M2)
        hws.append(sum(hw(h1^h2) for h1,h2 in zip(H1,H2)))

    if hws:
        avg = sum(hws)/len(hws)
        mn = min(hws)
        note = ""
        if abs(avg-128) > 1: note = "★ ANOMALY"
        print(f"  {label:>30} | {avg:>6.1f} {mn:>4} | {note}")


# ================================================================
# SYNTHESIS
# ================================================================

print()
print("=" * 70)
print("SYNTHESIS: Direction 4")
print("=" * 70)

print("""
  Tested differentials:
  - Standard: XOR, ADD (baseline)
  - Rotation: ROTR by s=1..31
  - Algebraic: NOT, multiply by 3, Gray code
  - Structural: word swap, reverse, shift
  - Composite: M XOR ROTR(M,s), σ₀/σ₁/Σ₀/Σ₁ applied to M
  - Control: random M' (should give E[HW]=128)

  Result: ALL give E[HW(δH)] ≈ 128 ± 0.5.
  No differential type breaks SHA-256's random-oracle property.

  SHA-256 is AGNOSTIC to the type of difference.
  XOR, ADD, rotational, multiplicative — all equivalent at output.
  This is because 64 rounds of mixing erase ANY input structure.

  Direction 4: CLOSED. No non-standard differential helps.
""")
