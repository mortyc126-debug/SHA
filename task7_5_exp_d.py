#!/usr/bin/env python3
"""
ЗАДАНИЕ 7.5 Experiment D: Low-bit partial collisions
k=1,2,4,8: find partial collisions, measure residual HW
"""
import struct, os, time
from collections import defaultdict

MASK = 0xFFFFFFFF
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
H0 = (0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19)

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def Ch(e,f,g): return (e&f)^((~e)&g)&MASK
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def add(*args):
    s=0
    for x in args: s=(s+x)&MASK
    return s
def rw(n): return list(struct.unpack(f'>{n}I', os.urandom(4*n)))

def sha256_hash(M):
    W = list(M[:16])
    for i in range(16,64):
        W.append(add(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    a,b,c,d,e,f,g,h = H0
    for r in range(64):
        T1=add(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add(Sig0(a),Maj(a,b,c))
        h,g,f,e,d,c,b,a = g,f,e,add(d,T1),c,b,a,add(T1,T2)
    return tuple(add(s,iv) for s,iv in zip((a,b,c,d,e,f,g,h), H0))

def hw256(ha, hb):
    """Hamming weight of ha XOR hb (256-bit)."""
    total = 0
    for i in range(8):
        total += bin((ha[i] ^ hb[i]) & MASK).count('1')
    return total

def low_k_key(h, k):
    """Extract lower k bits of each of 8 words as a key."""
    mask = (1 << k) - 1
    return tuple(w & mask for w in h)

def main():
    print("="*70)
    print("EXPERIMENT D: Low-bit partial collisions")
    print("="*70)
    t0 = time.time()

    for k in [1, 2, 4, 8]:
        print(f"\n  k={k}: match lower {k} bits of each word ({8*k} bits total)")
        print(f"    Birthday bound: 2^{4*k} = {2**(4*k)}")

        # Build hash table
        table = defaultdict(list)  # key → list of (M, H)
        n_needed = min(2**(4*k+1), 2_000_000)  # 2x birthday bound
        collisions = []
        batch = 1000

        generated = 0
        while generated < n_needed and len(collisions) < 200:
            for _ in range(batch):
                M = tuple(rw(16))
                H = sha256_hash(M)
                key = low_k_key(H, k)
                if key in table:
                    # Found partial collision
                    for M_prev, H_prev in table[key]:
                        if M != M_prev:
                            collisions.append((H_prev, H))
                            if len(collisions) >= 200:
                                break
                    if len(collisions) >= 200:
                        break
                table[key].append((M, H))
                generated += 1

        if not collisions:
            print(f"    No partial collisions found in {generated} hashes!")
            continue

        # Analyze collisions
        hws = [hw256(ha, hb) for ha, hb in collisions]
        expected_remaining = 256 - 8*k
        # Remaining bits: (256 - 8k) bits, each with P=0.5 of differing
        expected_hw = expected_remaining / 2
        mean_hw = sum(hws) / len(hws)
        min_hw = min(hws)
        max_hw = max(hws)

        print(f"    Found {len(collisions)} partial collisions in {generated} hashes")
        print(f"    Matched bits: {8*k}, remaining: {expected_remaining}")
        print(f"    Expected HW of remaining XOR: {expected_hw:.1f}")
        print(f"    Actual mean HW(H_a XOR H_b): {mean_hw:.1f}")
        print(f"    Min HW: {min_hw}, Max HW: {max_hw}")
        deviation = mean_hw - expected_hw
        print(f"    Deviation from expected: {deviation:+.2f}")
        if abs(deviation) > 2:
            print(f"    *** SIGNIFICANT DEVIATION — needs investigation ***")
        else:
            print(f"    Consistent with random (no correlation low→high bits)")

        # Extra: check if bits just above k are more likely to match
        if k <= 4:
            above_match = 0
            total_checked = 0
            for ha, hb in collisions:
                for i in range(8):
                    # Check bit k (just above matched region)
                    bit_a = (ha[i] >> k) & 1
                    bit_b = (hb[i] >> k) & 1
                    if bit_a == bit_b:
                        above_match += 1
                    total_checked += 1
            p_above = above_match / total_checked if total_checked else 0
            print(f"    P(bit {k} matches | bits 0..{k-1} match) = "
                  f"{above_match}/{total_checked} = {p_above:.4f} (expected 0.5)")

    elapsed = time.time() - t0
    print(f"\n  Total time: {elapsed:.1f}s")

if __name__ == "__main__":
    main()
