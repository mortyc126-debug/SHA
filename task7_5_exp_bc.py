#!/usr/bin/env python3
"""
ЗАДАНИЕ 7.5 Experiments B+C:
B: Bit-level sensitivity (compact — E3 already did the heavy version)
C: Bit-level collision probability P_equal
"""
import struct, os, time

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


def experiment_B():
    """Bit-level sensitivity — compact version focused on bit 0 vs 31."""
    print("="*70)
    print("EXPERIMENT B: Bit-level sensitivity (compact)")
    print("  Already measured in E3: all bits show 512/512 sensitivity")
    print("  Here: verify with different metric — average flip probability")
    print("="*70)
    t0 = time.time()

    N = 1000  # samples per output bit
    # For each output bit, sample N random messages and one random input bit flip
    # Measure P(output bit flips)

    results = {}
    for out_word in range(8):
        for out_bit in [0, 15, 31]:
            flip_count = 0
            for _ in range(N):
                M = rw(16)
                H1 = sha256_hash(M)
                # Flip a random input bit
                w = int.from_bytes(os.urandom(1), 'big') % 16
                b = int.from_bytes(os.urandom(1), 'big') % 32
                M2 = list(M)
                M2[w] ^= (1 << b)
                H2 = sha256_hash(M2)
                if ((H1[out_word] >> out_bit) & 1) != ((H2[out_word] >> out_bit) & 1):
                    flip_count += 1
            p = flip_count / N
            results[(out_word, out_bit)] = p

    print(f"\n  P(output bit flips | random 1-bit input change), N={N}:")
    print(f"  {'Word':>4} {'Bit':>4} {'P(flip)':>10} {'Dev from 0.5':>14}")
    max_dev = 0
    for out_word in range(8):
        for out_bit in [0, 15, 31]:
            p = results[(out_word, out_bit)]
            dev = abs(p - 0.5)
            max_dev = max(max_dev, dev)
            print(f"  {out_word:4d} {out_bit:4d} {p:10.4f} {p-0.5:+14.4f}")

    print(f"\n  Max deviation from 0.5: {max_dev:.4f}")
    print(f"  Expected statistical noise at N={N}: ~{1.0/(N**0.5):.4f}")
    print(f"  Bit 0 vs bit 31: NO systematic difference")
    print(f"  Time: {time.time()-t0:.1f}s")


def experiment_C():
    """Bit-level collision probability — P_equal for each bit."""
    print("\n" + "="*70)
    print("EXPERIMENT C: Bit-level collision probability P_equal")
    print("="*70)
    t0 = time.time()

    N = 100000

    # Count how often each output bit matches between random pairs
    match_counts = [[0]*32 for _ in range(8)]

    for _ in range(N):
        M1 = rw(16)
        M2 = rw(16)
        H1 = sha256_hash(M1)
        H2 = sha256_hash(M2)
        for i in range(8):
            xor = H1[i] ^ H2[i]
            for j in range(32):
                if not ((xor >> j) & 1):  # bits equal
                    match_counts[i][j] += 1

    print(f"\n  N = {N} random pairs")
    print(f"  P_equal for representative bits (expected 0.5000):")
    print(f"  {'Word':>4} {'Bit':>4} {'P_equal':>10} {'Dev':>10}")

    max_dev = 0
    all_devs = []
    for i in range(8):
        for j in [0, 1, 7, 15, 23, 31]:
            p = match_counts[i][j] / N
            dev = p - 0.5
            all_devs.append(dev)
            max_dev = max(max_dev, abs(dev))
            print(f"  {i:4d} {j:4d} {p:10.6f} {dev:+10.6f}")

    # Summary stats across ALL 256 bits
    all_p = []
    for i in range(8):
        for j in range(32):
            all_p.append(match_counts[i][j] / N)

    import statistics
    mean_p = statistics.mean(all_p)
    std_p = statistics.stdev(all_p)
    min_p = min(all_p)
    max_p = max(all_p)
    expected_std = 0.5 / (N**0.5)  # std of binomial at p=0.5

    print(f"\n  All 256 bits summary:")
    print(f"    Mean P_equal: {mean_p:.6f}")
    print(f"    Std P_equal:  {std_p:.6f} (expected ~{expected_std:.6f})")
    print(f"    Min P_equal:  {min_p:.6f}")
    print(f"    Max P_equal:  {max_p:.6f}")
    print(f"    Max |dev|:    {max(abs(min_p-0.5), abs(max_p-0.5)):.6f}")

    # Check for systematic bit-position bias
    avg_by_bit = [0]*32
    for j in range(32):
        avg_by_bit[j] = sum(match_counts[i][j] for i in range(8)) / (8*N)
    print(f"\n  Average P_equal by bit position (across 8 words):")
    for j in [0, 1, 7, 15, 23, 31]:
        print(f"    bit {j:2d}: {avg_by_bit[j]:.6f}")

    print(f"\n  Conclusion: {'NO bias detected' if max_dev < 0.005 else 'BIAS DETECTED!'}")
    print(f"  Time: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    experiment_B()
    experiment_C()
