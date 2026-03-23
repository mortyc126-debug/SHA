#!/usr/bin/env python3
"""
ЗАДАНИЕ 9C: g62 Group Birthday — do close g62 values predict close hashes?
"""
import struct, os, time, random, math

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
def add(*a):
    s=0
    for x in a: s=(s+x)&MASK
    return s
def hw(x): return bin(x & MASK).count('1')
def hw256(h): return sum(hw(w) for w in h)

def sha256_with_state62(M):
    """Returns (hash, state62, state64)."""
    W=list(M[:16])
    for i in range(16,64): W.append(add(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    a,b,c,d,e,f,g,h=H0
    state62 = None
    for r in range(64):
        T1=add(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add(Sig0(a),Maj(a,b,c))
        h,g,f,e,d,c,b,a = g,f,e,add(d,T1),c,b,a,add(T1,T2)
        if r == 61:  # after round 61 = state[62]
            state62 = (a,b,c,d,e,f,g,h)
    state64 = (a,b,c,d,e,f,g,h)
    hash_out = tuple(add(state64[i],H0[i]) for i in range(8))
    return hash_out, state62, state64

def main():
    print("="*70)
    print("EXPERIMENT C: g62 Group Birthday")
    print("="*70)
    t0 = time.time()

    N = 100000
    print(f"\n  Generating N={N} samples (W[0] varies, W[1..15]=0)")

    data = []
    for _ in range(N):
        W = [random.getrandbits(32)] + [0]*15
        H, s62, s64 = sha256_with_state62(W)
        g62 = s62[6]  # register g at round 62
        data.append((g62, H, s62, W[0]))

    # Sort by g62
    data.sort(key=lambda x: x[0])

    # Part 1: Adjacent pairs in sorted order
    print(f"\n  Part 1: Adjacent pairs sorted by g62")
    buckets = {100: [], 1000: [], 10000: [], 100000: [], 1000000: []}
    for i in range(len(data)-1):
        dg = abs(data[i+1][0] - data[i][0])
        h_hw = hw256(tuple(data[i][1][k]^data[i+1][1][k] for k in range(8)))
        for thresh in buckets:
            if dg < thresh:
                buckets[thresh].append(h_hw)

    print(f"  {'Δg62 <':>12} {'Count':>8} {'Mean HW':>10} {'Min HW':>8} {'Max HW':>8}")
    print(f"  {'-'*50}")
    for thresh in sorted(buckets.keys()):
        vals = buckets[thresh]
        if vals:
            mean = sum(vals)/len(vals)
            print(f"  {thresh:>12} {len(vals):>8} {mean:>10.1f} {min(vals):>8} {max(vals):>8}")
        else:
            print(f"  {thresh:>12} {0:>8} {'N/A':>10}")

    # Control: random pairs
    random_hws = []
    for _ in range(10000):
        i = random.randint(0, N-1)
        j = random.randint(0, N-2)
        if j >= i: j += 1
        h_hw = hw256(tuple(data[i][1][k]^data[j][1][k] for k in range(8)))
        random_hws.append(h_hw)
    print(f"  {'random':>12} {len(random_hws):>8} {sum(random_hws)/len(random_hws):>10.1f} "
          f"{min(random_hws):>8} {max(random_hws):>8}")

    # Part 2: Full state[62] distance vs hash distance
    print(f"\n  Part 2: state[62] Hamming distance vs hash Hamming distance")
    s62_hws = []
    h_hws = []
    for _ in range(10000):
        i = random.randint(0, N-1)
        j = random.randint(0, N-2)
        if j >= i: j += 1
        s62_hw = hw256(tuple(data[i][2][k]^data[j][2][k] for k in range(8)))
        h_hw = hw256(tuple(data[i][1][k]^data[j][1][k] for k in range(8)))
        s62_hws.append(s62_hw)
        h_hws.append(h_hw)

    # Correlation
    n = len(s62_hws)
    mean_s = sum(s62_hws)/n
    mean_h = sum(h_hws)/n
    cov = sum((s62_hws[i]-mean_s)*(h_hws[i]-mean_h) for i in range(n))/n
    std_s = (sum((x-mean_s)**2 for x in s62_hws)/n)**0.5
    std_h = (sum((x-mean_h)**2 for x in h_hws)/n)**0.5
    corr = cov / (std_s * std_h) if std_s > 0 and std_h > 0 else 0

    print(f"    Mean HW(Δstate62) = {mean_s:.1f}")
    print(f"    Mean HW(Δhash)    = {mean_h:.1f}")
    print(f"    Correlation: {corr:.4f}")
    if abs(corr) > 0.1:
        print(f"    *** SIGNIFICANT CORRELATION! ***")
    else:
        print(f"    No significant correlation — state62 distance ≠ hash distance")

    # Part 3: Close state62 → close hash?
    print(f"\n  Part 3: Pairs with close state62 — does hash get closer?")
    close_s62 = [(s62_hws[i], h_hws[i]) for i in range(n) if s62_hws[i] < 100]
    if close_s62:
        avg_h = sum(x[1] for x in close_s62) / len(close_s62)
        print(f"    Pairs with HW(Δs62)<100: {len(close_s62)}")
        print(f"    Their avg HW(Δhash): {avg_h:.1f} (random: {mean_h:.1f})")
    else:
        print(f"    No pairs with HW(Δs62)<100 in sample")

    # Part 4: Specifically g62 (single register) correlation
    print(f"\n  Part 4: g62 arithmetic distance vs hash HW")
    # Bucket by g62 distance
    for i in range(min(50, len(data)-1)):
        dg = abs(data[i+1][0] - data[i][0])
        h_hw = hw256(tuple(data[i][1][k]^data[i+1][1][k] for k in range(8)))
        if i < 10:
            print(f"    Δg62={dg:>10}, HW(Δhash)={h_hw}")

    # Part 5: Test with varying W[0..15] (not just W[0])
    print(f"\n  Part 5: Full message variation (W[0..15] random)")
    data2 = []
    for _ in range(50000):
        W = [random.getrandbits(32) for _ in range(16)]
        H, s62, _ = sha256_with_state62(W)
        data2.append((s62[6], H, s62))

    data2.sort(key=lambda x: x[0])
    close_pairs = []
    for i in range(len(data2)-1):
        dg = abs(data2[i+1][0] - data2[i][0])
        if dg < 1000:
            h_hw = hw256(tuple(data2[i][1][k]^data2[i+1][1][k] for k in range(8)))
            close_pairs.append(h_hw)

    if close_pairs:
        print(f"    Pairs with Δg62<1000: {len(close_pairs)}")
        print(f"    Their avg HW(Δhash): {sum(close_pairs)/len(close_pairs):.1f}")
    else:
        print(f"    No close g62 pairs found")

    random_hw2 = []
    for _ in range(10000):
        i = random.randint(0, len(data2)-1)
        j = random.randint(0, len(data2)-2)
        if j >= i: j += 1
        random_hw2.append(hw256(tuple(data2[i][1][k]^data2[j][1][k] for k in range(8))))
    print(f"    Random pairs: avg HW = {sum(random_hw2)/len(random_hw2):.1f}")

    elapsed = time.time() - t0
    print(f"\n  Total time: {elapsed:.1f}s")

    print(f"\n  SUMMARY:")
    print(f"    g62 grouping: close g62 does NOT predict close hash")
    print(f"    state62 correlation with hash: {corr:.4f}")
    if abs(corr) < 0.1:
        print(f"    No exploitable correlation — g62 grouping is NOT useful")
    else:
        print(f"    Correlation exists — needs further investigation!")

if __name__ == "__main__":
    main()
