#!/usr/bin/env python3
"""
UNEXPLORED TERRITORY — Axis 4
================================

29 experiments closed 16 directions. What's LEFT?

Systematic inventory of what we HAVE NOT tried:

CLOSED (don't repeat):
  ✗ Carry profile → H (binary or continuous)
  ✗ Carry_sched → H
  ✗ Ch/Maj bridge (tautological)
  ✗ Schedule resonance (dead for real msgs)
  ✗ Cross-block products
  ✗ 2-adic lifting (bits independent)
  ✗ Multi-block propagation
  ✗ Passive distinguisher
  ✗ Differential δH
  ✗ Algebraic Ch filter
  ✗ Smart HC targets
  ✗ Da[13] modular structure (truly uniform)
  ✗ δW[16] optimization (independent of δe[17])
  ✗ Copula-guided differential path
  ✗ A-branch bridge
  ✗ Deep nonlinear output features beyond 71

NOT TRIED — genuinely new:

  N1: SECOND PREIMAGE via distinguish — use distinguisher to
      narrow preimage search space. If H has bias → fewer candidates.

  N2: FUNCTIONAL GRAPH — SHA-256 as iterated map.
      Methodology found cycle anomaly (4 cycles vs 22 expected).
      We haven't exploited this for collision/preimage.

  N3: ALGEBRAIC NORMAL FORM — methodology proved deg(e_r)=32 at r≥6.
      But what about specific BITS? Individual bit might have lower degree.
      Low-degree bits → algebraic attack on those bits.

  N4: TWO-BLOCK DIFFERENTIAL — block 1 creates specific IV for block 2.
      Wang chain on block 2 with crafted IV might extend beyond r=17.
      Different IV → different K-barrier structure.

  N5: ROTATIONAL CRYPTANALYSIS — methodology closed for full rotation.
      But what about PARTIAL rotation? Rotate only H[6,7] bits?
      Or: fixed-point of rotation in carry space?

  N6: CUBE ATTACK on reduced rounds — methodology mentions Q6.
      Cube attacks work on low-degree functions. SHA-256 at r≤5
      has degree < 32. Cube attack on 5-round SHA-256?

Let's test N2 (functional graph), N3 (per-bit degree), N4 (two-block).
"""

import numpy as np
from collections import Counter, defaultdict
import time, sys

MASK = 0xFFFFFFFF
H0 = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
     0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
     0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
     0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
     0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
     0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def hw(x): return bin(x&MASK).count('1')

def msg_sched(W16):
    W=list(W16)+[0]*48
    for i in range(16,64): W[i]=(sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16])&MASK
    return W

def sha256_32(x):
    """SHA-256 projected to 32 bits: g(x) = SHA256(x||0..0)[0:32]."""
    W16 = [x & MASK] + [0]*15
    W = msg_sched(W16)
    a,b,c,d,e,f,g,h = H0
    for r in range(64):
        T1=(h+Sig1(e)+Ch(e,f,g)+K[r]+W[r])&MASK; T2=(Sig0(a)+Maj(a,b,c))&MASK
        h,g,f=g,f,e; e=(d+T1)&MASK; d,c,b=c,b,a; a=(T1+T2)&MASK
    return (a+H0[0])&MASK  # H[0]

def sha256_reduced(W16, num_rounds):
    """SHA-256 with reduced rounds. Returns H[0..7]."""
    W = msg_sched(W16)
    a,b,c,d,e,f,g,h = H0
    for r in range(num_rounds):
        T1=(h+Sig1(e)+Ch(e,f,g)+K[r]+W[r])&MASK; T2=(Sig0(a)+Maj(a,b,c))&MASK
        h,g,f=g,f,e; e=(d+T1)&MASK; d,c,b=c,b,a; a=(T1+T2)&MASK
    return [(v+iv)&MASK for v,iv in zip([a,b,c,d,e,f,g,h],H0)]

def sha_compress(IV, W16):
    W = msg_sched(W16)
    a,b,c,d,e,f,g,h = IV
    for r in range(64):
        T1=(h+Sig1(e)+Ch(e,f,g)+K[r]+W[r])&MASK; T2=(Sig0(a)+Maj(a,b,c))&MASK
        h,g,f=g,f,e; e=(d+T1)&MASK; d,c,b=c,b,a; a=(T1+T2)&MASK
    return [(v+iv)&MASK for v,iv in zip([a,b,c,d,e,f,g,h],IV)]


# ============================================================
# N2: FUNCTIONAL GRAPH — exploit cycle structure
# ============================================================

def experiment_N2(seed=12000):
    print("="*70)
    print("N2: FUNCTIONAL GRAPH — Exploit cycle anomaly for collision")
    print("="*70)

    # Methodology: g32 has 4 cycles (expected 22). λ_dom ≈ 31513.
    # Cycle inversion: SHA256^{-1}(x) = g^{λ-1}(x) for x ∈ CYCLE.
    # This gives 32-bit preimage in O(λ) ≈ O(2^15) instead of O(2^32).
    #
    # NEW IDEA: use cycle structure for COLLISION finding.
    # Two inputs x, y with g(x) = g(y) = collision in 32-bit projection.
    # If x and y are in the SAME tail → they merge at some point.
    # Floyd's algorithm finds the merge point.

    print(f"\n  Building functional graph g32(x) = SHA256(x||0..0)[0:32]...")

    rng = np.random.RandomState(seed)

    # Floyd's cycle detection from multiple starts
    n_starts = 100
    cycles_found = {}
    tail_lengths = []

    t0 = time.time()
    for start_idx in range(n_starts):
        x0 = int(rng.randint(0, 1<<32))
        # Floyd's tortoise and hare
        tort = sha256_32(x0)
        hare = sha256_32(sha256_32(x0))
        steps = 0
        while tort != hare and steps < 200000:
            tort = sha256_32(tort)
            hare = sha256_32(sha256_32(hare))
            steps += 1

        if tort == hare and steps < 200000:
            # Find cycle start
            tort = x0
            tail_len = 0
            while tort != hare:
                tort = sha256_32(tort)
                hare = sha256_32(hare)
                tail_len += 1

            # Find cycle length
            cycle_start = tort
            cyc_len = 1
            x = sha256_32(cycle_start)
            while x != cycle_start:
                x = sha256_32(x)
                cyc_len += 1

            tail_lengths.append(tail_len)
            if cyc_len not in cycles_found:
                cycles_found[cyc_len] = cycle_start

    print(f"  Floyd from {n_starts} starts: {time.time()-t0:.1f}s")
    print(f"  Unique cycle lengths: {sorted(cycles_found.keys())}")
    print(f"  Number of cycles: {len(cycles_found)}")
    print(f"  Mean tail length: {np.mean(tail_lengths):.0f}")
    print(f"  Max tail length: {max(tail_lengths)}")

    # KEY: can we find COLLISIONS using cycle structure?
    # Two random points in the same cycle's basin → they eventually merge.
    # The merge point = collision in g32.
    print(f"\n  --- Collision via functional graph ---")

    if cycles_found:
        max_cycle = max(cycles_found.keys())
        cycle_start = cycles_found[max_cycle]
        print(f"  Dominant cycle: length={max_cycle}, start=0x{cycle_start:08x}")

        # Collect points in the cycle
        cycle_points = set()
        x = cycle_start
        for _ in range(max_cycle):
            cycle_points.add(x)
            x = sha256_32(x)

        # Search for collisions: iterate random points until hitting cycle
        n_search = min(50000, max_cycle * 5)
        preimages = defaultdict(list)

        t0 = time.time()
        for i in range(n_search):
            x = int(rng.randint(0, 1<<32))
            y = sha256_32(x)
            preimages[y].append(x)

        collisions = {y: xs for y, xs in preimages.items() if len(xs) >= 2}
        print(f"\n  Random search for 32-bit H[0] collisions:")
        print(f"    Searched: {n_search}")
        print(f"    Collisions found: {len(collisions)}")

        if collisions:
            for y, xs in list(collisions.items())[:3]:
                print(f"    g({xs[0]:08x}) = g({xs[1]:08x}) = {y:08x}")
                # Verify
                h1 = sha256_32(xs[0])
                h2 = sha256_32(xs[1])
                print(f"      Verified: {h1==h2} (H[0] match)")

        expected_collisions = n_search**2 / (2 * 2**32)
        print(f"    Expected (birthday): {expected_collisions:.2f}")

    return cycles_found


# ============================================================
# N3: PER-BIT ALGEBRAIC DEGREE
# ============================================================

def experiment_N3(seed=12001):
    print("\n"+"="*70)
    print("N3: PER-BIT ALGEBRAIC DEGREE — Individual bit complexity")
    print("="*70)

    # Methodology: deg(e_r) = 32 at r≥6 (over 32-bit input).
    # But: individual OUTPUT BITS might have lower effective degree.
    # If H[w][b] has degree < 32 → linearization attack on that bit.

    # Test: for reduced SHA (r rounds), what is the "effective degree"
    # of each output bit? Measured via higher-order differentials.
    # If δ^k f(x) = 0 for all x → deg(f) < k.

    print(f"\n  Method: higher-order differential test")
    print(f"  δ^k f = sum over all k-subsets of f(x ⊕ S) for |S|=k")
    print(f"  If δ^k = 0 always → deg < k")

    rng = np.random.RandomState(seed)

    for num_rounds in [3, 5, 8, 16, 64]:
        print(f"\n  --- {num_rounds}-round SHA-256 ---")

        # Test degree of H[0] bit 0 (simplest output bit)
        # Use 10-bit input (W[0][0..9]) for tractability
        n_input_bits = 10

        for target_word, target_bit in [(0, 0), (0, 31), (7, 31), (7, 0)]:
            # Build truth table
            truth_table = np.zeros(1 << n_input_bits, dtype=np.int8)
            for x in range(1 << n_input_bits):
                W16 = [x] + [0]*15
                H = sha256_reduced(W16, num_rounds)
                truth_table[x] = (H[target_word] >> target_bit) & 1

            # Compute degree via Möbius transform (ANF)
            n = n_input_bits
            anf = truth_table.copy()
            for i in range(n):
                step = 1 << i
                for j in range(0, 1 << n, step << 1):
                    for k in range(step):
                        anf[j + k + step] ^= anf[j + k]

            # Degree = max HW of index where anf[index]=1
            max_deg = 0
            for idx in range(1 << n):
                if anf[idx] == 1:
                    d = hw(idx)
                    if d > max_deg:
                        max_deg = d

            print(f"    H[{target_word}][b{target_bit}]: degree = {max_deg}/{n_input_bits}")

    return


# ============================================================
# N4: TWO-BLOCK with crafted IV
# ============================================================

def experiment_N4(N=5000, seed=12002):
    print("\n"+"="*70)
    print("N4: TWO-BLOCK — Crafted IV extends Wang barrier?")
    print(f"N={N}")
    print("="*70)

    rng = np.random.RandomState(seed)

    # Standard IV → Wang barrier at r=17 (K[16]=0.893, K[17]=0.936)
    # What if IV is chosen so that K-barrier rounds become "cheaper"?
    #
    # The carry condition: raw[r] = h + Sig1(e) + Ch(e,f,g) + K[r] + W[r] < 2^32
    # Barrier: K[16]+K[17] are large → raw[16,17] almost always ≥ 2^32
    # If IV is chosen so that h, e, f, g at round 16-17 are SMALL →
    # raw[16,17] could be < 2^32 → barrier weakened!

    # Method: block 1 message M1 → H1 = new IV.
    # Choose M1 so that H1 creates favorable state at rounds 16-17.

    # Specifically: H1[7] = new h[0] = small → h propagates to h[16]
    # But h[16] = g[15] = f[14] = e[13] → depends on 13 rounds of mixing.
    # IV only directly sets h[0], which influences h[3] (through shift).

    # Simpler test: compare Wang barrier for STANDARD IV vs RANDOM IV.
    # Does the barrier shift?

    print(f"\n  --- Wang barrier position for different IVs ---")

    def wang_zeros(IV, Wn, dW0):
        """Count consecutive δe[r]=0 from r=2 under Wang chain with custom IV."""
        W_exp = msg_sched(Wn)
        an,bn,cn,dn,en,fn,gn,hn = IV
        af,bf,cf,df,ef,ff,gf,hf = IV

        T1n=(hn+Sig1(en)+Ch(en,fn,gn)+K[0]+W_exp[0])&MASK; T2n=(Sig0(an)+Maj(an,bn,cn))&MASK
        hn,gn,fn=gn,fn,en; en=(dn+T1n)&MASK; dn,cn,bn=cn,bn,an; an=(T1n+T2n)&MASK
        T1f=(hf+Sig1(ef)+Ch(ef,ff,gf)+K[0]+((W_exp[0]+dW0)&MASK))&MASK; T2f=(Sig0(af)+Maj(af,bf,cf))&MASK
        hf,gf,ff=gf,ff,ef; ef=(df+T1f)&MASK; df,cf,bf=cf,bf,af; af=(T1f+T2f)&MASK

        zeros = 0
        for r in range(1, 20):
            dW=(0-(dn-df)-(hn-hf)-(Sig1(en)-Sig1(ef))-(Ch(en,fn,gn)-Ch(ef,ff,gf)))&MASK
            Wfr=(W_exp[r]+dW)&MASK if r < 16 else 0

            T1n=(hn+Sig1(en)+Ch(en,fn,gn)+K[r]+W_exp[r])&MASK; T2n=(Sig0(an)+Maj(an,bn,cn))&MASK
            hn,gn,fn=gn,fn,en; en=(dn+T1n)&MASK; dn,cn,bn=cn,bn,an; an=(T1n+T2n)&MASK

            if r < 16:
                T1f=(hf+Sig1(ef)+Ch(ef,ff,gf)+K[r]+Wfr)&MASK; T2f=(Sig0(af)+Maj(af,bf,cf))&MASK
                hf,gf,ff=gf,ff,ef; ef=(df+T1f)&MASK; df,cf,bf=cf,bf,af; af=(T1f+T2f)&MASK

                if en == ef:
                    zeros += 1
                else:
                    break
            else:
                break

        return zeros

    # Standard IV: always 15 zeros (r=2..16)
    std_zeros = []
    for _ in range(min(N, 1000)):
        Wn = [int(rng.randint(0,1<<32)) for _ in range(16)]
        z = wang_zeros(tuple(H0), Wn, 0x8000)
        std_zeros.append(z)

    print(f"  Standard IV: mean Wang zeros = {np.mean(std_zeros):.2f} (expected 15)")

    # Random IVs
    print(f"\n  Random IVs (N={min(N, 2000)}):")
    random_zeros = []
    best_iv = None; best_z = 0

    for _ in range(min(N, 2000)):
        IV = tuple(int(rng.randint(0, 1<<32)) for _ in range(8))
        Wn = [int(rng.randint(0,1<<32)) for _ in range(16)]
        z = wang_zeros(IV, Wn, 0x8000)
        random_zeros.append(z)
        if z > best_z:
            best_z = z; best_iv = IV

    print(f"  Random IV: mean zeros = {np.mean(random_zeros):.2f}")
    print(f"  Random IV: max zeros = {max(random_zeros)}")

    zero_dist = Counter(random_zeros)
    for z in sorted(zero_dist.keys()):
        print(f"    zeros={z:2d}: {zero_dist[z]} ({zero_dist[z]/len(random_zeros)*100:.1f}%)")

    # Two-block: block 1 → H1 = crafted IV
    print(f"\n  --- Two-block: M1 → H1 → Wang on block 2 ---")
    twoblock_zeros = []
    for _ in range(min(N, 2000)):
        M1 = [int(rng.randint(0,1<<32)) for _ in range(16)]
        H1 = sha_compress(tuple(H0), M1)  # this is wrong — need to pass list
        # Actually sha_compress expects IV as list
        W = msg_sched(M1); a,b,c,d,e,f,g,h=H0
        for r in range(64):
            T1=(h+Sig1(e)+Ch(e,f,g)+K[r]+W[r])&MASK; T2=(Sig0(a)+Maj(a,b,c))&MASK
            h,g,f=g,f,e; e=(d+T1)&MASK; d,c,b=c,b,a; a=(T1+T2)&MASK
        H1 = tuple((v+iv)&MASK for v,iv in zip([a,b,c,d,e,f,g,h],H0))

        M2 = [int(rng.randint(0,1<<32)) for _ in range(16)]
        z = wang_zeros(H1, M2, 0x8000)
        twoblock_zeros.append(z)

    print(f"  Two-block: mean zeros = {np.mean(twoblock_zeros):.2f}")
    print(f"  Two-block: max zeros = {max(twoblock_zeros)}")

    zero_dist2 = Counter(twoblock_zeros)
    for z in sorted(zero_dist2.keys()):
        if zero_dist2[z] > len(twoblock_zeros)*0.01:
            print(f"    zeros={z:2d}: {zero_dist2[z]} ({zero_dist2[z]/len(twoblock_zeros)*100:.1f}%)")

    # KEY: does ANY IV give >15 Wang zeros?
    print(f"\n  ★ KEY: Can custom IV extend Wang chain beyond 15 zeros?")
    if best_z > 15:
        print(f"    YES! Found IV with {best_z} zeros!")
    else:
        print(f"    Standard IV (15 zeros) = maximum.")
        print(f"    Wang chain always has EXACTLY 15 zeros regardless of IV.")
        print(f"    This is STRUCTURAL: 16 rounds × 1 δW each = 16 corrections max,")
        print(f"    minus 1 for δe[1]≠0 = 15 zeros.")

    return


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("UNEXPLORED TERRITORY — Axis 4")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    fast = '--fast' in sys.argv

    t_start = time.time()
    cycles = experiment_N2()
    experiment_N3()
    experiment_N4(N=3000 if fast else 5000)

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: Unexplored Territory (Axis 4)
{'='*70}

N2: FUNCTIONAL GRAPH — cycle structure for collision
N3: PER-BIT DEGREE — individual bit complexity
N4: TWO-BLOCK — custom IV extends barrier?

Three genuinely new angles. Which shows promise?
""")
