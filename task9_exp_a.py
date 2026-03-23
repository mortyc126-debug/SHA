#!/usr/bin/env python3
"""
ЗАДАНИЕ 9A: SA Descent Beyond hw63=87
Minimize HW(SHA256(Wn) XOR SHA256(Wf)) where Wf = Wang pair for Wn.
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

def sha256_hash(M):
    W=list(M[:16])
    for i in range(16,64): W.append(add(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    a,b,c,d,e,f,g,h=H0
    for r in range(64):
        T1=add(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add(Sig0(a),Maj(a,b,c))
        h,g,f,e,d,c,b,a = g,f,e,add(d,T1),c,b,a,add(T1,T2)
    return tuple(add(s,iv) for s,iv in zip((a,b,c,d,e,f,g,h),H0))

def wang_pair(Wn):
    """Generate Wang pair with δW[0]=0x8000, adaptive δW[1..15]."""
    Wf = list(Wn)
    Wf[0] = Wn[0] ^ 0x8000
    sn = list(H0); sf = list(H0)
    T1n=add(sn[7],Sig1(sn[4]),Ch(sn[4],sn[5],sn[6]),K[0],Wn[0])
    T2n=add(Sig0(sn[0]),Maj(sn[0],sn[1],sn[2]))
    T1f=add(sf[7],Sig1(sf[4]),Ch(sf[4],sf[5],sf[6]),K[0],Wf[0])
    T2f=add(Sig0(sf[0]),Maj(sf[0],sf[1],sf[2]))
    sn=[add(T1n,T2n),sn[0],sn[1],sn[2],add(sn[3],T1n),sn[4],sn[5],sn[6]]
    sf=[add(T1f,T2f),sf[0],sf[1],sf[2],add(sf[3],T1f),sf[4],sf[5],sf[6]]
    for r in range(1,16):
        dd=(sf[3]-sn[3])&MASK; dh=(sf[7]-sn[7])&MASK
        dS=(Sig1(sf[4])-Sig1(sn[4]))&MASK; dC=(Ch(sf[4],sf[5],sf[6])-Ch(sn[4],sn[5],sn[6]))&MASK
        dWr=(-((dd+dh+dS+dC)&MASK))&MASK
        Wf[r]=add(Wn[r],dWr)
        T1n=add(sn[7],Sig1(sn[4]),Ch(sn[4],sn[5],sn[6]),K[r],Wn[r])
        T2n=add(Sig0(sn[0]),Maj(sn[0],sn[1],sn[2]))
        T1f=add(sf[7],Sig1(sf[4]),Ch(sf[4],sf[5],sf[6]),K[r],Wf[r])
        T2f=add(Sig0(sf[0]),Maj(sf[0],sf[1],sf[2]))
        sn=[add(T1n,T2n),sn[0],sn[1],sn[2],add(sn[3],T1n),sn[4],sn[5],sn[6]]
        sf=[add(T1f,T2f),sf[0],sf[1],sf[2],add(sf[3],T1f),sf[4],sf[5],sf[6]]
    return Wf

def evaluate(Wn):
    """Compute hw63 = HW(SHA256(Wn) XOR SHA256(Wf))."""
    Wf = wang_pair(Wn)
    Hn = sha256_hash(Wn)
    Hf = sha256_hash(Wf)
    return hw256(tuple(Hn[i]^Hf[i] for i in range(8))), Hn, Hf, Wf

def sa_search(steps=500000, restarts=100):
    print("="*70)
    print("EXPERIMENT A: SA Descent — minimize HW(δhash) via Wang pairs")
    print(f"  Steps={steps}, restarts={restarts}")
    print("="*70)
    t0 = time.time()

    global_best_hw = 999
    global_best_Wn = None
    all_bests = []

    # Population for crossover
    population = []  # (hw, Wn)

    for restart in range(restarts):
        # Initialize: random or crossover from population
        if len(population) >= 2 and random.random() < 0.3:
            p1 = random.choice(population[:16])
            p2 = random.choice(population[:16])
            Wn = list(p1[1][:8]) + list(p2[1][8:16])
            # Small perturbation
            for _ in range(3):
                idx = random.randint(0, 15)
                Wn[idx] ^= (1 << random.randint(0, 31))
        else:
            Wn = [random.getrandbits(32) for _ in range(16)]

        cur_hw, _, _, _ = evaluate(Wn)
        best_hw = cur_hw
        best_Wn = list(Wn)

        T = 50.0
        T_end = 0.01
        alpha = (T_end / T) ** (1.0 / steps)

        for step in range(steps):
            # Mutate: flip 1-3 random bits
            Wn_new = list(Wn)
            n_flips = 1 if random.random() < 0.7 else random.randint(2, 3)
            for _ in range(n_flips):
                widx = random.randint(0, 15)
                bidx = random.randint(0, 31)
                Wn_new[widx] ^= (1 << bidx)

            new_hw, _, _, _ = evaluate(Wn_new)
            delta = new_hw - cur_hw

            if delta < 0 or random.random() < math.exp(-delta / max(T, 0.001)):
                Wn = Wn_new
                cur_hw = new_hw
                if cur_hw < best_hw:
                    best_hw = cur_hw
                    best_Wn = list(Wn)

            T *= alpha

        all_bests.append(best_hw)
        population.append((best_hw, best_Wn))
        population.sort(key=lambda x: x[0])
        if len(population) > 32:
            population = population[:32]

        if best_hw < global_best_hw:
            global_best_hw = best_hw
            global_best_Wn = list(best_Wn)
            print(f"  Restart {restart:3d}: hw={best_hw} *** NEW GLOBAL BEST ***")
        elif restart % 20 == 0:
            print(f"  Restart {restart:3d}: hw={best_hw} (global best={global_best_hw})")

    elapsed = time.time() - t0
    print(f"\n  Total time: {elapsed:.1f}s")

    # Results
    print(f"\n  RESULTS:")
    print(f"    Global best hw63 = {global_best_hw}")
    print(f"    Reference (manual): hw63 = 87")
    print(f"    Distribution of best per restart:")
    from collections import Counter
    c = Counter(all_bests)
    for h in sorted(c.keys())[:20]:
        print(f"      hw={h}: {c[h]} restarts")
    print(f"    Mean best: {sum(all_bests)/len(all_bests):.1f}")
    print(f"    Min: {min(all_bests)}, Max: {max(all_bests)}")

    # Best pair details
    if global_best_Wn:
        hw_val, Hn, Hf, Wf = evaluate(global_best_Wn)
        dhash = tuple(Hn[i]^Hf[i] for i in range(8))
        print(f"\n  BEST PAIR (hw={hw_val}):")
        print(f"    Wn = {' '.join(f'{x:08x}' for x in global_best_Wn)}")
        print(f"    Wf = {' '.join(f'{x:08x}' for x in Wf)}")
        print(f"    Hn = {' '.join(f'{x:08x}' for x in Hn)}")
        print(f"    Hf = {' '.join(f'{x:08x}' for x in Hf)}")
        print(f"    δH = {' '.join(f'{x:08x}' for x in dhash)}")
        print(f"    Per-word HW:")
        for i in range(8):
            print(f"      H[{i}]: HW={hw(dhash[i])}/32")

    if global_best_hw < 87:
        print(f"\n    *** IMPROVED over manual's hw=87! ***")
        print(f"    Triple verification needed (Rule 14)")
        # Triple verify
        for trial in range(3):
            h, _, _, _ = evaluate(global_best_Wn)
            print(f"    Verify {trial+1}: hw={h}")
    else:
        print(f"\n    Result consistent with manual's hw=87 basin")

if __name__ == "__main__":
    sa_search(steps=200000, restarts=50)
