#!/usr/bin/env python3
"""
SCF: GREEDY BIT-FLIP SOLVER — построен на открытиях Anatomy.

Открытия:
1. ~15% бит помогают, 85% вредят → flip ТОЛЬКО полезные
2. Positive interference → некоторые биты помогают 20+ раундам
3. Раунды некоррелированы → можно оптимизировать подгруппами

Алгоритм:
1. Вычислить influence(word, bit) = Δ(total_mismatch) для ВСЕХ 512 бит
2. Отсортировать: от самых полезных к самым вредным
3. Жадно flip'ать полезные биты, проверяя улучшение
4. Если mismatch перестал уменьшаться → пересчитать influence (адаптивно)
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

def sha_all_states(W16):
    W=expand_real(W16); s=list(IV); states=[list(s)]
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        states.append(list(s))
    return states, W

def sha_compress(W16):
    st, _ = sha_all_states(W16)
    return [add32(IV[i], st[64][i]) for i in range(8)]

def compute_mismatch_total(W1, W2):
    states1, Wexp1 = sha_all_states(W1)
    states2, Wexp2 = sha_all_states(W2)
    total = 0
    for r in range(17, 52):
        s1r = states1[r]; s2r = states2[r]
        T1_1 = sub32(states1[r+1][4], s1r[3])
        T1_2_req = add32(T1_1, sub32(s1r[3], s2r[3]))
        W2_req = sub32(sub32(sub32(sub32(T1_2_req,
                  s2r[7]), Sig1(s2r[4])), Ch(s2r[4],s2r[5],s2r[6])), K[r])
        total += hw(W2_req ^ Wexp2[r])
    return total

def compute_influence_map(W1, W2):
    """Compute Δ(mismatch) for each of 512 input bit flips."""
    base = compute_mismatch_total(W1, W2)
    influences = []
    for word in range(16):
        for bit in range(32):
            W2_flip = list(W2)
            W2_flip[word] ^= (1 << bit)
            flip_score = compute_mismatch_total(W1, W2_flip)
            influences.append((flip_score - base, word, bit))
    return influences, base


def greedy_solver(W1, W2_init, max_rounds=20, recompute_every=50):
    """Greedy bit-flip solver with adaptive influence recomputation."""
    W2 = list(W2_init)
    best_W2 = list(W2)
    best_mismatch = compute_mismatch_total(W1, W2)
    initial = best_mismatch

    total_flips = 0

    for round_num in range(max_rounds):
        # Compute full influence map
        influences, current = compute_influence_map(W1, W2)

        # Sort: most helpful first
        influences.sort()

        # Greedy: flip all helpful bits (Δ < 0)
        improved = False
        flips_this_round = 0

        for delta, word, bit in influences:
            if delta >= 0:
                break  # No more helpful bits

            # Tentatively flip
            W2_try = list(W2)
            W2_try[word] ^= (1 << bit)
            actual_score = compute_mismatch_total(W1, W2_try)

            if actual_score < current:
                W2 = W2_try
                current = actual_score
                improved = True
                flips_this_round += 1
                total_flips += 1

                if current < best_mismatch:
                    best_mismatch = current
                    best_W2 = list(W2)

                # Don't flip too many without recomputing
                if flips_this_round >= recompute_every:
                    break

        print(f"  Round {round_num:2d}: mismatch={current:4d} "
              f"(flipped {flips_this_round}, total={total_flips})")

        if not improved:
            break

    return best_W2, best_mismatch, initial


def greedy_with_hash(W1, W2_init, max_rounds=15):
    """Greedy solver targeting HASH diff, not mismatch."""
    W2 = list(W2_init)
    H1 = sha_compress(W1)

    def hash_hw(W):
        H = sha_compress(W)
        return sum(hw(H1[i]^H[i]) for i in range(8))

    best_hw = hash_hw(W2)
    best_W2 = list(W2)

    for round_num in range(max_rounds):
        # Compute influence for hash
        base = hash_hw(W2)
        influences = []
        for word in range(16):
            for bit in range(32):
                W2f = list(W2); W2f[word] ^= (1<<bit)
                influences.append((hash_hw(W2f) - base, word, bit))
        influences.sort()

        improved = False
        for delta, word, bit in influences:
            if delta >= 0: break
            W2t = list(W2); W2t[word] ^= (1<<bit)
            s = hash_hw(W2t)
            if s < base:
                W2 = W2t; base = s; improved = True
                if s < best_hw: best_hw = s; best_W2 = list(W2)
                break  # Recompute after each flip

        if round_num < 5 or round_num % 5 == 0:
            print(f"  Round {round_num:2d}: hash HW={base}")
        if not improved: break

    return best_W2, best_hw


if __name__ == '__main__':
    print("="*70)
    print("SCF: GREEDY BIT-FLIP SOLVER")
    print("="*70)

    results = []

    for trial in range(5):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        print(f"\n{'='*50}")
        print(f"Trial {trial}: Mismatch solver")
        W2_opt, final_mm, init_mm = greedy_solver(W1, W2, max_rounds=15)

        H1 = sha_compress(W1); H2 = sha_compress(W2_opt)
        hash_hw_result = sum(hw(H1[i]^H2[i]) for i in range(8))
        n_diff = sum(1 for i in range(16) if W1[i] != W2_opt[i])

        print(f"  Mismatch: {init_mm} → {final_mm} ({init_mm-final_mm:+d})")
        print(f"  Hash HW: {hash_hw_result}")
        print(f"  Words changed: {n_diff}/16")

        # Now hash solver from mismatch result
        print(f"\n  Hash solver:")
        W2_hash, hash_best = greedy_with_hash(W1, W2_opt, max_rounds=15)
        n_diff2 = sum(1 for i in range(16) if W1[i] != W2_hash[i])
        print(f"  Hash HW: {hash_best}")

        results.append((final_mm, hash_hw_result, hash_best))

    # Random baseline
    print(f"\n{'='*50}")
    print("Random baseline (SA, 5000 iters):")
    for trial in range(3):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000
        H1 = sha_compress(W1)
        best = 256
        cur = list(W2)
        for it in range(5000):
            t = list(cur)
            t[int.from_bytes(os.urandom(1),'big')%16] ^= (1<<(int.from_bytes(os.urandom(1),'big')%32))
            H2 = sha_compress(t)
            s = sum(hw(H1[i]^H2[i]) for i in range(8))
            if s < best or (math.exp(-(s-best)/max(0.1,3-3*it/5000)) > int.from_bytes(os.urandom(4),'big')/(1<<32)):
                cur = t; best = min(best, s)
        print(f"  SA trial {trial}: hash HW = {best}")

    print(f"\n{'='*70}")
    print("ИТОГ")
    print(f"{'='*70}")
    avg_mm = sum(r[0] for r in results)/len(results)
    avg_hash_mm = sum(r[1] for r in results)/len(results)
    avg_hash_direct = sum(r[2] for r in results)/len(results)
    print(f"  Greedy mismatch: avg final = {avg_mm:.0f}")
    print(f"  Greedy hash (from mismatch): avg = {avg_hash_mm:.0f}")
    print(f"  Greedy hash (direct): avg = {avg_hash_direct:.0f}")
