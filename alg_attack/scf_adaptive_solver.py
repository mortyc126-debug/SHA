#!/usr/bin/env python3
"""
SCF: ADAPTIVE SOLVER — наш улучшенный Greedy.

Проблема старого Greedy: видит только 1-bit flips → застревает на HW≈99.
Улучшения:
1. Multi-bit: пробуем пары бит (top-K полезных × top-K)
2. Word-level: additive коррекция целых слов (±δ, ±2δ, ±carry_correction)
3. Escape: если застряли, пробуем K случайных 2-bit flips
4. History: запоминаем "хорошие" направления и пробуем их снова
"""
import os, sys, math

MASK32 = 0xFFFFFFFF
K_CONST = [
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

def sha_compress(W16):
    W=expand_real(W16); s=list(IV)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K_CONST[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return [add32(IV[i],s[i]) for i in range(8)]

def hash_diff(W1, W2):
    H1 = sha_compress(W1); H2 = sha_compress(W2)
    return sum(hw(H1[i]^H2[i]) for i in range(8))


def adaptive_solver(W1, fixed_word=0, fixed_delta=0x80000000, budget=3000):
    """
    Adaptive solver with multiple move types.

    Move types:
    A. Single bit flip (classic greedy)
    B. Double bit flip (pair of bits)
    C. Additive word correction (+small δ)
    D. Random escape (multi-bit perturbation)
    """
    W2 = list(W1)
    W2[fixed_word] ^= fixed_delta

    current_hw = hash_diff(W1, W2)
    best_hw = current_hw
    best_W2 = list(W2)

    # Track good moves
    good_moves = []  # (delta_hw, move_description)

    evals = 0
    phase = "A"  # Start with single-bit

    stale_count = 0

    while evals < budget:
        improved = False

        if phase == "A":
            # Move A: single bit flip (top-K scan)
            candidates = []
            for word in range(16):
                if word == fixed_word:
                    # Only skip the exact fixed bit
                    pass
                for bit in range(32):
                    if word == fixed_word and (fixed_delta >> bit) & 1:
                        continue
                    W2t = list(W2); W2t[word] ^= (1<<bit)
                    if W2t == W1: continue
                    s = hash_diff(W1, W2t); evals += 1
                    candidates.append((s, word, bit))
                    if evals >= budget: break
                if evals >= budget: break

            candidates.sort()
            if candidates and candidates[0][0] < current_hw:
                s, w, b = candidates[0]
                W2[w] ^= (1<<b)
                current_hw = s
                improved = True
                good_moves.append((-1, f"bit W[{w}][{b}]"))

            if not improved:
                phase = "B"
                continue

        elif phase == "B":
            # Move B: double bit flip from top-10 single-bit influences
            # Compute top-10 most helpful single bits
            singles = []
            for word in range(16):
                for bit in range(32):
                    if word == fixed_word and (fixed_delta >> bit) & 1: continue
                    W2t = list(W2); W2t[word] ^= (1<<bit)
                    if W2t == W1: continue
                    s = hash_diff(W1, W2t); evals += 1
                    singles.append((s - current_hw, word, bit))
                    if evals >= budget: break
                if evals >= budget: break

            singles.sort()
            top_k = min(15, len(singles))

            # Try pairs from top-K
            best_pair = None
            best_pair_hw = current_hw

            for i in range(top_k):
                for j in range(i+1, top_k):
                    w1, b1 = singles[i][1], singles[i][2]
                    w2, b2 = singles[j][1], singles[j][2]
                    W2t = list(W2)
                    W2t[w1] ^= (1<<b1)
                    W2t[w2] ^= (1<<b2)
                    if W2t == W1: continue
                    s = hash_diff(W1, W2t); evals += 1
                    if s < best_pair_hw:
                        best_pair_hw = s
                        best_pair = (w1, b1, w2, b2)
                    if evals >= budget: break
                if evals >= budget: break

            if best_pair and best_pair_hw < current_hw:
                w1,b1,w2,b2 = best_pair
                W2[w1] ^= (1<<b1)
                W2[w2] ^= (1<<b2)
                current_hw = best_pair_hw
                improved = True
                good_moves.append((-2, f"pair W[{w1}][{b1}]+W[{w2}][{b2}]"))

            if not improved:
                phase = "C"
                continue

        elif phase == "C":
            # Move C: additive word correction
            for word in range(16):
                for delta in [1, 2, 3, 7, 15, 31, 255, 0xFFFF, 0x80000000]:
                    for sign in [1, -1]:
                        W2t = list(W2)
                        W2t[word] = (W2[word] + sign * delta) & MASK32
                        if W2t == W1: continue
                        s = hash_diff(W1, W2t); evals += 1
                        if s < current_hw:
                            W2 = W2t
                            current_hw = s
                            improved = True
                            good_moves.append((-3, f"add W[{word}]+={sign*delta}"))
                            break
                        if evals >= budget: break
                    if improved or evals >= budget: break
                if improved or evals >= budget: break

            if not improved:
                phase = "D"
                continue

        elif phase == "D":
            # Move D: random escape — flip 2-4 random bits
            for attempt in range(50):
                W2t = list(W2)
                n_flips = 2 + (attempt % 3)
                for _ in range(n_flips):
                    w = int.from_bytes(os.urandom(1),'big') % 16
                    b = int.from_bytes(os.urandom(1),'big') % 32
                    W2t[w] ^= (1<<b)
                if W2t == W1: continue
                s = hash_diff(W1, W2t); evals += 1
                # Accept if better OR with probability (escape)
                if s < current_hw:
                    W2 = W2t; current_hw = s; improved = True
                    break
                elif s < current_hw + 5 and int.from_bytes(os.urandom(1),'big') < 30:
                    W2 = W2t; current_hw = s; improved = True
                    break
                if evals >= budget: break

            phase = "A"  # Reset

        if current_hw < best_hw:
            best_hw = current_hw
            best_W2 = list(W2)
            stale_count = 0
        else:
            stale_count += 1

        if improved:
            phase = "A"  # Go back to best move type after improvement

        if stale_count > 5:
            phase = "D"  # Force escape
            stale_count = 0

    return best_W2, best_hw, evals, good_moves


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    print("="*70)
    print("SCF: ADAPTIVE SOLVER")
    print("  Multi-bit + additive + escape + history")
    print("="*70)

    greedy_results = []
    adaptive_results = []
    sa_results = []

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        init_hw = hash_diff(W1, [W1[0]^0x80000000]+W1[1:])

        # Adaptive solver
        W2_adp, hw_adp, evals, moves = adaptive_solver(W1, budget=5000)
        adaptive_results.append(hw_adp)

        # SA comparison (same budget)
        W2_sa = list(W1); W2_sa[0] ^= 0x80000000
        best_sa = hash_diff(W1, W2_sa)
        cur_sa = list(W2_sa)
        for it in range(5000):
            t = list(cur_sa)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            if w==0 and b==31: continue
            t[w] ^= (1<<b)
            if t==W1: continue
            s = hash_diff(W1, t)
            T = max(0.01, 1-it/5000)
            if s < best_sa or math.exp(-(s-best_sa)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur_sa = t
                if s < best_sa: best_sa = s
        sa_results.append(best_sa)

        move_types = {}
        for _, desc in moves:
            key = desc.split()[0]
            move_types[key] = move_types.get(key, 0) + 1

        print(f"  Trial {trial:2d}: Adaptive={hw_adp:3d} SA={best_sa:3d} "
              f"(init={init_hw}) moves={move_types}")

    avg_adp = sum(adaptive_results)/len(adaptive_results)
    avg_sa = sum(sa_results)/len(sa_results)
    min_adp = min(adaptive_results)
    min_sa = min(sa_results)

    print(f"\n{'='*70}")
    print(f"SUMMARY ({N} trials, budget=5000 each)")
    print(f"  Adaptive: avg={avg_adp:.1f}, min={min_adp}")
    print(f"  SA:       avg={avg_sa:.1f}, min={min_sa}")
    print(f"  Advantage: avg {avg_sa-avg_adp:+.1f}, min {min_sa-min_adp:+d}")
