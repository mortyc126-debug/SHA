#!/usr/bin/env python3
"""
НОВАЯ МАТЕМАТИКА — с нуля, без ссылок на существующую.

АКСИОМЫ (то что мы ЗНАЕМ, без теорий):
A1. f: {0,1}^512 → {0,1}^256 — детерминированная функция
A2. f вычислима (~3000 операций)
A3. Хотим: x ≠ y, f(x) = f(y)

ОПРЕДЕЛЕНИЯ (наши, новые):
D1. ЛАНДШАФТ: L(x) = f(x). Точка в 512-мерном пространстве
    имеет "высоту" — 256-бит значение.
D2. УРОВЕНЬ: S(h) = {x : f(x) = h}. Множество прообразов.
    |S(h)| ≈ 2^{512}/2^{256} = 2^{256} (в среднем).
D3. СВЯЗНОСТЬ: x,y ∈ S(h) — пара на одном уровне.
    Вопрос: как найти пару?

Новое понятие: ТРАЕКТОРИЯ.
T(x, d) = {f(x), f(x⊕d), f(x⊕2d), f(x⊕3d), ...}
Это значения f вдоль НАПРАВЛЕНИЯ d в пространстве входов.

Если траектория ВОЗВРАЩАЕТСЯ (f(x⊕kd) = f(x)) → коллизия!

Вопрос 1: каков ПЕРИОД траектории?
Вопрос 2: зависит ли период от направления d?
Вопрос 3: можно ли выбрать d с КОРОТКИМ периодом?
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
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK32
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def add32(*args):
    s=0
    for x in args: s=(s+x)&MASK32
    return s
def hw(x): return bin(x).count('1')

def expand_real(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def f(x):
    """THE FUNCTION. 16 words in, 8 words out. No names, no theory."""
    W=expand_real(x); s=list(IV)
    for r in range(64):
        a,b,c,d,e,f_,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f_,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f_,g]
    return tuple(add32(IV[i],s[i]) for i in range(8))

def dist(h1, h2):
    """Distance between two outputs. Our own metric."""
    return sum(hw(h1[i]^h2[i]) for i in range(8))

def walk(x, d, steps):
    """Walk along direction d from x. Return trajectory of outputs."""
    trajectory = []
    current = list(x)
    for step in range(steps):
        h = f(current)
        trajectory.append(h)
        # Step: add d to current (word-by-word, wrapping)
        for i in range(16):
            current[i] = (current[i] + d[i]) & MASK32
    return trajectory


# ============================================================
# EXPLORATION 1: TRAJECTORY — walk along d, look for returns
# ============================================================
def explore1_trajectory(N_directions, N_steps):
    print("="*70)
    print("EXPLORATION 1: TRAJECTORY PERIODS")
    print("  Walk x, x+d, x+2d, ... and track f(x+kd)")
    print("  If f(x+kd) = f(x) for some k → collision!")
    print("="*70)

    x = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    h_start = f(x)

    best_return_dist = 256
    best_direction = None

    for di in range(N_directions):
        # Generate direction
        d = [0]*16
        if di < 16:
            # Single-word direction
            d[di] = 1
        elif di < 32:
            d[di-16] = 0x80000000
        else:
            # Random direction
            d = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        trajectory = walk(x, d, N_steps)

        # Find closest return to start
        min_dist = 256
        min_step = -1
        for step, h in enumerate(trajectory[1:], 1):
            dd = dist(h_start, h)
            if dd < min_dist:
                min_dist = dd
                min_step = step

        if min_dist < best_return_dist:
            best_return_dist = min_dist
            best_direction = (di, d, min_step)

        if di < 5 or min_dist < 100:
            d_desc = f"word[{di}]=1" if di < 16 else (f"word[{di-16}]=MSB" if di < 32 else "random")
            print(f"  d={d_desc}: closest return at step {min_step}, dist={min_dist}")

    print(f"\n  Best: direction #{best_direction[0]}, step={best_direction[2]}, "
          f"dist={best_return_dist}")


# ============================================================
# EXPLORATION 2: LEVEL SETS — how many x give same f(x)?
# ============================================================
def explore2_level_sets(N):
    print("\n" + "="*70)
    print("EXPLORATION 2: LEVEL SETS")
    print("  How clustered are f-values? Do some outputs repeat?")
    print("  If we find f(x1)=f(x2) with x1≠x2 → done!")
    print("="*70)

    # Truncated: match on first word of hash
    seen_trunc = {}
    closest_full = 256
    n_partial = 0

    for trial in range(N):
        x = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        h = f(x)

        # Check first word match
        key = h[0]
        if key in seen_trunc:
            old_x, old_h = seen_trunc[key]
            full_dist = dist(h, old_h)
            if full_dist < closest_full:
                closest_full = full_dist
                n_partial += 1
                if full_dist < 200:
                    print(f"  Partial match at trial {trial}: dist={full_dist}")
        else:
            seen_trunc[key] = (x, h)

    print(f"\n  {N} random x evaluated")
    print(f"  Partial matches (first word): {n_partial}")
    print(f"  Closest full distance: {closest_full}")
    print(f"  Expected first-word collisions at N={N}: ~{N**2/(2*2**32):.1f}")


# ============================================================
# EXPLORATION 3: NEIGHBORHOODS — structure AROUND each point
# ============================================================
def explore3_neighborhoods(N):
    print("\n" + "="*70)
    print("EXPLORATION 3: NEIGHBORHOODS")
    print("  For a point x: what does f look like NEARBY?")
    print("  Nearby = x with 1 bit flipped")
    print("  Is there STRUCTURE in the neighborhood distances?")
    print("="*70)

    for trial in range(min(N, 5)):
        x = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        h_center = f(x)

        # 512 neighbors (1-bit Hamming distance)
        neighbor_dists = []
        for word in range(16):
            for bit in range(32):
                x_neighbor = list(x)
                x_neighbor[word] ^= (1 << bit)
                h_neighbor = f(x_neighbor)
                d = dist(h_center, h_neighbor)
                neighbor_dists.append((d, word, bit))

        neighbor_dists.sort()

        avg_d = sum(d for d,_,_ in neighbor_dists)/len(neighbor_dists)
        min_d = neighbor_dists[0][0]
        max_d = neighbor_dists[-1][0]
        std_d = (sum((d-avg_d)**2 for d,_,_ in neighbor_dists)/len(neighbor_dists))**0.5

        print(f"\n  Point {trial}: avg_dist={avg_d:.1f}, min={min_d}, max={max_d}, std={std_d:.1f}")
        print(f"    5 closest neighbors:")
        for d, word, bit in neighbor_dists[:5]:
            print(f"      W[{word}][{bit:2d}]: dist={d}")
        print(f"    5 farthest:")
        for d, word, bit in neighbor_dists[-5:]:
            print(f"      W[{word}][{bit:2d}]: dist={d}")

        # KEY: is the neighborhood SYMMETRIC?
        # If not → there's a PREFERRED direction to move
        below_avg = sum(1 for d,_,_ in neighbor_dists if d < avg_d)
        above_avg = sum(1 for d,_,_ in neighbor_dists if d >= avg_d)
        print(f"    Below avg: {below_avg}/512, Above avg: {above_avg}/512")

        # Distribution shape
        buckets = [0]*20
        for d,_,_ in neighbor_dists:
            bucket = min(int(d / 256 * 20), 19)
            buckets[bucket] += 1
        print(f"    Distance distribution: ", end="")
        for b in buckets:
            print(f"{'█' * (b//5)}", end="")
        print()


# ============================================================
# EXPLORATION 4: MULTI-STEP DESCENT — greedy walk toward target
# ============================================================
def explore4_descent(N_starts, N_steps):
    print("\n" + "="*70)
    print("EXPLORATION 4: DESCENT")
    print("  Start at x. Target: find y with f(y) closest to f(x).")
    print("  At each step: move to neighbor with SMALLEST dist to target.")
    print("  NEW: target = f(random_other), not f(x).")
    print("="*70)

    for start in range(N_starts):
        x = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        # Target = hash of a DIFFERENT random message
        y = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        target = f(y)

        current = list(x)
        current_dist = dist(f(current), target)

        path = [current_dist]

        for step in range(N_steps):
            # Find best neighbor
            best_d = current_dist
            best_move = None

            for word in range(16):
                for bit in range(32):
                    neighbor = list(current)
                    neighbor[word] ^= (1 << bit)
                    d = dist(f(neighbor), target)
                    if d < best_d:
                        best_d = d
                        best_move = (word, bit)

            if best_move:
                current[best_move[0]] ^= (1 << best_move[1])
                current_dist = best_d
                path.append(current_dist)
            else:
                break

        final_dist = current_dist
        n_steps_taken = len(path) - 1

        if start < 5:
            print(f"\n  Start {start}: initial dist={path[0]}")
            for i, d in enumerate(path[1:6], 1):
                print(f"    Step {i}: dist={d}")
            if n_steps_taken > 5:
                print(f"    ... ({n_steps_taken} steps total)")
            print(f"    Final: dist={final_dist}")

        if final_dist == 0:
            print(f"\n  ★★★ COLLISION FOUND! Start {start}, {n_steps_taken} steps!")
            print(f"    x = {[hex(w) for w in x[:4]]}...")
            print(f"    y = {[hex(w) for w in y[:4]]}...")
            print(f"    current = {[hex(w) for w in current[:4]]}...")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    explore1_trajectory(40, min(N*10, 200))
    explore2_level_sets(min(N*1000, 50000))
    explore3_neighborhoods(3)
    explore4_descent(min(N, 5), 20)
