"""
MEETING COMPOSITION: два meeting window → full collision?

Window A (front): r=12-17 (a-repair, бесплатно)
Window B (back):  r=61-64 (backward asymmetry, бесплатно)
Gap: r=18-60 (43 раунда)

Вопрос: при ФИКСИРОВАННЫХ window A и B,
какова вероятность что gap "случайно сходится"?

Это НОВОЕ определение: collision = composition двух meeting.
Не парное равенство выходов, а совместное событие двух встреч.

Если P(gap сходится) > 1/C^4: это ДЕШЕВЛЕ стандартного!
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF
K_const = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def add32(x, y): return (x + y) & MASK32
def sub32(x, y): return (x - y) & MASK32
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def R(state, W_r, r_idx):
    a, b, c, d, e, f, g, h = state
    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e, f, g)), K_const[r_idx]), W_r)
    T2 = add32(Sigma0(a), Maj(a, b, c))
    return (add32(T1, T2), a, b, c, add32(d, T1), e, f, g)

def R_inv(sn, W_r, r_idx):
    ap, bp, cp, dp, ep, fp, gp, hp = sn
    a, b, c, e, f, g = bp, cp, dp, fp, gp, hp
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1 = sub32(ap, T2)
    d = sub32(ep, T1)
    h = sub32(sub32(sub32(sub32(T1, Sigma1(e)), Ch(e, f, g)), K_const[r_idx]), W_r)
    return (a, b, c, d, e, f, g, h)

def expand_schedule(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    return W

def a_repair_W(state2, target_a, r_idx):
    a, b, c, d, e, f, g, h = state2
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1_needed = sub32(target_a, T2)
    return sub32(sub32(sub32(sub32(T1_needed, h), Sigma1(e)), Ch(e, f, g)), K_const[r_idx])

def state_diff(s1, s2):
    return sum(hw(s1[i] ^ s2[i]) for i in range(8))

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())


def main():
    np.random.seed(42)

    print("=" * 70)
    print("MEETING COMPOSITION: front window + back window")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. ДВА WINDOW: forward(r=12-17) + backward(r=61-64)")
    print("=" * 70)

    N = 5000

    # Для каждой пары: forward window + backward window
    # Мера: δstate на "gap boundary" (r=18 и r=60)
    # Если оба малы → gap СЖАТ с двух сторон

    fwd_ds18 = []  # δstate[18] from forward
    bwd_ds60 = []  # δstate[60] from backward
    gap_min = []   # minimum δstate in gap
    gap_mid = []   # δstate at midpoint (r=39)

    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)

        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        # A-repair trace 2
        W2 = list(W1)
        W2[3] ^= (1 << 31)
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 16):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)

        # Forward trace 2
        s2_fwd = tuple(IV)
        fwd_states = [s2_fwd]
        for r in range(64):
            w = W2[r] if r < 16 else W2f[r]
            s2_fwd = R(s2_fwd, w, r)
            fwd_states.append(s2_fwd)

        # Backward trace from state1[64] with W2 schedule
        bwd_states = [None] * 65
        bwd_states[64] = states1[64]
        for r in range(63, -1, -1):
            bwd_states[r] = R_inv(bwd_states[r+1], W2f[r], r)

        # Forward δstate at gap entry (r=18)
        fwd_d18 = state_diff(states1[18], fwd_states[18])
        fwd_ds18.append(fwd_d18)

        # Backward δstate at gap exit (r=60)
        bwd_d60 = state_diff(states1[60], bwd_states[60])
        bwd_ds60.append(bwd_d60)

        # δstate through gap (forward trace vs trace1)
        min_gap = 256
        mid_gap = 0
        for r in range(18, 61):
            d = state_diff(states1[r], fwd_states[r])
            if d < min_gap:
                min_gap = d
            if r == 39:
                mid_gap = d
        gap_min.append(min_gap)
        gap_mid.append(mid_gap)

    print(f"\n  {N} пар с a-repair (front window r=12-17):")
    print(f"    Forward δstate[18]:  mean={np.mean(fwd_ds18):.1f}")
    print(f"    Backward δstate[60]: mean={np.mean(bwd_ds60):.1f}")
    print(f"    Gap midpoint [39]:   mean={np.mean(gap_mid):.1f}")
    print(f"    Gap minimum:         mean={np.mean(gap_min):.1f}, min={min(gap_min)}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. GAP COMPRESSION: forward squeeze + backward squeeze")
    print("=" * 70)

    # Forward window даёт δstate[17]=0, δstate[18]≈10
    # Backward window даёт δstate[64]=0, δstate[60]≈64
    # Между ними: gap r=18..60 = 43 раунда
    # Вопрос: forward и backward СЖИМАЮТ gap с двух сторон?

    # Измеряем: δstate по раундам для forward и backward ОДНОВРЕМЕННО
    combined_profile = np.zeros((N, 65))

    for trial in range(min(N, 1000)):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)
        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        W2 = list(W1)
        W2[3] ^= (1 << 31)
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 16):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)
        W2f = expand_schedule(W2)

        # Forward profile
        s2_f = tuple(IV)
        for r in range(64):
            w = W2[r] if r < 16 else W2f[r]
            s2_f = R(s2_f, w, r)
            combined_profile[trial, r+1] = state_diff(states1[r+1], s2_f)

    mean_combined = np.mean(combined_profile[:1000], axis=0)

    # backward δ profile
    bwd_profile = np.zeros(65)
    for trial in range(min(N, 500)):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)
        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        W2 = list(W1)
        W2[3] ^= (1 << 31)
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 16):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)
        W2f = expand_schedule(W2)

        sb = states1[64]
        for r in range(63, -1, -1):
            sb = R_inv(sb, W2f[r], r)
            bwd_profile[r] += state_diff(states1[r], sb)

    bwd_profile /= 500

    print(f"\n  Combined profile (forward + backward):")
    print(f"  {'r':>4} {'Fwd δ':>7} {'Bwd δ':>7} {'Min':>7}  Visual")
    for r in [0, 5, 10, 12, 16, 17, 18, 20, 25, 30, 35, 39, 40, 45, 50, 55, 58, 60, 61, 62, 63, 64]:
        f = mean_combined[r]
        b = bwd_profile[r]
        m = min(f, b)

        zone = "F" if f < 5 else "B" if b < 50 else "G"

        fbar = ">" * min(int(f/8), 16)
        bbar = "<" * min(int(b/8), 16)

        print(f"  {r:4d} {f:6.1f} {b:6.1f} {m:6.1f}  [{zone}] {fbar}|{bbar}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. ВСТРЕЧА В СЕРЕДИНЕ: δstate_fwd ≈ δstate_bwd?")
    print("=" * 70)

    # Forward: δstate растёт от 0 (r=17) до 128 (r=22+)
    # Backward: δstate растёт от 0 (r=64) до 128 (r=47-)
    # Встреча: где forward ≈ backward?

    # Forward reaches 128 at r≈22. Backward reaches 128 at r≈47.
    # Between r=22 and r=47: both ≈ 128. No squeeze.

    # НО: если мы ищем пары где forward МЕДЛЕННЕЕ
    # или backward МЕДЛЕННЕЕ → overlap zone шире

    # Pipe collision: δ на pipe-позициях
    # Forward: pipes зависят от a-chain (которая a-repaired)
    # Backward: pipes зависят от e-chain (backward asymmetry)

    # Ищем: пары где forward pipe-δ < 64 AND backward pipe-δ < 64
    # на каком-то раунде r → "pipe-meeting"

    pipe_meeting_count = 0
    any_meeting = 0

    for trial in range(min(N, 2000)):
        # Check if min(forward, backward) < 64 at any gap round
        min_combined = 256
        for r in range(20, 60):
            f = combined_profile[trial, r]
            b = bwd_profile[r]  # average, not per-trial
            m = min(f, b)
            if m < min_combined:
                min_combined = m

        if min_combined < 100:
            any_meeting += 1
        if min_combined < 64:
            pipe_meeting_count += 1

    print(f"\n  Gap meeting (min(fwd,bwd) < threshold):")
    print(f"    < 100: {any_meeting} / {min(N, 2000)}")
    print(f"    < 64:  {pipe_meeting_count} / {min(N, 2000)}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. НОВОЕ ОПРЕДЕЛЕНИЕ: степень collision")
    print("=" * 70)

    # Степень collision = число раундов где δstate = 0
    degree_dist = {}
    for trial in range(min(N, 2000)):
        degree = sum(1 for r in range(65) if combined_profile[trial, r] == 0)
        if degree not in degree_dist:
            degree_dist[degree] = 0
        degree_dist[degree] += 1

    print(f"\n  Распределение СТЕПЕНИ collision (число раундов с δ=0):")
    for deg in sorted(degree_dist.keys()):
        count = degree_dist[deg]
        pct = count / min(N, 2000) * 100
        cost_est = f"FREE" if pct > 99 else f"×{int(100/max(pct,0.1))}" if pct > 1 else "rare"
        bar = "█" * (count * 30 // min(N, 2000))
        print(f"    degree {deg:2d}: {count:5d} ({pct:5.1f}%)  cost≈{cost_est:>6}  {bar}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. PIPE-COLLISION: проекция на 6 pipe-позиций")
    print("=" * 70)

    # Pipe-collision: δH[1,2,3,5,6,7] = 0 (192 бита)
    # Стоимость в нашем мире: overlay на C^3 = (2^32)^3 = 2^96

    # Проверяем: pipe-collision сколько стоит?
    # И: при pipe-collision, чему равен δH[0,4] (node-позиции)?

    N_pipe = 200000
    pipe_dict = {}
    pipe_collisions = []

    for trial in range(N_pipe):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_words(W)
        pipe_key = (H[1], H[2], H[3], H[5], H[6], H[7])

        if pipe_key in pipe_dict:
            H_prev, W_prev = pipe_dict[pipe_key]
            pipe_collisions.append((H, H_prev))
        else:
            pipe_dict[pipe_key] = (H, W)

    print(f"\n  {N_pipe} хешей → pipe-collisions (6 word match): {len(pipe_collisions)}")
    print(f"  Expected (birthday 192 bit): {N_pipe**2 / 2**(192+1):.4f}")

    if pipe_collisions:
        for H1, H2 in pipe_collisions[:5]:
            node_diff = [hw(H1[i] ^ H2[i]) for i in [0, 4]]
            pipe_diff = [hw(H1[i] ^ H2[i]) for i in [1,2,3,5,6,7]]
            print(f"    Node δ: H[0]={node_diff[0]}, H[4]={node_diff[1]} | Pipe δ: {pipe_diff}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("6. НАТИВНАЯ СТОИМОСТЬ В НОВЫХ ОПРЕДЕЛЕНИЯХ")
    print("=" * 70)

    print(f"""
  СПЕКТР COLLISION В НАШЕМ МИРЕ:

  Тип                  Условие                     Стоимость
  ─────────────────── ─────────────────────────── ──────────
  Meeting window 5    δstate[12..16]=0             C^0 (FREE)
  Meeting window 7    δstate[12..18]=0             C^0 × 4
  Pipe-collision      δH[1,2,3,5,6,7]=0 (6 слов)  C^3 = 2^96
  Full collision      δH[0..7]=0 (8 слов)          C^4 = 2^128

  НОВЫЕ определения (НЕ парные):
  Degree-K collision  K раундов δstate=0            varies
  Composed meeting    front+back window overlap     ???
  Pipe-absorption     δ на pipes = 0                C^3

  Pipe-collision (C^3 = 2^96) — это ДЕШЕВЛЕ full collision!
  Но это 6/8 слов, не 8/8.

  Для ПРАКТИКИ: pipe-collision + brute node = 2^96 + 2^64 = 2^96.
  Если node-позиции зависят от pipe → может ещё дешевле.
""")


if __name__ == "__main__":
    main()
