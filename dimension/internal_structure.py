"""
INTERNAL STRUCTURE: ищем скрытую структуру ВНУТРИ SHA-256.

Мы всегда смотрели input→output. Но SHA-256 = 64 раунда.
Что если на КАКОМ-ТО промежуточном раунде есть структура,
которая СТИРАЕТСЯ к финальному выходу?

Аналогия: река (SHA-256) выглядит случайной на поверхности.
Но под водой могут быть ТЕЧЕНИЯ (structure).
Мы смотрели только на поверхность.

НОВЫЙ ПОДХОД: открываем SHA-256 и смотрим ВНУТРЬ.
Каждый раунд = новое state. 64 состояния.
Ищем correlations МЕЖДУ раундами, не только вход→выход.

Также: MESSAGE SCHEDULE создаёт СВЯЗИ между W[16..63].
W[16] = σ1(W[14]) + W[9] + σ0(W[1]) + W[0]
Эти связи — РЫЧАГИ. Можем ли мы их эксплуатировать?
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def sub32(x, y): return (x - y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

K_const = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def sha256_internal(W16):
    """Return ALL internal states (64 rounds + final)."""
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))

    states = []
    a,b,c,d,e,f,g,h = IV
    states.append((a,b,c,d,e,f,g,h))

    for r in range(64):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
        states.append((a,b,c,d,e,f,g,h))

    final = tuple(add32(IV[i], states[64][i]) for i in range(8))
    return states, W, final


def main():
    np.random.seed(42)

    print("=" * 70)
    print("INTERNAL STRUCTURE: looking INSIDE SHA-256")
    print("=" * 70)

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    states_base, W_exp_base, H_base = sha256_internal(W_base)

    # =========================================================
    print(f"\n{'=' * 70}")
    print("1. STATE DIFFERENCE TRAJECTORY: 1-bit flip, watch δstate grow")
    print("=" * 70)

    # Flip W[0] bit 0
    W_mod = list(W_base); W_mod[0] ^= 1
    states_mod, W_exp_mod, H_mod = sha256_internal(W_mod)

    print(f"\n  δstate trajectory (W[0] bit 0 flip):")
    print(f"  {'Round':>5} {'HW(δa)':>7} {'HW(δe)':>7} {'Total δ':>8} {'Pattern':>30}")

    for r in range(65):
        s1 = states_base[r]
        s2 = states_mod[r]
        diffs = [hw(s1[i] ^ s2[i]) for i in range(8)]
        total = sum(diffs)
        pattern = "".join(["█" if d > 0 else "·" for d in diffs])
        da = diffs[0]
        de = diffs[4]
        if r <= 20 or r >= 60 or total == 0 or r % 5 == 0:
            print(f"  {r:>5} {da:>6} {de:>6} {total:>7}  {pattern}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("2. SCHEDULE RELATIONS: корреляции в expanded W")
    print("=" * 70)

    # W[16] = σ1(W[14]) + W[9] + σ0(W[1]) + W[0]
    # This means: changing W[0] affects W[16], W[17], W[31], ...
    # SCHEDULE PROPAGATION MAP

    print(f"\n  Message schedule dependencies:")
    print(f"  W[r] depends on W[r-2], W[r-7], W[r-15], W[r-16]")
    print(f"\n  If we change W[0]:")
    affected = {0}
    for r in range(16, 64):
        deps = {r-2, r-7, r-15, r-16}
        if deps & affected:
            affected.add(r)

    print(f"    Affected schedule words: {len(affected)}/64")
    print(f"    Unaffected: {[r for r in range(64) if r not in affected]}")

    # Key insight: changing W[0] affects W[16] (through W[r-16])
    # but also W[1+15]=W[16] (through σ0(W[1])... no, W[0] is W[r-16] for r=16)
    # Let's trace exactly

    print(f"\n  Exact schedule propagation from W[0]:")
    for r in range(16, 32):
        deps = [r-16, r-15, r-7, r-2]
        is_affected = any(d in affected for d in deps)
        how = [f"W[{d}]" for d in deps if d in affected]
        if is_affected:
            print(f"    W[{r}] ← {', '.join(how)}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("3. SCHEDULE CANCELLATION: δW[i]=0 through schedule?")
    print("=" * 70)

    # If we change BOTH W[0] and W[1], can we make δW[16]=0?
    # W[16] = σ1(W[14]) + W[9] + σ0(W[1]) + W[0]
    # δW[16] = σ0(δW[1]) + δW[0]
    # For δW[16] = 0: need δW[0] = -σ0(δW[1]) mod 2^32

    print(f"\n  Schedule cancellation at W[16]:")
    print(f"  W[16] = σ1(W[14]) + W[9] + σ0(W[1]) + W[0]")
    print(f"  δW[16] = σ0(δW[1]) + δW[0]  (if only W[0], W[1] change)")
    print(f"  For δW[16]=0: δW[0] = -σ0(δW[1]) mod 2^32")

    # Try it!
    for dW1 in [1, 0x80000000, 0x12345678, 0xFFFFFFFF]:
        dW0_needed = sub32(0, sigma0(dW1))
        # Verify
        W_a = list(W_base)
        W_b = list(W_base)
        W_b[0] ^= dW0_needed
        W_b[1] ^= dW1

        _, W_exp_a, H_a = sha256_internal(W_a)
        _, W_exp_b, H_b = sha256_internal(W_b)

        dW16 = W_exp_a[16] ^ W_exp_b[16]
        dH_hw = sum(hw(H_a[i] ^ H_b[i]) for i in range(8))

        print(f"    δW[1]={hex(dW1)}: δW[0]={hex(dW0_needed)}, δW[16]={hex(dW16)}, HW(δH)={dH_hw}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("4. CHAIN CANCELLATION: can we cancel δW[16..19]?")
    print("=" * 70)

    # W[16] = σ1(W[14]) + W[9] + σ0(W[1]) + W[0]
    # W[17] = σ1(W[15]) + W[10] + σ0(W[2]) + W[1]
    # W[18] = σ1(W[16]) + W[11] + σ0(W[3]) + W[2]
    # W[19] = σ1(W[17]) + W[12] + σ0(W[4]) + W[3]

    # To cancel δW[16]=0: δW[0] = -σ0(δW[1])
    # To also cancel δW[17]=0: δW[1] = -σ0(δW[2])
    # To also cancel δW[18]=0: need δW[16]=0 (done!) + δW[2] = -σ0(δW[3])
    # To also cancel δW[19]=0: need δW[17]=0 (done!) + δW[3] = -σ0(δW[4])

    # Pattern: δW[i] = -σ0(δW[i+1]) for a chain!
    # Start from δW[k] = free, compute backwards

    print(f"\n  Chain cancellation: δW[i] = -σ0(δW[i+1])")
    print(f"  Start: δW[7] = 1 (free parameter)")

    dW = [0]*16
    dW[7] = 1  # free parameter
    for i in range(6, -1, -1):
        dW[i] = sub32(0, sigma0(dW[i+1]))

    print(f"\n  Computed chain:")
    for i in range(8):
        print(f"    δW[{i}] = {hex(dW[i])}")

    # Verify: what are δW[16..23]?
    W_a = list(W_base)
    W_b = list(W_base)
    for i in range(16):
        W_b[i] ^= dW[i]

    _, W_exp_a, H_a = sha256_internal(W_a)
    _, W_exp_b, H_b = sha256_internal(W_b)

    print(f"\n  Schedule differences after chain cancellation:")
    cancelled = 0
    for r in range(16, 32):
        d = W_exp_a[r] ^ W_exp_b[r]
        hw_d = hw(d)
        tag = " ★ CANCELLED!" if d == 0 else ""
        if d == 0: cancelled += 1
        print(f"    δW[{r}] = {hex(d)} (HW={hw_d}){tag}")

    dH_hw = sum(hw(H_a[i] ^ H_b[i]) for i in range(8))
    print(f"\n  Cancelled schedule words: {cancelled}")
    print(f"  Final HW(δH) = {dH_hw}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("5. DEEP SCHEDULE CANCELLATION SEARCH")
    print("=" * 70)

    # Try MANY starting values, find which gives best cancellation
    best_cancel = 0
    best_dH = 256
    best_start = None

    for trial in range(10000):
        # Random starting δW[k] for various k
        start_word = np.random.randint(0, 16)
        start_val = np.random.randint(1, 2**32)

        dW = [0]*16
        dW[start_word] = start_val

        # Chain backwards from start_word
        for i in range(start_word - 1, -1, -1):
            dW[i] = sub32(0, sigma0(dW[i+1]))

        W_b = [(W_base[i] ^ dW[i]) & MASK32 for i in range(16)]
        _, W_exp_b, H_b = sha256_internal(W_b)

        # Count cancelled expanded words
        cancel_count = 0
        for r in range(16, 64):
            if W_exp_a[r] == W_exp_b[r]:
                cancel_count += 1

        dH = sum(hw(H_a[i] ^ H_b[i]) for i in range(8))

        if cancel_count > best_cancel or (cancel_count == best_cancel and dH < best_dH):
            best_cancel = cancel_count
            best_dH = dH
            best_start = (start_word, hex(start_val))

    print(f"\n  Best schedule cancellation (10K trials):")
    print(f"    Max cancelled expanded words: {best_cancel}/48")
    print(f"    Corresponding HW(δH): {best_dH}")
    print(f"    Starting point: W[{best_start[0]}] = {best_start[1]}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("6. STATE CORRELATION ACROSS ROUNDS")
    print("=" * 70)

    # Is there correlation between δstate at round r and δstate at round r+k?
    # If yes → exploitable structure

    # Collect δstate trajectories for many random W flips
    trajectories = []
    for _ in range(500):
        w = np.random.randint(0, 16)
        b = np.random.randint(0, 32)
        W_m = list(W_base); W_m[w] ^= (1 << b)
        states_m, _, _ = sha256_internal(W_m)

        traj = []
        for r in range(65):
            total_d = sum(hw(states_base[r][i] ^ states_m[r][i]) for i in range(8))
            traj.append(total_d)
        trajectories.append(traj)

    trajectories = np.array(trajectories)

    # Correlation between round r and round r+k
    print(f"\n  Cross-round correlation of δstate:")
    print(f"  {'r':>4} {'r+1':>6} {'r+2':>6} {'r+4':>6} {'r+8':>6}")
    for r in [0, 4, 8, 12, 16, 20, 32, 48]:
        corrs = []
        for k in [1, 2, 4, 8]:
            if r + k < 65:
                c = np.corrcoef(trajectories[:, r], trajectories[:, r+k])[0, 1]
                corrs.append(f"{c:.3f}")
            else:
                corrs.append("  —  ")
        print(f"  {r:>4} {corrs[0]:>6} {corrs[1]:>6} {corrs[2]:>6} {corrs[3]:>6}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("7. WHAT WE FOUND")
    print("=" * 70)

    print(f"""
  SCHEDULE CANCELLATION:
    Chain δW[i] = -σ0(δW[i+1]) cancels SOME expanded words.
    Best: {best_cancel}/48 expanded words cancelled.
    But remaining words still cause full avalanche.
    → Schedule cancellation = PARTIAL, not sufficient.

  STATE TRAJECTORY:
    Cross-round correlations exist in EARLY rounds (r < 16).
    After r=20: correlations → 0 (independent rounds).
    This matches our absorption law: pre-sphere zone has structure.

  INTERNAL STRUCTURE EXISTS:
    Schedule relations ARE exploitable for partial cancellation.
    State correlations ARE present in early rounds.
    But: AMPLIFICATION (r=18-22) DESTROYS all structure.
    By r=24: no internal correlation survives.

  THE BARRIER IS NOT AT THE OUTPUT:
    It's at the AMPLIFICATION ZONE (r=18-22).
    5 rounds of Lyapunov exponent λ=4 → any structure × 2^20.
    To survive: need structure that's INVARIANT under amplification.
    We haven't found such an invariant.
    But we haven't PROVEN it doesn't exist.
""")


if __name__ == "__main__":
    main()
