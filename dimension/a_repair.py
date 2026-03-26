"""
A-REPAIR в нашем измерении.

Из потокового измерения (Flow Mathematics):
  Wang чинит СИМПТОМ (δe=0).
  a-repair чинит ПРИЧИНУ (δa=0).

  δa=0 → shift register: δb[+1]=0, δc[+2]=0, δd[+3]=0
  → через 3 раунда ВСЕ входы T2 нулевые → δT2=0
  → a-ветвь стабилизирована
  → e-ветвь ТОЖЕ стабилизируется (d=0 → e_new = 0 + T1)

Формула a-repair:
  Цель: a₂[r+1] = a₁[r+1]
  a_new = T1 + T2
  T2₂ = Σ0(a₂) + Maj(a₂,b₂,c₂)        — вычислимо
  T1_needed = a₁[r+1] - T2₂             — знаем цель
  W₂[r] = T1_needed - h₂ - Σ1(e₂) - Ch - K[r]  — одна формула

Break на r=3 оптимален:
  δW[0]=δW[1]=δW[2]=0 → δW[16]=0 ТОЧНО (все deps нулевые)

Тестируем: a-repair в нашем измерении с backward asymmetry.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF
K = [
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
    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e, f, g)), K[r_idx]), W_r)
    T2 = add32(Sigma0(a), Maj(a, b, c))
    return (add32(T1, T2), a, b, c, add32(d, T1), e, f, g)

def expand_schedule(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    return W

def state_diff(s1, s2):
    return sum(hw(s1[i] ^ s2[i]) for i in range(8))

reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']


def a_repair_W(state2, target_a, r_idx):
    """
    a-repair: вычисляем W₂[r] чтобы a₂[r+1] = target_a.

    a_new = T1 + T2
    T2 = Σ0(a₂) + Maj(a₂,b₂,c₂)    — известно из state2
    T1_needed = target_a - T2
    W₂[r] = T1_needed - h₂ - Σ1(e₂) - Ch(e₂,f₂,g₂) - K[r]
    """
    a, b, c, d, e, f, g, h = state2
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1_needed = sub32(target_a, T2)
    W = sub32(sub32(sub32(sub32(T1_needed, h), Sigma1(e)), Ch(e, f, g)), K[r_idx])
    return W


def test_a_repair():
    np.random.seed(42)

    print("=" * 70)
    print("A-REPAIR: починка ПРИЧИНЫ (δa=0) вместо симптома (δe=0)")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. A-REPAIR С BREAK НА r=3")
    print("=" * 70)

    n_trials = 500
    all_profiles = []
    all_both_zero = []
    all_dw16 = []

    for trial in range(n_trials):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)

        # Trace 1: полный forward
        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        # Trace 2: a-repair
        # δW[0]=δW[1]=δW[2]=0 (до break), break на r=3
        W2 = list(W1)  # W2[0..2] = W1[0..2]
        W2[3] = W1[3] ^ (1 << (trial % 32))  # break: 1 бит в W[3]

        s2 = tuple(IV)
        # Раунды 0-2: одинаковые W → одинаковые states
        for r in range(3):
            s2 = R(s2, W2[r], r)

        # Раунд 3: break (δW[3] ≠ 0)
        s2 = R(s2, W2[3], 3)

        # Раунды 4-15: a-repair (W₂[r] подобран чтобы a₂=a₁)
        for r in range(4, 16):
            target_a = states1[r + 1][0]  # a₁[r+1]
            W2[r] = a_repair_W(s2, target_a, r)
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)

        # Раунды 16-63: schedule
        for r in range(16, 64):
            s2 = R(s2, W2f[r], r)

        # Profile
        profile = []
        s2_check = tuple(IV)
        for r in range(64):
            w = W2[r] if r < 16 else W2f[r]
            s2_check = R(s2_check, w, r)

            ds = state_diff(states1[r + 1], s2_check)
            da = hw(states1[r + 1][0] ^ s2_check[0])
            de = hw(states1[r + 1][4] ^ s2_check[4])
            both = 1 if ds == 0 else 0
            profile.append({'r': r + 1, 'ds': ds, 'da': da, 'de': de, 'both': both})

        all_profiles.append(profile)

        # Count BOTH=0 rounds
        both_count = sum(p['both'] for p in profile)
        all_both_zero.append(both_count)

        # δW[16]
        dw16 = hw(W1f[16] ^ W2f[16])
        all_dw16.append(dw16)

    # Средний profile
    mean_ds = [np.mean([p[r]['ds'] for p in all_profiles]) for r in range(64)]
    mean_da = [np.mean([p[r]['da'] for p in all_profiles]) for r in range(64)]
    mean_de = [np.mean([p[r]['de'] for p in all_profiles]) for r in range(64)]
    both_pct = [np.mean([p[r]['both'] for p in all_profiles]) * 100 for r in range(64)]

    print(f"\n  {'r':>4} {'δstate':>7} {'δa':>5} {'δe':>5} {'BOTH=0%':>8}  {'Status':>12}")
    print(f"  {'-'*4} {'-'*7} {'-'*5} {'-'*5} {'-'*8}  {'-'*12}")

    for r in range(25):
        ds = mean_ds[r]
        da = mean_da[r]
        de = mean_de[r]
        bp = both_pct[r]

        status = ""
        if bp > 99: status = "★ BOTH=0"
        elif bp > 50: status = "~ partial"
        elif r < 3: status = "pre-break"
        elif r == 3: status = "BREAK"
        elif r < 12 and da < 1: status = "a-repaired"
        elif r >= 12 and bp > 90: status = "REBOOT!"

        print(f"  {r+1:4d} {ds:6.1f} {da:4.1f} {de:4.1f} {bp:7.1f}%  {status:>12}")

    # Раунды после schedule
    print(f"\n  Schedule zone:")
    for r in [16, 17, 18, 19, 20, 24, 32, 48, 63]:
        print(f"  {r+1:4d} {mean_ds[r]:6.1f} {mean_da[r]:4.1f} {mean_de[r]:4.1f} {both_pct[r]:7.1f}%")

    print(f"\n  BOTH=0 раундов: mean={np.mean(all_both_zero):.1f}, min={min(all_both_zero)}, max={max(all_both_zero)}")
    print(f"  δW[16]: mean HW={np.mean(all_dw16):.1f}")

    # HW(δH)
    dh_list = []
    for trial in range(n_trials):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)
        s1 = tuple(IV)
        for r in range(64):
            s1 = R(s1, W1f[r], r)
        H1 = tuple(add32(IV[i], s1[i]) for i in range(8))

        W2 = list(W1)
        W2[3] = W1[3] ^ (1 << (trial % 32))
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)

        s1_states = [tuple(IV)]
        s1t = tuple(IV)
        for r in range(64):
            s1t = R(s1t, W1f[r], r)
            s1_states.append(s1t)

        for r in range(4, 16):
            W2[r] = a_repair_W(s2, s1_states[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)
        for r in range(16, 64):
            s2 = R(s2, W2f[r], r)

        H2 = tuple(add32(IV[i], s2[i]) for i in range(8))
        dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))
        dh_list.append(dh)

    print(f"\n  HW(δH): mean={np.mean(dh_list):.1f}, min={min(dh_list)}, std={np.std(dh_list):.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. СРАВНЕНИЕ: a-repair vs forward fix (Wang)")
    print("=" * 70)

    # Forward fix (Wang): δW[0]=1бит, fix W[1..15]
    wang_both = []
    wang_dh = []
    for trial in range(500):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)
        s1 = tuple(IV)
        s1_states = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            s1_states.append(s1)

        W2 = list(W1)
        W2[0] ^= (1 << (trial % 32))
        s2 = tuple(IV)
        s2 = R(s2, W2[0], 0)
        for r in range(1, 16):
            # Wang: fix a' to match (same as forward fix)
            a_target = s1_states[r + 1][0]
            T2 = add32(Sigma0(s2[0]), Maj(s2[0], s2[1], s2[2]))
            T1_need = sub32(a_target, T2)
            W2[r] = sub32(sub32(sub32(sub32(T1_need, s2[7]), Sigma1(s2[4])),
                           Ch(s2[4], s2[5], s2[6])), K[r])
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)
        for r in range(16, 64):
            s2 = R(s2, W2f[r], r)

        both = sum(1 for r in range(64) if state_diff(s1_states[r+1],
            (lambda: (s_tmp := tuple(IV), [s_tmp := R(s_tmp, W2[rr] if rr < 16 else W2f[rr], rr) for rr in range(r+1)], s_tmp)[-1])()) == 0)

        # Simpler: just compute δH
        H1 = tuple(add32(IV[i], s1_states[64][i]) for i in range(8))
        H2 = tuple(add32(IV[i], s2[i]) for i in range(8))
        dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))
        wang_dh.append(dh)

    print(f"\n  {'Method':<20} {'BOTH=0 rounds':>15} {'HW(δH) mean':>12} {'HW(δH) min':>12}")
    print(f"  {'-'*20} {'-'*15} {'-'*12} {'-'*12}")
    print(f"  {'a-repair (break=3)':<20} {np.mean(all_both_zero):>14.1f} {np.mean(dh_list):>11.1f} {min(dh_list):>11}")
    print(f"  {'Wang (forward fix)':<20} {'~17':>14} {np.mean(wang_dh):>11.1f} {min(wang_dh):>11}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. REBOOT: проверяем феномен самосхождения")
    print("=" * 70)

    # a-repair фиксирует δa=0 на r=4..15
    # Вопрос: δe сходится к 0 САМО? На каком раунде?

    reboot_round = []
    for trial in range(n_trials):
        profile = all_profiles[trial]
        rebooted = False
        for r in range(4, 20):
            if profile[r]['both'] == 1 and profile[r]['ds'] == 0:
                if r > 5:  # не тривиальный
                    reboot_round.append(r + 1)
                    rebooted = True
                    break
        if not rebooted:
            reboot_round.append(-1)

    valid_reboots = [r for r in reboot_round if r > 0]
    print(f"\n  Reboot (δstate→0 после break):")
    print(f"    Occurred: {len(valid_reboots)}/{n_trials} ({len(valid_reboots)/n_trials*100:.1f}%)")
    if valid_reboots:
        print(f"    Mean round: {np.mean(valid_reboots):.1f}")
        print(f"    Min round:  {min(valid_reboots)}")
        print(f"    Max round:  {max(valid_reboots)}")

        # Distribution
        from collections import Counter
        cnt = Counter(valid_reboots)
        for r in sorted(cnt):
            bar = "█" * (cnt[r] * 50 // n_trials)
            print(f"    r={r:2d}: {cnt[r]:4d} ({cnt[r]/n_trials*100:5.1f}%)  {bar}")


if __name__ == "__main__":
    test_a_repair()
