#!/usr/bin/env python3
"""
SCF: BRIDGE ALGEBRA — строим алгебру с нуля для THE GAP.

Идея: p-адическая фильтрация. SHA mod 2^k для k=1,2,...,32.
На уровне k: carry ограничен k-1 позициями.
Вопрос: на каком уровне k информация ВЫЖИВАЕТ через gap (r=17..51)?

Если найдём k* где collision-information сохраняется:
→ решаем collision в Z/2^{k*} (пространство 2^{8k*})
→ поднимаем Hensel-лифтом на следующий уровень
→ поэтапный подъём до полных 32 бит

Это НОВАЯ конструкция: не атакуем SHA целиком, а решаем ПОСЛОЙНО.
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

def sha_states(W16):
    W=expand_real(W16); s=list(IV); states=[list(s)]
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        states.append(list(s))
    return states


# ============================================================
# LEVEL 1: Информация на уровне mod 2^k
# ============================================================
def exp1_information_survival(N):
    """
    Для пары (W, W⊕δ): на каком уровне mod 2^k
    δstate[r] содержит exploitable информацию?

    Мерим: P(δstate[r] ≡ 0 mod 2^k) для k=1,2,...,16
    Если P >> 2^{-8k} → информация о collision СОХРАНЯЕТСЯ на уровне k.
    """
    print("="*70)
    print("LEVEL 1: INFORMATION SURVIVAL — P(δstate ≡ 0 mod 2^k) по раундам")
    print("="*70)

    results = {}

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= 0x80000000

        st1 = sha_states(W1)
        st2 = sha_states(W2)

        for r in range(65):
            if r not in results:
                results[r] = {k: 0 for k in [1,2,3,4,6,8,12,16]}

            for k in results[r]:
                mask = (1 << k) - 1
                if all((st1[r][i] & mask) == (st2[r][i] & mask) for i in range(8)):
                    results[r][k] += 1

    print(f"\n  P(δstate ≡ 0 mod 2^k) — higher = more info survives")
    print(f"  Expected random: P = 2^{{-8k}}")
    print(f"  {'r':>3} |", end="")
    for k in [1,2,3,4,6,8,12,16]:
        print(f" {'mod'+str(2**k):>8}", end="")
    print()
    print("  " + "-"*75)

    for r in [0,1,2,5,10,15,16,17,18,20,25,30,35,40,45,50,55,60,64]:
        print(f"  {r:3d} |", end="")
        for k in [1,2,3,4,6,8,12,16]:
            p = results[r][k] / N
            expected = 2**(-8*k) if r > 5 else 1.0

            if p > 0:
                marker = "★" if p > expected * 10 else " "
            else:
                marker = " "

            print(f" {p:7.4f}{marker}", end="")
        print()


# ============================================================
# LEVEL 2: Hensel compatibility — can we lift mod 2^k collision?
# ============================================================
def exp2_hensel_lift_gap(N):
    """
    Start: find (W, W') with δstate[17] ≡ 0 mod 2^k (Wang-like).
    Question: does this mod-2^k collision SURVIVE through the gap to r=51?

    If yes at level k → Hensel approach works for the gap!
    """
    print("\n" + "="*70)
    print("LEVEL 2: HENSEL LIFT THROUGH THE GAP")
    print("  Start from δstate[17] ≡ 0 mod 2^k")
    print("  Measure: P(δstate[r] ≡ 0 mod 2^k) for r=17..64")
    print("="*70)

    for k in [1, 2, 3, 4]:
        mask = (1 << k) - 1
        print(f"\n  --- mod 2^{k} (mask = 0x{mask:x}) ---")

        # Find pairs where δstate[17] ≡ 0 mod 2^k
        surviving = {r: 0 for r in range(17, 65)}
        n_valid = 0

        for trial in range(N * 10):
            W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            W2 = list(W1)
            W2[0] ^= (trial + 1)

            st1 = sha_states(W1)
            st2 = sha_states(W2)

            # Check if δstate[17] ≡ 0 mod 2^k
            if all((st1[17][i] & mask) == (st2[17][i] & mask) for i in range(8)):
                n_valid += 1
                for r in range(17, 65):
                    if all((st1[r][i] & mask) == (st2[r][i] & mask) for i in range(8)):
                        surviving[r] += 1

                if n_valid >= N:
                    break

        if n_valid == 0:
            print(f"  No pairs found with δstate[17] ≡ 0 mod 2^{k}")
            continue

        print(f"  Found {n_valid} pairs with δstate[17] ≡ 0 mod 2^{k}")
        print(f"  {'r':>4} | P(survives) | bar")
        print("  " + "-"*40)

        for r in range(17, 65, 2):
            p = surviving[r] / n_valid
            bar = "█" * int(p * 40)
            marker = " ★" if p > 0.1 and r > 20 else ""
            print(f"  {r:4d} | {p:11.4f} | {bar}{marker}")


# ============================================================
# LEVEL 3: Truncated SHA collision — solve mod 2^k
# ============================================================
def exp3_truncated_collision(N):
    """
    SHA mod 2^k: find collision for low k.
    Cost should be O(2^{4k}) by birthday (output = 8k bits).
    If we find it cheaper → structure helps!
    """
    print("\n" + "="*70)
    print("LEVEL 3: TRUNCATED COLLISION — SHA mod 2^k")
    print("  Output size = 8k bits. Birthday = O(2^{4k})")
    print("="*70)

    for k in [1, 2, 3, 4]:
        mask = (1 << k) - 1
        output_bits = 8 * k
        birthday = 1 << (4 * k)

        # Search for truncated collision
        seen = {}
        found = False
        tries = 0

        for trial in range(min(N * birthday, 500000)):
            W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            st = sha_states(W)
            # Final state mod 2^k
            final_trunc = tuple(st[64][i] & mask for i in range(8))

            tries += 1

            if final_trunc in seen:
                W_prev = seen[final_trunc]
                if W_prev != W:
                    # Verify it's a real truncated collision
                    st_prev = sha_states(W_prev)
                    if tuple(st_prev[64][i] & mask for i in range(8)) == final_trunc:
                        # Full diff
                        full_diff = sum(hw(st[64][i] ^ st_prev[64][i]) for i in range(8))
                        print(f"  mod 2^{k}: COLLISION at try {tries} "
                              f"(birthday={birthday}, ratio={tries/birthday:.2f}x) "
                              f"full_diff_HW={full_diff}")
                        found = True
                        break
            else:
                seen[final_trunc] = W

        if not found:
            print(f"  mod 2^{k}: No collision in {tries} tries (birthday={birthday})")


# ============================================================
# LEVEL 4: Conditional collision propagation
# ============================================================
def exp4_conditional_propagation(N):
    """
    Key question: if δstate[r] ≡ 0 mod 2^k, what is the
    probability that δstate[r+1] ≡ 0 mod 2^k?

    This is the ONE-STEP survival probability for level k.
    If P_survive > 2^{-8k} → information persists!
    """
    print("\n" + "="*70)
    print("LEVEL 4: ONE-STEP SURVIVAL PROBABILITY")
    print("  P(δstate[r+1] ≡ 0 mod 2^k | δstate[r] ≡ 0 mod 2^k)")
    print("="*70)

    for k in [1, 2, 3, 4]:
        mask = (1 << k) - 1
        random_prob = 2**(-8*k)

        # For each round in the gap, measure conditional probability
        cond_match = {r: 0 for r in range(17, 64)}
        cond_total = {r: 0 for r in range(17, 64)}

        for trial in range(N * 5):
            W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            W2 = list(W1)
            W2[0] ^= int.from_bytes(os.urandom(4),'big')

            st1 = sha_states(W1)
            st2 = sha_states(W2)

            for r in range(17, 64):
                # Check if δstate[r] ≡ 0 mod 2^k
                if all((st1[r][i] & mask) == (st2[r][i] & mask) for i in range(8)):
                    cond_total[r] += 1
                    # Check if δstate[r+1] ≡ 0 mod 2^k too
                    if all((st1[r+1][i] & mask) == (st2[r+1][i] & mask) for i in range(8)):
                        cond_match[r] += 1

        print(f"\n  --- mod 2^{k} (random P = {random_prob:.2e}) ---")
        print(f"  {'r':>4} | P(surv) | ratio_vs_random | samples")
        print("  " + "-"*50)

        for r in range(17, 64, 3):
            if cond_total[r] > 0:
                p = cond_match[r] / cond_total[r]
                ratio = p / random_prob if random_prob > 0 else 0
                marker = " ★★" if ratio > 5 else (" ★" if ratio > 2 else "")
                print(f"  {r:4d} | {p:.4f} | {ratio:>14.1f}x | {cond_total[r]:5d}{marker}")
            else:
                print(f"  {r:4d} | --- no samples ---")


# ============================================================
# LEVEL 5: Multi-step propagation — does mod-2^k collision chain?
# ============================================================
def exp5_chain_length(N):
    """
    Starting from δstate[17] ≡ 0 mod 2^k:
    How many consecutive rounds does the mod-2^k match SURVIVE?
    This is the CHAIN LENGTH — the key metric for bridge algebra.
    """
    print("\n" + "="*70)
    print("LEVEL 5: CHAIN LENGTH — consecutive mod-2^k survival")
    print("  Start: δstate[17] ≡ 0 mod 2^k")
    print("  Count: how many rounds r=18,19,... continue to match")
    print("="*70)

    for k in [1, 2, 3, 4]:
        mask = (1 << k) - 1
        chain_lengths = []

        for trial in range(N * 20):
            W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            W2 = list(W1)
            W2[0] ^= (trial + 1)

            st1 = sha_states(W1)
            st2 = sha_states(W2)

            # Check entry condition
            if not all((st1[17][i] & mask) == (st2[17][i] & mask) for i in range(8)):
                continue

            # Count chain
            chain = 0
            for r in range(18, 65):
                if all((st1[r][i] & mask) == (st2[r][i] & mask) for i in range(8)):
                    chain += 1
                else:
                    break

            chain_lengths.append(chain)
            if len(chain_lengths) >= N:
                break

        if not chain_lengths:
            print(f"  mod 2^{k}: no valid starting pairs found")
            continue

        avg = sum(chain_lengths) / len(chain_lengths)
        mx = max(chain_lengths)
        print(f"\n  mod 2^{k} ({len(chain_lengths)} chains):")
        print(f"    Avg chain: {avg:.1f} rounds")
        print(f"    Max chain: {mx} rounds")
        print(f"    Chain ≥ 34 (bridges gap): {sum(1 for c in chain_lengths if c >= 34)}/{len(chain_lengths)}")

        # Distribution
        from collections import Counter
        dist = Counter(chain_lengths)
        print(f"    Distribution: ", end="")
        for length in sorted(dist.keys())[:10]:
            print(f"{length}:{dist[length]} ", end="")
        if len(dist) > 10:
            print("...", end="")
        print()


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200

    print("="*70)
    print("SCF: BRIDGE ALGEBRA — p-ADIC FILTRATION ЧЕРЕЗ THE GAP")
    print("="*70)

    exp1_information_survival(min(N, 500))
    exp4_conditional_propagation(min(N, 200))
    exp5_chain_length(N)
    exp3_truncated_collision(N)
    exp2_hensel_lift_gap(min(N, 100))

    print("\n" + "="*70)
    print("ИТОГ: BRIDGE ALGEBRA")
    print("="*70)
    print("""
  Фильтрационная алгебра: SHA mod 2^k для k = 1, 2, ..., 32.

  Для каждого уровня k:
    - Collision space: 8k бит → birthday O(2^{4k})
    - Carry complexity: k-1 бит carry (vs 31 для полного SHA)
    - Chain length: сколько раундов gap mod-2^k collision выживает

  ЕСЛИ chain_length(k*) ≥ 34:
    → mod 2^{k*} collision пересекает весь gap
    → Решаем collision mod 2^{k*}: стоимость O(2^{4k*})
    → Hensel-lift до mod 2^{k*+1}: стоимость O(2^{4(k*+1)})
    → ... до mod 2^32: стоимость O(2^{128})
    → НО: если лифт работает, промежуточные шаги O(2^{4k})
    → Total: O(32 × 2^{4k*}) вместо O(2^128)
    → Для k*=4: O(2^{21}) — ПРОРЫВ!
    → Для k*=8: O(2^{37}) — всё ещё прорыв!
    """)
