#!/usr/bin/env python3
"""
SCF: y-BASIS GAP TEST — выживает ли информация через THE GAP в y-базисе?

Факт: в y-базисе (y=x+1, R=GF(2)[y]/(y^32)):
  - ROTR линеен (умножение на (1+y)^r)
  - v(carry)>0 в 55% случаев → y^0 компонента защищена
  - Σ-операторы: v(S0)=v(S1)=0 → сохраняют y^0

Вопрос: если δstate[17] ≡ 0 по y^0 компоненте,
сколько раундов gap это выживает?

Если chain_length >> 1 в y-базисе (vs chain=1 в mod-2 башне)
→ y-фильтрация = правильная алгебра для bridging the gap!
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

def to_y_basis(word):
    """x-basis → y-basis (y = x+1)."""
    result = 0
    for k in range(32):
        coeff = 0
        for i in range(32):
            if (word >> i) & 1:
                c = 1
                ii, kk = i, k
                while kk > 0:
                    if (kk & 1) > (ii & 1):
                        c = 0; break
                    ii >>= 1; kk >>= 1
                coeff ^= c
        result |= (coeff << k)
    return result

def y_component(word, level):
    """Extract y^level coefficient from y-basis representation."""
    yb = to_y_basis(word)
    return (yb >> level) & 1

def y_low_components(word, n_levels):
    """Extract y^0, y^1, ..., y^{n_levels-1} as an n_levels-bit value."""
    yb = to_y_basis(word)
    return yb & ((1 << n_levels) - 1)

def y_state_component(state, level):
    """y^level component of full 8-register state → 8-bit value."""
    return sum(y_component(state[i], level) << i for i in range(8))

def y_state_low(state, n_levels):
    """Low y-components of full state → 8*n_levels bits."""
    result = 0
    for reg in range(8):
        low = y_low_components(state[reg], n_levels)
        result |= low << (reg * n_levels)
    return result


# ============================================================
# EXP 1: y^0 chain length through THE GAP
# ============================================================
def exp1_y0_chain(N):
    print("="*70)
    print("EXP 1: y^0 CHAIN LENGTH THROUGH THE GAP")
    print("  y^0(word) = parity = XOR of all bits")
    print("  y^0(state) = 8-bit parity vector")
    print("  Chain = consecutive rounds where δ(y^0(state)) = 0")
    print("="*70)

    chains_from_17 = []
    chains_from_0 = []

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= (trial + 1) | 1  # Ensure nonzero delta

        st1 = sha_states(W1)
        st2 = sha_states(W2)

        # y^0 component at each round
        y0_match = []
        for r in range(65):
            c1 = y_state_component(st1[r], 0)
            c2 = y_state_component(st2[r], 0)
            y0_match.append(c1 == c2)

        # Chain from r=17
        chain = 0
        if y0_match[17]:
            for r in range(18, 65):
                if y0_match[r]:
                    chain += 1
                else:
                    break
        chains_from_17.append(chain)

        # Chain from r=0
        chain0 = 0
        for r in range(0, 65):
            if y0_match[r]:
                chain0 += 1
            else:
                break
        chains_from_0.append(chain0)

    # Filter: only cases where y^0 matches at r=17
    valid = [c for i, c in enumerate(chains_from_17)
             if sum(y_state_component(sha_states([int.from_bytes(os.urandom(4),'big') for _ in range(16)])[17], 0) == 0 for _ in [0])]

    # Simple stats
    n_start = sum(1 for i in range(N) if chains_from_17[i] >= 0)
    # Count how many had y0 match at r=17
    n_y0_17 = 0
    chain_when_matched = []

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= (trial * 7 + 3)
        if W2[0] == W1[0]: W2[0] ^= 1

        st1 = sha_states(W1)
        st2 = sha_states(W2)

        c1_17 = y_state_component(st1[17], 0)
        c2_17 = y_state_component(st2[17], 0)

        if c1_17 == c2_17:
            n_y0_17 += 1
            # Count chain
            chain = 0
            for r in range(18, 65):
                if y_state_component(st1[r], 0) == y_state_component(st2[r], 0):
                    chain += 1
                else:
                    break
            chain_when_matched.append(chain)

    print(f"\n  y^0 match at r=17: {n_y0_17}/{N} = {n_y0_17/N:.3f}")
    print(f"  (expected if random: {1/256:.4f})")

    if chain_when_matched:
        avg = sum(chain_when_matched) / len(chain_when_matched)
        mx = max(chain_when_matched)
        bridges = sum(1 for c in chain_when_matched if c >= 34)
        print(f"\n  Chain length from r=17 (when y^0 matches):")
        print(f"    Avg: {avg:.1f}")
        print(f"    Max: {mx}")
        print(f"    Bridges gap (≥34): {bridges}/{len(chain_when_matched)}")

        from collections import Counter
        dist = Counter(chain_when_matched)
        print(f"    Distribution: ", end="")
        for k in sorted(dist.keys())[:15]:
            print(f"{k}:{dist[k]} ", end="")
        print()
    else:
        print("  No y^0 matches at r=17 found")

    return chain_when_matched


# ============================================================
# EXP 2: Multi-level y chain (y^0 through y^k)
# ============================================================
def exp2_multilevel(N):
    print("\n" + "="*70)
    print("EXP 2: MULTI-LEVEL y-CHAIN (y^0 through y^k)")
    print("  Match first k y-components simultaneously")
    print("="*70)

    for n_levels in [1, 2, 3, 4]:
        output_bits = 8 * n_levels
        random_match = 2**(-output_bits)

        n_match_17 = 0
        chains = []

        for trial in range(N * 5):
            W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            W2 = list(W1)
            W2[0] ^= int.from_bytes(os.urandom(4),'big')
            if W2[0] == W1[0]: continue

            st1 = sha_states(W1)
            st2 = sha_states(W2)

            # Check y-low match at r=17
            yl1 = y_state_low(st1[17], n_levels)
            yl2 = y_state_low(st2[17], n_levels)

            if yl1 == yl2:
                n_match_17 += 1
                chain = 0
                for r in range(18, 65):
                    if y_state_low(st1[r], n_levels) == y_state_low(st2[r], n_levels):
                        chain += 1
                    else:
                        break
                chains.append(chain)

                if len(chains) >= min(N, 50):
                    break

        tested = N * 5
        print(f"\n  y-levels 0..{n_levels-1} ({output_bits} bits):")
        print(f"    Matches at r=17: {n_match_17}/{tested} = {n_match_17/tested:.6f}")
        print(f"    Expected random: {random_match:.6f}")

        if n_match_17 > 0:
            ratio = (n_match_17/tested) / random_match
            print(f"    Ratio vs random: {ratio:.1f}x")

        if chains:
            avg = sum(chains)/len(chains)
            mx = max(chains)
            print(f"    Chain: avg={avg:.1f}, max={mx} ({len(chains)} samples)")
            bridges = sum(1 for c in chains if c >= 34)
            print(f"    Bridges gap: {bridges}/{len(chains)}")
        else:
            print(f"    No matches found (need more samples)")


# ============================================================
# EXP 3: y^0 propagation — ONE STEP survival
# ============================================================
def exp3_one_step(N):
    print("\n" + "="*70)
    print("EXP 3: y^0 ONE-STEP SURVIVAL IN THE GAP")
    print("  P(y^0 match at r+1 | y^0 match at r)")
    print("="*70)

    per_round = {r: [0, 0] for r in range(17, 64)}  # [total, survived]

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= int.from_bytes(os.urandom(4),'big')
        if W2[0] == W1[0]: continue

        st1 = sha_states(W1)
        st2 = sha_states(W2)

        for r in range(17, 64):
            c1r = y_state_component(st1[r], 0)
            c2r = y_state_component(st2[r], 0)
            if c1r == c2r:
                per_round[r][0] += 1
                c1r1 = y_state_component(st1[r+1], 0)
                c2r1 = y_state_component(st2[r+1], 0)
                if c1r1 == c2r1:
                    per_round[r][1] += 1

    random_p = 1/256
    print(f"\n  {'r':>4} | P(surv) | total | ratio_vs_random")
    print("  " + "-"*50)
    for r in range(17, 64, 3):
        total = per_round[r][0]
        surv = per_round[r][1]
        if total > 5:
            p = surv / total
            ratio = p / random_p
            marker = " ★★" if ratio > 5 else (" ★" if ratio > 2 else "")
            print(f"  {r:4d} | {p:.4f} | {total:5d} | {ratio:6.1f}x{marker}")
        else:
            print(f"  {r:4d} | - (only {total} samples)")


# ============================================================
# EXP 4: Compare y-basis vs x-basis (mod 2) chain
# ============================================================
def exp4_compare(N):
    print("\n" + "="*70)
    print("EXP 4: y-BASIS vs x-BASIS (mod 2) — DIRECT COMPARISON")
    print("  Same pairs, measure chain in both bases")
    print("="*70)

    y_chains = []
    x_chains = []

    for trial in range(N * 10):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= (trial + 1)

        st1 = sha_states(W1)
        st2 = sha_states(W2)

        # y^0 chain from r=0
        yc = 0
        for r in range(65):
            if y_state_component(st1[r], 0) == y_state_component(st2[r], 0):
                yc += 1
            else:
                break

        # x mod 2 chain from r=0
        xc = 0
        for r in range(65):
            if all((st1[r][i] & 1) == (st2[r][i] & 1) for i in range(8)):
                xc += 1
            else:
                break

        y_chains.append(yc)
        x_chains.append(xc)

    avg_y = sum(y_chains) / len(y_chains)
    avg_x = sum(x_chains) / len(x_chains)
    max_y = max(y_chains)
    max_x = max(x_chains)

    print(f"\n  y-basis (y^0 = parity):")
    print(f"    Avg chain from r=0: {avg_y:.2f}")
    print(f"    Max chain: {max_y}")

    print(f"\n  x-basis (mod 2 = LSB):")
    print(f"    Avg chain from r=0: {avg_x:.2f}")
    print(f"    Max chain: {max_x}")

    print(f"\n  y-basis advantage: {avg_y - avg_x:+.2f} rounds avg, {max_y - max_x:+d} max")

    if avg_y > avg_x + 0.5:
        print(f"  ★ y-basis is BETTER than x-basis!")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 300

    print("="*70)
    print("SCF: y-BASIS GAP BRIDGE TEST")
    print("="*70)

    exp4_compare(N)
    exp3_one_step(min(N, 500))
    exp1_y0_chain(N)
    exp2_multilevel(min(N, 200))

    print("\n" + "="*70)
    print("ИТОГ")
    print("="*70)
