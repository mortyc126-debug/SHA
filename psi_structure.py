"""
Ψ STRUCTURE: Properties of the carry difference function.

Ψ(M) = Φ(M) ⊕ Φ(M⊕δM) = ⊕_r ΔE_r(M, δM)

Each ΔE_r = carry_correction(round r, M) ⊕ carry_correction(round r, M⊕δM)
         = how the carry at round r CHANGES when message is perturbed.

Questions:
1. How does ΔE_r behave? Random or structured?
2. Are ΔE_r at different rounds independent or correlated?
3. Does Ψ have algebraic degree < 256 as function of M?
4. What is the IMAGE SIZE of Ψ? (If < 2^256 → structured)
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def per_round_carry_diff(M, delta_word, delta_bit, R=64):
    """
    Compute ΔE_r for each round r.
    Returns list of 64 values, each = XOR of carry corrections at that round.
    """
    M2 = list(M)
    M2[delta_word] ^= (1 << delta_bit)

    states1, W1 = sha256_round_trace(M, rounds=R)
    states2, W2 = sha256_round_trace(M2, rounds=R)

    delta_E = []
    for r in range(R):
        # Carry correction at round r for M:
        # E_r = (state_out_real XOR state_out_xor) for a and e
        a1,b1,c1,d1,e1,f1,g1,h1 = states1[r]
        T1_1 = (h1 + Sig1(e1) + Ch(e1,f1,g1) + K[r] + W1[r]) & MASK
        T2_1 = (Sig0(a1) + Maj(a1,b1,c1)) & MASK
        a_real1 = (T1_1 + T2_1) & MASK
        a_xor1 = T1_1 ^ T2_1
        e_real1 = (d1 + T1_1) & MASK
        e_xor1 = d1 ^ T1_1
        E_a1 = a_real1 ^ a_xor1
        E_e1 = e_real1 ^ e_xor1

        # Same for M2
        a2,b2,c2,d2,e2,f2,g2,h2 = states2[r]
        T1_2 = (h2 + Sig1(e2) + Ch(e2,f2,g2) + K[r] + W2[r]) & MASK
        T2_2 = (Sig0(a2) + Maj(a2,b2,c2)) & MASK
        a_real2 = (T1_2 + T2_2) & MASK
        a_xor2 = T1_2 ^ T2_2
        e_real2 = (d2 + T1_2) & MASK
        e_xor2 = d2 ^ T1_2
        E_a2 = a_real2 ^ a_xor2
        E_e2 = e_real2 ^ e_xor2

        dE_a = E_a1 ^ E_a2
        dE_e = E_e1 ^ E_e2
        delta_E.append((dE_a, dE_e))

    return delta_E


def experiment_deltaE_profile():
    """How does ΔE_r behave across rounds?"""
    print("=" * 80)
    print("ΔE_r PROFILE: Per-round carry difference")
    print("=" * 80)

    N = 300
    delta_word, delta_bit = 0, 0

    # Average HW of ΔE_r at each round
    hw_per_round_a = [0.0] * 64
    hw_per_round_e = [0.0] * 64

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        dE = per_round_carry_diff(M, delta_word, delta_bit)

        for r in range(64):
            hw_per_round_a[r] += bin(dE[r][0]).count('1')
            hw_per_round_e[r] += bin(dE[r][1]).count('1')

    print(f"\n  E[HW(ΔE_a[r])] and E[HW(ΔE_e[r])] per round (N={N}):")
    for r in [0, 1, 2, 4, 8, 12, 16, 24, 32, 48, 63]:
        a_hw = hw_per_round_a[r] / N
        e_hw = hw_per_round_e[r] / N
        print(f"    r={r:>2}: ΔE_a={a_hw:>5.2f}  ΔE_e={e_hw:>5.2f}  (random=16)")


def experiment_deltaE_correlation():
    """Are ΔE_r at different rounds correlated?"""
    print("\n" + "=" * 80)
    print("ΔE CORRELATION: Are per-round carry diffs independent?")
    print("=" * 80)

    N = 500
    delta_word, delta_bit = 0, 0

    # XOR consecutive ΔE_r: if independent, HW ≈ 16
    xor_hw = [0.0] * 63

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        dE = per_round_carry_diff(M, delta_word, delta_bit)

        for r in range(63):
            xor_val = dE[r][0] ^ dE[r+1][0]
            xor_hw[r] += bin(xor_val).count('1')

    print(f"\n  E[HW(ΔE_a[r] ⊕ ΔE_a[r+1])] per round (N={N}):")
    for r in [0, 1, 4, 8, 16, 32, 62]:
        avg = xor_hw[r] / N
        print(f"    r={r:>2}→{r+1:>2}: {avg:.2f} (random=16)")


def experiment_psi_degree():
    """
    What is the algebraic degree of Ψ as function of M?

    Ψ = ⊕ ΔE_r. Each ΔE_r comes from carry differences.
    If Ψ has degree d < full → algebraic structure.

    Test: degree-1 (linearity) and degree-2 of Ψ.
    """
    print("\n" + "=" * 80)
    print("Ψ DEGREE: Algebraic degree of carry difference function")
    print("=" * 80)

    delta_word, delta_bit = 0, 0
    rng = random.Random(42)
    M_base = [rng.randint(0, MASK) for _ in range(16)]

    def psi_bit0(M):
        """Bit 0 of Ψ(M) = bit 0 of [Φ(M) ⊕ Φ(M⊕δ)]."""
        M2 = list(M); M2[delta_word] ^= (1 << delta_bit)
        states1, W1 = sha256_round_trace(M)
        states2, W2 = sha256_round_trace(M2)
        H1 = [(states1[64][i] + H0[i]) & MASK for i in range(8)]
        H2 = [(states2[64][i] + H0[i]) & MASK for i in range(8)]

        def xor_sha(msg):
            W = list(msg)
            for i in range(16, 64):
                s1 = rotr(W[i-2],17)^rotr(W[i-2],19)^(W[i-2]>>10)
                s0 = rotr(W[i-15],7)^rotr(W[i-15],18)^(W[i-15]>>3)
                W.append(s1^W[i-7]^s0^W[i-16])
            a,b,c,d,e,f,g,h = H0
            for r in range(64):
                T1 = h^Sig1(e)^Ch(e,f,g)^K[r]^W[r]
                T2 = Sig0(a)^Maj(a,b,c)
                h,g,f = g,f,e; e=d^T1; d,c,b=c,b,a; a=T1^T2
            return [a^H0[0],b^H0[1],c^H0[2],d^H0[3],e^H0[4],f^H0[5],g^H0[6],h^H0[7]]

        L1 = xor_sha(M)
        L2 = xor_sha(M2)
        phi1 = H1[0] ^ L1[0]
        phi2 = H2[0] ^ L2[0]
        return (phi1 ^ phi2) & 1

    # Degree-1 test
    N = 200
    d1_viol = 0
    d2_viol = 0
    total = 0

    for trial in range(N):
        rng2 = random.Random(trial * 1000)
        w1,b1 = rng2.randint(0,15), rng2.randint(0,31)
        w2,b2 = rng2.randint(0,15), rng2.randint(0,31)
        w3,b3 = rng2.randint(0,15), rng2.randint(0,31)
        if (w1,b1)==(w2,b2) or (w1,b1)==(w3,b3) or (w2,b2)==(w3,b3):
            continue

        def flip(M, w, b):
            M2 = list(M); M2[w] ^= (1<<b); return M2

        f0 = psi_bit0(M_base)
        f1 = psi_bit0(flip(M_base, w1, b1))
        f2 = psi_bit0(flip(M_base, w2, b2))
        f12 = psi_bit0(flip(flip(M_base, w1, b1), w2, b2))

        if f0 ^ f1 ^ f2 ^ f12:
            d1_viol += 1

        f3 = psi_bit0(flip(M_base, w3, b3))
        f13 = psi_bit0(flip(flip(M_base, w1, b1), w3, b3))
        f23 = psi_bit0(flip(flip(M_base, w2, b2), w3, b3))
        f123 = psi_bit0(flip(flip(flip(M_base, w1, b1), w2, b2), w3, b3))

        if f0^f1^f2^f3^f12^f13^f23^f123:
            d2_viol += 1

        total += 1

    print(f"  Ψ bit 0 degree test (N={total}):")
    print(f"    Degree-1 violations: {d1_viol}/{total} = {d1_viol/total:.3f}")
    print(f"    Degree-2 violations: {d2_viol}/{total} = {d2_viol/total:.3f}")

    if d1_viol/total > 0.4:
        print(f"    Ψ is HIGH DEGREE (≥ 3)")
    elif d2_viol/total < 0.1:
        print(f"    Ψ is nearly degree 2!")
    elif d1_viol/total < 0.1:
        print(f"    Ψ is nearly linear!")


if __name__ == "__main__":
    experiment_deltaE_profile()
    experiment_deltaE_correlation()
    experiment_psi_degree()
