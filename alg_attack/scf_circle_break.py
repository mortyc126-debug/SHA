#!/usr/bin/env python3
"""
AI-MATH: РАЗРЫВ КРУГА — фиксируем state[16], решаем две задачи.

Я вижу: круг (state → required_W → schedule → state).
Разрываю: фиксирую state[16]. Тогда:
  Задача A: state[0]→state[16] через W[0..15] (W свободны, 1:1)
  Задача B: state[16]→state[64] через W[16..63] (schedule из W[0..15])

Для collision: state₁[16] и state₂[16] — РАЗНЫЕ.
Но W₁[0..15] и W₂[0..15] оба порождают schedule.
Нужно: schedule₂[16..63] ведёт state₂ к тому же hash что и state₁.

Задача A решается Wang chain (sens=1, 15 нулей).
Задача B — это "gap". Но теперь я вижу gap ПО-ДРУГОМУ:

Gap = 48 слов schedule порождённых 16 словами.
Schedule ЛИНЕЕН по W[0..15] (над GF(2), без carry).
Значит: gap = 48 ЛИНЕЙНЫХ функций от 16 переменных.
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
def sub32(a,b): return (a-b)&MASK32
def hw(x): return bin(x).count('1')

def expand_real(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def sha_from_state(state, W_expanded, start_round, end_round):
    """Run SHA from given state for rounds [start_round, end_round)."""
    s = list(state)
    for r in range(start_round, end_round):
        a,b,c,d,e,f,g,h = s
        T1 = add32(h,Sig1(e),Ch(e,f,g),K[r],W_expanded[r])
        T2 = add32(Sig0(a),Maj(a,b,c))
        s = [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return s

def sha_full(W16):
    W = expand_real(W16)
    s = sha_from_state(IV, W, 0, 64)
    return [add32(IV[i],s[i]) for i in range(8)]


def inverse_round(state_after, state_before_partial, W_r, K_r):
    """Given state after round r and W[r], verify consistency."""
    a_new = state_after[0]
    # T2 = Sig0(state_before[0]) + Maj(...)
    # T1 = a_new - T2
    # But we need state_before which we might not have fully.
    pass


# ============================================================
# CORE: Split at round 16
# ============================================================
def split_at_16(W1):
    """Split SHA into two halves at round 16.
    Returns: state[16], W_expanded"""
    W = expand_real(W1)
    state16 = sha_from_state(IV, W, 0, 16)
    return state16, W


def score_from_split(W1, W2):
    """For two messages: compute hash diff using split perspective."""
    Wexp1 = expand_real(W1); Wexp2 = expand_real(W2)

    # First half: rounds 0-15 (W[0..15] direct)
    state1_16 = sha_from_state(IV, Wexp1, 0, 16)
    state2_16 = sha_from_state(IV, Wexp2, 0, 16)

    # δstate at round 16
    ds16 = sum(hw(state1_16[i] ^ state2_16[i]) for i in range(8))

    # Second half: rounds 16-63 (schedule from W[0..15])
    state1_64 = sha_from_state(state1_16, Wexp1, 16, 64)
    state2_64 = sha_from_state(state2_16, Wexp2, 16, 64)

    H1 = [add32(IV[i], state1_64[i]) for i in range(8)]
    H2 = [add32(IV[i], state2_64[i]) for i in range(8)]
    dH = sum(hw(H1[i]^H2[i]) for i in range(8))

    return ds16, dH


# ============================================================
# EXP 1: What does δstate[16] look like? (Wang gives δ≈0!)
# ============================================================
def exp1_state16():
    print("="*70)
    print("EXP 1: δstate[16] — how different are states at the split?")
    print("="*70)

    for trial in range(5):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        ds16, dH = score_from_split(W1, W2)

        # Wang chain: set W2[1..15] to cancel δe
        W2_wang = list(W1); W2_wang[0] ^= 0x80000000
        Wexp1 = expand_real(W1)
        s1 = list(IV); s2 = list(IV)
        for r in range(16):
            a1,b1,c1,d1,e1,f1,g1,h1 = s1
            T1_1 = add32(h1,Sig1(e1),Ch(e1,f1,g1),K[r],Wexp1[r])
            T2_1 = add32(Sig0(a1),Maj(a1,b1,c1))
            s1 = [add32(T1_1,T2_1),a1,b1,c1,add32(d1,T1_1),e1,f1,g1]

            a2,b2,c2,d2,e2,f2,g2,h2 = s2
            # Force T1₂ = T1₁
            W2_wang[r] = sub32(sub32(sub32(sub32(T1_1,h2),Sig1(e2)),Ch(e2,f2,g2)),K[r])
            T1_2 = T1_1  # Forced
            T2_2 = add32(Sig0(a2),Maj(a2,b2,c2))
            s2 = [add32(T1_2,T2_2),a2,b2,c2,add32(d2,T1_2),e2,f2,g2]

        ds16_wang = sum(hw(s1[i]^s2[i]) for i in range(8))

        # Now: second half with Wang-modified W
        Wexp2_wang = expand_real(W2_wang)
        state1_64 = sha_from_state(s1, Wexp1, 16, 64)
        state2_64 = sha_from_state(s2, Wexp2_wang, 16, 64)
        H1 = [add32(IV[i],state1_64[i]) for i in range(8)]
        H2 = [add32(IV[i],state2_64[i]) for i in range(8)]
        dH_wang = sum(hw(H1[i]^H2[i]) for i in range(8))

        print(f"\n  Trial {trial}:")
        print(f"    Raw MSB flip: δstate[16]={ds16}, δH={dH}")
        print(f"    T1-forced:    δstate[16]={ds16_wang}, δH={dH_wang}")

        # KEY: after T1-forcing, W2[0..15] are KNOWN.
        # The schedule Wexp2_wang[16..63] is DETERMINED.
        # Does this schedule match what's needed for collision?

        # Required W for collision at each round 16..63:
        total_sched_gap = 0
        for r in range(16, 64):
            # Required: W that makes T1₂ = T1₁ at round r
            s1_r = sha_from_state(s1 if r == 16 else sha_from_state(s1, Wexp1, 16, r),
                                  Wexp1, r, r+1) if r > 16 else s1
            # This is getting complicated. Simplified: just measure gap.
            pass

        # Instead: directly measure schedule gap
        sched_gap = sum(hw(Wexp1[r] ^ Wexp2_wang[r]) for r in range(16, 64))
        print(f"    Schedule δ[16..63]: {sched_gap} bits ({sched_gap/48:.1f}/round)")

        if ds16_wang == 0:
            print(f"    ★ T1-forcing gives δstate[16]=0! Potentials solvable!")


# ============================================================
# EXP 2: If δstate[16]=0, what determines δH?
# ============================================================
def exp2_same_state_diff_schedule():
    print("\n" + "="*70)
    print("EXP 2: SAME STATE, DIFFERENT SCHEDULE")
    print("  If state₁[16] = state₂[16] but schedule₁ ≠ schedule₂")
    print("  → δH comes ENTIRELY from schedule difference!")
    print("="*70)

    for trial in range(5):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        # Create W2 with same state[16] but different schedule
        # T1-forcing gives δstate=0 at every round 0..15
        # So state₂[16] = state₁[16] (T1-forced matching from Vision 2!)
        # But W₂[0..15] ≠ W₁[0..15] → schedule₂[16..63] ≠ schedule₁[16..63]

        W2 = list(W1); W2[0] ^= 0x80000000
        Wexp1 = expand_real(W1)

        # T1-force W2 for rounds 0..15
        s1 = list(IV); s2 = list(IV)
        for r in range(16):
            a1,b1,c1,d1,e1,f1,g1,h1 = s1
            T1_1 = add32(h1,Sig1(e1),Ch(e1,f1,g1),K[r],Wexp1[r])
            T2_1 = add32(Sig0(a1),Maj(a1,b1,c1))
            s1 = [add32(T1_1,T2_1),a1,b1,c1,add32(d1,T1_1),e1,f1,g1]

            a2,b2,c2,d2,e2,f2,g2,h2 = s2
            W2[r] = sub32(sub32(sub32(sub32(T1_1,h2),Sig1(e2)),Ch(e2,f2,g2)),K[r])
            T2_2 = add32(Sig0(a2),Maj(a2,b2,c2))
            s2 = [add32(T1_1,T2_2),a2,b2,c2,add32(d2,T1_1),e2,f2,g2]

        # Check: states should be identical!
        ds = sum(hw(s1[i]^s2[i]) for i in range(8))

        # Now expand schedules from different W[0..15]
        Wexp2 = expand_real(W2)
        sched_diff = [hw(Wexp1[r] ^ Wexp2[r]) for r in range(16, 64)]
        total_sched = sum(sched_diff)

        # Run second half
        state1_64 = sha_from_state(s1, Wexp1, 16, 64)
        state2_64 = sha_from_state(s2, Wexp2, 16, 64)
        H1 = [add32(IV[i],state1_64[i]) for i in range(8)]
        H2 = [add32(IV[i],state2_64[i]) for i in range(8)]
        dH = sum(hw(H1[i]^H2[i]) for i in range(8))

        print(f"\n  Trial {trial}:")
        print(f"    δstate[16] = {ds} {'← ZERO!' if ds==0 else ''}")
        print(f"    δschedule[16..63] = {total_sched} bits ({total_sched/48:.1f}/round)")
        print(f"    δH = {dH}")

        if ds == 0:
            print(f"    ★ PURE SCHEDULE COLLISION: same state, different schedule!")
            print(f"    Collision ⟺ find W₂[0..15] where schedule₂[16..63]")
            print(f"    produces SAME hash from SAME state₁[16].")
            print(f"    This is 48 words of freedom (schedule) vs 256 bits (hash).")

        # Show schedule diff profile
        print(f"    Schedule diff by phase:")
        for phase, s_r, e_r in [("r=16..31", 16, 32), ("r=32..47", 32, 48), ("r=48..63", 48, 64)]:
            phase_diff = sum(sched_diff[r-16] for r in range(s_r, e_r))
            print(f"      {phase}: {phase_diff} bits ({phase_diff/16:.1f}/round)")


if __name__ == '__main__':
    exp1_state16()
    exp2_same_state_diff_schedule()
