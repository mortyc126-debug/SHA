#!/usr/bin/env python3
"""
SCF: UNIFIED SHA SOLVER — построен на 6 знаниях о SHA-256.

НЕ использует: SA, random search, greedy bit-flip, birthday.
ИСПОЛЬЗУЕТ: exact round inversion, Newton per-round, schedule
shaping, shift register tracking, Wang chain propagation.

Архитектура:
  Phase 1: Wang Chain (r=0..16) — δe=0 для 15 раундов, FREE
  Phase 2: Per-Round Newton through gap — для каждого r, вычислить
           ТРЕБУЕМЫЙ W2[r] и ФАКТИЧЕСКИЙ W2[r], скорректировать
           через schedule shaping (изменить W2[0..15] чтобы
           schedule[r] стал ближе к required[r])
  Phase 3: 2D Tracking — используем shift register для получения
           "бесплатных" нулей от каждого solved раунда

Ключевое отличие от старых методов: мы НЕ ИЩЕМ — мы ВЫЧИСЛЯЕМ.
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

def sha_round_single(state, W_r, K_r):
    a,b,c,d,e,f,g,h = state
    T1 = add32(h,Sig1(e),Ch(e,f,g),K_r,W_r)
    T2 = add32(Sig0(a),Maj(a,b,c))
    return [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

def sha_all_states(W16):
    W=expand_real(W16); s=list(IV); states=[list(s)]
    for r in range(64):
        s = sha_round_single(s, W[r], K[r])
        states.append(list(s))
    return states, W

def sha_compress(W16):
    st, _ = sha_all_states(W16)
    return [add32(IV[i], st[64][i]) for i in range(8)]


# ============================================================
# Tool 1: EXACT ROUND INVERTER
# Given states at round r for both messages, compute required W2[r]
# ============================================================
def invert_round_e(state1_r, state2_r, e1_target, r):
    """Compute W2[r] needed to make e2[r+1] = e1_target."""
    a2,b2,c2,d2,e2,f2,g2,h2 = state2_r
    # e_new = d + T1, so T1_needed = e1_target - d2
    T1_needed = sub32(e1_target, d2)
    # T1 = h + Σ1(e) + Ch(e,f,g) + K[r] + W[r]
    W2_needed = sub32(sub32(sub32(sub32(T1_needed, h2), Sig1(e2)), Ch(e2,f2,g2)), K[r])
    return W2_needed


# ============================================================
# Tool 2: SCHEDULE SHAPER
# Given target W_expanded[r], find best W[0..15] adjustment
# ============================================================
def schedule_sensitivity(W16, target_r):
    """For each word w, compute ∂W_expanded[target_r]/∂W[w]."""
    Wexp = expand_real(W16)
    sens = []
    for w in range(16):
        W_test = list(W16)
        W_test[w] = add32(W16[w], 1)
        Wexp_test = expand_real(W_test)
        s = sub32(Wexp_test[target_r], Wexp[target_r])
        sens.append(s)
    return sens


def shape_schedule(W16, target_r, target_val):
    """Adjust W16 so that expanded_schedule[target_r] ≈ target_val.
    Uses schedule sensitivity to pick best word and correction."""
    Wexp = expand_real(W16)
    current_val = Wexp[target_r]
    mismatch = sub32(target_val, current_val)

    if mismatch == 0:
        return W16, 0  # Already correct

    sens = schedule_sensitivity(W16, target_r)

    # Find word with invertible sensitivity
    best_word = -1
    best_residual_hw = 32

    for w in range(16):
        if sens[w] & 1 == 0:
            continue  # Not invertible

        inv = pow(sens[w], -1, 1 << 32)
        correction = (mismatch * inv) & MASK32

        W_corrected = list(W16)
        W_corrected[w] = add32(W16[w], correction)
        Wexp_c = expand_real(W_corrected)
        residual = sub32(Wexp_c[target_r], target_val)

        if hw(residual) < best_residual_hw:
            best_residual_hw = hw(residual)
            best_word = w
            best_W = list(W_corrected)

    if best_word >= 0:
        return best_W, best_residual_hw
    return W16, hw(mismatch)


# ============================================================
# Tool 3: WANG CHAIN PROPAGATOR
# Compute Wang chain δW for first 16 rounds
# ============================================================
def wang_chain(W1, delta_w0):
    """Build W2 with Wang chain: δe[2..16]=0."""
    W2 = list(W1)
    W2[0] ^= delta_w0

    states1, Wexp1 = sha_all_states(W1)
    states2, Wexp2 = sha_all_states(W2)

    # Wang: for r=1..15, set δW[r] = -δe_natural[r+1]
    # sensitivity = 1 when δe[r-1..r-3]=0
    for r in range(1, 16):
        s1 = sha_all_states(W1)[0]
        s2_current, Wexp2_current = sha_all_states(W2)

        # δe_natural at r+1 (what happens if we don't correct)
        de_nat = sub32(s2_current[r+1][4], s1[r+1][4])

        if de_nat != 0:
            # Correct: δW[r] -= de_nat (additive Wang)
            W2[r] = sub32(W2[r], de_nat)

    return W2


# ============================================================
# Tool 4: 2D RECURRENCE TRACKER
# Track which rounds have δe=0 and propagate through shift register
# ============================================================
def track_zeros(W1, W2):
    """Return list of rounds where δe=0 and δa=0."""
    st1, _ = sha_all_states(W1)
    st2, _ = sha_all_states(W2)

    de_zeros = []
    da_zeros = []
    for r in range(65):
        if st1[r][4] == st2[r][4]: de_zeros.append(r)
        if st1[r][0] == st2[r][0]: da_zeros.append(r)

    return de_zeros, da_zeros


# ============================================================
# UNIFIED SOLVER: Phase 1 (Wang) + Phase 2 (Gap) + Phase 3 (Track)
# ============================================================
def unified_solve(W1, delta_w0=0x80000000, gap_iterations=5):
    """
    Full pipeline:
    1. Wang chain → δe[2..16]=0
    2. For each gap round: compute required W2[r], shape schedule
    3. Track zeros, iterate
    """
    print(f"  Phase 1: Wang Chain")
    W2 = wang_chain(W1, delta_w0)

    de_z, da_z = track_zeros(W1, W2)
    print(f"    δe=0 rounds: {len(de_z)} → {de_z[:20]}{'...' if len(de_z)>20 else ''}")
    print(f"    δa=0 rounds: {len(da_z)} → {da_z[:20]}{'...' if len(da_z)>20 else ''}")

    H1 = sha_compress(W1); H2 = sha_compress(W2)
    hw_init = sum(hw(H1[i]^H2[i]) for i in range(8))
    print(f"    Hash diff after Wang: HW = {hw_init}")

    # Phase 2: Gap solving
    print(f"\n  Phase 2: Gap Round Solving")
    best_hw = hw_init
    best_W2 = list(W2)

    for iteration in range(gap_iterations):
        st1, Wexp1 = sha_all_states(W1)
        st2, Wexp2 = sha_all_states(W2)

        # For each gap round, compute required vs actual
        # Pick the round with SMALLEST mismatch → shape it
        mismatches = []
        for r in range(17, 52):
            e1_target = st1[r+1][4]  # Want e2[r+1] = e1[r+1]
            W2_needed = invert_round_e(st1[r], st2[r], e1_target, r)
            actual = Wexp2[r]
            mm_hw = hw(W2_needed ^ actual)
            mismatches.append((mm_hw, r, W2_needed))

        mismatches.sort()

        if not mismatches:
            break

        # Try to shape the easiest round
        shaped = False
        for mm_hw, target_r, W2_needed in mismatches[:5]:
            W2_trial, residual = shape_schedule(W2, target_r, W2_needed)

            if W2_trial == W1:
                continue  # Skip trivial

            H2_trial = sha_compress(W2_trial)
            trial_hw = sum(hw(H1[i]^H2_trial[i]) for i in range(8))

            if trial_hw < best_hw:
                best_hw = trial_hw
                best_W2 = list(W2_trial)
                W2 = W2_trial
                shaped = True

                de_z, da_z = track_zeros(W1, W2)
                print(f"    Iter {iteration}: shaped r={target_r} "
                      f"(mm={mm_hw}→{residual}), HW(δH)={trial_hw}, "
                      f"δe_zeros={len(de_z)}")
                break

        if not shaped:
            # Try shaping for δa instead
            for r in range(17, 52):
                a1_target = st1[r+1][0]
                # a_new = T1 + T2, more complex to invert
                # Skip for now — δe is the primary target
                pass
            print(f"    Iter {iteration}: no improvement found")
            break

    # Phase 3: Final tracking
    print(f"\n  Phase 3: Final State")
    de_z, da_z = track_zeros(W1, best_W2)
    H2_final = sha_compress(best_W2)
    final_hw = sum(hw(H1[i]^H2_final[i]) for i in range(8))
    n_diff = sum(1 for i in range(16) if W1[i] != best_W2[i])

    print(f"    δe=0 rounds: {len(de_z)}")
    print(f"    δa=0 rounds: {len(da_z)}")
    print(f"    Hash diff: HW = {final_hw}")
    print(f"    Words changed: {n_diff}/16")

    # Show per-register hash diff
    per_reg = [hw(H1[i]^H2_final[i]) for i in range(8)]
    regs = "abcdefgh"
    print(f"    Per-register: " + " ".join(f"{regs[i]}={per_reg[i]}" for i in range(8)))

    return best_W2, final_hw


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    print("="*70)
    print("SCF: UNIFIED SHA SOLVER")
    print("  Phase 1: Wang Chain (exact, P=1)")
    print("  Phase 2: Per-Round Exact Inversion + Schedule Shaping")
    print("  Phase 3: 2D Zero Tracking")
    print("="*70)

    results = []

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        print(f"\n{'='*50}")
        print(f"Trial {trial}")

        W2, final_hw = unified_solve(W1, gap_iterations=10)
        results.append(final_hw)

    print(f"\n{'='*70}")
    print(f"SUMMARY ({N} trials)")
    print(f"  Avg HW(δH): {sum(results)/len(results):.1f}")
    print(f"  Min HW(δH): {min(results)}")
    print(f"  Max HW(δH): {max(results)}")
