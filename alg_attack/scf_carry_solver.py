#!/usr/bin/env python3
"""
SCF: CARRY SOLVER — инструмент, отменяющий carry-mismatch.

Exact Solver показал: GF(2) коррекция отменяет XOR-часть.
Остаток = carry ≈ 16 бит/раунд.
Carry ВЫЧИСЛИМ: Ψ(a,b) = (a+b) ⊕ (a⊕b).

Новый инструмент: iterative carry correction.
1. GF(2) solve → отменяет XOR-часть mismatch
2. Остаток = carry → вычисляем ТОЧНО
3. Корректируем carry → создаём новый XOR-mismatch
4. Снова GF(2) solve → отменяем новый XOR
5. Повторяем пока carry → 0

Это p-адический подъём (Hensel) но НЕ по mod 2^k,
а по РАУНДАМ: solve GF(2), correct carry, repeat.
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

def sha_all_states(W16):
    W=expand_real(W16); s=list(IV); states=[list(s)]
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        states.append(list(s))
    return states, W

def sha_compress_full(W16):
    st, _ = sha_all_states(W16)
    return [add32(IV[i], st[64][i]) for i in range(8)]


def iterative_carry_correct(W1, W2_init, target_r, max_iter=50):
    """
    Iterative carry correction for a SINGLE round.

    Goal: find W2 (close to W2_init) with δe[target_r+1] = 0.

    Method:
    1. Compute states for both messages
    2. Compute required W2[target_r] for δe=0 (exact inverse)
    3. Required ≠ actual → mismatch
    4. Mismatch = XOR part + carry part
    5. Correct W2[0..15] to match required W2[target_r]
       by adjusting the SOURCE words that feed into schedule at target_r

    The schedule: W[r] = σ1(W[r-2]) + W[r-7] + σ0(W[r-15]) + W[r-16]
    For r in gap: W[r] depends on specific W[0..15] words.
    Changing one source word changes W[r] by a COMPUTABLE amount.
    """
    W2 = list(W2_init)

    for iteration in range(max_iter):
        states1, Wexp1 = sha_all_states(W1)
        states2, Wexp2 = sha_all_states(W2)

        # Current δe
        de = states1[target_r+1][4] ^ states2[target_r+1][4]
        if de == 0:
            return W2, iteration, True

        # Required W2[target_r]
        s1r = states1[target_r]; s2r = states2[target_r]
        T1_1 = sub32(states1[target_r+1][4], s1r[3])  # e_new = d + T1
        T1_2_req = add32(T1_1, sub32(s1r[3], s2r[3]))
        W2_req = sub32(sub32(sub32(sub32(T1_2_req,
                  s2r[7]), Sig1(s2r[4])), Ch(s2r[4],s2r[5],s2r[6])), K[target_r])

        mismatch = sub32(W2_req, Wexp2[target_r])  # additive mismatch

        if mismatch == 0:
            return W2, iteration, True

        # Which W[0..15] word affects W[target_r] most directly?
        # Schedule dependencies: W[r] depends on W[r-2], W[r-7], W[r-15], W[r-16]
        # For r=target_r: deps = {target_r-2, target_r-7, target_r-15, target_r-16}
        # These may be in 0..15 (direct) or 16+ (indirect)

        # Strategy: try correcting each W[0..15] word, pick best
        best_word = -1
        best_de_hw = hw(de)

        for word in range(16):
            # How much does W2[word]+1 change Wexp2[target_r]?
            W2_test = list(W2)
            W2_test[word] = add32(W2[word], 1)
            Wexp_test = expand_real(W2_test)
            effect = sub32(Wexp_test[target_r], Wexp2[target_r])

            if effect == 0:
                continue

            # Ideal correction: W2[word] += mismatch / effect (mod 2^32)
            if effect & 1:  # Invertible
                inv = pow(effect, -1, 1 << 32)
                correction = (mismatch * inv) & MASK32

                W2_corrected = list(W2)
                W2_corrected[word] = add32(W2[word], correction)

                # Recompute states and check
                sc, Wc = sha_all_states(W2_corrected)
                de_new = states1[target_r+1][4] ^ sc[target_r+1][4]

                if hw(de_new) < best_de_hw:
                    best_de_hw = hw(de_new)
                    best_word = word
                    best_W2 = list(W2_corrected)

                    if de_new == 0:
                        break

        if best_word >= 0:
            W2 = best_W2
            if iteration < 5 or best_de_hw == 0:
                print(f"    Iter {iteration}: HW(δe)={best_de_hw} (corrected W[{best_word}])")
            if best_de_hw == 0:
                return W2, iteration, True
        else:
            # Random perturbation
            w = int.from_bytes(os.urandom(1),'big') % 16
            W2[w] ^= (1 << (int.from_bytes(os.urandom(1),'big') % 32))

    final_de = states1[target_r+1][4] ^ sha_all_states(W2)[0][target_r+1][4]
    return W2, max_iter, False


# ============================================================
# EXP 1: Single-round carry solver
# ============================================================
def exp1_single_round():
    print("="*70)
    print("EXP 1: ITERATIVE CARRY SOLVER — single round")
    print("="*70)

    for trial in range(5):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        for target_r in [20, 30, 40, 50]:
            W2_result, iters, success = iterative_carry_correct(W1, W2, target_r, max_iter=30)
            marker = "★★★" if success else ""
            if success:
                # Verify
                s1, _ = sha_all_states(W1)
                s2, _ = sha_all_states(W2_result)
                de = s1[target_r+1][4] ^ s2[target_r+1][4]
                n_diff = sum(1 for i in range(16) if W1[i] != W2_result[i])
                print(f"  Trial {trial}, r={target_r}: SOLVED in {iters} iters, "
                      f"δe=0 ✓, {n_diff} words changed {marker}")
            else:
                print(f"  Trial {trial}, r={target_r}: not converged in {iters} iters")


# ============================================================
# EXP 2: Chain — solve multiple rounds sequentially
# ============================================================
def exp2_chain():
    print("\n" + "="*70)
    print("EXP 2: CHAIN CARRY SOLVER — sequential rounds")
    print("  Solve r=20, then r=24, then r=28, ... keeping previous fixes")
    print("="*70)

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    W2 = list(W1); W2[0] ^= 0x80000000

    target_rounds = list(range(20, 52, 4))
    solved = 0

    for tr in target_rounds:
        W2, iters, success = iterative_carry_correct(W1, W2, tr, max_iter=20)
        if success:
            solved += 1
            # Check: did previous solutions survive?
            s1, _ = sha_all_states(W1)
            s2, _ = sha_all_states(W2)
            still_zero = sum(1 for prev_r in target_rounds[:solved]
                           if s1[prev_r+1][4] == s2[prev_r+1][4])
            print(f"  r={tr}: SOLVED ✓ (previous still zero: {still_zero}/{solved})")
        else:
            print(f"  r={tr}: FAILED")

    # Final hash diff
    H1 = sha_compress_full(W1)
    H2 = sha_compress_full(W2)
    hash_hw = sum(hw(H1[i]^H2[i]) for i in range(8))
    n_diff = sum(1 for i in range(16) if W1[i] != W2[i])
    print(f"\n  Chain result: {solved}/{len(target_rounds)} rounds solved")
    print(f"  Words changed: {n_diff}/16")
    print(f"  Hash diff: HW = {hash_hw}")


# ============================================================
# EXP 3: Full pipeline — Wang + carry solver through gap + hash
# ============================================================
def exp3_full_pipeline():
    print("\n" + "="*70)
    print("EXP 3: FULL PIPELINE — carry solver on every gap round")
    print("="*70)

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    W2 = list(W1); W2[0] ^= 0x80000000

    # Phase 1: try to zero δe at every round from 17 to 51
    zeros_achieved = 0
    for tr in range(17, 52):
        W2, iters, success = iterative_carry_correct(W1, W2, tr, max_iter=15)
        if success:
            zeros_achieved += 1

    print(f"  Rounds with δe=0: {zeros_achieved}/35")

    # Check full state
    s1, _ = sha_all_states(W1)
    s2, _ = sha_all_states(W2)

    actual_de_zeros = sum(1 for r in range(17, 52) if s1[r+1][4] == s2[r+1][4])
    actual_da_zeros = sum(1 for r in range(17, 52) if s1[r+1][0] == s2[r+1][0])

    print(f"  Verified δe=0 rounds: {actual_de_zeros}/35")
    print(f"  Verified δa=0 rounds: {actual_da_zeros}/35")

    # Hash
    H1 = sha_compress_full(W1)
    H2 = sha_compress_full(W2)
    hash_hw = sum(hw(H1[i]^H2[i]) for i in range(8))
    print(f"  Hash diff: HW = {hash_hw}")
    print(f"  Words changed: {sum(1 for i in range(16) if W1[i]!=W2[i])}/16")

    if hash_hw < 100:
        print(f"  ★ Below random!")


if __name__ == '__main__':
    exp1_single_round()
    exp2_chain()
    exp3_full_pipeline()
