"""
Round-by-round breakthrough solver.

Phase 1 (R=1-16): Wang chain, de[2..16]=0, FREE
Phase 2 (R=17):   Birthday on W[0], cost O(2^32)
Phase 3 (R=18+):  Each extra round = birthday on neutral bits

Goal: go as far as possible, measuring cost at each round.
"""

import random
import time

from qt_solver.sha256_traced import (
    MASK32, IV256, K256, sha256_compress, sha256_compress_traced,
    get_all_carry_chains, add32, big_sigma0, big_sigma1, ch, maj,
    sigma0, sigma1,
)


def full_sha_pair(W_n, W_f, num_rounds):
    """Compute SHA on both messages, return de at each round."""
    a_n, b_n, c_n, d_n = IV256[0], IV256[1], IV256[2], IV256[3]
    e_n, f_n, g_n, h_n = IV256[4], IV256[5], IV256[6], IV256[7]
    a_f, b_f, c_f, d_f = IV256[0], IV256[1], IV256[2], IV256[3]
    e_f, f_f, g_f, h_f = IV256[4], IV256[5], IV256[6], IV256[7]

    # Expand schedule for both
    wn = list(W_n)
    wf = list(W_f)
    for i in range(16, num_rounds):
        wn.append(add32(add32(add32(sigma0(wn[i-15]), wn[i-16]), wn[i-7]), sigma1(wn[i-2])))
        wf.append(add32(add32(add32(sigma0(wf[i-15]), wf[i-16]), wf[i-7]), sigma1(wf[i-2])))

    de_list = []
    state_n = []
    state_f = []

    for r in range(num_rounds):
        T1_n = add32(add32(add32(add32(h_n, big_sigma1(e_n)), ch(e_n, f_n, g_n)), K256[r]), wn[r])
        T2_n = add32(big_sigma0(a_n), maj(a_n, b_n, c_n))
        T1_f = add32(add32(add32(add32(h_f, big_sigma1(e_f)), ch(e_f, f_f, g_f)), K256[r]), wf[r])
        T2_f = add32(big_sigma0(a_f), maj(a_f, b_f, c_f))

        e_new_n = add32(d_n, T1_n)
        a_new_n = add32(T1_n, T2_n)
        e_new_f = add32(d_f, T1_f)
        a_new_f = add32(T1_f, T2_f)

        de = (e_new_f - e_new_n) & MASK32
        de_list.append(de)

        h_n, g_n, f_n, e_n = g_n, f_n, e_n, e_new_n
        d_n, c_n, b_n, a_n = c_n, b_n, a_n, a_new_n
        h_f, g_f, f_f, e_f = g_f, f_f, e_f, e_new_f
        d_f, c_f, b_f, a_f = c_f, b_f, a_f, a_new_f

    # Final hash
    H_n = [add32(x, iv) for x, iv in zip([a_n,b_n,c_n,d_n,e_n,f_n,g_n,h_n], IV256)]
    H_f = [add32(x, iv) for x, iv in zip([a_f,b_f,c_f,d_f,e_f,f_f,g_f,h_f], IV256)]

    return de_list, H_n, H_f


def wang_adapt(W_n, delta_W0=0x8000):
    """Wang chain: adapt W_f[0..15] to force de[2..16]=0."""
    W_f = list(W_n)
    W_f[0] = add32(W_n[0], delta_W0)

    a_n, b_n, c_n, d_n = IV256[0], IV256[1], IV256[2], IV256[3]
    e_n, f_n, g_n, h_n = IV256[4], IV256[5], IV256[6], IV256[7]
    a_f, b_f, c_f, d_f = IV256[0], IV256[1], IV256[2], IV256[3]
    e_f, f_f, g_f, h_f = IV256[4], IV256[5], IV256[6], IV256[7]

    for r in range(16):
        T1_n = add32(add32(add32(add32(h_n, big_sigma1(e_n)), ch(e_n, f_n, g_n)), K256[r]), W_n[r])
        T2_n = add32(big_sigma0(a_n), maj(a_n, b_n, c_n))
        T1_f = add32(add32(add32(add32(h_f, big_sigma1(e_f)), ch(e_f, f_f, g_f)), K256[r]), W_f[r])
        T2_f = add32(big_sigma0(a_f), maj(a_f, b_f, c_f))

        e_new_n = add32(d_n, T1_n)
        a_new_n = add32(T1_n, T2_n)
        e_new_f = add32(d_f, T1_f)
        a_new_f = add32(T1_f, T2_f)

        h_n, g_n, f_n, e_n = g_n, f_n, e_n, e_new_n
        d_n, c_n, b_n, a_n = c_n, b_n, a_n, a_new_n
        h_f, g_f, f_f, e_f = g_f, f_f, e_f, e_new_f
        d_f, c_f, b_f, a_f = c_f, b_f, a_f, a_new_f

        # Adapt next round
        if r < 15:
            delta_d = (d_f - d_n) & MASK32
            delta_h = (h_f - h_n) & MASK32
            delta_S1 = (big_sigma1(e_f) - big_sigma1(e_n)) & MASK32
            delta_Ch = (ch(e_f, f_f, g_f) - ch(e_n, f_n, g_n)) & MASK32
            dw = (0 - delta_d - delta_h - delta_S1 - delta_Ch) & MASK32
            W_f[r + 1] = add32(W_n[r + 1], dw)

    return W_f


def breakthrough(target_rounds=64, max_birthday=2**25, seed=42, verbose=True):
    """
    Push through round by round.

    Phase 1: Wang chain (rounds 1-16, free)
    Phase 2: Birthday for round 17 (search over W[0])
    Phase 3: Birthday for rounds 18+ (search over neutral bits)
    """
    rng = random.Random(seed)

    if verbose:
        print(f"{'═'*60}")
        print(f"BREAKTHROUGH: pushing to R={target_rounds}")
        print(f"{'═'*60}")

    # Base message
    W_base = [rng.randint(0, MASK32) for _ in range(16)]

    # Phase 1: Wang chain (FREE, rounds 1-16)
    t0 = time.time()
    W_f = wang_adapt(W_base)
    de_list, H_n, H_f = full_sha_pair(W_base, W_f, min(target_rounds, 16))

    zeros = sum(1 for de in de_list[1:] if de == 0)
    if verbose:
        print(f"\nPhase 1 (Wang chain): de[2..{min(16,target_rounds)}]=0")
        print(f"  {zeros} consecutive zeros (rounds 2-{zeros+1})")
        print(f"  Cost: O(1), Time: {time.time()-t0:.3f}s")

    if target_rounds <= 16:
        # Check hash at target rounds
        Hn = sha256_compress(W_base, target_rounds)
        Hf = sha256_compress(W_f, target_rounds)
        match = Hn == Hf
        if verbose:
            print(f"  Hash match at R={target_rounds}: {match}")
            if match and W_base != W_f:
                print(f"  *** COLLISION FOR R={target_rounds}! ***")
        return {'rounds_solved': zeros + 1, 'cost': 0,
                'W_n': W_base, 'W_f': W_f, 'hash_match': match}

    # Phase 2: Birthday for round 17
    # de[17] depends on W[0] (through Da[13] and DW[16])
    # Search: vary W[0], check de[17]
    if verbose:
        print(f"\nPhase 2 (Birthday R=17): searching over W[0]...")

    t1 = time.time()
    best_hw = 32
    attempts = 0

    for attempt in range(max_birthday):
        W_try = list(W_base)
        W_try[0] = rng.randint(0, MASK32)
        W_f_try = wang_adapt(W_try)

        de_list, _, _ = full_sha_pair(W_try, W_f_try, 17)

        # de[1] != 0 (always), de[2..16] = 0 (Wang), check de[17]
        de17 = de_list[16]  # 0-indexed: round 17 = index 16
        hw = bin(de17).count('1')

        if hw < best_hw:
            best_hw = hw
            best_W = list(W_try)
            best_Wf = list(W_f_try)

        if de17 == 0:
            elapsed = time.time() - t1
            if verbose:
                print(f"  FOUND de[17]=0 at attempt {attempt}! ({elapsed:.1f}s)")
                print(f"  W[0] = 0x{W_try[0]:08x}")

            # Continue to round 18+
            if target_rounds <= 17:
                Hn = sha256_compress(W_try, 17)
                Hf = sha256_compress(W_f_try, 17)
                if verbose:
                    print(f"  Hash match R=17: {Hn == Hf}")
                return {'rounds_solved': 17, 'cost': attempt,
                        'W_n': W_try, 'W_f': W_f_try}

            # Phase 3: push further
            return _push_further(W_try, W_f_try, 17, target_rounds,
                                  max_birthday, rng, verbose)

        attempts = attempt

        if verbose and attempt > 0 and attempt % 5000000 == 0:
            elapsed = time.time() - t1
            print(f"  {attempt} attempts, best HW(de17)={best_hw}, "
                  f"{attempt/elapsed:.0f}/s")

    if verbose:
        print(f"  R=17 not found in {max_birthday} attempts. "
              f"Best HW={best_hw}")
    return {'rounds_solved': 16, 'cost': max_birthday, 'best_hw17': best_hw}


def _push_further(W_n, W_f, solved_rounds, target_rounds, max_per_round,
                   rng, verbose):
    """
    Push past solved_rounds using neutral bit search.
    For each round R: vary neutral message words, check de[R+1].

    Neutral words for de[R+1]: any W[j] where j >= R-1 and j not in
    the schedule dependency of W64[R].
    """
    current_Wn = list(W_n)
    current_Wf = list(W_f)

    for R in range(solved_rounds + 1, min(target_rounds + 1, 65)):
        if verbose:
            print(f"\nPhase 3: pushing to R={R}...")

        # Which message words are "neutral" for this round?
        # W64[R] depends on W[R-16], sigma0(W[R-15]), W[R-7], sigma1(W[R-2])
        # For R=18: depends on W[2], W[3], W[11], W[16]
        # W[16] itself depends on W[0,1,9,14]
        # Neutral: words NOT in the dependency chain
        # Simplification: vary W[0] (affects Da and DW through schedule)

        t0 = time.time()
        found = False

        for attempt in range(max_per_round):
            W_try = list(current_Wn)
            # Vary W[0] (changes everything through Wang chain)
            W_try[0] = rng.randint(0, MASK32)
            W_f_try = wang_adapt(W_try)

            de_list, _, _ = full_sha_pair(W_try, W_f_try, R)

            # Check all de[2..R]
            all_zero = all(de == 0 for de in de_list[1:])

            if all_zero:
                elapsed = time.time() - t0
                if verbose:
                    print(f"  R={R}: FOUND at attempt {attempt}! ({elapsed:.1f}s)")
                current_Wn = W_try
                current_Wf = W_f_try
                found = True
                break

            if attempt % 5000000 == 0 and attempt > 0 and verbose:
                print(f"  R={R}: {attempt} attempts, {attempt/(time.time()-t0):.0f}/s")

        if not found:
            if verbose:
                print(f"  R={R}: not found in {max_per_round} attempts")
            return {'rounds_solved': R - 1, 'cost': max_per_round,
                    'W_n': current_Wn, 'W_f': current_Wf}

    Hn = sha256_compress(current_Wn, target_rounds)
    Hf = sha256_compress(current_Wf, target_rounds)
    if verbose:
        print(f"\nFinal: R={target_rounds}, hash match: {Hn == Hf}")
    return {'rounds_solved': target_rounds,
            'W_n': current_Wn, 'W_f': current_Wf}


if __name__ == '__main__':
    # Phase 1: free rounds
    breakthrough(target_rounds=16, verbose=True)
    print()
    # Phase 2: birthday for R=17
    breakthrough(target_rounds=17, max_birthday=2**25, verbose=True)
