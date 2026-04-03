"""
OMEGA Replay: compute ALL intermediate variables from a SHA-256 trace.

For each round: replay T1 chain (4 adds), T2 (1 add), e_new (1 add), a_new (1 add).
For each addition: compute result bits and carry chain bits.
Map everything to OMEGA variable IDs.
"""

from qt_solver.sha256_traced import (
    MASK32, IV256, K256, get_bit, add32_traced,
    big_sigma0, big_sigma1, ch, maj, sigma0, sigma1,
    sha256_compress, sha256_compress_traced,
)


def omega_replay(sys, msg, R):
    """
    Replay SHA-256 in OMEGA coordinates.
    Returns complete assignment dict: var_id → bit_value.
    """
    assignment = {}
    trace = sha256_compress_traced(msg, R)

    # 1. Message words
    for w in range(16):
        for b in range(32):
            assignment[sys.W[w][b]] = get_bit(msg[w], b)

    # 2. Schedule words W[16..R-1]
    for w in range(16, R):
        for b in range(32):
            assignment[sys.W[w][b]] = get_bit(trace.schedule[w], b)

    # 3. State registers
    for r in range(1, R + 1):
        for reg in range(8):
            for b in range(32):
                assignment[sys.state[(r, reg)][b]] = get_bit(trace.states[r][reg], b)

    # 4. ALL intermediates: replay round by round
    # We need to walk through sys.V.names to find variable IDs
    # and compute their values from the SHA trace.

    # Build name→var_id lookup
    name_to_var = {name: var for var, name in sys.V.names.items()}

    # For each round: replay the computation
    for r in range(R):
        state = trace.states[r]
        a, b_reg, c_reg, d = state[0], state[1], state[2], state[3]
        e, f_reg, g, h = state[4], state[5], state[6], state[7]
        w_val = trace.schedule[r]

        # S1 = Sigma1(e)
        S1_val = big_sigma1(e)
        _assign_word(assignment, name_to_var, f"S1_{r}", S1_val)

        # ch = Ch(e, f, g)
        ch_val = ch(e, f_reg, g)
        _assign_word(assignment, name_to_var, f"ch_{r}", ch_val)

        # T1 chain: h + S1 + ch + K + W
        # add1: h + S1
        t1a, c1, _ = add32_traced(h, S1_val)
        _assign_word(assignment, name_to_var, f"t1a_{r}_r", t1a)
        _assign_word(assignment, name_to_var, f"t1a_{r}_c", c1)

        # add2: t1a + ch
        t1b, c2, _ = add32_traced(t1a, ch_val)
        _assign_word(assignment, name_to_var, f"t1b_{r}_r", t1b)
        _assign_word(assignment, name_to_var, f"t1b_{r}_c", c2)

        # add3: t1b + K
        t1c, c3, _ = add32_traced(t1b, K256[r])
        _assign_word(assignment, name_to_var, f"t1c_{r}_r", t1c)
        _assign_word(assignment, name_to_var, f"t1c_{r}_c", c3)

        # add4: t1c + W
        T1, c4, _ = add32_traced(t1c, w_val)
        _assign_word(assignment, name_to_var, f"t1d_{r}_r", T1)
        _assign_word(assignment, name_to_var, f"t1d_{r}_c", c4)

        # S0 = Sigma0(a)
        S0_val = big_sigma0(a)
        _assign_word(assignment, name_to_var, f"S0_{r}", S0_val)

        # mj = Maj(a, b, c)
        mj_val = maj(a, b_reg, c_reg)
        _assign_word(assignment, name_to_var, f"mj_{r}", mj_val)

        # T2 = S0 + mj
        T2, c5, _ = add32_traced(S0_val, mj_val)
        _assign_word(assignment, name_to_var, f"t2_{r}_r", T2)
        _assign_word(assignment, name_to_var, f"t2_{r}_c", c5)

        # e_new = d + T1
        e_new, c6, _ = add32_traced(d, T1)
        _assign_word(assignment, name_to_var, f"enew_{r}_r", e_new)
        _assign_word(assignment, name_to_var, f"enew_{r}_c", c6)

        # a_new = T1 + T2
        a_new, c7, _ = add32_traced(T1, T2)
        _assign_word(assignment, name_to_var, f"anew_{r}_r", a_new)
        _assign_word(assignment, name_to_var, f"anew_{r}_c", c7)

    # 5. Schedule intermediates (for R > 16)
    for w in range(16, R):
        s1_val = sigma1(trace.schedule[w - 2])
        s0_val = sigma0(trace.schedule[w - 15])

        _assign_word(assignment, name_to_var, f"ss1_{w}", s1_val)
        _assign_word(assignment, name_to_var, f"ss0_{w}", s0_val)

        # 3 schedule additions
        t1, c1, _ = add32_traced(trace.schedule[w - 16], s0_val)
        _assign_word(assignment, name_to_var, f"sw1_{w}_r", t1)
        _assign_word(assignment, name_to_var, f"sw1_{w}_c", c1)

        t2, c2, _ = add32_traced(t1, trace.schedule[w - 7])
        _assign_word(assignment, name_to_var, f"sw2_{w}_r", t2)
        _assign_word(assignment, name_to_var, f"sw2_{w}_c", c2)

        t3, c3, _ = add32_traced(t2, s1_val)
        _assign_word(assignment, name_to_var, f"sw3_{w}_r", t3)
        _assign_word(assignment, name_to_var, f"sw3_{w}_c", c3)

    # 6. Final addition carries
    final_state = trace.states[R]
    for reg in range(8):
        _, fc, _ = add32_traced(final_state[reg], IV256[reg])
        _assign_word(assignment, name_to_var, f"cf{reg}", fc)

    return assignment


def _assign_word(assignment, name_to_var, prefix, value):
    """Assign 32-bit word to variables with given prefix."""
    for b in range(32):
        name = f"{prefix}[{b}]"
        if name in name_to_var:
            assignment[name_to_var[name]] = get_bit(value, b)
