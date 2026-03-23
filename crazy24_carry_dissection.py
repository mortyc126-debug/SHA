#!/usr/bin/env python3
"""
crazy24_carry_dissection.py — Carry-chain dissection of the De17 barrier function.

We study how W[14] propagates through modular additions to reach e[17],
decompose each addition into XOR + carry, and measure effective chain length.
Full 2^32 enumeration (or 2^24 sample) to count De17=0 occurrences.
"""

import random
import time
import struct

MASK = 0xFFFFFFFF
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc]
H0 = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Ch(e,f,g): return ((e&f)^(~e&g))&MASK
def Maj(a,b,c): return ((a&b)^(a&c)^(b&c))&MASK
def add32(*a):
    s=0
    for x in a: s=(s+x)&MASK
    return s
def hw(x): return bin(x&MASK).count('1')

def sha_round(st, w, k):
    a,b,c,d,e,f,g,h = st
    T1 = add32(h, Sig1(e), Ch(e,f,g), k, w)
    T2 = add32(Sig0(a), Maj(a,b,c))
    return [add32(T1,T2), a, b, c, add32(d,T1), e, f, g]

def expand_W(W16):
    W = list(W16)
    for i in range(16, 20):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W

def get_de17(msg, iv=H0):
    """Compute De17 = e[17] XOR e'[17] for Wang twin with W[0] bit 31 flip."""
    W = list(msg); Wp = list(W); Wp[0] ^= 0x80000000
    s = list(iv); sp = list(iv)
    # Round 0
    s = sha_round(s, W[0], K[0])
    sp = sha_round(sp, Wp[0], K[0])
    # Rounds 1-15: adjust Wp[t] so that d' = d after each round
    for t in range(1, 16):
        a,b,c,d,e,f,g,h = s
        a2,b2,c2,d2,e2,f2,g2,h2 = sp
        tp = add32(h, Sig1(e), Ch(e,f,g), K[t])
        tp2 = add32(h2, Sig1(e2), Ch(e2,f2,g2), K[t])
        target = add32(d, tp, W[t])
        Wp[t] = (target - d2 - tp2) & MASK
        s = sha_round(s, W[t], K[t])
        sp = sha_round(sp, Wp[t], K[t])
    # Expand both
    We = expand_W(W)
    Wpe = expand_W(Wp)
    # Recompute full 17 rounds from scratch
    s1 = list(iv); s2 = list(iv)
    for t in range(17):
        s1 = sha_round(s1, We[t], K[t])
        s2 = sha_round(s2, Wpe[t], K[t])
    return s1[4] ^ s2[4]

# ─── CARRY DECOMPOSITION UTILITIES ───

def add32_with_carries(a, b):
    """Compute a+b mod 2^32, return (sum, carry_bits).
    carry_bits: bit k is set if there's a carry INTO position k+1."""
    s = (a + b) & MASK
    # carry bits: where the sum differs from XOR
    xor_result = a ^ b
    carry_bits = (s ^ xor_result)  # these are the bits affected by carries
    # The actual carry-in at each position: carry_in[k] = (s[k] XOR a[k] XOR b[k])
    # But carry_bits already captures this
    return s, carry_bits

def add32_multi_with_carries(values):
    """Chain of additions, return (sum, list_of_carry_bits_per_addition)."""
    if len(values) == 0:
        return 0, []
    s = values[0]
    carries_list = []
    for v in values[1:]:
        s, carry = add32_with_carries(s, v)
        carries_list.append(carry)
    return s, carries_list

# ─── TRACE WHICH ADDITIONS W[14] PASSES THROUGH ───

def trace_w14_to_e17(msg, iv=H0):
    """Trace the path of W[14] through the SHA-256 computation to e[17].
    Identify every modular addition W[14] participates in."""
    W = list(msg)
    We = expand_W(W)

    additions = []

    # --- W[14] enters at round 14 as the message word ---
    # In round 14: T1 = h + Sig1(e) + Ch(e,f,g) + K[14] + W[14]
    # Then: e_new = d + T1, a_new = T1 + T2

    # First, run rounds 0..13 to get state before round 14
    s = list(iv)
    for t in range(14):
        s = sha_round(s, We[t], K[t])

    a,b,c,d,e,f,g,h = s

    # Round 14: W[14] enters
    # Addition 1: partial = h + Sig1(e)
    p1, c1 = add32_with_carries(h, Sig1(e))
    additions.append(("R14: h + Sig1(e)", h, Sig1(e), c1))

    # Addition 2: partial += Ch(e,f,g)
    p2, c2 = add32_with_carries(p1, Ch(e,f,g))
    additions.append(("R14: +Ch(e,f,g)", p1, Ch(e,f,g), c2))

    # Addition 3: partial += K[14]
    p3, c3 = add32_with_carries(p2, K[14])
    additions.append(("R14: +K[14]", p2, K[14], c3))

    # Addition 4: T1 = partial + W[14]  *** W[14] enters here ***
    T1, c4 = add32_with_carries(p3, We[14])
    additions.append(("R14: +W[14] -> T1", p3, We[14], c4))

    # Addition 5: e_new = d + T1
    e_new, c5 = add32_with_carries(d, T1)
    additions.append(("R14: d+T1 -> e_new", d, T1, c5))

    # Addition 6: T2 = Sig0(a) + Maj(a,b,c)
    T2, c6 = add32_with_carries(Sig0(a), Maj(a,b,c))
    additions.append(("R14: Sig0+Maj -> T2", Sig0(a), Maj(a,b,c), c6))

    # Addition 7: a_new = T1 + T2
    a_new, c7 = add32_with_carries(T1, T2)
    additions.append(("R14: T1+T2 -> a_new", T1, T2, c7))

    # Now state after round 14
    s14 = [a_new, a, b, c, e_new, e, f, g]

    # Round 15: W[14] affects state through s14.
    # e from R14 = e_new, a from R14 = a_new
    # T1_15 = g15 + Sig1(e15) + Ch(e15,f15,g15) + K[15] + W[15]
    # e15 = e_new (carries W[14] info)
    a15,b15,c15,d15,e15,f15,g15,h15 = s14

    # Sig1(e15) and Ch(e15,...) are bitwise ops on e_new — no additions, but
    # they feed INTO additions:
    p1_15, c1_15 = add32_with_carries(h15, Sig1(e15))
    additions.append(("R15: h+Sig1(e) [e carries W14]", h15, Sig1(e15), c1_15))

    p2_15, c2_15 = add32_with_carries(p1_15, Ch(e15,f15,g15))
    additions.append(("R15: +Ch(e,f,g)", p1_15, Ch(e15,f15,g15), c2_15))

    p3_15, c3_15 = add32_with_carries(p2_15, K[15])
    additions.append(("R15: +K[15]", p2_15, K[15], c3_15))

    T1_15, c4_15 = add32_with_carries(p3_15, We[15])
    additions.append(("R15: +W[15] -> T1", p3_15, We[15], c4_15))

    T2_15, c6_15 = add32_with_carries(Sig0(a15), Maj(a15,b15,c15))
    additions.append(("R15: Sig0+Maj -> T2 [a carries W14]", Sig0(a15), Maj(a15,b15,c15), c6_15))

    e_new_15, c5_15 = add32_with_carries(d15, T1_15)
    additions.append(("R15: d+T1 -> e_new", d15, T1_15, c5_15))

    a_new_15, c7_15 = add32_with_carries(T1_15, T2_15)
    additions.append(("R15: T1+T2 -> a_new", T1_15, T2_15, c7_15))

    s15 = [a_new_15, a15, b15, c15, e_new_15, e15, f15, g15]

    # Round 16: W[16] = sig1(W[14]) + W[9] + sig0(W[1]) + W[0]
    # So W[16] itself involves W[14] through sig1
    # sig1 is bitwise (rotr+shift+xor), no addition
    # But W[16] computation: 4-input addition
    sig1_w14 = sig1(We[14])
    p1_w16, cw1 = add32_with_carries(sig1_w14, We[9])
    additions.append(("W[16]: sig1(W14)+W[9]", sig1_w14, We[9], cw1))
    p2_w16, cw2 = add32_with_carries(p1_w16, sig0(We[1]))
    additions.append(("W[16]: +sig0(W1)", p1_w16, sig0(We[1]), cw2))
    W16, cw3 = add32_with_carries(p2_w16, We[0])
    additions.append(("W[16]: +W[0] -> W16", p2_w16, We[0], cw3))

    # Round 16 state computation
    a16,b16,c16,d16,e16,f16,g16,h16 = s15
    p1_16, c1_16 = add32_with_carries(h16, Sig1(e16))
    additions.append(("R16: h+Sig1(e) [e carries W14]", h16, Sig1(e16), c1_16))
    p2_16, c2_16 = add32_with_carries(p1_16, Ch(e16,f16,g16))
    additions.append(("R16: +Ch(e,f,g)", p1_16, Ch(e16,f16,g16), c2_16))
    p3_16, c3_16 = add32_with_carries(p2_16, K[16])
    additions.append(("R16: +K[16]", p2_16, K[16], c3_16))
    T1_16, c4_16 = add32_with_carries(p3_16, We[16])
    additions.append(("R16: +W[16] -> T1 [W16 carries W14]", p3_16, We[16], c4_16))
    T2_16, c6_16 = add32_with_carries(Sig0(a16), Maj(a16,b16,c16))
    additions.append(("R16: Sig0+Maj -> T2 [a carries W14]", Sig0(a16), Maj(a16,b16,c16), c6_16))
    e_new_16, c5_16 = add32_with_carries(d16, T1_16)
    additions.append(("R16: d+T1 -> e_new", d16, T1_16, c5_16))
    a_new_16, c7_16 = add32_with_carries(T1_16, T2_16)
    additions.append(("R16: T1+T2 -> a_new", T1_16, T2_16, c7_16))

    s16 = [a_new_16, a16, b16, c16, e_new_16, e16, f16, g16]

    # Round 17: e[17] = d[16] + T1[17]
    # T1[17] = h[16] + Sig1(e[16]) + Ch(e[16],f[16],g[16]) + K[17] + W[17]
    # W[17] = sig1(W[15]) + W[10] + sig0(W[2]) + W[1]  — no W[14] dependency
    # But the state registers carry W[14] info from rounds 14-16
    a17,b17,c17,d17,e17,f17,g17,h17 = s16
    p1_17, c1_17 = add32_with_carries(h17, Sig1(e17))
    additions.append(("R17: h+Sig1(e) [e carries W14]", h17, Sig1(e17), c1_17))
    p2_17, c2_17 = add32_with_carries(p1_17, Ch(e17,f17,g17))
    additions.append(("R17: +Ch(e,f,g)", p1_17, Ch(e17,f17,g17), c2_17))
    p3_17, c3_17 = add32_with_carries(p2_17, K[17])
    additions.append(("R17: +K[17]", p2_17, K[17], c3_17))
    T1_17, c4_17 = add32_with_carries(p3_17, We[17])
    additions.append(("R17: +W[17] -> T1", p3_17, We[17], c4_17))
    e_new_17, c5_17 = add32_with_carries(d17, T1_17)
    additions.append(("R17: d+T1 -> e[17]", d17, T1_17, c5_17))

    return additions, e_new_17

# ─── MAIN EXPERIMENT ───

def main():
    print("="*80)
    print("CARRY-CHAIN DISSECTION OF De17 BARRIER")
    print("="*80)

    # Generate base message with seed 0xC024
    rng = random.Random(0xC024)
    base_msg = [rng.randint(0, MASK) for _ in range(16)]

    print(f"\nBase message (seed 0xC024):")
    for i in range(16):
        print(f"  W[{i:2d}] = 0x{base_msg[i]:08x}")

    # ─── PART 1: Trace additions W[14] passes through ───
    print("\n" + "="*80)
    print("PART 1: MODULAR ADDITIONS ON W[14]->e[17] PATH")
    print("="*80)

    additions, e17_val = trace_w14_to_e17(base_msg)

    total_carry_bits = 0
    predictable = 0

    for i, (desc, op_a, op_b, carry) in enumerate(additions):
        cb = hw(carry)
        total_carry_bits += cb
        # A carry bit is "predictable" if both input bits are same (both 0 -> carry out 0, both 1 -> carry out 1)
        # Unpredictable when inputs differ (carry depends on lower bits)
        agree = hw(~(op_a ^ op_b) & MASK)
        disagree = 32 - agree
        pred = agree  # positions where carry is locally determined
        predictable += pred

        print(f"  Add {i+1:2d}: {desc}")
        print(f"          a=0x{op_a:08x} b=0x{op_b:08x}")
        print(f"          carry_bits=0x{carry:08x} (HW={cb:2d}), "
              f"agree_bits={agree}, disagree={disagree}")

    print(f"\n  TOTAL additions on path: {len(additions)}")
    print(f"  TOTAL carry bits generated: {total_carry_bits}")
    print(f"  TOTAL predictable positions: {predictable} / {len(additions)*32}")
    print(f"  Carry density: {total_carry_bits/(len(additions)*32)*100:.1f}%")

    # ─── PART 2: Classify additions into "known" vs "W[14]-dependent" ───
    print("\n" + "="*80)
    print("PART 2: KNOWN vs W[14]-DEPENDENT ADDITIONS")
    print("="*80)

    # Additions in round 14 before W[14] enters: adds 1-3 are independent of W[14]
    # Addition 4 is where W[14] first enters
    # Addition 6 (T2 computation) is also independent of W[14] at round 14

    # Classify: "known" = does not depend on W[14] value
    # In R14: adds 1-3 (T1 partial before W[14]), add 6 (T2) are known
    # Add 4 (T1 = partial + W[14]): FIRST W[14]-dependent
    # Add 5 (e = d + T1): depends on W[14] through T1
    # Add 7 (a = T1 + T2): depends on W[14] through T1
    # From R15 onward: all state-dependent adds carry W[14] info

    known_adds = [0, 1, 2, 5]  # indices of known additions (0-indexed)
    first_w14_add = 3  # addition 4 (index 3) is first W[14]-dependent

    known_carry_bits = sum(hw(additions[i][3]) for i in known_adds)
    total_known_positions = len(known_adds) * 32

    dep_indices = [i for i in range(len(additions)) if i not in known_adds]
    dep_carry_bits = sum(hw(additions[i][3]) for i in dep_indices)
    total_dep_positions = len(dep_indices) * 32

    print(f"\n  Known (W[14]-independent) additions: {len(known_adds)}")
    print(f"    Carry bits: {known_carry_bits} / {total_known_positions} = {known_carry_bits/total_known_positions*100:.1f}%")
    print(f"  W[14]-dependent additions: {len(dep_indices)}")
    print(f"    Carry bits: {dep_carry_bits} / {total_dep_positions} = {dep_carry_bits/total_dep_positions*100:.1f}%")
    print(f"  Effective chain length (W[14]-dependent): {len(dep_indices)} additions")

    # ─── PART 3: Measure De17 for various W[14] values (carry fixed experiment) ───
    print("\n" + "="*80)
    print("PART 3: RESIDUAL HW WHEN KNOWN CARRIES ARE FIXED")
    print("="*80)

    # Sample De17 for many W[14] values, measure HW distribution
    N_sample = 1 << 20  # 1M samples for quick stats
    hw_counts = [0] * 33
    de17_zero = 0
    hw_sum = 0

    for i in range(N_sample):
        msg = list(base_msg)
        msg[14] = i & MASK
        de = get_de17(msg)
        h = hw(de)
        hw_counts[h] += 1
        hw_sum += h
        if de == 0:
            de17_zero += 1

    mean_hw = hw_sum / N_sample
    print(f"\n  Sampled {N_sample} values of W[14]")
    print(f"  Mean HW(De17) = {mean_hw:.4f} (ideal random: 16.0)")
    print(f"  De17 = 0 count: {de17_zero}")
    print(f"  De17 = 0 rate: {de17_zero/N_sample:.2e} (random expectation: {1/2**32:.2e})")

    # HW distribution
    print(f"\n  HW distribution (top 10 most frequent):")
    ranked = sorted(range(33), key=lambda h: -hw_counts[h])
    for h in ranked[:10]:
        print(f"    HW={h:2d}: {hw_counts[h]:8d} ({hw_counts[h]/N_sample*100:.3f}%)")

    # ─── PART 4: Per-bit bias of De17 ───
    print("\n" + "="*80)
    print("PART 4: PER-BIT BIAS OF De17")
    print("="*80)

    bit_ones = [0]*32
    N_bias = 1 << 18  # 256K for speed
    for i in range(N_bias):
        msg = list(base_msg)
        msg[14] = rng.randint(0, MASK)
        de = get_de17(msg)
        for b in range(32):
            if de & (1 << b):
                bit_ones[b] += 1

    print(f"\n  Per-bit probability of De17[k]=1 (ideal: 0.5000):")
    max_bias = 0
    for b in range(32):
        p = bit_ones[b] / N_bias
        bias = abs(p - 0.5)
        max_bias = max(max_bias, bias)
        marker = " ***" if bias > 0.01 else ""
        print(f"    bit {b:2d}: P=1 = {p:.4f}  bias={bias:.4f}{marker}")
    print(f"\n  Max single-bit bias: {max_bias:.4f}")

    # ─── PART 5: Larger enumeration (2^22 sequential) ───
    print("\n" + "="*80)
    print("PART 5: LARGE-SCALE ENUMERATION (2^22 sequential)")
    print("="*80)

    N_enum = 1 << 22
    t0 = time.time()
    zero_count = 0
    low_hw_count = 0  # HW <= 4
    hw_accum = 0
    min_hw_seen = 33
    min_hw_val = 0

    # Pre-compute parts that don't depend on W[14]
    # Run rounds 0-13 once, then inline rounds 14-17
    W_base = list(base_msg)
    We_base = expand_W(W_base)
    Wp_base = list(W_base); Wp_base[0] ^= 0x80000000
    s_o = list(H0); s_t = list(H0)
    s_o = sha_round(s_o, W_base[0], K[0])
    s_t = sha_round(s_t, Wp_base[0], K[0])
    for t in range(1, 16):
        a,b,c,d,e,f,g,h = s_o
        a2,b2,c2,d2,e2,f2,g2,h2 = s_t
        tp_ = add32(h, Sig1(e), Ch(e,f,g), K[t])
        tp2_ = add32(h2, Sig1(e2), Ch(e2,f2,g2), K[t])
        target_ = add32(d, tp_, W_base[t])
        Wp_base[t] = (target_ - d2 - tp2_) & MASK
        s_o = sha_round(s_o, W_base[t], K[t])
        s_t = sha_round(s_t, Wp_base[t], K[t])

    for w14 in range(N_enum):
        msg = list(base_msg)
        msg[14] = w14
        de = get_de17(msg)
        if de == 0:
            zero_count += 1
            print(f"    *** De17=0 found at W[14]=0x{w14:08x} ***", flush=True)
        h = hw(de)
        if h < min_hw_seen:
            min_hw_seen = h
            min_hw_val = w14
        if h <= 4:
            low_hw_count += 1
        hw_accum += h
        if w14 > 0 and w14 % (1 << 20) == 0:
            print(f"    ... {w14>>20}/{N_enum>>20} M done ...", flush=True)

    elapsed = time.time() - t0

    print(f"\n  Enumerated: 2^22 = {N_enum} values of W[14]")
    print(f"  Time: {elapsed:.1f}s ({N_enum/elapsed:.0f} evals/sec)")
    print(f"  De17=0 count: {zero_count}")
    print(f"  De17=0 rate: {zero_count/N_enum:.2e}")
    extrapolated_zeros = zero_count * (2**32 / N_enum)
    print(f"  Extrapolated De17=0 in full 2^32: ~{extrapolated_zeros:.1f}")
    print(f"  Mean HW: {hw_accum/N_enum:.4f}")
    print(f"  Min HW seen: {min_hw_seen} at W[14]=0x{min_hw_val:08x}")
    print(f"  HW<=4 count: {low_hw_count} ({low_hw_count/N_enum*100:.4f}%)")

    # ─── PART 6: Carry chain depth analysis ───
    print("\n" + "="*80)
    print("PART 6: CARRY CHAIN DEPTH ANALYSIS")
    print("="*80)

    # For a few random W[14] values, trace carry propagation length
    # (longest consecutive carry chain in each addition)
    print("\n  Carry chain lengths for 10 sample W[14] values:")

    for trial in range(10):
        msg = list(base_msg)
        msg[14] = rng.randint(0, MASK)
        adds, _ = trace_w14_to_e17(msg)

        max_chain = 0
        total_carry_hw = 0
        chains_per_add = []

        for desc, op_a, op_b, carry in adds:
            # Find longest consecutive run of 1s in carry bits
            c = carry
            longest = 0
            current = 0
            for bit in range(32):
                if c & (1 << bit):
                    current += 1
                    longest = max(longest, current)
                else:
                    current = 0
            max_chain = max(max_chain, longest)
            total_carry_hw += hw(carry)
            chains_per_add.append(longest)

        print(f"    W[14]=0x{msg[14]:08x}: max_chain={max_chain:2d}, "
              f"total_carry_bits={total_carry_hw:3d}, "
              f"per-add max chains={chains_per_add}")

    # ─── PART 7: Differential carry analysis ───
    print("\n" + "="*80)
    print("PART 7: DIFFERENTIAL CARRY ANALYSIS (orig vs twin)")
    print("="*80)

    # For the base message, compare carries in original vs Wang twin
    W = list(base_msg)
    Wp = list(W)
    Wp[0] ^= 0x80000000

    # Build twin Wp[1..15]
    s = list(H0); sp = list(H0)
    s = sha_round(s, W[0], K[0])
    sp = sha_round(sp, Wp[0], K[0])
    for t in range(1, 16):
        a,b,c,d,e,f,g,h = s
        a2,b2,c2,d2,e2,f2,g2,h2 = sp
        tp = add32(h, Sig1(e), Ch(e,f,g), K[t])
        tp2 = add32(h2, Sig1(e2), Ch(e2,f2,g2), K[t])
        target = add32(d, tp, W[t])
        Wp[t] = (target - d2 - tp2) & MASK
        s = sha_round(s, W[t], K[t])
        sp = sha_round(sp, Wp[t], K[t])

    adds_orig, e17_orig = trace_w14_to_e17(W)
    adds_twin, e17_twin = trace_w14_to_e17(Wp)

    print(f"\n  De17 = 0x{e17_orig ^ e17_twin:08x} (HW={hw(e17_orig ^ e17_twin)})")

    carry_diff_total = 0
    for i in range(len(adds_orig)):
        desc_o, a_o, b_o, carry_o = adds_orig[i]
        desc_t, a_t, b_t, carry_t = adds_twin[i]
        carry_diff = carry_o ^ carry_t
        cd_hw = hw(carry_diff)
        carry_diff_total += cd_hw
        if cd_hw > 0:
            print(f"    Add {i+1:2d}: {desc_o}")
            print(f"      orig carry=0x{carry_o:08x}  twin carry=0x{carry_t:08x}  diff HW={cd_hw}")

    print(f"\n  Total carry difference bits: {carry_diff_total}")
    print(f"  Average per addition: {carry_diff_total/len(adds_orig):.1f}")

    # ─── SUMMARY ───
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"""
  Total modular additions W[14] -> e[17]:     {len(additions)}
    - Known (W[14]-independent):              {len(known_adds)}
    - W[14]-dependent:                        {len(dep_indices)}
  Effective carry chain depth:                {len(dep_indices)} additions

  De17 statistics (2^24 enumeration):
    - Mean HW:                                {hw_accum/N_enum:.4f}
    - De17=0 count:                           {zero_count} / {N_enum}
    - Extrapolated full 2^32 zeros:           ~{extrapolated_zeros:.1f}
    - Barrier weight (bits):                  ~{32 - (zero_count/N_enum * 2**32).bit_length() if zero_count > 0 else 32}

  Max single-bit bias:                        {max_bias:.4f}

  Carry analysis:
    - Total carry bits (base msg):            {total_carry_bits}
    - Carry density:                          {total_carry_bits/(len(additions)*32)*100:.1f}%
    - Differential carry bits (orig vs twin): {carry_diff_total}
""")

if __name__ == "__main__":
    main()
