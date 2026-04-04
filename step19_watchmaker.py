"""
Step 19: Watchmaker — disassemble SHA-256 into individual mechanisms

Not trying to SOLVE. Trying to UNDERSTAND.
Remove each component one by one. See what breaks and what doesn't.

Components to test:
1. Ch(e,f,g) → replace with 0 or with e⊕f⊕g (linear version)
2. Maj(a,b,c) → replace with 0 or a⊕b⊕c
3. Σ₀ (rotation) → replace with identity
4. Σ₁ (rotation) → replace with identity
5. Carry (mod addition) → replace with XOR
6. Schedule (W[16+]) → set to 0 or keep W[r]=W[r%16]
7. Constants K → set to 0
8. IV → set to 0
9. Feedforward (H = IV + state) → remove

For each variant: measure collision fiber structure (from Step 18).
The component whose removal CREATES fiber structure = the "lock".
The component whose removal DESTROYS security = the "key mechanism".
"""

import numpy as np
from collections import Counter

N = 4
MASK = (1 << N) - 1
N_MSG = 4
N_INPUT = N * N_MSG
N_TOTAL = 1 << N_INPUT

IV_ORIG = [0x6, 0xB, 0x3, 0xA, 0x5, 0x9, 0x1, 0xF]
K_ORIG = [0x4, 0x2, 0xB, 0x7, 0xA, 0x3, 0xE, 0x5,
          0x9, 0x1, 0xD, 0x6, 0x0, 0x8, 0xC, 0xF]


def rotr(x, s):
    return ((x >> s) | (x << (N - s))) & MASK


def modular_sha(msg, R, config):
    """
    Configurable mini-SHA. Each component can be toggled on/off.

    config dict:
      ch: 'normal', 'zero', 'linear' (e⊕f⊕g)
      maj: 'normal', 'zero', 'linear' (a⊕b⊕c)
      sig0: 'normal', 'zero' (identity)
      sig1: 'normal', 'zero'
      add: 'normal', 'xor' (replace + with ⊕)
      schedule: 'normal', 'zero' (W[r]=0 for r≥N_MSG), 'repeat'
      K: 'normal', 'zero'
      iv: 'normal', 'zero'
    """
    ch_mode = config.get('ch', 'normal')
    maj_mode = config.get('maj', 'normal')
    sig0_mode = config.get('sig0', 'normal')
    sig1_mode = config.get('sig1', 'normal')
    add_mode = config.get('add', 'normal')
    k_mode = config.get('K', 'normal')
    iv_mode = config.get('iv', 'normal')

    def add_op(a, b):
        if add_mode == 'xor':
            return a ^ b
        return (a + b) & MASK

    def ch_op(e, f, g):
        if ch_mode == 'zero':
            return 0
        elif ch_mode == 'linear':
            return e ^ f ^ g
        return (e & f) ^ (~e & g) & MASK

    def maj_op(a, b, c):
        if maj_mode == 'zero':
            return 0
        elif maj_mode == 'linear':
            return a ^ b ^ c
        return (a & b) ^ (a & c) ^ (b & c)

    def sig0_op(x):
        if sig0_mode == 'zero':
            return x  # identity, not zero
        return rotr(x, 1) ^ rotr(x, 2) ^ rotr(x, 3)

    def sig1_op(x):
        if sig1_mode == 'zero':
            return x
        return rotr(x, 1) ^ rotr(x, 3) ^ (x >> 1)

    if iv_mode == 'zero':
        a, b, c, d, e, f, g, h = 0, 0, 0, 0, 0, 0, 0, 0
    else:
        a, b, c, d, e, f, g, h = IV_ORIG[:]

    W = list(msg) + [0] * max(0, R - N_MSG)

    for r in range(R):
        w_r = W[r] if r < len(W) else 0
        k_r = K_ORIG[r % len(K_ORIG)] if k_mode == 'normal' else 0

        T1 = add_op(add_op(add_op(add_op(h, sig1_op(e)), ch_op(e, f, g)), k_r), w_r)
        T2 = add_op(sig0_op(a), maj_op(a, b, c))
        a_new = add_op(T1, T2)
        e_new = add_op(d, T1)
        h, g, f = g, f, e
        e = e_new
        d, c, b = c, b, a
        a = a_new

    return (a, e)


def compute_fiber_mod(M, R, config):
    """Compute collision fiber for modular SHA."""
    base = modular_sha(M, R, config)
    fiber = set()
    for didx in range(1, N_TOTAL):
        dw = []
        tmp = didx
        for w in range(N_MSG):
            dw.append(tmp & MASK)
            tmp >>= N
        M2 = [(M[w] ^ dw[w]) for w in range(N_MSG)]
        out = modular_sha(M2, R, config)
        if out == base:
            fiber.add(didx)
    return fiber


def xor_closure_rate(fiber):
    """Quick XOR closure test."""
    if len(fiber) < 2:
        return 0, 0, 0
    fiber_z = fiber | {0}
    closed = 0
    total = 0
    sample = sorted(fiber)[:200]
    for i in range(len(sample)):
        for j in range(i+1, len(sample)):
            total += 1
            if (sample[i] ^ sample[j]) in fiber_z:
                closed += 1
    return closed, total, closed/max(total,1)


def analyze_component(name, config, M, R):
    """Analyze one SHA variant."""
    fiber = compute_fiber_mod(M, R, config)
    n_outputs = len(set(modular_sha(
        [(M[w] ^ ((didx >> (w*N)) & MASK)) for w in range(N_MSG)], R, config)
        for didx in range(N_TOTAL)))

    closed, total, rate = xor_closure_rate(fiber)
    expected_rate = len(fiber) / N_TOTAL

    return {
        'name': name,
        'fiber_size': len(fiber),
        'n_outputs': n_outputs,
        'xor_rate': rate,
        'expected_rate': expected_rate,
        'enrichment': rate / max(expected_rate, 1e-10),
    }


def main():
    import random
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(N_MSG)]
    R = 8

    configs = [
        ("FULL SHA",           {}),
        ("No Ch (=0)",         {'ch': 'zero'}),
        ("Ch linear (e⊕f⊕g)", {'ch': 'linear'}),
        ("No Maj (=0)",        {'maj': 'zero'}),
        ("Maj linear",         {'maj': 'linear'}),
        ("No Σ₀ (identity)",   {'sig0': 'zero'}),
        ("No Σ₁ (identity)",   {'sig1': 'zero'}),
        ("No rotations",       {'sig0': 'zero', 'sig1': 'zero'}),
        ("XOR only (no carry)",{'add': 'xor'}),
        ("No K constants",     {'K': 'zero'}),
        ("Zero IV",            {'iv': 'zero'}),
        ("No Ch, no Maj",      {'ch': 'zero', 'maj': 'zero'}),
        ("No Ch,Maj + XOR add",{'ch': 'zero', 'maj': 'zero', 'add': 'xor'}),
        ("Linear Ch+Maj",      {'ch': 'linear', 'maj': 'linear'}),
        ("Linear Ch+Maj+XOR",  {'ch': 'linear', 'maj': 'linear', 'add': 'xor'}),
        ("No rot, no carry",   {'sig0': 'zero', 'sig1': 'zero', 'add': 'xor'}),
        ("Only carry (no Ch,Maj,rot)", {'ch': 'zero', 'maj': 'zero', 'sig0': 'zero', 'sig1': 'zero'}),
    ]

    print(f"WATCHMAKER DECOMPOSITION — R={R}, M={[hex(w) for w in M]}")
    print(f"{'='*90}")
    print(f"{'Component':>25} {'|Fiber|':>8} {'Outputs':>8} {'XOR%':>7} {'Expect%':>8} {'Enrich':>7}")
    print(f"{'-'*90}")

    results = []
    for name, config in configs:
        r = analyze_component(name, config, M, R)
        results.append(r)
        print(f"{r['name']:>25} {r['fiber_size']:>8} {r['n_outputs']:>8} "
              f"{r['xor_rate']*100:>6.1f}% {r['expected_rate']*100:>7.2f}% "
              f"{r['enrichment']:>6.1f}×")

    # Find the components that MOST affect security
    print(f"\n{'='*90}")
    print(f"ANALYSIS: Which component contributes WHAT?")
    print(f"{'='*90}")

    full = results[0]
    for r in results[1:]:
        fiber_change = r['fiber_size'] / max(full['fiber_size'], 1)
        output_change = r['n_outputs'] / max(full['n_outputs'], 1)
        enrich_change = r['enrichment'] / max(full['enrichment'], 0.01)

        marker = ""
        if r['fiber_size'] > full['fiber_size'] * 5:
            marker = " ★ SECURITY BROKEN"
        elif r['fiber_size'] < full['fiber_size'] / 5:
            marker = " ★ HARDER"
        if r['enrichment'] > 10:
            marker += " ★★ STRUCTURE APPEARS"
        if r['n_outputs'] < 100:
            marker += " ★ FEW OUTPUTS"

        print(f"  {r['name']:>25}: fiber {fiber_change:>5.1f}×, "
              f"outputs {output_change:.2f}×, "
              f"XOR enrich {r['enrichment']:.1f}×{marker}")


if __name__ == "__main__":
    main()
