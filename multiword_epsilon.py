#!/usr/bin/env python3
"""
Multi-Word ε — Axis 8.2: Where exactly does ε collapse?
==========================================================

Axis 8: ε=0.73 at R=64 for W[1..15]=0 (1-word message).
Question: how many active W words kill ε?

  W[0] only:     ε ≈ 0.73 (confirmed)
  W[0..1]:       ε = ?
  W[0..3]:       ε = ?
  W[0..7]:       ε = ?
  W[0..15]:      ε = 0 (expected from RAYON)

The transition point = maximum message length for RAYON-style attack.
"""

import numpy as np
import time, sys
from math import log2

MASK = 0xFFFFFFFF
H0 = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
     0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
     0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
     0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
     0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
     0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)

def msg_sched(W16):
    W=list(W16)+[0]*48
    for i in range(16,64): W[i]=(sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16])&MASK
    return W

def sha256_bit(W16, target_w, target_b):
    """Compute single output bit."""
    W=msg_sched(W16); a,b,c,d,e,f,g,h=H0
    for r in range(64):
        T1=(h+Sig1(e)+Ch(e,f,g)+K[r]+W[r])&MASK; T2=(Sig0(a)+Maj(a,b,c))&MASK
        h,g,f=g,f,e; e=(d+T1)&MASK; d,c,b=c,b,a; a=(T1+T2)&MASK
    out = [(v+iv)&MASK for v,iv in zip([a,b,c,d,e,f,g,h],H0)]
    return (out[target_w] >> target_b) & 1


def dfs_multiword(n_bits_per_word, n_words, target_w, target_b, target_val, rng):
    """
    DFS with propagation on n_words active W words, n_bits_per_word bits each.
    Other W words = 0.
    Returns number of nodes explored.
    """
    n_total = n_bits_per_word * n_words
    brute = 1 << n_total
    nodes = 0
    max_nodes = min(brute, 100000)  # safety cap

    def make_W16(partial_bits):
        W16 = [0]*16
        for w in range(n_words):
            val = 0
            for b in range(n_bits_per_word):
                bit_idx = w * n_bits_per_word + b
                val |= ((partial_bits >> bit_idx) & 1) << b
            W16[w] = val
        return W16

    def search(depth, partial):
        nonlocal nodes
        nodes += 1
        if nodes > max_nodes:
            return False

        if depth == n_total:
            W16 = make_W16(partial)
            return sha256_bit(W16, target_w, target_b) == target_val

        # Try both values for current bit
        for v in [0, 1]:
            new_partial = partial | (v << depth)

            # Propagation: if few free bits left, check exhaustively
            n_free = n_total - depth - 1
            if n_free <= 7:
                vals = set()
                for rem in range(1 << n_free):
                    full = new_partial | (rem << (depth + 1))
                    W16 = make_W16(full)
                    ov = sha256_bit(W16, target_w, target_b)
                    vals.add(ov)
                    if len(vals) > 1:
                        break
                if len(vals) == 1:
                    if target_val in vals:
                        return True
                    else:
                        continue  # prune
            if search(depth + 1, new_partial):
                return True
        return False

    search(0, 0)
    return min(nodes, max_nodes)


def experiment_multiword(seed=25000):
    print("="*70)
    print("MULTI-WORD ε — Where does DFS ε collapse?")
    print("="*70)

    rng = np.random.RandomState(seed)
    target_w, target_b = 4, 0
    N_trials = 15

    # Test configurations: (n_bits_per_word, n_words)
    configs = [
        (10, 1, "1 word, 10 bits"),
        (8,  2, "2 words, 8 bits each"),
        (6,  3, "3 words, 6 bits each"),
        (5,  4, "4 words, 5 bits each"),
        (4,  5, "5 words, 4 bits each"),
        (4,  4, "4 words, 4 bits"),
        (3,  6, "6 words, 3 bits each"),
        (3,  5, "5 words, 3 bits"),
        (2,  8, "8 words, 2 bits each"),
    ]

    print(f"\n  Target: H[{target_w}][b{target_b}], {N_trials} trials each")
    print(f"\n  {'config':>25} | {'n_total':>8} | {'brute':>8} | {'mean nodes':>11} | {'ε':>7} | {'vs 1-word':>10}")
    print(f"  {'-'*25}-+-{'-'*8}-+-{'-'*8}-+-{'-'*11}-+-{'-'*7}-+-{'-'*10}")

    baseline_eps = None

    for n_bpw, n_w, desc in configs:
        n_total = n_bpw * n_w
        brute = 1 << n_total
        node_counts = []

        t0 = time.time()
        for trial in range(N_trials):
            # Random target
            target_input = [int(rng.randint(0, 1<<n_bpw)) for _ in range(n_w)] + [0]*(16-n_w)
            target_val = sha256_bit(target_input, target_w, target_b)
            nodes = dfs_multiword(n_bpw, n_w, target_w, target_b, target_val, rng)
            node_counts.append(nodes)

        mean_nodes = np.mean(node_counts)
        eps = 1 - log2(max(mean_nodes, 1)) / max(n_total, 1)
        elapsed = time.time() - t0

        if baseline_eps is None:
            baseline_eps = eps
        ratio = eps / max(baseline_eps, 0.001)

        capped = sum(1 for n in node_counts if n >= min(brute, 100000))
        marker = " ★" if eps > 0.5 else (" ✗" if eps < 0.1 or capped > N_trials//2 else "")
        print(f"  {desc:>25} | {n_total:8d} | {brute:8d} | {mean_nodes:11.0f} | {eps:7.3f} | {ratio:10.2f}×{marker} ({elapsed:.1f}s)")
        if capped > 0:
            print(f"  {'':>25}   ({capped}/{N_trials} capped at {min(brute,100000)})")

    # Schedule zero analysis for multi-word
    print(f"\n  --- Schedule zeros for different n_words ---")
    for n_w in [1, 2, 4, 8, 16]:
        W16 = [42]*n_w + [0]*(16-n_w)
        W = msg_sched(W16)
        zeros = [r for r in range(16, 64) if W[r] == 0]
        print(f"    {n_w:2d} words active: {len(zeros)} zero schedule words: {zeros[:10]}{'...' if len(zeros)>10 else ''}")

    return


if __name__ == '__main__':
    print("Multi-Word ε — Axis 8.2")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    experiment_multiword()
