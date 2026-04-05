#!/usr/bin/env python3
"""
Local symmetries of SHA-256.

Global symmetry: Aut(SHA-256) = {id} (known).
But LOCAL symmetry: do some rounds behave SIMILARLY?

If rounds r1 and r2 have similar structure → differential path
through r1 can be "reused" at r2. This is NOT periodicity
but a GROUP-THEORETIC property: the round functions at r1 and r2
are conjugate (related by a transformation).

Questions:
Q1: Do any rounds have similar GPK profiles?
Q2: Do any rounds have similar sensitivity patterns?
Q3: Does the round function at r depend on K[r] in a way that
    creates "families" of similar rounds?
Q4: Can we define a DISTANCE between rounds and find clusters?
"""

import random, math
from collections import defaultdict

MASK32 = 0xFFFFFFFF
K=[0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV=[0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
def rotr(x,n):return((x>>n)|(x<<(32-n)))&MASK32
def ssig0(x):return rotr(x,7)^rotr(x,18)^(x>>3)
def ssig1(x):return rotr(x,17)^rotr(x,19)^(x>>10)
def sig0(x):return rotr(x,2)^rotr(x,13)^rotr(x,22)
def sig1(x):return rotr(x,6)^rotr(x,11)^rotr(x,25)
def ch(e,f,g):return(e&f)^(~e&g)&MASK32
def maj(a,b,c):return(a&b)^(a&c)^(b&c)
def hw(x):return bin(x&MASK32).count('1')


def one_round(state, W_r, K_r):
    """Execute one SHA-256 round, return new state."""
    a,b,c,d,e,f,g,h = state
    T1=(h+sig1(e)+ch(e,f,g)+K_r+W_r)&MASK32
    T2=(sig0(a)+maj(a,b,c))&MASK32
    return ((T1+T2)&MASK32, a, b, c, (d+T1)&MASK32, e, f, g)


# ================================================================
# Q1: Round fingerprint via sensitivity
# For each round r: flip each bit of state, measure output change.
# The 256×256 sensitivity matrix IS the round's fingerprint.
# Two rounds with similar matrices = similar behavior.
# ================================================================

print("=" * 70)
print("Q1: Round sensitivity fingerprints")
print("=" * 70)

N = 500
random.seed(42)

# For each round: compute average HW of output change per input bit flip
# This gives a 256-element vector = "sensitivity fingerprint"

def round_fingerprint(r, n_samples):
    """Compute sensitivity fingerprint of round r."""
    # Average over random states and W values
    fingerprint = [0.0] * 256  # 8 regs × 32 bits input

    for _ in range(n_samples):
        # Random state
        state = tuple(random.randint(0, MASK32) for _ in range(8))
        W_r = random.randint(0, MASK32)

        out_base = one_round(state, W_r, K[r])

        for reg in range(8):
            for bit in range(32):
                # Flip bit (reg, bit) of input state
                state_flip = list(state)
                state_flip[reg] ^= (1 << bit)
                state_flip = tuple(state_flip)

                out_flip = one_round(state_flip, W_r, K[r])

                # Total HW of output change
                total_change = sum(hw(out_base[i] ^ out_flip[i]) for i in range(8))
                fingerprint[reg * 32 + bit] += total_change / n_samples

    return fingerprint

# Compute for a sample of rounds (full 64 too slow)
print(f"\n  Computing fingerprints for selected rounds (N={N} samples each)...")

rounds_to_test = [0, 1, 7, 8, 15, 16, 17, 30, 31, 32, 47, 48, 55, 56, 62, 63]
fingerprints = {}

for r in rounds_to_test:
    fingerprints[r] = round_fingerprint(r, N)

# Distance between fingerprints
def fp_distance(fp1, fp2):
    return math.sqrt(sum((a-b)**2 for a,b in zip(fp1, fp2)))

print(f"\n  Pairwise distances between round fingerprints:")
print(f"  {'':>4}", end="")
for r2 in rounds_to_test[:8]:
    print(f"  r={r2:<3}", end="")
print()

for r1 in rounds_to_test:
    print(f"  r={r1:<2}", end="")
    for r2 in rounds_to_test[:8]:
        d = fp_distance(fingerprints[r1], fingerprints[r2])
        if r1 == r2:
            print(f"    -- ", end="")
        elif d < 1.0:
            print(f"  {d:.2f}★", end="")
        else:
            print(f"  {d:.2f} ", end="")
    print()

# Find most similar pair
min_dist = 9999
min_pair = (0, 0)
for i, r1 in enumerate(rounds_to_test):
    for j, r2 in enumerate(rounds_to_test):
        if r1 >= r2: continue
        d = fp_distance(fingerprints[r1], fingerprints[r2])
        if d < min_dist:
            min_dist = d
            min_pair = (r1, r2)

print(f"\n  Most similar rounds: r={min_pair[0]} and r={min_pair[1]}, distance={min_dist:.3f}")


# ================================================================
# Q2: K-constant structure
# K[r] is the ONLY thing that differs between rounds
# (besides state, which is input-dependent).
# Group rounds by K[r] properties.
# ================================================================

print()
print("=" * 70)
print("Q2: K-constant structure — round families")
print("=" * 70)

# Group by HW(K[r])
hw_groups = defaultdict(list)
for r in range(64):
    hw_groups[hw(K[r])].append(r)

print(f"\n  Rounds grouped by HW(K[r]):")
for h in sorted(hw_groups.keys()):
    rounds = hw_groups[h]
    print(f"    HW={h:>2}: {len(rounds)} rounds: {rounds[:10]}{'...' if len(rounds)>10 else ''}")

# K[r] mod small numbers — any patterns?
print(f"\n  K[r] mod 8 distribution:")
mod8_groups = defaultdict(list)
for r in range(64):
    mod8_groups[K[r] % 8].append(r)
for m in range(8):
    print(f"    K mod 8 = {m}: {len(mod8_groups[m])} rounds")

# LSB pattern of K
lsb_pattern = ''.join(str(K[r] & 1) for r in range(64))
print(f"\n  K[r] LSB pattern: {lsb_pattern}")
print(f"  Ones: {lsb_pattern.count('1')}/64")

# Do similar K values give similar round behavior?
# Compare rounds with similar K (small |K[r1]-K[r2]|)
k_diffs = []
for r1 in range(64):
    for r2 in range(r1+1, 64):
        k_diff = abs(K[r1] - K[r2]) / 2**32
        k_diffs.append((k_diff, r1, r2))

k_diffs.sort()
print(f"\n  Most similar K-values:")
for kd, r1, r2 in k_diffs[:5]:
    print(f"    K[{r1}] ≈ K[{r2}]: |diff| = {kd:.6f} × 2^32")


# ================================================================
# Q3: Does round function commute with any transformation?
# Find T such that: round_r(T(state)) = T(round_r(state))
# for some nontrivial T.
#
# If T exists → T is a SYMMETRY of round r.
# If same T works for multiple rounds → exploitable.
# ================================================================

print()
print("=" * 70)
print("Q3: Round symmetries — does any transformation commute?")
print("=" * 70)

N3 = 1000
random.seed(123)

# Test candidate transformations:
# T1: bit rotation of each word by s
# T2: word permutation
# T3: XOR with constant
# T4: negation (NOT)

def apply_rotr_all(state, s):
    return tuple(rotr(x, s) for x in state)

def apply_not_all(state):
    return tuple((~x) & MASK32 for x in state)

def apply_xor_const(state, c):
    return tuple((x ^ c) & MASK32 for x in state)

def apply_swap_ae(state):
    """Swap a-chain and e-chain registers."""
    a,b,c,d,e,f,g,h = state
    return (e,f,g,h,a,b,c,d)

transforms = [
    ("ROTR all by 1", lambda s: apply_rotr_all(s, 1)),
    ("ROTR all by 16", lambda s: apply_rotr_all(s, 16)),
    ("NOT all", apply_not_all),
    ("XOR 0x80000000", lambda s: apply_xor_const(s, 0x80000000)),
    ("XOR 0xFFFFFFFF", lambda s: apply_xor_const(s, MASK32)),
    ("Swap a↔e chains", apply_swap_ae),
]

print(f"\n  Testing {len(transforms)} candidate symmetries on rounds 0,16,32,48,63")
print(f"  N={N3} random states per test\n")

for t_name, t_fn in transforms:
    commutes = defaultdict(int)

    for r in [0, 16, 32, 48, 63]:
        count = 0
        for _ in range(N3):
            state = tuple(random.randint(0, MASK32) for _ in range(8))
            W_r = random.randint(0, MASK32)

            # round(T(state)) vs T(round(state))
            lhs = one_round(t_fn(state), W_r, K[r])
            rhs = t_fn(one_round(state, W_r, K[r]))

            if lhs == rhs:
                count += 1

        commutes[r] = count / N3

    avg = sum(commutes.values()) / len(commutes)
    flag = " ★ SIGNAL" if avg > 0.01 else ""
    print(f"  {t_name:>25}: avg P(commute)={avg:.4f} per-round={dict(commutes)}{flag}")


# ================================================================
# Q4: Functional distance between rounds
# Define: d(r1, r2) = E[HW(round_r1(state, W) XOR round_r2(state, W))]
# Averaged over random state and W.
# If d(r1,r2) ≈ 0 → rounds behave identically.
# If d(r1,r2) < d(random) → rounds partially similar.
# ================================================================

print()
print("=" * 70)
print("Q4: Functional distance between rounds")
print("=" * 70)

N4 = 2000
random.seed(456)

# Pre-compute states and W values
test_states = [tuple(random.randint(0,MASK32) for _ in range(8)) for _ in range(N4)]
test_W = [random.randint(0, MASK32) for _ in range(N4)]

# Compute all round outputs
round_outputs = {}
for r in range(64):
    outputs = []
    for i in range(N4):
        out = one_round(test_states[i], test_W[i], K[r])
        outputs.append(out)
    round_outputs[r] = outputs

# Pairwise distance
print(f"\n  N={N4}, computing 64×64 distance matrix...")

dist_matrix = [[0.0]*64 for _ in range(64)]

for r1 in range(64):
    for r2 in range(r1+1, 64):
        total_hw = 0
        for i in range(N4):
            o1 = round_outputs[r1][i]
            o2 = round_outputs[r2][i]
            total_hw += sum(hw(o1[j] ^ o2[j]) for j in range(8))
        avg_hw = total_hw / N4
        dist_matrix[r1][r2] = avg_hw
        dist_matrix[r2][r1] = avg_hw

# Find closest pairs
pairs = []
for r1 in range(64):
    for r2 in range(r1+1, 64):
        pairs.append((dist_matrix[r1][r2], r1, r2))
pairs.sort()

print(f"\n  Top 10 CLOSEST round pairs:")
for d, r1, r2 in pairs[:10]:
    k_diff = abs(K[r1] - K[r2]) / 2**32
    print(f"    r={r1:>2} ↔ r={r2:>2}: dist={d:.2f} (K diff={k_diff:.4f})")

print(f"\n  Bottom 5 FARTHEST:")
for d, r1, r2 in pairs[-5:]:
    print(f"    r={r1:>2} ↔ r={r2:>2}: dist={d:.2f}")

# Overall distribution
all_dists = [d for d, _, _ in pairs]
print(f"\n  Distance stats: min={min(all_dists):.2f} max={max(all_dists):.2f} avg={sum(all_dists)/len(all_dists):.2f} std={(sum((x-sum(all_dists)/len(all_dists))**2 for x in all_dists)/len(all_dists))**0.5:.2f}")

# Correlation: round distance vs K distance
k_diffs_list = []
d_list = []
for d, r1, r2 in pairs:
    kd = abs(K[r1] - K[r2]) / 2**32
    k_diffs_list.append(kd)
    d_list.append(d)

def pearson(x,y):
    n=len(x);mx=sum(x)/n;my=sum(y)/n
    sx=math.sqrt(sum((a-mx)**2 for a in x)/n)
    sy=math.sqrt(sum((a-my)**2 for a in y)/n)
    if sx*sy==0:return 0
    return sum((a-mx)*(b-my) for a,b in zip(x,y))/(n*sx*sy)

corr_kd = pearson(k_diffs_list, d_list)
print(f"\n  Correlation(K-distance, round-distance) = {corr_kd:.4f}")
print(f"  {'K determines round similarity!' if abs(corr_kd) > 0.3 else 'K does NOT determine round similarity.'}")


# ================================================================
# Q5: The actual symmetry question for attack
# If rounds r1 and r2 are "similar": can a diff path at r1
# be REUSED at r2 with small cost?
#
# Concretely: take Wang pair, look at δstate[r1].
# Apply same δstate at r2. How much does it cost?
# ================================================================

print()
print("=" * 70)
print("Q5: Can differential structure be reused across rounds?")
print("=" * 70)

# For Wang pair: δstate[r] for each r.
# Take δstate[r1], inject it at round r2.
# Measure: does the diff propagate similarly?

N5 = 200
random.seed(789)

def get_wang_states():
    """Get one Wang pair and all intermediate states."""
    M = [random.randint(0,MASK32) for _ in range(16)]
    W1 = list(M); DW=[0]*16; DW[0]=1

    def states_full(M_in):
        W=list(M_in)+[0]*(64-len(M_in))
        for i in range(16,64):W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
        a,b,c,d,e,f,g,h=IV
        ss=[(a,b,c,d,e,f,g,h)]
        for r in range(64):
            T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32;T2=(sig0(a)+maj(a,b,c))&MASK32
            h,g,f=g,f,e;e=(d+T1)&MASK32;d,c,b=c,b,a;a=(T1+T2)&MASK32
            ss.append((a,b,c,d,e,f,g,h))
        return ss, W

    # Build Wang cascade
    for step in range(15):
        wi=step+1
        W2=[(W1[i]+DW[i])&MASK32 for i in range(16)]
        s1,_=states_full(W1);s2,_=states_full(W2)
        de=(s2[step+2][4]-s1[step+2][4])&MASK32;DW[wi]=(-de)&MASK32

    W2=[(W1[i]+DW[i])&MASK32 for i in range(16)]
    ss1, Wexp1 = states_full(W1)
    ss2, Wexp2 = states_full(W2)
    return ss1, ss2, Wexp1, Wexp2

# Measure: δstate at rounds 16-20 for Wang pairs
# Then: if we had the SAME δstate but at a different round, what happens?

random.seed(789)
ss1, ss2, W1, W2 = get_wang_states()

# δstate at round 17 (first post-barrier)
dstate17 = tuple((ss2[17][i] - ss1[17][i]) & MASK32 for i in range(8))
hw17 = sum(hw(d) for d in dstate17)

print(f"\n  Wang pair: δstate[17] HW = {hw17}")
print(f"  δstate[17] = {tuple(hex(d) for d in dstate17)}")

# Now: inject this same δstate at different rounds
# Take fresh random state, add δstate, run one round, measure output δ
print(f"\n  Injecting δstate[17] at different rounds:")
print(f"  {'r':>3} | {'HW(δ_out)':>10} | {'HW(δ_in)':>9} | {'amplification':>13}")

for r_inject in [0, 4, 8, 16, 17, 20, 32, 48, 63]:
    total_hw_out = 0

    for _ in range(N5):
        state = tuple(random.randint(0, MASK32) for _ in range(8))
        W_r = random.randint(0, MASK32)

        state2 = tuple((state[i] + dstate17[i]) & MASK32 for i in range(8))

        out1 = one_round(state, W_r, K[r_inject])
        out2 = one_round(state2, W_r, K[r_inject])

        hw_out = sum(hw(out1[i] ^ out2[i]) for i in range(8))
        total_hw_out += hw_out

    avg_out = total_hw_out / N5
    amp = avg_out / max(hw17, 1)
    print(f"  {r_inject:>3} | {avg_out:>10.1f} | {hw17:>9} | {amp:>13.3f}")


# ================================================================
# SYNTHESIS
# ================================================================

print()
print("=" * 70)
print("SYNTHESIS: Local symmetries")
print("=" * 70)

print("""
  Q1: Round sensitivity fingerprints — how similar are rounds?
  Q2: K-constant families — do K-values create round clusters?
  Q3: Commuting transformations — does any T commute with rounds?
  Q4: Functional distance — which rounds are closest?
  Q5: Differential reuse — can δstate at r1 work at r2?

  For attack: need rounds r1, r2 where:
  - differential structure at r1 "works" at r2 (low-cost transfer)
  - this creates a "shortcut" through the 44-round gap

  If all rounds equally distant → no symmetry to exploit.
  If clusters exist → potential for path reuse.
""")
