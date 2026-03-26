#!/usr/bin/env python3
"""
SCF: DECOUPLED ROUND PARAMETERIZATION — четвёртый фреймворк.

Старая параметризация: W[0..15] → schedule → 64 раунда (СВЯЗАНЫ).
Новая: a-sequence a[0..64] → required W[r] для каждого раунда.

Свобода: a[0..8] фиксированы (IV), a[9..64] = 56 свободных значений.
Из a-sequence: W[r] = inverse_round(a[r], a[r+1], a[r-1..r-7]).
Schedule constraint: W[r] = σ1(W[r-2])+W[r-7]+σ0(W[r-15])+W[r-16].

COLLISION = a1-sequence и a2-sequence совпадают на r=53..64.
CONSTRAINT = W[r] из обеих a-sequences удовлетворяют schedule.

Новая задача: найти a2-sequence (56 свободных 32-бит значений),
такую что:
1. a2[53..64] = a1[53..64] (collision condition)
2. W2[r] = inverse(a2-seq) совместимы с schedule (16 уравнений)

56 свободных - 12 collision - 16 schedule = 28 свободных!
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


def D_val(a_seq, r):
    """D[r] = Σ0(a[r-1]) + Maj(a[r-1],a[r-2],a[r-3]) - a[r-4]"""
    return sub32(add32(Sig0(a_seq[r-1]), Maj(a_seq[r-1],a_seq[r-2],a_seq[r-3])),
                 a_seq[r-4])

def e_val(a_seq, r):
    """e[r] = a[r] - D[r]"""
    if r < 4:
        return IV[4]  # Simplified for early rounds
    return sub32(a_seq[r], D_val(a_seq, r))

def W_from_a_seq(a_seq, r):
    """Compute W[r] from a-sequence using inverse round.
    a[r+1] = T1 + T2 → T1 = a[r+1] - T2
    T2 = Σ0(a[r]) + Maj(a[r],a[r-1],a[r-2])
    W[r] = T1 - h - Σ1(e) - Ch(e,f,g) - K[r]
    where h=e[r-3], e=e[r], f=e[r-1], g=e[r-2]
    """
    if r < 4 or r+1 >= len(a_seq):
        return 0

    T2 = add32(Sig0(a_seq[r]), Maj(a_seq[r], a_seq[r-1], a_seq[r-2]))
    T1 = sub32(a_seq[r+1], T2)

    e_r = e_val(a_seq, r)
    e_r1 = e_val(a_seq, r-1) if r >= 5 else IV[5]
    e_r2 = e_val(a_seq, r-2) if r >= 6 else IV[6]
    h_r = e_val(a_seq, r-3) if r >= 7 else IV[7]

    return sub32(sub32(sub32(sub32(T1, h_r), Sig1(e_r)), Ch(e_r, e_r1, e_r2)), K[r])


def schedule_consistency(W_seq):
    """Check: how well does W[0..63] satisfy schedule constraints?
    W[r] should equal σ1(W[r-2])+W[r-7]+σ0(W[r-15])+W[r-16] for r≥16."""
    total_mm = 0
    for r in range(16, 64):
        expected = add32(sig1(W_seq[r-2]), W_seq[r-7], sig0(W_seq[r-15]), W_seq[r-16])
        total_mm += hw(W_seq[r] ^ expected)
    return total_mm


# ============================================================
# EXP 1: A-sequence → W-sequence → schedule consistency
# ============================================================
def exp1_a_to_W_consistency(N):
    print("="*70)
    print("EXP 1: A-SEQUENCE → W-SEQUENCE → SCHEDULE CHECK")
    print("  Random a-sequence → compute W[r] for each r → check schedule")
    print("="*70)

    # First: verify with REAL SHA (a-sequence from actual message)
    print(f"\n  Verification with real SHA:")
    for trial in range(3):
        W_real = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W_exp = expand_real(W_real)
        s = list(IV); a_seq = [s[0]]
        for r in range(64):
            a,b,c,d,e,f,g,h = s
            T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W_exp[r])
            T2=add32(Sig0(a),Maj(a,b,c))
            s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
            a_seq.append(s[0])

        # Reconstruct W from a-sequence
        W_reconstructed = [0]*64
        for r in range(4, 64):
            W_reconstructed[r] = W_from_a_seq(a_seq, r)

        # Compare with actual
        match = sum(1 for r in range(4, 64) if W_reconstructed[r] == W_exp[r])
        print(f"    Trial {trial}: W match {match}/60 rounds")

        # Schedule consistency of reconstructed W
        sc = schedule_consistency(W_reconstructed)
        sc_real = schedule_consistency(W_exp)
        print(f"    Schedule consistency: reconstructed={sc}, real={sc_real}")

    # Now: RANDOM a-sequence (not from any real message)
    print(f"\n  Random a-sequence (not from real SHA):")
    for trial in range(3):
        a_seq = [IV[0]]  # Start from IV[0]
        for r in range(1, 65):
            a_seq.append(int.from_bytes(os.urandom(4),'big'))

        W_from_random_a = [0]*64
        for r in range(4, 64):
            W_from_random_a[r] = W_from_a_seq(a_seq, r)

        sc = schedule_consistency(W_from_random_a)
        print(f"    Trial {trial}: schedule consistency = {sc} bits mismatch")
        print(f"    (expected if random: ~{48*16} = {48*16})")


# ============================================================
# EXP 2: Collision a-sequence → schedule constraint count
# ============================================================
def exp2_collision_constraints():
    print("\n" + "="*70)
    print("EXP 2: COLLISION CONSTRAINTS IN A-SEQUENCE SPACE")
    print("="*70)

    # Real message gives a-sequence with schedule consistency = 0.
    # Collision requires a2[53..64] = a1[53..64] AND schedule consistency = 0.

    # Freedom count:
    # a[0..8]: fixed by IV (9 values, but a[0]=IV[0] only)
    # Actually: a[0] = IV[0] fixed.
    # a[1..64]: 64 free values × 32 bits = 2048 bits

    # BUT: not all a-sequences are reachable!
    # Only those from valid W[0..15] (512 bits).
    # SHA maps 512 bits → 2048 bits (a[1..64]).
    # Image has dimension ≤ 512. So only 512/2048 of a-sequences are valid.

    print(f"""
  Dimension counting:
    Total a-sequence space: 64 values × 32 bits = 2048 bits
    Reachable from W[0..15]: ≤ 512 bits (= dim of input)

    Collision conditions: a2[53..64] = a1[53..64] → 12×32 = 384 bits

    IN INPUT SPACE (W[0..15]):
      Freedom: 512 bits
      Collision: 384 bits (projected from a-space)
      Free: 512 - 384 = 128 bits → O(2^64) birthday

    IN A-SEQUENCE SPACE:
      Freedom: 2048 bits (but only 512 reachable)
      Collision: 384 bits
      Schedule: 48×32 = 1536 bits of constraints
      512 reachable - 384 collision = 128 free → SAME answer

    The constraint count is CONSISTENT regardless of parameterization.
    The 128-bit kernel is a TOPOLOGICAL invariant of the collision variety.
    """)

    # Verify: GF(2) rank of a-sequence Jacobian
    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    s=list(IV); a1=[s[0]]
    W_exp = expand_real(W1)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W_exp[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        a1.append(s[0])

    # Jacobian: ∂a[53..64]/∂W[0..15] — 384×512 over GF(2)
    rows = []
    for target_r in range(53, 65):
        for bit in range(32):
            row = 0
            for word in range(16):
                for wbit in range(32):
                    W2 = list(W1); W2[word] ^= (1<<wbit)
                    s2=list(IV); W2_exp=expand_real(W2)
                    for r in range(target_r+1):
                        a,b,c,d,e,f,g,h=s2
                        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W2_exp[r])
                        T2=add32(Sig0(a),Maj(a,b,c))
                        s2=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
                    a2_r = s2[0]
                    if (a1[target_r] ^ a2_r) >> bit & 1:
                        row |= 1 << (word*32+wbit)
            rows.append(row)

    # Rank
    m = list(rows); rank = 0
    for bp in range(511, -1, -1):
        mask = 1 << bp
        pv = -1
        for i in range(rank, len(m)):
            if m[i] & mask: pv=i; break
        if pv == -1: continue
        m[rank],m[pv] = m[pv],m[rank]
        for i in range(len(m)):
            if i!=rank and m[i]&mask: m[i]^=m[rank]
        rank += 1

    print(f"  GF(2) Jacobian ∂a[53..64]/∂W[0..15]:")
    print(f"    Size: {len(rows)} × 512")
    print(f"    Rank: {rank}")
    print(f"    Kernel dim: {512-rank}")
    print(f"    → Birthday in kernel: O(2^{(512-rank)//2})")


# ============================================================
# EXP 3: Can we decouple at all?
# ============================================================
def exp3_decouple_test():
    print("\n" + "="*70)
    print("EXP 3: DECOUPLING TEST")
    print("  Idea: split W[0..15] into W_early[0..7] and W_late[8..15]")
    print("  If a[53..64] depends mainly on W_late → semi-decoupled!")
    print("="*70)

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    s1=list(IV); W1_exp=expand_real(W1)
    for r in range(65):
        if r < 64:
            a,b,c,d,e,f,g,h=s1
            T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W1_exp[r])
            T2=add32(Sig0(a),Maj(a,b,c))
            s1=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

    a1_seq = []
    s=list(IV)
    for r in range(65):
        a1_seq.append(s[0])
        if r < 64:
            a,b,c,d,e,f,g,h=s
            T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W1_exp[r])
            T2=add32(Sig0(a),Maj(a,b,c))
            s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

    # Sensitivity of a[53..64] to early (W[0..7]) vs late (W[8..15])
    early_sens = 0; late_sens = 0

    for word in range(16):
        total_hw_change = 0
        for bit in [0, 8, 16, 24, 31]:
            W2=list(W1); W2[word]^=(1<<bit)
            s2=list(IV); W2_exp=expand_real(W2)
            for r in range(65):
                if r < 64:
                    a,b,c,d,e,f,g,h=s2
                    T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W2_exp[r])
                    T2=add32(Sig0(a),Maj(a,b,c))
                    s2=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

            a2_tail = s2[0]  # This is a[64]
            # Get full a-seq for W2 would be expensive, just check hash
            H2 = [add32(IV[i],s2[i]) for i in range(8)]

            for r in range(53, 65):
                # Approximate: check a[64] only
                pass

            total_hw_change += hw(a1_seq[64] ^ a2_tail) if len(a1_seq) > 64 else 0

        if word < 8:
            early_sens += total_hw_change
        else:
            late_sens += total_hw_change

    print(f"  Sensitivity of a[64] to:")
    print(f"    W[0..7] (early):  {early_sens} total HW change")
    print(f"    W[8..15] (late):  {late_sens} total HW change")
    print(f"    Ratio late/early: {late_sens/(early_sens+1):.2f}")

    if late_sens > early_sens * 2:
        print(f"  ★ Late words dominate! Semi-decoupling possible!")
    else:
        print(f"  Both halves contribute equally. Full coupling.")


if __name__ == '__main__':
    exp1_a_to_W_consistency(1)
    exp2_collision_constraints()
    exp3_decouple_test()
