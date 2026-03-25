#!/usr/bin/env python3
"""
SCF: Bitslice Attack — атакуем SHA-256 побитово снизу вверх.

Нелогичная идея: никто не атакует хеш-функцию по 1 биту за раз.
Но бит 0 CARRY-FREE (F3). Это чистый GF(2) канал.

LSB-SHA: 8-bit state, 16-bit input, degree 2 (Ch/Maj).
Total space: 2^24 — можно перебрать!

Потом: Hensel lifting bit 0 → bit 1 → ... → bit 31.
Бит k зависит от бит 0..k-1 через carry. Это ПОСЛЕДОВАТЕЛЬНО решаемо!
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
def hw(x): return bin(x).count('1')
def bit(x,i): return (x>>i)&1

def expand_schedule(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W


# ============================================================
# LSB-SHA: SHA-256 restricted to bit position k
# ============================================================
def sha_bitslice(W16, target_bit):
    """
    Run SHA-256 but track only bit position `target_bit`.
    For target_bit=0: NO CARRY (pure GF(2)).
    For target_bit=k: carry depends on bits 0..k-1 (which we've already solved).

    Returns: 8-tuple of single bits (LSB of each register).
    """
    W = expand_schedule(W16)
    s = list(IV)

    for r in range(64):
        a,b,c,d,e,f,g,h = s
        T1 = add32(h, Sig1(e), Ch(e,f,g), K[r], W[r])
        T2 = add32(Sig0(a), Maj(a,b,c))
        s = [add32(T1,T2), a, b, c, add32(d,T1), e, f, g]

    return tuple(bit(s[i], target_bit) for i in range(8))


def sha_bitslice_xor(W16, target_bit):
    """SHA-256 with XOR-only (no carry), restricted to bit target_bit.
    For bit 0 this should equal sha_bitslice (since carry-free)."""
    W = list(W16)
    for i in range(16,64):
        W.append(sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16])
    s = list(IV)
    for r in range(64):
        a,b,c,d,e,f,g,h = s
        T1 = h^Sig1(e)^Ch(e,f,g)^K[r]^W[r]
        T2 = Sig0(a)^Maj(a,b,c)
        s = [T1^T2,a,b,c,d^T1,e,f,g]
    return tuple(bit(s[i], target_bit) for i in range(8))


# ============================================================
# EXP 1: Verify bit 0 is carry-free
# ============================================================
def exp1_verify_lsb(N):
    print("="*70)
    print("EXP 1: VERIFY LSB IS CARRY-FREE")
    print("  sha_bitslice(W, bit=0) should equal sha_bitslice_xor(W, bit=0)")
    print("="*70)

    match = 0
    for _ in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        real = sha_bitslice(W, 0)
        xor = sha_bitslice_xor(W, 0)
        if real == xor:
            match += 1

    print(f"\n  {match}/{N} matches → {'✓ CARRY-FREE' if match==N else '✗ HAS CARRY'}")

    # Also check bit 1 (should have carry)
    match1 = 0
    for _ in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        real = sha_bitslice(W, 1)
        xor = sha_bitslice_xor(W, 1)
        if real == xor:
            match1 += 1

    print(f"  Bit 1: {match1}/{N} matches → {'carry-free' if match1==N else f'has carry ({N-match1} differ)'}")

    # Bit 2
    match2 = 0
    for _ in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        real = sha_bitslice(W, 2)
        xor = sha_bitslice_xor(W, 2)
        if real == xor:
            match2 += 1
    print(f"  Bit 2: {match2}/{N} matches → carry differs in {N-match2}")


# ============================================================
# EXP 2: LSB collision search (bit 0 only)
# ============================================================
def exp2_lsb_collision(N):
    print("\n" + "="*70)
    print("EXP 2: LSB COLLISION SEARCH")
    print("  Find (W1, W2) where sha(W1)[bit0] = sha(W2)[bit0] for all 8 regs")
    print("  This is an 8-bit target → expect collision in ~2^4 = 16 tries")
    print("="*70)

    # Fix W1, search for W2 differing only in W[0]
    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    target = sha_bitslice(W1, 0)

    found = 0
    tries = 0
    first_found = None

    for delta in range(1, min(N*100, 100000)):
        W2 = list(W1)
        W2[0] ^= delta
        result = sha_bitslice(W2, 0)
        tries += 1

        if result == target:
            found += 1
            if first_found is None:
                first_found = delta
            if found <= 3:
                print(f"  LSB collision #{found}: δW[0] = 0x{delta:08x} (HW={hw(delta)}), tries={tries}")

    rate = found/tries if tries > 0 else 0
    expected = 1/256  # 2^8 possible outputs
    print(f"\n  Total: {found} collisions in {tries} tries")
    print(f"  Rate: {rate:.6f} (expected ~{expected:.6f} = 1/256)")
    print(f"  Ratio: {rate/expected:.2f}x expected")

    return first_found


# ============================================================
# EXP 3: Multi-bit collision — bits 0..k
# ============================================================
def exp3_multibit_collision(N):
    print("\n" + "="*70)
    print("EXP 3: MULTI-BIT COLLISION — bits 0..k simultaneously")
    print("  How many low bits can we match?")
    print("="*70)

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    # Get full SHA output for W1
    s1 = sha_rounds_full(W1)

    for k_bits in [1, 2, 3, 4, 5, 6, 8, 10, 12, 16]:
        mask = (1 << k_bits) - 1  # Low k_bits mask
        target = tuple(s1[i] & mask for i in range(8))
        target_space = 1 << (8 * k_bits)

        found = 0
        tries = 0
        max_tries = min(N * 500, target_space * 4, 500000)

        for trial in range(max_tries):
            W2 = list(W1)
            # Modify only low bits of W[0] for small k, more words for larger k
            if k_bits <= 8:
                W2[0] ^= (trial + 1) & ((1 << min(k_bits * 2, 32)) - 1)
                if W2[0] == W1[0]: continue
            else:
                W2[trial % 16] ^= int.from_bytes(os.urandom(4),'big')

            s2 = sha_rounds_full(W2)
            result = tuple(s2[i] & mask for i in range(8))
            tries += 1

            if result == target:
                found += 1
                if found == 1:
                    print(f"  k={k_bits:2d}: First collision at try {tries} "
                          f"(expected ~{target_space**0.5:.0f})")
                if found >= 3:
                    break

        if found == 0:
            print(f"  k={k_bits:2d}: No collision in {tries} tries "
                  f"(space=2^{8*k_bits}, need ~2^{4*k_bits})")
        else:
            print(f"  k={k_bits:2d}: {found} collisions in {tries} tries "
                  f"(space=2^{8*k_bits})")


def sha_rounds_full(W16):
    W=expand_schedule(W16); s=list(IV)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return s


# ============================================================
# EXP 4: Hensel lifting — extend LSB collision to higher bits
# ============================================================
def exp4_hensel_lift(N):
    print("\n" + "="*70)
    print("EXP 4: HENSEL LIFTING — LSB collision → full collision")
    print("  Start with bit-0 collision. Try to extend to bits 0..k.")
    print("  At each level, carry from lower bits is KNOWN.")
    print("="*70)

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    s1_full = sha_rounds_full(W1)

    # Step 1: Find bit-0 collision
    print("\n  Step 1: Find LSB collision")
    target_0 = tuple(bit(s1_full[i], 0) for i in range(8))

    lsb_deltas = []
    for delta in range(1, 100000):
        W2 = list(W1)
        W2[0] ^= delta
        s2 = sha_rounds_full(W2)
        if tuple(bit(s2[i], 0) for i in range(8)) == target_0:
            lsb_deltas.append(delta)
            if len(lsb_deltas) >= 20:
                break

    print(f"  Found {len(lsb_deltas)} LSB collisions")
    if not lsb_deltas:
        print("  No LSB collision found!")
        return

    # Step 2: Among LSB collisions, find ones that also match bit 1
    print("\n  Step 2: Hensel lift — which LSB collisions extend to bit 1?")
    bit1_matches = []
    for delta in lsb_deltas:
        W2 = list(W1)
        W2[0] ^= delta
        s2 = sha_rounds_full(W2)
        match_bits = 0
        for b in range(32):
            if all(bit(s1_full[i], b) == bit(s2[i], b) for i in range(8)):
                match_bits += 1
            else:
                break
        if match_bits >= 2:
            bit1_matches.append((delta, match_bits))

    print(f"  {len(bit1_matches)}/{len(lsb_deltas)} extend to bit 1+")
    for delta, mb in bit1_matches[:5]:
        print(f"    δ=0x{delta:08x}: matches {mb} consecutive low bits")

    # Step 3: General Hensel — how many consecutive bits can we match?
    print("\n  Step 3: Best consecutive bit match from LSB collisions")
    best_match = 0
    best_delta = 0

    for delta in lsb_deltas:
        W2 = list(W1)
        W2[0] ^= delta
        s2 = sha_rounds_full(W2)
        consec = 0
        for b in range(32):
            if all(bit(s1_full[i], b) == bit(s2[i], b) for i in range(8)):
                consec += 1
            else:
                break
        if consec > best_match:
            best_match = consec
            best_delta = delta

    print(f"  Best: {best_match} consecutive bits from LSB (δ=0x{best_delta:08x})")
    expected = 1  # Random: ~1 bit match on average
    print(f"  Expected random: ~1 bit")
    if best_match > 2:
        print(f"  ★ {best_match} > 2: Hensel lifting shows promise!")

    # Step 4: Broader search — multi-word delta, max consecutive low-bit match
    print("\n  Step 4: Broader search — multi-word δ")
    best_broad = 0
    best_broad_hw = 256

    for trial in range(min(N * 50, 50000)):
        W2 = list(W1)
        # Small perturbation
        word = trial % 16
        W2[word] ^= (trial // 16 + 1)

        s2 = sha_rounds_full(W2)
        consec = 0
        for b in range(32):
            if all(bit(s1_full[i], b) == bit(s2[i], b) for i in range(8)):
                consec += 1
            else:
                break

        total_hw = sum(hw(s1_full[i] ^ s2[i]) for i in range(8))

        if consec > best_broad or (consec == best_broad and total_hw < best_broad_hw):
            best_broad = consec
            best_broad_hw = total_hw
            if consec >= 3:
                print(f"    {consec} bits match (HW diff={total_hw}), δW[{word}]=0x{W1[word]^W2[word]:08x}")

    print(f"\n  Best broad: {best_broad} consecutive low bits, HW={best_broad_hw}")


# ============================================================
# EXP 5: Bit-plane independence — are low bits of output correlated?
# ============================================================
def exp5_bitplane_correlation(N):
    print("\n" + "="*70)
    print("EXP 5: BIT-PLANE CORRELATION")
    print("  Are output bits 0,1,2,... of SHA correlated?")
    print("  If bit k of output depends weakly on bit k of input → structure!")
    print("="*70)

    # For each input bit flip in W[0], which output bits change?
    bit_response = [[0]*32 for _ in range(32)]  # input_bit × output_bit

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        s1 = sha_rounds_full(W1)

        for in_bit in range(32):
            W2 = list(W1)
            W2[0] ^= (1 << in_bit)
            s2 = sha_rounds_full(W2)

            # Check which bits of a[64] changed
            da = s1[0] ^ s2[0]
            for out_bit in range(32):
                if bit(da, out_bit):
                    bit_response[in_bit][out_bit] += 1

    # Normalize
    for i in range(32):
        for j in range(32):
            bit_response[i][j] /= N

    # Show diagonal and off-diagonal
    print(f"\n  P(output bit j flips | input bit i flips) for register a:")
    print(f"  Diagonal = P(out[i] flips | in[i] flips)")
    print(f"  Expected: 0.5 for all (random)")

    hdr = 'in\\out'
    print(f"\n  {hdr:>6}", end="")
    for j in [0,1,2,3,4,8,16,24,31]:
        print(f" {j:>5}", end="")
    print()

    for i in [0,1,2,3,4,8,16,24,31]:
        print(f"  {i:>6}", end="")
        for j in [0,1,2,3,4,8,16,24,31]:
            v = bit_response[i][j]
            marker = "*" if abs(v-0.5) > 0.05 else " "
            print(f" {v:.3f}{marker}" if abs(v-0.5)>0.03 else f" {v:.3f} ", end="")
        print()

    # Diagonal analysis
    diag = [bit_response[i][i] for i in range(32)]
    avg_diag = sum(diag)/32
    print(f"\n  Diagonal average: {avg_diag:.4f} (0.5 = no bit-plane structure)")
    print(f"  Off-diagonal average: {sum(bit_response[i][j] for i in range(32) for j in range(32) if i!=j)/(32*31):.4f}")


# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200

    print("="*70)
    print("SCF: BITSLICE ATTACK — побитовая атака снизу вверх")
    print("  Нелогичный подход: атакуем 1 бит за раз")
    print("="*70)

    exp1_verify_lsb(min(N, 500))
    exp2_lsb_collision(N)
    exp3_multibit_collision(N)
    exp4_hensel_lift(N)
    exp5_bitplane_correlation(min(N, 100))

    print("\n" + "="*70)
    print("ИТОГ: BITSLICE ATTACK")
    print("="*70)
    print("""
  Подход: атаковать SHA-256 побитово, начиная с carry-free LSB.

  Бит 0: carry-free → чистый GF(2) → collision за O(2^4) ← ТРИВИАЛЬНО
  Бит 1: carry из бита 0 ИЗВЕСТЕН → quasi-GF(2)
  Бит k: carry из бит 0..k-1 ПОСЛЕДОВАТЕЛЬНО вычислим

  Вопрос: можно ли Hensel-lift LSB collision → full collision?
  Стоимость: O(2^4) на бит × 32 бита = O(2^4 × 32) = O(2^9)?

  НО: carry на бите k зависит от ВСЕХ нижних бит ВСЕГО state,
  а state зависит от message через ВСЕ 64 раунда.
  Это создаёт экспоненциальный blowup при подъёме.

  Реальная стоимость: O(2^k) на подъём к биту k → суммарно O(2^32).
  Это лучше birthday (2^128)? Нет — это стоимость ОДНОГО
  НЕПОЛНОГО matching, не collision.
    """)
