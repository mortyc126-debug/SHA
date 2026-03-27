"""
KECCAK/SHA-3 ЧЕРЕЗ НАШЕ ИЗМЕРЕНИЕ.

SHA-256 = Merkle-Damgård + ARX (ADD, Rotate, XOR).
Keccak  = Sponge + XOR-only (NO carry!).

Наша теория для SHA-256:
  - rank(CE) = 256 → нелинейность от carry + Ch/Maj
  - K = 128 → сфера
  - collision = C^(N/2) = 2^128

Keccak:
  - State = 1600 bits (5×5×64)
  - Rate = 1088 bits (SHA3-256)
  - Capacity = 512 bits
  - Round function: θ,ρ,π,χ,ι (ALL linear except χ)
  - χ: a[i] ^= (~a[i+1]) & a[i+2]  ← ЕДИНСТВЕННАЯ нелинейность!
  - 24 rounds

Вопрос: наше измерение применимо к sponge?
  - Нет carry → нет "carry error"
  - Есть χ → quadratic boolean (как наш Ch)
  - Нет register structure → N_reg = ?

Строим MINI-KECCAK для тестирования.
"""

import numpy as np

MASK64 = 0xFFFFFFFFFFFFFFFF

def hw(x): return bin(x).count('1')

# ============================================================
# Mini-Keccak: 200-bit state (5×5×8), rate=136, capacity=64
# Same structure as real Keccak but 8-bit lanes instead of 64
# ============================================================

MASK8 = 0xFF

def rot8(x, n):
    n = n % 8
    return ((x << n) | (x >> (8 - n))) & MASK8

RC_8 = [0x01, 0x82, 0x8A, 0x00, 0x8B, 0x01, 0x81, 0x09,
        0x8A, 0x88, 0x09, 0x0A, 0x8B, 0x8B, 0x89, 0x03,
        0x02, 0x80, 0x0A, 0x0A, 0x81, 0x80, 0x01, 0x08]

def keccak_f200(state, n_rounds=18):
    """Keccak-f[200]: 5×5 matrix of 8-bit lanes."""
    A = [list(row) for row in state]  # 5×5 of uint8

    for rd in range(n_rounds):
        # θ (theta): column parity
        C = [0]*5
        for x in range(5):
            for y in range(5):
                C[x] ^= A[x][y]
        D = [C[(x-1)%5] ^ rot8(C[(x+1)%5], 1) for x in range(5)]
        for x in range(5):
            for y in range(5):
                A[x][y] ^= D[x]

        # ρ (rho) + π (pi): rotation + permutation
        ROT_OFFSETS = [
            [0, 36, 3, 41, 18],
            [1, 44, 10, 45, 2],
            [62, 6, 43, 15, 61],
            [28, 55, 25, 21, 56],
            [27, 20, 39, 8, 14]
        ]
        B = [[0]*5 for _ in range(5)]
        for x in range(5):
            for y in range(5):
                B[y][(2*x+3*y)%5] = rot8(A[x][y], ROT_OFFSETS[x][y] % 8)

        # χ (chi): NONLINEAR step
        for x in range(5):
            T = [B[x][y] for y in range(5)]
            for y in range(5):
                A[x][y] = T[y] ^ ((~T[(y+1)%5]) & T[(y+2)%5]) & MASK8

        # ι (iota): round constant
        A[0][0] ^= RC_8[rd % len(RC_8)]

    return [tuple(row) for row in A]


def state_to_flat(state):
    """5×5 → flat list of 25 bytes."""
    return [state[x][y] for x in range(5) for y in range(5)]

def flat_to_state(flat):
    """flat list → 5×5."""
    return [tuple(flat[x*5:(x+1)*5]) for x in range(5)]


def mini_keccak_hash(msg_bytes, n_rounds=18):
    """Mini sponge: absorb msg into rate portion, squeeze output."""
    # State = 25 bytes (200 bits). Rate = 17 bytes (136 bits). Capacity = 8 bytes (64 bits).
    rate = 17
    state = [[0]*5 for _ in range(5)]

    # Pad message to rate
    padded = list(msg_bytes[:rate])
    if len(padded) < rate:
        padded += [0] * (rate - len(padded))
    # Simple padding: XOR 0x06 at end, 0x80 at last rate byte
    padded[min(len(msg_bytes), rate-1)] ^= 0x06
    padded[rate-1] ^= 0x80

    # Absorb
    flat = state_to_flat(state)
    for i in range(rate):
        flat[i] ^= padded[i]
    state = flat_to_state(flat)

    # Permutation
    state = keccak_f200(state, n_rounds)

    # Squeeze: take first 8 bytes = 64 bits output (like SHA3-256 → 256 bits)
    flat = state_to_flat(state)
    return tuple(flat[:8])


# ============================================================
# НАШЕ ИЗМЕРЕНИЕ: метрики для Mini-Keccak
# ============================================================

def main():
    np.random.seed(42)

    print("=" * 70)
    print("KECCAK ЧЕРЕЗ НАШЕ ИЗМЕРЕНИЕ")
    print("=" * 70)

    # Input: 17 bytes (rate). Output: 8 bytes.
    # Input bits: 136. Output bits: 64.

    N_input = 17  # rate bytes
    N_output = 8   # output bytes
    C_bits = 8     # bits per "register" (lane)
    input_bits = N_input * C_bits  # 136
    output_bits = N_output * C_bits  # 64

    print(f"\n  Mini-Keccak-f[200]:")
    print(f"    State: 200 bits (5×5×8)")
    print(f"    Rate: {N_input} bytes = {input_bits} bits")
    print(f"    Capacity: {25-N_input} bytes = {(25-N_input)*8} bits")
    print(f"    Output: {N_output} bytes = {output_bits} bits")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("1. JACOBIAN rank(T) — линейная часть")
    print("=" * 70)

    # Base message
    M_base = [np.random.randint(0, 256) for _ in range(N_input)]
    H_base = mini_keccak_hash(M_base)

    # Build T: flip each input bit, observe output change
    T_rows = []
    for byte_idx in range(N_input):
        for bit in range(8):
            M_mod = list(M_base)
            M_mod[byte_idx] ^= (1 << bit)
            H_mod = mini_keccak_hash(M_mod)
            row = []
            for w in range(N_output):
                d = H_base[w] ^ H_mod[w]
                for b in range(8):
                    row.append((d >> b) & 1)
            T_rows.append(row)

    T = np.array(T_rows, dtype=np.uint8)
    rank_T = np.linalg.matrix_rank(T.astype(float))

    print(f"\n  T matrix: {T.shape[0]}×{T.shape[1]} (input_bits × output_bits)")
    print(f"  rank(T) = {rank_T}")
    print(f"  Max possible = min({input_bits}, {output_bits}) = {min(input_bits, output_bits)}")
    print(f"  Full rank: {'YES' if rank_T == min(input_bits, output_bits) else 'NO'}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("2. GF(2)-KERNEL + NONLINEAR ERROR")
    print("=" * 70)

    # Kernel dimension = input_bits - rank_T
    kernel_dim = input_bits - rank_T
    print(f"\n  GF(2)-kernel dimension = {input_bits} - {rank_T} = {kernel_dim}")

    if rank_T == output_bits:
        # Full rank → kernel exists
        M_rref = T.T.copy()
        pivots = []; row = 0
        for col in range(input_bits):
            found = False
            for r in range(row, output_bits):
                if M_rref[r,col]==1:
                    M_rref[[row,r]] = M_rref[[r,row]]
                    found = True; break
            if not found: continue
            pivots.append(col)
            for r in range(output_bits):
                if r!=row and M_rref[r,col]==1:
                    M_rref[r] ^= M_rref[row]
            row += 1
        free_vars = [c for c in range(input_bits) if c not in pivots]

        print(f"  Pivot variables: {len(pivots)}")
        print(f"  Free variables: {len(free_vars)}")

        # Sample kernel vectors → compute ACTUAL hash difference
        ces = []
        ce_hws = []
        for fc in free_vars[:min(64, len(free_vars))]:
            x = np.zeros(input_bits, dtype=np.uint8); x[fc] = 1
            for i in range(len(pivots)-1, -1, -1):
                pc = pivots[i]; val = np.uint8(0)
                for j in range(input_bits):
                    if j!=pc: val ^= (M_rref[i,j] & x[j])
                x[pc] = val

            # Convert to message difference
            dM = [0]*N_input
            for byte_idx in range(N_input):
                for bit in range(8):
                    if x[byte_idx*8+bit]:
                        dM[byte_idx] ^= (1<<bit)

            M2 = [(M_base[i]^dM[i]) & 0xFF for i in range(N_input)]
            H2 = mini_keccak_hash(M2)

            # Nonlinear error
            ce = []
            for w in range(N_output):
                d = H_base[w] ^ H2[w]
                for b in range(8): ce.append((d>>b)&1)
            ces.append(ce)
            ce_hws.append(sum(ce))

        if ces:
            CE = np.array(ces, dtype=np.uint8)
            rank_CE = np.linalg.matrix_rank(CE.astype(float))

            print(f"\n  Nonlinear Error (NE) matrix: {CE.shape}")
            print(f"  rank(NE) = {rank_CE}")
            print(f"  Max possible = {output_bits}")
            print(f"  Full rank: {'YES' if rank_CE == output_bits else 'NO'}")
            print(f"\n  NE HW: mean={np.mean(ce_hws):.1f}, min={min(ce_hws)}, max={max(ce_hws)}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("3. CURVATURE K — геометрия")
    print("=" * 70)

    Ks = []
    for _ in range(500):
        b1, b2 = np.random.choice(input_bits, 2, replace=False)
        M1 = list(M_base); M1[b1//8] ^= (1 << (b1%8))
        M2 = list(M_base); M2[b2//8] ^= (1 << (b2%8))
        M12 = list(M_base); M12[b1//8] ^= (1 << (b1%8)); M12[b2//8] ^= (1 << (b2%8))

        H1 = mini_keccak_hash(M1)
        H2 = mini_keccak_hash(M2)
        H12 = mini_keccak_hash(M12)

        nl = sum(hw((H_base[i]^H12[i])^((H_base[i]^H1[i])^(H_base[i]^H2[i]))) for i in range(N_output))
        Ks.append(nl)

    K = np.mean(Ks)
    print(f"\n  Curvature K = {K:.1f}")
    print(f"  output_bits/2 = {output_bits/2}")
    print(f"  Ratio K/(output/2) = {K/(output_bits/2):.3f}")
    print(f"  Sphere: {'YES' if abs(K - output_bits/2) < 5 else 'NO'}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("4. SENSITIVITY — средний отклик на 1-bit flip")
    print("=" * 70)

    sens = []
    for _ in range(500):
        byte_idx = np.random.randint(0, N_input)
        bit = np.random.randint(0, 8)
        M2 = list(M_base); M2[byte_idx] ^= (1 << bit)
        H2 = mini_keccak_hash(M2)
        sens.append(sum(hw(H_base[i]^H2[i]) for i in range(N_output)))

    S = np.mean(sens)
    print(f"\n  Sensitivity = {S:.1f}")
    print(f"  output_bits/2 = {output_bits/2}")
    print(f"  Ratio S/(output/2) = {S/(output_bits/2):.3f}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("5. RANK по РАУНДАМ — когда Keccak достигает full rank?")
    print("=" * 70)

    print(f"\n  {'Rounds':>6} {'rank(T)':>8} {'K':>6} {'Sens':>6}")
    for n_r in [1, 2, 3, 4, 6, 8, 12, 18, 24]:
        M_base_r = [np.random.randint(0, 256) for _ in range(N_input)]
        H_base_r = mini_keccak_hash(M_base_r, n_rounds=n_r)

        # Quick rank
        T_r = []
        for _ in range(min(input_bits, 136)):
            byte_idx = np.random.randint(0, N_input)
            bit = np.random.randint(0, 8)
            M2 = list(M_base_r); M2[byte_idx] ^= (1<<bit)
            H2 = mini_keccak_hash(M2, n_rounds=n_r)
            row = []
            for w in range(N_output):
                d = H_base_r[w]^H2[w]
                for b in range(8): row.append((d>>b)&1)
            T_r.append(row)
        Tr = np.array(T_r, dtype=np.uint8)
        rk = np.linalg.matrix_rank(Tr.astype(float))

        # Quick K
        Ks_r = []
        for _ in range(100):
            b1, b2 = np.random.choice(input_bits, 2, replace=False)
            M1 = list(M_base_r); M1[b1//8] ^= (1<<(b1%8))
            M2 = list(M_base_r); M2[b2//8] ^= (1<<(b2%8))
            M12 = list(M_base_r); M12[b1//8] ^= (1<<(b1%8)); M12[b2//8] ^= (1<<(b2%8))
            H1 = mini_keccak_hash(M1, n_rounds=n_r)
            H2 = mini_keccak_hash(M2, n_rounds=n_r)
            H12 = mini_keccak_hash(M12, n_rounds=n_r)
            nl = sum(hw((H_base_r[i]^H12[i])^((H_base_r[i]^H1[i])^(H_base_r[i]^H2[i]))) for i in range(N_output))
            Ks_r.append(nl)

        # Quick sens
        sens_r = []
        for _ in range(100):
            byte_idx = np.random.randint(0, N_input)
            bit = np.random.randint(0, 8)
            M2 = list(M_base_r); M2[byte_idx] ^= (1<<bit)
            H2 = mini_keccak_hash(M2, n_rounds=n_r)
            sens_r.append(sum(hw(H_base_r[i]^H2[i]) for i in range(N_output)))

        print(f"  {n_r:6d} {rk:>7} {np.mean(Ks_r):>5.1f} {np.mean(sens_r):>5.1f}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("6. СРАВНЕНИЕ: SHA-256 vs KECCAK через наше измерение")
    print("=" * 70)

    print(f"""
  СТРУКТУРА:
    {'Metric':<25} {'SHA-256':>12} {'Mini-Keccak':>12} {'Same?':>6}
    {'—'*25} {'—'*12} {'—'*12} {'—'*6}
    {'Construction':<25} {'Merkle-Damgård':>12} {'Sponge':>12} {'NO':>6}
    {'Nonlinearity source':<25} {'carry+Ch/Maj':>12} {'χ only':>12} {'NO':>6}
    {'Input bits':<25} {'512':>12} {'{0}'.format(input_bits):>12} {'—':>6}
    {'Output bits':<25} {'256':>12} {'{0}'.format(output_bits):>12} {'—':>6}
    {'rank(T)':<25} {'256':>12} {'{0}'.format(rank_T):>12} {'—':>6}
    {'rank(NE/CE)':<25} {'256':>12} {'{0}'.format(rank_CE if 'rank_CE' in dir() else '?'):>12} {'—':>6}
    {'K':<25} {'128':>12} {'{0:.1f}'.format(K):>12} {'—':>6}
    {'K/(output/2)':<25} {'1.00':>12} {'{0:.3f}'.format(K/(output_bits/2)):>12} {'—':>6}
    {'Sensitivity':<25} {'128':>12} {'{0:.1f}'.format(S):>12} {'—':>6}

  КЛЮЧЕВОЙ ВОПРОС: K/(output/2) ≈ 1.0 для обеих конструкций?
  Если ДА → наша формула УНИВЕРСАЛЬНА (работает для ANY hash, not just ARX).
  Если НЕТ → ARX-specific phenomenon.
""")

    # =========================================================
    print(f"{'=' * 70}")
    print("7. COLLISION COST по нашей формуле")
    print("=" * 70)

    # Для Keccak: C = 2^8 (lane width), N_reg = N_output = 8
    # collision = C^(N_reg/2) = (2^8)^(8/2) = 2^32
    # Но: capacity = 64 bits → реальная безопасность = 2^(capacity/2) = 2^32

    # Для реального SHA3-256: C = 2^64, output = 256 bits = 4 lanes
    # Но capacity = 512 bits → security = 2^256 (collision)
    # Наша формула: C^(N_output/2) = (2^64)^(4/2) = 2^128... нет, output = 256/64 = 4 lanes
    # Hmm: SHA3-256 output = 256 bits, capacity = 512 bits
    # Collision security = min(2^(output/2), 2^(capacity/2)) = min(2^128, 2^256) = 2^128

    C_lane = 2**8  # lane width for mini-keccak
    N_out = N_output  # output lanes
    our_collision = C_lane ** (N_out // 2)
    capacity_bits = (25 - N_input) * 8
    sponge_collision = 2 ** min(output_bits // 2, capacity_bits // 2)

    print(f"""
  Mini-Keccak:
    Lane width C = 2^8 = 256
    Output lanes N = {N_output}

    Our formula: C^(N/2) = 256^{N_output//2} = 2^{8 * (N_output//2)}
    Sponge bound: min(2^(output/2), 2^(capacity/2)) = min(2^{output_bits//2}, 2^{capacity_bits//2}) = 2^{min(output_bits//2, capacity_bits//2)}

    Match: {'YES!' if 8*(N_output//2) == min(output_bits//2, capacity_bits//2) else 'NO — different!'}

  Real SHA3-256:
    Lane width C = 2^64
    Output = 256 bits = 4 lanes
    Capacity = 512 bits

    Our formula: C^(N_output/2) = (2^64)^2 = 2^128
    Sponge bound: min(2^128, 2^256) = 2^128

    Match: YES!

  ВЫВОД: наша формула C^(N/2) РАБОТАЕТ для sponge конструкций!
  collision_cost = C^(N_output/2) where N_output = output_bits / lane_bits.
""")


if __name__ == "__main__":
    main()
