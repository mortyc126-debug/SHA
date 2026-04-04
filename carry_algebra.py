"""
P3.3: CARRY OPERATOR ALGEBRA

Define: for fixed y ∈ {0,1}^n, the carry operator C_y: {0,1}^n → {0,1}^n
  C_y(x)[k] = carry into position k when computing x + y.

  C_y(x)[0] = 0 always
  C_y(x)[k] = MAJ(x[k-1], y[k-1], C_y(x)[k-1])

Properties to investigate:
  1. Is C_y linear? (NO — MAJ is degree 2)
  2. Is C_y ∘ C_y = C_y? (idempotent?)
  3. Is there k such that C_y^k = 0? (nilpotent?)
  4. What is the GF(2) rank of C_y's Jacobian?
  5. Does C_y have fixed points?
  6. C_x ∘ C_y vs C_y ∘ C_x — commutative?

These properties characterize the "carry algebra" and might
reveal structure that our linear/quadratic analysis missed.
"""

import random


def carry_operator(x, y, n):
    """Compute carry vector for x + y in n bits."""
    c = 0
    result = []
    for k in range(n):
        result.append(c)
        xk = (x >> k) & 1
        yk = (y >> k) & 1
        c = (xk & yk) | (xk & c) | (yk & c)
    return result  # list of n bits


def carry_as_int(x, y, n):
    """Carry vector as integer."""
    bits = carry_operator(x, y, n)
    return sum(b << k for k, b in enumerate(bits))


def experiment_basic_properties():
    """Test basic algebraic properties of C_y."""
    print("=" * 80)
    print("CARRY OPERATOR: Basic algebraic properties")
    print("=" * 80)

    n = 8
    N_total = 2**n

    # Property 1: Linearity
    # C_y(x1 XOR x2) = C_y(x1) XOR C_y(x2)?
    print(f"\n  n={n}")
    linear_violations = 0
    total_tests = 0
    for y in range(0, N_total, 17):  # Sample y values
        for x1 in range(N_total):
            for x2 in range(x1+1, min(x1+10, N_total)):
                c1 = carry_as_int(x1, y, n)
                c2 = carry_as_int(x2, y, n)
                c12 = carry_as_int(x1 ^ x2, y, n)
                if c12 != c1 ^ c2:
                    linear_violations += 1
                total_tests += 1
    print(f"  Linearity: {linear_violations}/{total_tests} violations → {'NONLINEAR' if linear_violations > 0 else 'LINEAR'}")

    # Property 2: Idempotent (C_y ∘ C_y = C_y)
    # C_y(C_y(x)) = C_y(x)?
    # Note: C_y maps {0,1}^n → {0,1}^n. Can compose.
    idempotent_count = 0
    total = 0
    for y in range(0, N_total, 7):
        for x in range(N_total):
            c1 = carry_as_int(x, y, n)
            c2 = carry_as_int(c1, y, n)
            if c1 == c2:
                idempotent_count += 1
            total += 1
    print(f"  Idempotent: C_y(C_y(x)) = C_y(x) for {idempotent_count}/{total} = {idempotent_count/total:.3f}")

    # Property 3: Nilpotent (∃k: C_y^k = 0)
    # Check: C_y^k(x) → 0?
    nilpotent_count = 0
    max_depth = 0
    total = 0
    for y in range(0, N_total, 7):
        for x in range(N_total):
            v = x
            for depth in range(n + 5):
                v = carry_as_int(v, y, n)
                if v == 0:
                    if depth + 1 > max_depth:
                        max_depth = depth + 1
                    nilpotent_count += 1
                    break
            total += 1
    print(f"  Nilpotent: C_y^k(x) → 0 for {nilpotent_count}/{total} = {nilpotent_count/total:.3f}")
    print(f"  Max depth to 0: {max_depth}")

    # Property 4: Fixed points (C_y(x) = x)
    fixed_count = 0
    total = 0
    for y in range(0, N_total, 7):
        n_fixed = 0
        for x in range(N_total):
            if carry_as_int(x, y, n) == x:
                n_fixed += 1
        fixed_count += n_fixed
        total += N_total
    avg_fixed = fixed_count / (N_total // 7 + 1)
    print(f"  Fixed points: avg {avg_fixed:.1f} per y (out of {N_total})")

    # Property 5: Image size (how many distinct C_y(x) values?)
    image_sizes = []
    for y in range(N_total):
        image = set()
        for x in range(N_total):
            image.add(carry_as_int(x, y, n))
        image_sizes.append(len(image))
    avg_image = sum(image_sizes) / len(image_sizes)
    min_image = min(image_sizes)
    max_image = max(image_sizes)
    print(f"  Image size: min={min_image}, max={max_image}, avg={avg_image:.1f} (out of {N_total})")

    # Property 6: Commutativity C_x(C_y(z)) vs C_y(C_x(z))
    commute_violations = 0
    total = 0
    for x in range(0, N_total, 11):
        for y in range(0, N_total, 13):
            for z in range(0, N_total, 7):
                cyx = carry_as_int(carry_as_int(z, y, n), x, n)
                cxy = carry_as_int(carry_as_int(z, x, n), y, n)
                if cyx != cxy:
                    commute_violations += 1
                total += 1
    print(f"  Commutativity: {commute_violations}/{total} violations → {'NON-COMMUTATIVE' if commute_violations > 0 else 'COMMUTATIVE'}")


def experiment_jacobian_rank():
    """GF(2) Jacobian rank of C_y for each y."""
    print("\n" + "=" * 80)
    print("CARRY OPERATOR: Jacobian rank")
    print("=" * 80)

    n = 8
    N_total = 2**n

    ranks = []
    for y in range(N_total):
        # Jacobian of C_y: for each input bit j, which output bits change?
        J = []
        base = carry_as_int(0, y, n)  # C_y(0)
        for j in range(n):
            x_flip = 1 << j
            c_flip = carry_as_int(x_flip, y, n)
            row = [(base >> k) & 1 ^ (c_flip >> k) & 1 for k in range(n)]
            J.append(row)

        # GF(2) rank
        m = [list(row) for row in J]
        rank = 0
        for col in range(n):
            pivot = None
            for row in range(rank, n):
                if m[row][col] == 1:
                    pivot = row
                    break
            if pivot is None:
                continue
            m[rank], m[pivot] = m[pivot], m[rank]
            for row in range(n):
                if row != rank and m[row][col] == 1:
                    for c in range(n):
                        m[row][c] ^= m[rank][c]
            rank += 1
        ranks.append(rank)

    from collections import Counter
    rank_dist = Counter(ranks)
    print(f"  n={n}: Jacobian rank distribution of C_y:")
    for r, count in sorted(rank_dist.items()):
        print(f"    rank={r}: {count}/{N_total} ({count/N_total:.1%})")


def experiment_carry_composition():
    """
    Key: how does carry compose with itself?

    In SHA-256, one round has ~6 additions. Each addition has its own carry.
    The COMPOSITION of carries across additions within one round:
      result1 = x + y (carry C1)
      result2 = result1 + z (carry C2, depends on result1 which depends on C1)

    So: C2 ∘ (x + y) + z — the carry of the second addition depends on
    the RESULT of the first, not its carry directly.

    What's the effective "composed carry" of two sequential additions?
    """
    print("\n" + "=" * 80)
    print("CARRY COMPOSITION: Sequential additions")
    print("=" * 80)

    n = 8
    mask = (1 << n) - 1

    # For x + y + z (two additions):
    # Step 1: s1 = x + y, carry1 = C(x,y)
    # Step 2: s2 = s1 + z, carry2 = C(s1, z) where s1 = (x+y) mod 2^n
    # Total carry effect on result: result = x XOR y XOR z XOR total_carry
    # total_carry = carry1 XOR carry2 XOR correction
    # Actually: result = (x + y + z) mod 2^n
    #          = x XOR y XOR z XOR carry_of_3way_add(x,y,z)

    # The 3-way carry: more complex than 2-way.
    # Is it the XOR of carry(x,y) and carry(x+y, z)? NO — because carry(x+y, z)
    # depends on x+y which includes carry(x,y).

    # Let's measure: carry_3way vs carry_2way compositions
    compositions_match = 0
    total = 0

    for x in range(0, 2**n, 3):
        for y in range(0, 2**n, 5):
            for z in range(0, 2**n, 7):
                # 3-way result
                r3 = (x + y + z) & mask
                xor3 = x ^ y ^ z
                carry3 = r3 ^ xor3  # total carry effect

                # 2-way composed
                s1 = (x + y) & mask
                c1 = carry_as_int(x, y, n)  # As integer
                c1_effect = s1 ^ (x ^ y)  # = carry contribution as bits

                c2 = carry_as_int(s1, z, n)
                c2_effect = r3 ^ (s1 ^ z)

                # Is carry3 = c1_effect XOR c2_effect?
                if carry3 == (c1_effect ^ c2_effect):
                    compositions_match += 1
                total += 1

    print(f"  carry(x+y+z) = carry(x,y) XOR carry(x+y, z): "
          f"{compositions_match}/{total} = {compositions_match/total:.4f}")

    if compositions_match == total:
        print(f"  ✓ CARRIES COMPOSE BY XOR!")
        print(f"  Total carry = XOR of individual carries. This is ADDITIVE.")
        print(f"  This means: carry is a COCYCLE in some cohomology!")
    else:
        print(f"  ✗ Carries do NOT compose by XOR.")


if __name__ == "__main__":
    experiment_basic_properties()
    experiment_jacobian_rank()
    experiment_carry_composition()
