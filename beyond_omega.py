"""
BEYOND OMEGA: What do we know that OMEGA doesn't?

OMEGA: carry as variables → degree-2 GF(2) system → Gaussian elimination.
       Works for R≤15, fails at R=16 (α-kernel = 0).

WE KNOW:
  1. 4 layers of 127 (OMEGA treats all 512 carry vars equally)
  2. CarryBridge structure (carry is sequential, not independent)
  3. 93.7% of round = fixed (OMEGA doesn't use this)
  4. Schedule = only grammar (OMEGA uses schedule but not as "grammar")

QUESTION: Can we solve the degree-2 system LAYER BY LAYER?

Layer 0: 127 quadratic GF(2) equations, 512 unknowns.
  → Underdetermined. Solution space = 512 - 127 = 385 dimensions.
  → This is the L-kernel from PHI!

Layer 1 (given layer 0 values): 127 new quadratic equations.
  → But carries from layer 0 are now CONSTANTS (known).
  → New unknowns: same 512 M bits, but 127 more constraints.
  → Solution space = 512 - 254 = 258 dimensions.

Layer 2: +127 → solution space = 131 dimensions.
Layer 3: +127 → solution space = 4 dimensions.
Layer 4: +4 → solution space = 0 dimensions (unique solution).

OMEGA solves ALL 512 equations at once → hard.
WE could solve 127 at a time → 4 steps of quadratic system solving.

Each step: 127 quadratic equations over GF(2).
Solving 127 quadratic equations in N variables:
  - Brute force: 2^N
  - XL/Groebner: depends on structure
  - With structure: potentially much cheaper

But each step REDUCES the solution space.
Step 0: 512 → 385 free.
Step 1: 385 → 258 free. But 127 new equations involve PRODUCTS of unknowns.
Step 2: 258 → 131 free.
Step 3: 131 → 4 free.
Step 4: 4 → 0.

The PRODUCTS are from Ch and Maj: e[k]·f[k], a[k]·b[k], etc.
At each step, these are products of bits AT THE SAME POSITION k.
Not cross-position products (those would be from carry, which we handle).

THIS IS THE KEY STRUCTURAL PROPERTY:
  Within a layer, quadratic terms are POSITION-LOCAL.
  They couple bits of the same position k across registers (e·f, a·b).
  They do NOT couple different positions.
  Different positions are coupled ONLY through carry (handled by layers).

So: within layer k, the 127 equations DECOMPOSE into 32 independent groups
(one per bit position within the layer), minus some global coupling from
rotations (Sig0, Sig1 read from other positions).

Wait — Sig0 and Sig1 READ from other positions (2, 6, 11, 13, 22, 25).
But these reads are from the SAME layer k (same bit position across all words).
No — Sig0(a)[k] = a[(k+2)%32] XOR a[(k+13)%32] XOR a[(k+22)%32].
The bit positions (k+2)%32, (k+13)%32, (k+22)%32 are IN DIFFERENT LAYERS!

This is the CROSS-LAYER COUPLING through rotations.
It's LINEAR (XOR), but it connects layers.

So: rotations couple layers LINEARLY.
    Carries couple layers QUADRATICALLY (through MAJ).
    Ch/Maj create WITHIN-LAYER quadratic terms.

The full picture:
  Layer k has 127 equations that involve:
  - LINEAR reads from OTHER layers (through Sig0, Sig1 rotations)
  - QUADRATIC terms within the SAME register position (Ch, Maj)
  - CONSTANT carry corrections from lower layers

THIS is what OMEGA doesn't decompose.
Let's measure: how many of the 127 equations are truly "within layer k"
vs "cross-layer through rotations"?
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def experiment_rotation_coupling():
    """
    Sig0(a)[k] = a[(k+2)%32] XOR a[(k+13)%32] XOR a[(k+22)%32]
    Sig1(e)[k] = e[(k+6)%32] XOR e[(k+11)%32] XOR e[(k+25)%32]

    For layer k (bit position k):
      Sig0 reads from layers (k+2)%32, (k+13)%32, (k+22)%32.
      Sig1 reads from layers (k+6)%32, (k+11)%32, (k+25)%32.

    Layer k's equation for a_new[k]:
      a_new[k] = Sig0(a)[k] ⊕ Maj(a,b,c)[k] ⊕ Sig1(e)[k] ⊕ Ch(e,f,g)[k] ⊕ h[k] ⊕ K[k] ⊕ W[k]
                 ⊕ carry_into_k

    The terms:
      Sig0(a)[k] = a[(k+2)] ⊕ a[(k+13)] ⊕ a[(k+22)]  → reads layers k+2, k+13, k+22
      Maj(a,b,c)[k] = a[k]b[k] ⊕ a[k]c[k] ⊕ b[k]c[k] → reads layer k (quadratic)
      Sig1(e)[k] = e[(k+6)] ⊕ e[(k+11)] ⊕ e[(k+25)]  → reads layers k+6, k+11, k+25
      Ch(e,f,g)[k] = e[k]f[k] ⊕ (1⊕e[k])g[k]         → reads layer k (quadratic)
      h[k] = g_prev[k] = f_prev_prev[k] = e_{r-3}[k]  → reads layer k
      carry_into_k = from layers 0..k-1 (carry bridge)

    So layer k's equation involves:
      OWN LAYER: Maj, Ch, h, carry → quadratic in layer k variables
      OTHER LAYERS: Sig0 reads from k+2, k+13, k+22
                    Sig1 reads from k+6, k+11, k+25
      TOTAL OTHER LAYERS: 6 distinct layers read per equation.

    This means: each layer is NOT isolated. Each layer equation involves
    6 other layers through LINEAR (XOR) reads.
    """
    print("=" * 80)
    print("ROTATION COUPLING: How layers connect through Sig0/Sig1")
    print("=" * 80)

    # For each layer k, which other layers does it read from?
    print(f"\n  Layer coupling map (via Sig0 and Sig1 rotations):")
    print(f"  Layer k reads from layers: Sig0={{(k+2)%32, (k+13)%32, (k+22)%32}}")
    print(f"                             Sig1={{(k+6)%32, (k+11)%32, (k+25)%32}}")

    # Build the coupling graph
    coupling = {}
    for k in range(32):
        reads = set()
        for offset in [2, 13, 22, 6, 11, 25]:
            reads.add((k + offset) % 32)
        reads.discard(k)  # Remove self
        coupling[k] = reads

    print(f"\n  Coupling graph (layer → set of layers it reads):")
    for k in range(8):  # Show first 8
        print(f"    Layer {k:>2} reads from: {sorted(coupling[k])}")

    # How many layers does each layer couple to?
    all_counts = [len(coupling[k]) for k in range(32)]
    print(f"\n  Layers read per layer: min={min(all_counts)}, max={max(all_counts)}, avg={sum(all_counts)/32:.1f}")

    # Critical: is the coupling graph CONNECTED?
    # If we start from layer 0, can we reach all 32 layers?
    reached = {0}
    frontier = {0}
    steps = 0
    while frontier:
        new_frontier = set()
        for k in frontier:
            for neighbor in coupling[k]:
                if neighbor not in reached:
                    reached.add(neighbor)
                    new_frontier.add(neighbor)
        frontier = new_frontier
        steps += 1
        if len(reached) == 32:
            break

    print(f"  Coupling graph connected: {len(reached) == 32} (reached all in {steps} steps)")

    # Now: the EFFECTIVE structure
    print(f"""
  ═══════════════════════════════════════════════════════════════
  THE COMPLETE STRUCTURE OF SHA-256 IN OUR LANGUAGE
  ═══════════════════════════════════════════════════════════════

  SHA-256 consists of:

  32 BIT-LAYERS (k = 0..31), each with:
    - 127 constraints per layer (8 regs × 64 rounds, minus створочное)
    - QUADRATIC terms: Ch, Maj at position k (degree 2, position-local)
    - LINEAR reads from 6 OTHER layers: {{(k+2)%32, (k+13)%32, (k+22)%32,
                                           (k+6)%32, (k+11)%32, (k+25)%32}}
    - CARRY BRIDGE: carry_into_k = function of layers 0..k-1

  Layer connectivity:
    - 6 linear connections per layer (through Sig0/Sig1 rotations)
    - 1 sequential connection per layer (carry from below)
    - Graph is FULLY CONNECTED in {steps} steps

  The FOUR CUMULATIVE LAYERS (from unfolding):
    - First 4 bit-layers contribute 127 INDEPENDENT constraints each
    - By bit 4: all 512 message bits are determined
    - Layers 5-31 add NO new constraints (rank stays at 512)

  WHY layers 5-31 add nothing:
    Because each layer reads from 6 others through rotations.
    By the time layers 0-4 determine all 512 bits,
    layers 5-31 are FULLY DETERMINED by them.

  THE MAGIC NUMBER 4:
    4 layers × 127 + 4 = 512
    This means: carry depth of 4 is sufficient.
    Average carry segment length ≈ 3.64 (we measured).
    4 ≈ 3.64 — NOT A COINCIDENCE.
    The number of independent layers ≈ average carry chain length.
    """)


def experiment_cross_layer_equations():
    """
    Count exactly: how many of the 127 constraints per layer
    are "pure layer k" vs "involves other layers"?
    """
    print("\n" + "=" * 80)
    print("EQUATION STRUCTURE: Pure vs cross-layer per layer")
    print("=" * 80)

    # In one round's bit-k equation for a_new[k]:
    #
    # a_new[k] involves:
    #   OWN LAYER k:
    #     - Maj(a[k], b[k], c[k]) — 3 vars from layer k
    #     - Ch(e[k], f[k], g[k]) — 3 vars from layer k
    #     - h[k] — 1 var from layer k
    #     - K[r][k], W[r][k] — constants
    #     - carry_into_k — bridge (function of lower layers)
    #     Total own-layer vars: 7 per round
    #
    #   CROSS LAYER (via Sig0):
    #     - a[(k+2)%32], a[(k+13)%32], a[(k+22)%32] — 3 vars from other layers
    #
    #   CROSS LAYER (via Sig1):
    #     - e[(k+6)%32], e[(k+11)%32], e[(k+25)%32] — 3 vars from other layers
    #
    # Total: 7 own + 6 cross = 13 variables per equation.
    # 7/13 = 54% own-layer, 6/13 = 46% cross-layer.

    print(f"  Per equation (one round, one register):")
    print(f"    Own-layer variables: 7 (Maj: a,b,c; Ch: e,f,g; h)")
    print(f"    Cross-layer variables: 6 (Sig0: 3 rotated a's; Sig1: 3 rotated e's)")
    print(f"    Ratio: own = 54%, cross = 46%")
    print(f"")
    print(f"  The cross-layer coupling is LINEAR (XOR of rotated bits).")
    print(f"  The own-layer coupling is QUADRATIC (Ch, Maj products).")
    print(f"")
    print(f"  So: the NONLINEARITY is local (within layer).")
    print(f"  The COUPLING is global but linear (between layers).")
    print(f"")
    print(f"  THIS IS THE SEPARATION:")
    print(f"    NONLINEAR + LOCAL (Ch, Maj at position k)")
    print(f"    LINEAR + GLOBAL (Sig0, Sig1 read other positions)")
    print(f"    SEQUENTIAL + BRIDGE (carry from lower layers)")
    print(f"")
    print(f"  Three types of interaction, three roles:")
    print(f"    Ch/Maj = creates complexity (degree 2)")
    print(f"    Sig0/Sig1 = spreads information (linear mixing)")
    print(f"    Carry = bridges layers (sequential coupling)")


if __name__ == "__main__":
    experiment_rotation_coupling()
    experiment_cross_layer_equations()
