/*
 * CARRY×NLF: поиск обратного оператора.
 *
 * Оператор F: (state, W) → (a_new, e_new)
 *   a_new = (h + Σ₁(e) + Ch(e,f,g) + K + W) + (Σ₀(a) + Maj(a,b,c))
 *   e_new = d + (h + Σ₁(e) + Ch(e,f,g) + K + W)
 *
 * F = Linear_part + Carry×NLF_correction
 *
 * Вопрос: существует ли F⁻¹?
 *
 * Подходы:
 * 1. ПОБИТОВАЯ ИНВЕРСИЯ: bit 0 → bit 1 → ... (triangular)
 * 2. ALGEBRAIC INVERSION: Ch⊕Maj=g&~f identity
 * 3. FEEDBACK LOOP: 4-round period → algebraic equation
 * 4. EXHAUSTIVE TEST: for small n, check if F bijective
 *
 * gcc -O3 -march=native -o cnlf_inv cnlf_inverse.c
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))
static inline uint32_t Ch(uint32_t e, uint32_t f, uint32_t g) {
    return (e&f)^(~e&g);
}
static inline uint32_t Maj(uint32_t a, uint32_t b, uint32_t c) {
    return (a&b)^(a&c)^(b&c);
}

/* ══════════════════════════════════════════════════ */
/* TEST 1: Is carry×NLF BIJECTIVE for one addition?  */
/*                                                    */
/* f(x) = x + Ch(x, c1, c2) for fixed c1,c2.        */
/* Is f a PERMUTATION of Z/2^n?                       */
/* If yes: invertible! If no: information lost.       */
/* ══════════════════════════════════════════════════ */

void test_bijectivity_small() {
    printf("TEST 1: Is x → x + NLF(x, c1, c2) a BIJECTION?\n");
    printf("─────────────────────────────────────────────\n");

    /* For small n, exhaustive check */
    for(int n = 4; n <= 16; n += 4) {
        uint32_t size = 1u << n;
        uint32_t mask = size - 1;

        /* Try several (c1, c2) values */
        int bijective_count = 0;
        int total_tests = 0;

        srand(n * 42);
        for(int test = 0; test < 20; test++) {
            uint32_t c1 = rand() & mask;
            uint32_t c2 = rand() & mask;

            /* Check: f(x) = (x + Ch(x,c1,c2)) & mask. Is f bijective? */
            uint8_t *seen = calloc(size, 1);
            int is_bij = 1;
            for(uint32_t x = 0; x < size; x++) {
                uint32_t ch_val = ((x & c1) ^ (~x & c2)) & mask;
                uint32_t y = (x + ch_val) & mask;
                if(seen[y]) { is_bij = 0; break; }
                seen[y] = 1;
            }
            free(seen);
            if(is_bij) bijective_count++;
            total_tests++;
        }
        printf("  n=%2d: x→x+Ch(x,c1,c2) bijective in %d/%d cases\n",
               n, bijective_count, total_tests);

        /* Same for Maj */
        bijective_count = 0;
        for(int test = 0; test < 20; test++) {
            uint32_t c1 = rand() & mask;
            uint32_t c2 = rand() & mask;
            uint8_t *seen = calloc(size, 1);
            int is_bij = 1;
            for(uint32_t x = 0; x < size; x++) {
                uint32_t mj = ((x&c1)^(x&c2)^(c1&c2)) & mask;
                uint32_t y = (x + mj) & mask;
                if(seen[y]) { is_bij = 0; break; }
                seen[y] = 1;
            }
            free(seen);
            if(is_bij) bijective_count++;
            total_tests++;
        }
        printf("  n=%2d: x→x+Maj(x,c1,c2) bijective in %d/20 cases\n",
               n, bijective_count);
    }
}

/* ══════════════════════════════════════════════════ */
/* TEST 2: Is the FULL round function bijective?     */
/*                                                    */
/* Round: (a,b,c,d,e,f,g,h) → (a',b',c',d',e',f',g',h') */
/* Given W: is this a permutation of 256-bit state?  */
/* SHA-256 rounds ARE bijective (known). But is the  */
/* carry×NLF component SEPARATELY bijective?         */
/* ══════════════════════════════════════════════════ */

void test_round_bijectivity() {
    printf("\nTEST 2: Round function bijectivity components\n");
    printf("─────────────────────────────────────────────\n");

    /* The round function: known to be bijective given W.
     * Inverse: standard backward round.
     * But this doesn't mean carry×NLF is SEPARABLY invertible.
     *
     * Key question: can we invert carry×NLF WITHOUT knowing
     * the intermediate carry values?
     *
     * For one addition: f(x) = x + g(x) where g = Ch or Maj.
     * f is bijective iff no two x give same f(x).
     * g(x) depends on x QUADRATICALLY → f = linear + quadratic.
     * For linear part: always bijective (x + const).
     * Adding quadratic: may or may not be bijective. */

    printf("  SHA-256 round: BIJECTIVE (proven, standard inverse exists).\n");
    printf("  Inverse needs: state[r+1] and W[r].\n");
    printf("  WITHOUT W[r]: state[r] has 32 unknown bits (h[r]).\n");
    printf("  The 'lock' is not bijectivity but h[r]+W[r] coupling.\n");
}

/* ══════════════════════════════════════════════════ */
/* TEST 3: ПОБИТОВАЯ ИНВЕРСИЯ carry×NLF.              */
/*                                                    */
/* Given: y = x + Ch(x, c1, c2) (mod 2^n), find x.  */
/* Bit 0: y[0] = x[0] ⊕ Ch(x,c1,c2)[0] ⊕ carry[0] */
/*        carry[0] = 0.                               */
/*        Ch(x,c1,c2)[0] = x[0]&c1[0] ⊕ ~x[0]&c2[0]*/
/*                        = x[0]?(c1[0]):(c2[0])      */
/*        y[0] = x[0] ⊕ (x[0]?c1[0]:c2[0])          */
/*                                                    */
/* If c1[0]=c2[0]: y[0] = x[0] ⊕ c1[0] → x[0] = y[0] ⊕ c1[0]  */
/* If c1[0]≠c2[0]: y[0] = x[0] ⊕ x[0] = 0 or 1... need case split. */
/* ══════════════════════════════════════════════════ */

void test_bitwise_inversion() {
    printf("\nTEST 3: BITWISE INVERSION of x + Ch(x, c1, c2)\n");
    printf("─────────────────────────────────────────────\n");

    /* For n=32: given y, c1, c2, find x such that (x + Ch(x,c1,c2)) = y.
     *
     * Bit-by-bit from bit 0:
     *   Ch(x,c1,c2)[k] = x[k]?c1[k]:c2[k] = x[k]&(c1[k]⊕c2[k]) ⊕ c2[k]
     *   Let m[k] = c1[k]⊕c2[k] (mask: where x matters in Ch)
     *       v[k] = c2[k] (value when x[k]=0)
     *   Ch[k] = m[k]&x[k] ⊕ v[k]
     *
     *   y[k] = x[k] ⊕ Ch[k] ⊕ carry[k]
     *        = x[k] ⊕ (m[k]&x[k] ⊕ v[k]) ⊕ carry[k]
     *
     *   If m[k]=0: y[k] = x[k] ⊕ v[k] ⊕ carry[k]
     *              → x[k] = y[k] ⊕ v[k] ⊕ carry[k]  (UNIQUE!)
     *
     *   If m[k]=1: y[k] = x[k] ⊕ x[k] ⊕ v[k] ⊕ carry[k]
     *              Wait: x[k] ⊕ (1&x[k]) = x[k] ⊕ x[k] = 0? No!
     *              x[k] XOR (m[k] AND x[k]) = x[k] XOR x[k] = 0 when m=1.
     *              So: y[k] = 0 ⊕ v[k] ⊕ carry[k] = v[k] ⊕ carry[k]
     *              This is INDEPENDENT of x[k]!
     *              → x[k] is FREE (any value gives same y[k])
     *              → NOT INVERTIBLE at this bit!
     *
     * NUMBER OF FREE BITS = HW(m) = HW(c1 ⊕ c2).
     * Expected: HW ≈ 16 out of 32 → ~16 bits free per addition!
     */

    int N = 100000;
    double avg_free = 0;
    srand(42);
    for(int trial = 0; trial < N; trial++) {
        uint32_t c1 = rand() | (rand()<<16);
        uint32_t c2 = rand() | (rand()<<16);
        uint32_t m = c1 ^ c2;
        avg_free += __builtin_popcount(m);
    }
    avg_free /= N;

    printf("  For y = x + Ch(x, c1, c2):\n");
    printf("  At bit k where c1[k]=c2[k]: x[k] DETERMINED (unique inverse)\n");
    printf("  At bit k where c1[k]≠c2[k]: x[k] FREE (not in equation!)\n");
    printf("  Average free bits: %.1f/32 = %.0f%%\n", avg_free, avg_free/32*100);
    printf("\n  BUT: carry at bit k+1 DEPENDS on x[k]!\n");
    printf("  carry[k+1] = MAJ(x[k], Ch[k], carry[k])\n");
    printf("  When m[k]=1: Ch[k] = v[k] (known), but x[k] is free.\n");
    printf("  carry[k+1] = MAJ(x[k], v[k], carry[k])\n");
    printf("  If v[k]=carry[k]: carry[k+1]=v[k] regardless of x[k] → carry STABLE.\n");
    printf("  If v[k]≠carry[k]: carry[k+1]=x[k] → carry DEPENDS on free bit!\n");
    printf("  → The free bit PROPAGATES through carry to higher bits.\n\n");

    /* COUNT: how many free bits actually affect carry? */
    double avg_affecting_carry = 0;
    for(int trial = 0; trial < 10000; trial++) {
        uint32_t c1 = rand() | (rand()<<16);
        uint32_t c2 = rand() | (rand()<<16);
        uint32_t y = rand() | (rand()<<16);
        uint32_t m = c1 ^ c2;
        uint32_t v = c2; /* Ch when x[k]=0 */

        /* Simulate bitwise: find where free bits affect carry */
        int carry = 0;
        int affecting = 0;
        for(int k = 0; k < 32; k++) {
            int mk = (m >> k) & 1;
            int vk = (v >> k) & 1;
            if(mk) {
                /* x[k] free. Does it affect carry[k+1]? */
                if(vk != carry) {
                    /* carry[k+1] = x[k] → depends on free bit! */
                    affecting++;
                    carry = 0; /* unknown: could be 0 or 1 */
                } else {
                    carry = vk; /* stable regardless of x[k] */
                }
            } else {
                /* x[k] determined: x[k] = y[k] ⊕ v[k] ⊕ carry */
                int xk = ((y >> k) & 1) ^ vk ^ carry;
                int chk = vk; /* Ch[k] = v[k] when m[k]=0 */
                /* Actually Ch[k] = m[k]&x[k] ⊕ v[k] = 0&x[k] ⊕ v[k] = v[k] */
                /* carry[k+1] = MAJ(x[k], ch[k]+carry components...) */
                /* Simplified: operand1=x[k], operand2=ch[k] */
                int sum_k = xk + chk + carry;
                carry = sum_k >> 1;
            }
        }
        avg_affecting_carry += affecting;
    }
    avg_affecting_carry /= 10000;
    printf("  Free bits that PROPAGATE through carry: %.1f/%.0f\n",
           avg_affecting_carry, avg_free);
    printf("  (These are the bits that CREATE ambiguity in higher positions)\n");

    /* VERIFICATION: try to invert for real */
    printf("\n  INVERSION TEST: solve x from y = x + Ch(x, c1, c2)\n");
    int solved = 0, failed = 0, multiple = 0;
    for(int trial = 0; trial < 10000; trial++) {
        uint32_t c1 = rand() | (rand()<<16);
        uint32_t c2 = rand() | (rand()<<16);
        uint32_t x_actual = rand() | (rand()<<16);
        uint32_t ch_val = Ch(x_actual, c1, c2);
        uint32_t y = x_actual + ch_val;

        /* Try to recover x from y, c1, c2 */
        /* Bit-by-bit with branching on free bits */
        uint32_t m = c1 ^ c2;
        int n_free = __builtin_popcount(m);

        if(n_free > 20) { /* too many branches, skip */
            continue;
        }

        /* Enumerate all 2^n_free possibilities for free bits */
        int n_solutions = 0;
        uint32_t found_x = 0;

        for(uint32_t free_val = 0; free_val < (1u << n_free); free_val++) {
            uint32_t x_try = 0;
            int carry = 0;
            int free_idx = 0;
            int valid = 1;

            for(int k = 0; k < 32 && valid; k++) {
                int mk = (m >> k) & 1;
                int vk = (c2 >> k) & 1;
                int yk = (y >> k) & 1;
                int xk;

                if(mk) {
                    /* Free bit: take from free_val */
                    xk = (free_val >> free_idx) & 1;
                    free_idx++;
                } else {
                    /* Determined: xk = yk ⊕ vk ⊕ carry */
                    xk = yk ^ vk ^ carry;
                }
                x_try |= (xk << k);

                /* Compute carry for next bit */
                int chk = (mk & xk) ^ vk;
                int sum = xk + chk + carry;
                int result_k = sum & 1;
                carry = sum >> 1;

                /* Check: does result match y[k]? */
                if(result_k != yk) valid = 0;
            }

            if(valid && carry == 0) {
                /* Verify */
                uint32_t ch_check = Ch(x_try, c1, c2);
                if((x_try + ch_check) == y) {
                    n_solutions++;
                    found_x = x_try;
                }
            }
        }

        if(n_solutions == 1 && found_x == x_actual) solved++;
        else if(n_solutions > 1) multiple++;
        else failed++;
    }
    printf("  N=%d (skipped large n_free):\n", solved+multiple+failed);
    printf("    Unique solution (= correct x): %d\n", solved);
    printf("    Multiple solutions: %d\n", multiple);
    printf("    No solution found: %d\n", failed);

    if(multiple > 0) {
        printf("\n  ★ MULTIPLE SOLUTIONS exist for x + Ch(x,c1,c2) = y!\n");
        printf("  → Carry×NLF is NOT always bijective!\n");
        printf("  → Average solutions per case: %.2f\n",
               (double)(solved + multiple*2) / (solved + multiple + failed));
        printf("  → This IS the information loss mechanism.\n");
    } else if(solved > 0 && failed == 0) {
        printf("\n  ★ ALWAYS UNIQUE SOLUTION!\n");
        printf("  → Carry×NLF IS bijective!\n");
        printf("  → An inverse operator EXISTS!\n");
    }
}

int main() {
    printf("CARRY×NLF INVERSE: does it exist?\n");
    printf("==================================\n\n");

    test_bijectivity_small();
    test_round_bijectivity();
    test_bitwise_inversion();

    printf("\n══════════════════════════════════════════\n");
    printf("FINAL ANSWER\n");
    printf("══════════════════════════════════════════\n");

    return 0;
}
