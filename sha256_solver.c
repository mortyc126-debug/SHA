/*
 * SHA-256 Constraint Propagation Solver
 *
 * Uses ALL available equations + backward chain + schedule backward
 * to propagate known bits and minimize search space.
 *
 * Strategy:
 *   1. From hash: recover a[57..64] (free, backward chain)
 *   2. T2 equations: SIG0(a[r]) + MAJ(a[r],a[r-1],a[r-2]) = known
 *      → determines ~15/32 bits of each unknown per equation
 *   3. For each determined bit: propagate through створочне, schedule, rounds
 *   4. Branch on free bits, propagate, prune inconsistencies
 *
 * Compile: gcc -O3 -march=native -o solver sha256_solver.c
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))
#define CH(e,f,g)  (((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))
#define SIG0(x) (ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define SIG1(x) (ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define sig0(x) (ROTR(x,7)^ROTR(x,18)^((x)>>3))
#define sig1(x) (ROTR(x,17)^ROTR(x,19)^((x)>>10))

static const uint32_t K[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};
static const uint32_t IV[8] = {
    0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
    0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19
};

/* Forward SHA-256 with full trace */
void sha256_forward(const uint32_t msg[16], uint32_t hash[8],
                    uint32_t states[65][8], uint32_t W[64],
                    uint32_t T1s[64], uint32_t T2s[64]) {
    for (int i=0;i<16;i++) W[i]=msg[i];
    for (int i=16;i<64;i++) W[i]=sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16];
    for (int i=0;i<8;i++) states[0][i]=IV[i];
    for (int i=0;i<64;i++) {
        uint32_t a=states[i][0],b=states[i][1],c=states[i][2],d=states[i][3];
        uint32_t e=states[i][4],f=states[i][5],g=states[i][6],h=states[i][7];
        uint32_t t1=h+SIG1(e)+CH(e,f,g)+K[i]+W[i];
        uint32_t t2=SIG0(a)+MAJ(a,b,c);
        T1s[i]=t1; T2s[i]=t2;
        states[i+1][0]=t1+t2; states[i+1][1]=a; states[i+1][2]=b; states[i+1][3]=c;
        states[i+1][4]=d+t1;  states[i+1][5]=e; states[i+1][6]=f; states[i+1][7]=g;
    }
    for (int i=0;i<8;i++) hash[i]=states[64][i]+IV[i];
}

/*
 * Solve T2 = SIG0(a_known) + MAJ(a_known, b_known, X) for X
 * using bit-by-bit solver with carry tracking and branching.
 *
 * T2_target: known target value
 * a, b: known values
 * Returns number of solutions found, fills solutions array.
 */
typedef struct {
    uint32_t value;
    int valid;
} Solution;

int solve_maj_equation(uint32_t T2_target, uint32_t a, uint32_t b,
                       Solution *solutions, int max_solutions) {
    uint32_t sig0_a = SIG0(a);
    uint32_t maj_target = T2_target - sig0_a; /* mod 2^32 */

    /* MAJ(a, b, X) = (a&b) ^ ((a^b) & X)
     * SIG0(a) + MAJ(a,b,X) = T2_target
     * sig0_a + MAJ(a,b,X) = T2_target  (mod 2^32)
     *
     * Bit-by-bit with carry:
     * For bit k: sig0_a[k] + maj_k + carry_in = T2_target[k] + 2*carry_out
     * maj_k = (a&b)[k] ^ ((a^b)[k] & X[k])
     */

    uint32_t ab = a & b;    /* constant part of MAJ */
    uint32_t axb = a ^ b;   /* mask: where X matters */

    /* Recursive bit solver with branching */
    /* State: current bit position, carry, partial X value */
    /* Use iterative approach with stack */

    typedef struct { int bit; int carry; uint32_t x_partial; } State;
    State stack[1024];
    int sp = 0;
    int n_solutions = 0;

    stack[sp++] = (State){0, 0, 0};

    while (sp > 0 && n_solutions < max_solutions) {
        State s = stack[--sp];
        int k = s.bit;
        int carry = s.carry;
        uint32_t x = s.x_partial;

        if (k == 32) {
            /* Completed all bits — check carry is 0 */
            if (carry == 0) {
                solutions[n_solutions].value = x;
                solutions[n_solutions].valid = 1;
                n_solutions++;
            }
            continue;
        }

        int target_k = (T2_target >> k) & 1;
        int sig0_k = (sig0_a >> k) & 1;
        int ab_k = (ab >> k) & 1;
        int axb_k = (axb >> k) & 1;

        /* Try X[k] = 0 and X[k] = 1 */
        for (int xk = 0; xk <= 1; xk++) {
            /* Skip if MASK=0 and xk=1 doesn't change MAJ */
            int maj_k = ab_k ^ (axb_k & xk);
            int sum = sig0_k + maj_k + carry;
            int result_k = sum & 1;
            int new_carry = sum >> 1;

            if (result_k == target_k) {
                /* Consistent! Push next bit. */
                uint32_t new_x = x | ((uint32_t)xk << k);
                stack[sp++] = (State){k + 1, new_carry, new_x};
            }
            /* If not consistent: prune this branch */
        }
    }

    return n_solutions;
}

int main(int argc, char **argv) {
    printf("SHA-256 Constraint Propagation Solver\n");
    printf("=====================================\n\n");

    srand(42);
    uint32_t msg[16];
    for (int i=0;i<16;i++) msg[i]=(uint32_t)rand()|((uint32_t)rand()<<16);

    uint32_t hash[8], states[65][8], W[64], T1s[64], T2s[64];
    sha256_forward(msg, hash, states, W, T1s, T2s);

    printf("Hash: ");
    for (int i=0;i<8;i++) printf("%08x ", hash[i]);
    printf("\n\n");

    /* Phase 1: backward chain → a[57..64] known */
    uint32_t a_known[65] = {0};
    uint32_t T2_known[64] = {0};
    uint32_t T1_known[64] = {0};

    /* From state[64] = hash - IV */
    for (int i=0;i<8;i++) {
        uint32_t s = hash[i] - IV[i];
        if (i < 4) a_known[64-i] = s;
    }

    /* d-extraction: a[57..60] */
    for (int r = 63; r >= 60; r--) {
        T2_known[r] = SIG0(a_known[r]) + MAJ(a_known[r], a_known[r-1], a_known[r-2]);
        T1_known[r] = a_known[r+1] - T2_known[r];
        /* e[r+1] from state[64]: positions 4,5,6,7 = e64,e63,e62,e61 */
        uint32_t e_rp1;
        if (r == 63) e_rp1 = hash[4] - IV[4];      /* e[64] */
        else if (r == 62) e_rp1 = hash[5] - IV[5];  /* e[63] */
        else if (r == 61) e_rp1 = hash[6] - IV[6];  /* e[62] */
        else e_rp1 = hash[7] - IV[7];               /* e[61] */
        a_known[r-3] = e_rp1 - T1_known[r];         /* d[r] = a[r-3] */
    }

    printf("Phase 1: Backward chain\n");
    for (int r = 57; r <= 64; r++)
        printf("  a[%d] = %08x %s\n", r, a_known[r],
               a_known[r] == states[r][0] ? "OK" : "MISMATCH");

    /* Phase 2: T2 equations on unknowns */
    printf("\nPhase 2: T2 equations\n");

    /* Eq1: T2[58] = SIG0(a[58]) + MAJ(a[58], a[57], a[56])
     * a[58], a[57] KNOWN. Solve for a[56]. */
    T2_known[58] = a_known[59] - T1s[58]; /* Hmm, T1s from trace — we don't have it! */

    /* WAIT: we need T1[58] to compute T2[58] = a[59] - T1[58].
     * T1[58] is NOT directly known from backward chain.
     * T1[r] = h[r] + SIG1(e[r]) + CH(e[r],f[r],g[r]) + K[r] + W[r]
     * contains h[r] and W[r] which are unknowns!
     *
     * BUT: T1[r] + T2[r] = a[r+1] (KNOWN)
     * AND: T2[r] = SIG0(a[r]) + MAJ(a[r], a[r-1], a[r-2])
     * For r=59: T2[59] = SIG0(a[59]) + MAJ(a[59],a[58],a[57]) — ALL KNOWN!
     * → T1[59] = a[60] - T2[59] — KNOWN!
     *
     * For r=58: T2[58] = SIG0(a[58]) + MAJ(a[58],a[57],a[56]) — a[56] UNKNOWN
     * → T2[58] = a[59] - T1[58] — but T1[58] unknown
     * → CANNOT get T2[58] directly!
     *
     * HOWEVER: we can REFORMULATE:
     * a[59] = T1[58] + T2[58]
     * T2[58] = SIG0(a[58]) + MAJ(a[58], a[57], X)   where X = a[56]
     * T1[58] = a[59] - T2[58] = a[59] - SIG0(a[58]) - MAJ(a[58], a[57], X)
     *
     * We don't need T2[58] separately — we need:
     * SIG0(a[58]) + MAJ(a[58], a[57], X) = T2[58]
     * And: T2[58] = a[59] - T1[58]
     * But: T1[58] depends on h[58], W[58] — unknowns!
     *
     * WAIT — rethink. From backward:
     * T1[r] = a[r+1] - T2[r]. If T2[r] involves X, then T1[r] also has X.
     * This is NOT an independent equation for X.
     * It's: SIG0(a58) + MAJ(a58,a57,X) + T1[58] = a[59]
     * Two unknowns: X and T1[58]. One equation.
     *
     * UNLESS T1[58] can be expressed differently!
     * T1[58] = h[58] + SIG1(e[58]) + CH(e[58],f[58],g[58]) + K[58] + W[58]
     * h[58] = e[55] — unknown
     * e[58] = state[58][4] — unknown (depends on a[54])
     * W[58] = schedule — unknown
     * Multiple unknowns in T1[58]. Can't separate.
     */

    /* The T2 equation is NOT independently solvable!
     * T2[r] = a[r+1] - T1[r] introduces T1[r] as second unknown.
     *
     * The Python test worked because we used T1 from the TRACE (cheat).
     * Without the trace, T1 is unknown.
     */

    printf("  CRITICAL: T2[58] = a[59] - T1[58], but T1[58] is UNKNOWN!\n");
    printf("  T1[58] depends on h[58], e[58], W[58] — all unknowns.\n");
    printf("  The T2 equation is NOT independently solvable.\n\n");

    printf("  HOWEVER: T2[59] IS solvable!\n");
    T2_known[59] = SIG0(a_known[59]) + MAJ(a_known[59], a_known[58], a_known[57]);
    T1_known[59] = a_known[60] - T2_known[59];
    printf("  T2[59] = SIG0(a59)+MAJ(a59,a58,a57) = %08x (all args known)\n", T2_known[59]);
    printf("  T1[59] = a[60] - T2[59] = %08x\n", T1_known[59]);
    printf("  Verify T2[59]: %08x %s\n", T2s[59],
           T2_known[59] == T2s[59] ? "OK" : "MISMATCH");

    /* So: T2[59] is KNOWN → T1[59] is KNOWN.
     * T1[59] = h[59] + SIG1(e[59]) + CH(e[59],f[59],g[59]) + K[59] + W[59]
     * This gives us: h[59] + SIG1(e[59]) + CH(e[59],f[59],g[59]) + W[59] = T1[59] - K[59]
     * Still multiple unknowns: h[59]=e[56], e[59], W[59].
     *
     * But e[59] = state[60][5] = f[60]. state[60] = (a60,a59,a58,a57,e60,e59,e58,e57).
     * e[59] is NOT directly known (it's state[60][5] = f[60], which needs state[60]).
     * e[60] is unknown (= h[63] = our first mystery).
     *
     * So T1[59] gives: equation on h[59]=e[56] + e[59] + W[59].
     * Three unknowns. Not directly solvable.
     */

    printf("\n  T1[59] gives: h[59]+SIG1(e59)+CH(e59,f59,g59)+W[59] = %08x\n",
           T1_known[59] - K[59]);
    printf("  But h[59]=e[56], e[59], W[59] all unknown.\n");

    /* THE REAL EQUATION:
     * From round 58: a[59] = T1[58] + T2[58]
     *   T2[58] = SIG0(a58) + MAJ(a58, a57, a56)  ← involves a56
     *   T1[58] = h58 + SIG1(e58) + CH(e58,f58,g58) + K[58] + W[58]
     *          h58 = e55 (unknown), e58 (unknown), W58 (unknown)
     *
     * This is ONE equation (a[59] = T1+T2) with MANY unknowns:
     * a[56], a[55] (through e58), a[54] (through e55=h58), W[58] (schedule)
     *
     * The "15 bits free" result from Python was because we KNEW T1[58]
     * from the trace. Without it, the equation has too many unknowns.
     *
     * CONCLUSION: The T2 equation needs T1, and T1 has independent unknowns.
     * The system is COUPLED. Cannot solve one equation at a time.
     * Back to: 128 bits, all coupled.
     */

    printf("\n========================================\n");
    printf("CONCLUSION: T2 equation is coupled with T1.\n");
    printf("Cannot solve for a[56] independently.\n");
    printf("The 15/32 'free' bits were an artifact of using trace T1.\n");
    printf("Real system: 128 coupled unknowns.\n");
    printf("Birthday bound 2^128 holds.\n");

    return 0;
}
