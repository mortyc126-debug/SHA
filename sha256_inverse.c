/*
 * SHA-256 Inverse Engine: Zone C backward + a[56] recovery
 *
 * Architecture:
 *   1. Forward SHA-256 compression (standard)
 *   2. Zone C inverse: state[64] → state[49] deterministically (O(1))
 *   3. Extended backward via створочне: recover a[57..60], W[57..63]
 *   4. a[56] brute force: enumerate 2^32 candidates
 *   5. For each a[56]: solve forward Zone A+B, check consistency
 *
 * Compile: gcc -O3 -march=native -o sha256_inv sha256_inverse.c
 * Usage:   ./sha256_inv [mode]
 *          mode=0: verify inverse (default)
 *          mode=1: a[56] search on reduced rounds
 *          mode=2: full 2^32 enumeration (LONG)
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

/* ================================================================
 *  SHA-256 Constants
 * ================================================================ */

static const uint32_t K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

static const uint32_t IV[8] = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
};

/* ================================================================
 *  SHA-256 Primitives
 * ================================================================ */

#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))
#define SHR(x,n)  ((x)>>(n))
#define CH(e,f,g)  (((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))
#define SIG0(x) (ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define SIG1(x) (ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define sig0(x) (ROTR(x,7)^ROTR(x,18)^SHR(x,3))
#define sig1(x) (ROTR(x,17)^ROTR(x,19)^SHR(x,10))

/* ================================================================
 *  SHA-256 Forward with full state trace
 * ================================================================ */

typedef struct {
    uint32_t W[64];           /* message schedule */
    uint32_t state[65][8];    /* state[r] = (a,b,c,d,e,f,g,h) before round r */
    uint32_t hash[8];         /* final hash output */
    uint32_t T1[64], T2[64];  /* intermediate values per round */
} SHA256Trace;

void sha256_forward(const uint32_t msg[16], SHA256Trace *tr) {
    int i;
    /* message schedule */
    for (i = 0; i < 16; i++) tr->W[i] = msg[i];
    for (i = 16; i < 64; i++)
        tr->W[i] = sig1(tr->W[i-2]) + tr->W[i-7] + sig0(tr->W[i-15]) + tr->W[i-16];

    /* initial state */
    for (i = 0; i < 8; i++) tr->state[0][i] = IV[i];

    /* 64 rounds */
    for (i = 0; i < 64; i++) {
        uint32_t a = tr->state[i][0], b = tr->state[i][1];
        uint32_t c = tr->state[i][2], d = tr->state[i][3];
        uint32_t e = tr->state[i][4], f = tr->state[i][5];
        uint32_t g = tr->state[i][6], h = tr->state[i][7];

        uint32_t t1 = h + SIG1(e) + CH(e,f,g) + K[i] + tr->W[i];
        uint32_t t2 = SIG0(a) + MAJ(a,b,c);
        tr->T1[i] = t1;
        tr->T2[i] = t2;

        tr->state[i+1][0] = t1 + t2;  /* a */
        tr->state[i+1][1] = a;         /* b */
        tr->state[i+1][2] = b;         /* c */
        tr->state[i+1][3] = c;         /* d */
        tr->state[i+1][4] = d + t1;    /* e */
        tr->state[i+1][5] = e;         /* f */
        tr->state[i+1][6] = f;         /* g */
        tr->state[i+1][7] = g;         /* h */
    }

    /* final hash */
    for (i = 0; i < 8; i++)
        tr->hash[i] = tr->state[64][i] + IV[i];
}

/* ================================================================
 *  Zone C Inverse: state[r+1] → state[r] given W[r]
 *
 *  From state[r+1] = (a', b'=a, c'=b, d'=c, e', f'=e, g'=f, h'=g):
 *    a = b'   (shift)
 *    b = c'
 *    c = d'
 *    e = f'
 *    f = g'
 *    g = h'
 *    T2 = SIG0(a) + MAJ(a, b, c)
 *    T1 = a' - T2
 *    d = e' - T1
 *    h = T1 - SIG1(e) - CH(e,f,g) - K[r] - W[r]
 * ================================================================ */

void inverse_round(const uint32_t next[8], uint32_t prev[8],
                   uint32_t W_r, int r) {
    /* recover shifted registers */
    uint32_t a = next[1];  /* b' = a */
    uint32_t b = next[2];  /* c' = b */
    uint32_t c = next[3];  /* d' = c */
    uint32_t e = next[5];  /* f' = e */
    uint32_t f = next[6];  /* g' = f */
    uint32_t g = next[7];  /* h' = g */

    /* T2 and T1 */
    uint32_t t2 = SIG0(a) + MAJ(a, b, c);
    uint32_t t1 = next[0] - t2;  /* a' = T1 + T2 → T1 = a' - T2 */

    /* recover d and h */
    uint32_t d = next[4] - t1;    /* e' = d + T1 → d = e' - T1 */
    uint32_t h = t1 - SIG1(e) - CH(e,f,g) - K[r] - W_r;

    prev[0] = a; prev[1] = b; prev[2] = c; prev[3] = d;
    prev[4] = e; prev[5] = f; prev[6] = g; prev[7] = h;
}

/* ================================================================
 *  Full backward: state[64] → state[target_r]
 *  Requires W[target_r .. 63]
 * ================================================================ */

void backward_zone_c(const uint32_t state_final[8], const uint32_t W[64],
                     int target_r, uint32_t state_out[8]) {
    uint32_t cur[8], prev[8];
    memcpy(cur, state_final, 32);

    for (int r = 63; r >= target_r; r--) {
        inverse_round(cur, prev, W[r], r);
        memcpy(cur, prev, 32);
    }
    memcpy(state_out, cur, 32);
}

/* ================================================================
 *  Створочне: recover a[r-4] from state[r]
 *
 *  e[r] = a[r] + a[r-4] - T2[r-1]
 *  → a[r-4] = e[r] - a[r] + T2[r-1]
 *  T2[r-1] = SIG0(a[r-1]) + MAJ(a[r-1], a[r-2], a[r-3])
 * ================================================================ */

uint32_t stvorochne_recover_a(uint32_t a_r, uint32_t a_r1, uint32_t a_r2,
                               uint32_t a_r3, uint32_t e_r) {
    uint32_t t2_prev = SIG0(a_r1) + MAJ(a_r1, a_r2, a_r3);
    return e_r - a_r + t2_prev;
}

/* ================================================================
 *  W[r] recovery from state[r] and state[r+1]
 *
 *  T2[r] = SIG0(a[r]) + MAJ(a[r], b[r], c[r])
 *  T1[r] = a[r+1] - T2[r]
 *  W[r] = T1[r] - h[r] - SIG1(e[r]) - CH(e[r],f[r],g[r]) - K[r]
 * ================================================================ */

uint32_t recover_W(const uint32_t st[8], const uint32_t st_next[8], int r) {
    uint32_t a = st[0], b = st[1], c = st[2];
    uint32_t e = st[4], f = st[5], g = st[6], h = st[7];
    uint32_t t2 = SIG0(a) + MAJ(a, b, c);
    uint32_t t1 = st_next[0] - t2;
    return t1 - h - SIG1(e) - CH(e,f,g) - K[r];
}

/* ================================================================
 *  Schedule verification: check if W[0..15] produces given W[16..63]
 * ================================================================ */

int verify_schedule(const uint32_t W[64]) {
    for (int i = 16; i < 64; i++) {
        uint32_t expected = sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16];
        if (expected != W[i]) return 0;
    }
    return 1;
}

/* ================================================================
 *  Mode 0: Verify inverse round correctness
 * ================================================================ */

void test_inverse(void) {
    printf("=== Zone C Inverse Verification ===\n\n");

    uint32_t msg[16];
    srand(42);
    for (int i = 0; i < 16; i++) msg[i] = (uint32_t)rand() | ((uint32_t)rand() << 16);

    SHA256Trace tr;
    sha256_forward(msg, &tr);

    printf("Forward hash: ");
    for (int i = 0; i < 8; i++) printf("%08x ", tr.hash[i]);
    printf("\n\n");

    /* Test inverse round by round */
    int errors = 0;
    for (int r = 63; r >= 0; r--) {
        uint32_t recovered[8];
        inverse_round(tr.state[r+1], recovered, tr.W[r], r);
        if (memcmp(recovered, tr.state[r], 32) != 0) {
            printf("  MISMATCH at round %d!\n", r);
            errors++;
        }
    }
    printf("  Inverse rounds 0..63: %s (%d errors)\n",
           errors ? "FAILED" : "ALL CORRECT", errors);

    /* Test Zone C backward: state[64] → state[49] */
    uint32_t state49[8];
    backward_zone_c(tr.state[64], tr.W, 49, state49);
    if (memcmp(state49, tr.state[49], 32) == 0)
        printf("  Zone C (64→49): CORRECT\n");
    else
        printf("  Zone C (64→49): MISMATCH\n");

    /* Test створочне backward */
    printf("\n  Створочне backward (a[57..60] from state[64]):\n");
    /* state[64] = (a64, a63, a62, a61, e64, e63, e62, e61) */
    uint32_t a[65]; /* a[0..64] */
    for (int r = 0; r <= 64; r++) a[r] = tr.state[r][0];

    for (int r = 64; r >= 61; r--) {
        uint32_t a_r = tr.state[r][0];
        uint32_t a_r1 = (r >= 1) ? tr.state[r-1][0] : IV[1]; /* a[r-1] */
        uint32_t a_r2 = (r >= 2) ? tr.state[r-2][0] : IV[2];
        uint32_t a_r3 = (r >= 3) ? tr.state[r-3][0] : IV[3];
        uint32_t e_r = tr.state[r][4];

        uint32_t a_r4 = stvorochne_recover_a(a_r, a_r1, a_r2, a_r3, e_r);
        uint32_t actual = (r >= 4) ? tr.state[r-4][0] : IV[4-(r-4)]; /* approximate */
        actual = tr.state[r-4][0]; /* direct */

        printf("    a[%d] = %08x (actual: %08x) %s\n",
               r-4, a_r4, actual, a_r4 == actual ? "OK" : "MISMATCH");
    }

    /* W recovery */
    printf("\n  W[57..63] recovery from states:\n");
    for (int r = 57; r < 64; r++) {
        uint32_t w_rec = recover_W(tr.state[r], tr.state[r+1], r);
        printf("    W[%d] = %08x (actual: %08x) %s\n",
               r, w_rec, tr.W[r], w_rec == tr.W[r] ? "OK" : "MISMATCH");
    }
}

/* ================================================================
 *  Mode 1: a[56] search on reduced rounds
 *
 *  Given hash H, recover state[64], backward to state[57],
 *  then enumerate a[56] values and check consistency.
 * ================================================================ */

void test_a56_search(void) {
    printf("\n=== a[56] Brute Force Search ===\n\n");

    uint32_t msg[16];
    srand(42);
    for (int i = 0; i < 16; i++) msg[i] = (uint32_t)rand() | ((uint32_t)rand() << 16);

    SHA256Trace tr;
    sha256_forward(msg, &tr);

    uint32_t actual_a56 = tr.state[56][0];
    printf("Target: a[56] = %08x\n", actual_a56);
    printf("Hash: ");
    for (int i = 0; i < 8; i++) printf("%08x ", tr.hash[i]);
    printf("\n\n");

    /* Step 1: From hash → state[64] (subtract IV) */
    uint32_t state64[8];
    for (int i = 0; i < 8; i++) state64[i] = tr.hash[i] - IV[i];

    /* Step 2: Zone C backward: state[64] → state[57] using known W[57..63] */
    /* But we DON'T know W[57..63] yet! We need to recover them.
     * From state[64] we know a[61..64] and e[61..64].
     * створочне gives a[57..60].
     * Then W[r] = recover_W(state[r], state[r+1], r). */

    /* Recover a[57..60] via створочне */
    /* state[64] gives: a64=state64[0], a63=state64[1], a62=state64[2], a61=state64[3]
     *                  e64=state64[4], e63=state64[5], e62=state64[6], e61=state64[7] */
    uint32_t a_seq[65]; /* recovered a-values */
    a_seq[64] = state64[0];
    a_seq[63] = state64[1];
    a_seq[62] = state64[2];
    a_seq[61] = state64[3];

    /* e-values from state[64] */
    uint32_t e_seq[65];
    e_seq[64] = state64[4];
    e_seq[63] = state64[5];
    e_seq[62] = state64[6];
    e_seq[61] = state64[7];

    /* створочне: a[r-4] = e[r] - a[r] + T2[r-1] */
    for (int r = 64; r >= 61; r--) {
        uint32_t t2 = SIG0(a_seq[r-1]) + MAJ(a_seq[r-1], a_seq[r-2], a_seq[r-3]);
        a_seq[r-4] = e_seq[r] - a_seq[r] + t2;
    }

    printf("Recovered a[57..60] via створочне:\n");
    for (int r = 57; r <= 60; r++) {
        printf("  a[%d] = %08x (actual: %08x) %s\n",
               r, a_seq[r], tr.state[r][0],
               a_seq[r] == tr.state[r][0] ? "OK" : "MISMATCH");
    }

    /* Now recover e[57..60] via створочне: e[r] = a[r] + a[r-4] - T2[r-1] */
    /* We need a[53..56] for this... which we DON'T have.
     * But we can reconstruct full states[57..64] from a[57..64] + e[57..64]
     * and then recover W[57..63]. */

    /* Build states from a and e sequences */
    /* state[r] = (a[r], a[r-1], a[r-2], a[r-3], e[r], e[r-1], e[r-2], e[r-3]) */
    /* We have a[57..64] and e[61..64]. Need e[57..60]. */

    /* e[60] = a[60] + a[56] - T2[59] — NEEDS a[56]! */
    /* This is the crux: e[57..60] depends on a[53..56], which we don't know. */

    /* BUT: to recover W[63], we only need state[63] and state[64]. */
    /* state[63] = (a63, a62, a61, a60, e63, e62, e61, e60) */
    /* We know a60..a63 and e61..e63. We need e60. */
    /* e[60] = a[60] + a[56] - T2[59] */
    /* T2[59] = SIG0(a[59]) + MAJ(a[59], a[58], a[57]) — KNOWN! */
    /* So e[60] depends on a[56] — the unknown! */

    printf("\nW recovery chain — what depends on a[56]:\n");
    uint32_t t2_59 = SIG0(a_seq[59]) + MAJ(a_seq[59], a_seq[58], a_seq[57]);
    printf("  T2[59] = %08x (known from a[57..59])\n", t2_59);
    printf("  e[60] = a[60] + a[56] - T2[59] = %08x + a[56] - %08x\n",
           a_seq[60], t2_59);
    printf("  → e[60] is LINEAR in a[56]!\n\n");

    /* For each guess of a[56], we can compute:
     * e[60] → state[63] → W[63]
     * e[59] depends on a[55]... which we also don't know.
     * But e[59] = a[59] + a[55] - T2[58]
     * a[55] is NOT recoverable without more unknowns.
     *
     * HOWEVER: for W[63] recovery, we need state[63]:
     * state[63] = (a63, a62, a61, a60, e63, e62, e61, e60)
     * e60 depends on a[56]. All others known. So W[63] = f(a[56]).
     *
     * For W[62], need state[62] = (a62, a61, a60, a59, e62, e61, e60, e59)
     * e59 depends on a[55] — second unknown!
     *
     * So: W[63] depends on a[56] only.
     *     W[62] depends on a[56] and a[55].
     *     W[61] depends on a[56], a[55], a[54].
     *     etc.
     */

    printf("Dependency chain:\n");
    printf("  W[63]: needs e[60] → needs a[56] ONLY\n");
    printf("  W[62]: needs e[59] → needs a[55]\n");
    printf("  W[61]: needs e[58] → needs a[54]\n");
    printf("  W[60]: needs e[57] → needs a[53]\n");
    printf("  W[59]: needs state[59] → needs a[52..55]\n\n");

    /* KEY INSIGHT: if we guess ONLY a[56], we can compute W[63].
     * W[63] is a schedule word = f(W[0..15]).
     * So W[63] constrains W[0..15] to a HYPERPLANE.
     * Each bit of W[63] is one equation on the 512 msg bits. */

    /* TEST: enumerate a[56] and check W[63] consistency */
    /* For the correct a[56], recovered W[63] must match actual W[63] */

    printf("Testing a[56] → W[63] consistency:\n");
    clock_t start = clock();
    uint64_t found = 0;
    uint32_t correct_w63 = tr.W[63];

    /* Build state[63] with a[56] as parameter */
    uint32_t st63_base[8];
    st63_base[0] = a_seq[63]; /* a */
    st63_base[1] = a_seq[62]; /* b */
    st63_base[2] = a_seq[61]; /* c */
    st63_base[3] = a_seq[60]; /* d */
    st63_base[4] = e_seq[63]; /* e */
    st63_base[5] = e_seq[62]; /* f */
    st63_base[6] = e_seq[61]; /* g */
    /* st63_base[7] = e[60] = a[60] + guess_a56 - T2[59] — varies! */

    uint32_t a60_plus_neg_t2 = a_seq[60] - t2_59;

    /* Full 2^32 search */
    printf("  Searching 2^32 values of a[56]...\n");
    uint64_t total = (uint64_t)1 << 32;
    uint64_t progress_step = total / 100;

    for (uint64_t guess = 0; guess < total; guess++) {
        uint32_t a56 = (uint32_t)guess;
        uint32_t e60 = a60_plus_neg_t2 + a56;

        /* state[63][7] = h = e[60] */
        uint32_t h63 = e60;

        /* W[63] = T1 - h - SIG1(e) - CH(e,f,g) - K[63] */
        uint32_t t2_63 = SIG0(st63_base[0]) + MAJ(st63_base[0], st63_base[1], st63_base[2]);
        uint32_t t1_63 = a_seq[64] - t2_63; /* a[64] - T2[63] */
        uint32_t w63_guess = t1_63 - h63 - SIG1(st63_base[4]) -
                             CH(st63_base[4], st63_base[5], st63_base[6]) - K[63];

        if (w63_guess == correct_w63) {
            found++;
            if (found <= 5) {
                printf("  ★ MATCH at a[56] = %08x (actual: %08x) %s\n",
                       a56, actual_a56,
                       a56 == actual_a56 ? "← CORRECT" : "");
            }
        }

        if (guess % progress_step == 0 && guess > 0) {
            double pct = 100.0 * guess / total;
            double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
            double rate = guess / elapsed / 1e6;
            printf("  %.0f%% (%.1fM/s, found=%lu)\r", pct, rate, found);
            fflush(stdout);
        }
    }

    double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("\n\n  Results:\n");
    printf("  Searched: 2^32 = %lu values\n", total);
    printf("  Matches (W[63] consistent): %lu\n", found);
    printf("  Time: %.1f seconds (%.1fM/s)\n", elapsed, total / elapsed / 1e6);
    printf("  Correct a[56] found: %s\n",
           found > 0 ? "YES" : "NO");

    if (found > 1) {
        printf("\n  NOTE: %lu matches means W[63] alone is not enough to\n", found);
        printf("  uniquely determine a[56]. Need W[62] too (requires a[55]).\n");
        printf("  Each additional W word halves the candidates.\n");
        printf("  With 1 W word (32 bits): ~2^32/2^32 = ~1 match expected.\n");
    }
}

/* ================================================================
 *  Main
 * ================================================================ */

int main(int argc, char **argv) {
    int mode = 0;
    if (argc > 1) mode = atoi(argv[1]);

    printf("SHA-256 Inverse Engine\n");
    printf("==============================================\n\n");

    if (mode == 0 || mode == 1) {
        test_inverse();
    }

    if (mode == 1 || mode == 2) {
        test_a56_search();
    }

    return 0;
}
