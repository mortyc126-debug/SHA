/*
 * Step 17c: Reduced-round SHA-256 collision search.
 *
 * Strategy: DW[0] = +1 (additive), force De=0 for rounds 1-15
 * via computed DW[t], then check if hash matches for N rounds.
 *
 * For each trial:
 * 1. Random W[0..15]
 * 2. W'[0] = W[0] + 1
 * 3. For rounds 1-15: compute DW[t] to force De=0
 * 4. Expand both W and W' through schedule
 * 5. Continue compression and check hash diff
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define MASK32 0xFFFFFFFFU

static inline uint32_t rotr(uint32_t x, int n) {
    return (x >> n) | (x << (32 - n));
}

static inline uint32_t Sigma0(uint32_t x) {
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22);
}

static inline uint32_t Sigma1(uint32_t x) {
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25);
}

static inline uint32_t sigma0(uint32_t x) {
    return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3);
}

static inline uint32_t sigma1(uint32_t x) {
    return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10);
}

static inline uint32_t Ch(uint32_t e, uint32_t f, uint32_t g) {
    return (e & f) ^ (~e & g);
}

static inline uint32_t Maj(uint32_t a, uint32_t b, uint32_t c) {
    return (a & b) ^ (a & c) ^ (b & c);
}

static inline int hw(uint32_t x) {
    return __builtin_popcount(x);
}

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

/* Simple xorshift64 PRNG */
static uint64_t rng_state = 0x123456789ABCDEF0ULL;

static inline uint32_t rand32(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (uint32_t)(rng_state & MASK32);
}

static void expand_message(const uint32_t W16[16], uint32_t W64[64]) {
    for (int i = 0; i < 16; i++) W64[i] = W16[i];
    for (int t = 16; t < 64; t++) {
        W64[t] = sigma1(W64[t-2]) + W64[t-7] + sigma0(W64[t-15]) + W64[t-16];
    }
}

/* SHA-256 round function: updates state in-place */
static inline void sha256_round(uint32_t s[8], uint32_t Wt, uint32_t Kt) {
    uint32_t T1 = s[7] + Sigma1(s[4]) + Ch(s[4], s[5], s[6]) + Kt + Wt;
    uint32_t T2 = Sigma0(s[0]) + Maj(s[0], s[1], s[2]);
    s[7] = s[6];
    s[6] = s[5];
    s[5] = s[4];
    s[4] = s[3] + T1;
    s[3] = s[2];
    s[2] = s[1];
    s[1] = s[0];
    s[0] = T1 + T2;
}

int main(void) {
    rng_state = (uint64_t)time(NULL) ^ 0xDEADBEEFCAFEULL;

    printf("Step 17c: Reduced-round SHA-256 collision search\n");
    printf("Strategy: DW[0]=+1, De=0 forced for rounds 1-15\n\n");

    uint64_t trials = 0;
    int best_hw_17 = 999, best_hw_18 = 999, best_hw_19 = 999;
    int best_hw_20 = 999, best_hw_21 = 999;
    uint32_t best_W[16], best_W2[16];
    int best_round = 0;

    uint64_t max_trials = 100000000ULL;  /* 10^8 */

    while (trials < max_trials) {
        uint32_t W[16], W2[16];
        for (int i = 0; i < 16; i++) {
            W[i] = rand32();
            W2[i] = W[i];
        }
        W2[0] = W[0] + 1;  /* DW[0] = +1 */

        /* Compute rounds 0-15 with De=0 forcing */
        uint32_t s1[8], s2[8];
        for (int i = 0; i < 8; i++) { s1[i] = IV[i]; s2[i] = IV[i]; }

        /* Round 0: fixed DW */
        sha256_round(s1, W[0], K[0]);
        sha256_round(s2, W2[0], K[0]);

        /* Rounds 1-15: force De=0 */
        for (int t = 1; t < 16; t++) {
            uint32_t T1_p1 = s1[7] + Sigma1(s1[4]) + Ch(s1[4], s1[5], s1[6]) + K[t];
            uint32_t T1_p2 = s2[7] + Sigma1(s2[4]) + Ch(s2[4], s2[5], s2[6]) + K[t];
            uint32_t dw = (s1[3] + T1_p1) - (s2[3] + T1_p2);
            W2[t] = W[t] + dw;

            sha256_round(s1, W[t], K[t]);
            sha256_round(s2, W2[t], K[t]);
        }

        /* Expand messages */
        uint32_t Wexp[64], W2exp[64];
        expand_message(W, Wexp);
        expand_message(W2, W2exp);

        /* Continue rounds 16-20 */
        for (int t = 16; t < 21; t++) {
            sha256_round(s1, Wexp[t], K[t]);
            sha256_round(s2, W2exp[t], K[t]);

            /* Check hash diff at this round */
            int thw = 0;
            for (int j = 0; j < 8; j++) {
                thw += hw((IV[j] + s1[j]) ^ (IV[j] + s2[j]));
            }

            int nr = t + 1;
            int *best = NULL;
            if (nr == 17) best = &best_hw_17;
            else if (nr == 18) best = &best_hw_18;
            else if (nr == 19) best = &best_hw_19;
            else if (nr == 20) best = &best_hw_20;
            else if (nr == 21) best = &best_hw_21;

            if (best && thw < *best) {
                *best = thw;
                if (thw <= 10) {
                    printf("[trial %lu] %d-round hash diff HW = %d\n",
                           (unsigned long)trials, nr, thw);
                    if (nr > best_round) {
                        best_round = nr;
                        memcpy(best_W, W, sizeof(W));
                        memcpy(best_W2, W2, sizeof(W2));
                    }
                }
                if (thw == 0) {
                    printf("\n*** %d-ROUND COLLISION FOUND! ***\n", nr);
                    printf("W  = {");
                    for (int i = 0; i < 16; i++) printf("0x%08x%s", W[i], i<15?", ":"}\n");
                    printf("W' = {");
                    for (int i = 0; i < 16; i++) printf("0x%08x%s", W2[i], i<15?", ":"}\n");

                    /* Verify */
                    uint32_t vs1[8], vs2[8], vW[64], vW2_64[64];
                    expand_message(W, vW);
                    expand_message(W2, vW2_64);
                    for (int j = 0; j < 8; j++) { vs1[j] = IV[j]; vs2[j] = IV[j]; }
                    for (int j = 0; j < nr; j++) {
                        sha256_round(vs1, vW[j], K[j]);
                        sha256_round(vs2, vW2_64[j], K[j]);
                    }
                    printf("H  = {");
                    for (int j = 0; j < 8; j++) printf("0x%08x%s", IV[j]+vs1[j], j<7?", ":"}\n");
                    printf("H' = {");
                    for (int j = 0; j < 8; j++) printf("0x%08x%s", IV[j]+vs2[j], j<7?", ":"}\n");
                }
            }
        }

        trials++;
        if (trials % 10000000 == 0) {
            printf("... %luM trials | best: 17r=%d 18r=%d 19r=%d 20r=%d 21r=%d\n",
                   (unsigned long)(trials / 1000000),
                   best_hw_17, best_hw_18, best_hw_19, best_hw_20, best_hw_21);
        }
    }

    printf("\nFinal results after %lu trials:\n", (unsigned long)trials);
    printf("  17 rounds: best hash diff HW = %d\n", best_hw_17);
    printf("  18 rounds: best hash diff HW = %d\n", best_hw_18);
    printf("  19 rounds: best hash diff HW = %d\n", best_hw_19);
    printf("  20 rounds: best hash diff HW = %d\n", best_hw_20);
    printf("  21 rounds: best hash diff HW = %d\n", best_hw_21);

    return 0;
}
