/*
 * de17_search.c — Fast search for De17=0 in SHA-256 Wang chain.
 *
 * Uses multi-target strategy:
 *   - Generate M random base messages W[0..11]
 *   - For each, sweep W[14] through 2^32/M values
 *   - Check if De17(message, DW[0]=0x80000000, W14) = 0
 *
 * Compile: gcc -O3 -o de17_search de17_search.c
 * Run: ./de17_search [num_messages] [trials_per_message]
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define MASK 0xFFFFFFFFU

static inline uint32_t rotr(uint32_t x, int n) {
    return (x >> n) | (x << (32 - n));
}

static inline uint32_t Sig0(uint32_t x) { return rotr(x,2) ^ rotr(x,13) ^ rotr(x,22); }
static inline uint32_t Sig1(uint32_t x) { return rotr(x,6) ^ rotr(x,11) ^ rotr(x,25); }
static inline uint32_t sig0(uint32_t x) { return rotr(x,7) ^ rotr(x,18) ^ (x>>3); }
static inline uint32_t sig1(uint32_t x) { return rotr(x,17) ^ rotr(x,19) ^ (x>>10); }
static inline uint32_t Ch(uint32_t e, uint32_t f, uint32_t g) { return (e&f)^((~e)&g); }
static inline uint32_t Maj(uint32_t a, uint32_t b, uint32_t c) { return (a&b)^(a&c)^(b&c); }

static const uint32_t IV[8] = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
};

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

#define DW_BIT 0x80000000U

typedef struct {
    uint32_t s[8];
} State;

static inline void sha_round_fn(State *st, uint32_t W_r, uint32_t K_r) {
    uint32_t a=st->s[0], b=st->s[1], c=st->s[2], d=st->s[3];
    uint32_t e=st->s[4], f=st->s[5], g=st->s[6], h=st->s[7];
    uint32_t T1 = h + Sig1(e) + Ch(e,f,g) + K_r + W_r;
    uint32_t T2 = Sig0(a) + Maj(a,b,c);
    st->s[0]=T1+T2; st->s[1]=a; st->s[2]=b; st->s[3]=c;
    st->s[4]=d+T1; st->s[5]=e; st->s[6]=f; st->s[7]=g;
}

static uint32_t xrand(void) {
    return ((uint32_t)rand() << 16) ^ (uint32_t)rand() ^ ((uint32_t)rand() << 28);
}

static int popcount32(uint32_t x) {
    return __builtin_popcount(x);
}

int main(int argc, char **argv) {
    int M = 1024;      /* number of base messages */
    long trials = 1L << 22;  /* trials per message */

    if (argc > 1) M = atoi(argv[1]);
    if (argc > 2) trials = atol(argv[2]);

    srand(time(NULL));

    printf("De17=0 search: M=%d messages, %ld trials/msg, total=%ld\n",
           M, trials, (long)M * trials);

    int best_hw = 32;
    uint32_t best_de17 = 0xFFFFFFFF;
    int best_msg = -1;
    uint32_t best_w14 = 0;
    int found = 0;

    for (int msg = 0; msg < M && !found; msg++) {
        /* Generate random base message W[0..15] */
        uint32_t W[64];
        for (int i = 0; i < 16; i++) W[i] = xrand();

        /* Precompute schedule entries that don't depend on W[14] */
        /* W[16] = sig1(W[14]) + W[9] + sig0(W[1]) + W[0] */
        /* W[14] only affects W[16] (and beyond, but we only go to r=17) */
        uint32_t W16_base = W[9] + sig0(W[1]) + W[0];
        /* For perturbed: W'[0] = W[0] ^ DW_BIT */
        uint32_t W16_base_f = W[9] + sig0(W[1]) + (W[0] ^ DW_BIT);

        /* Precompute states through round 13 (uses W[0..13]) */
        /* (need full schedule through W[13]) */
        /* Actually: let's just precompute state through round 11 */
        /* Rounds 0-11 use W[0..11], W[12] enters at round 12 */
        /* But W[12..15] are free, so we precompute through round 11 */

        /* Full schedule for W[0..15] */
        for (int i = 16; i < 17; i++) {
            W[i] = sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16];
        }

        /* Also compute W' schedule */
        uint32_t Wf[64];
        memcpy(Wf, W, sizeof(uint32_t)*64);
        Wf[0] = W[0] ^ DW_BIT;
        /* DW only in W[0], so Wf[i] = W[i] for i=1..15 */
        /* Wf[16] = sig1(Wf[14]) + Wf[9] + sig0(Wf[1]) + Wf[0] */
        /* = sig1(W[14]) + W[9] + sig0(W[1]) + W[0]^DW_BIT */
        /* = W[16] - W[0] + (W[0]^DW_BIT) */

        /* Precompute state through round 13 */
        State sn, sf;
        memcpy(sn.s, IV, sizeof(IV));
        memcpy(sf.s, IV, sizeof(IV));

        for (int r = 0; r < 14; r++) {
            sha_round_fn(&sn, W[r], K[r]);
            sha_round_fn(&sf, (r==0 ? Wf[0] : W[r]), K[r]);
        }
        /* sn, sf are states after round 13 (used W[0..13]) */
        /* Rounds 14, 15, 16 will use W[14], W[15], W[16] */

        /* Now sweep W[14] */
        for (long t = 0; t < trials; t++) {
            uint32_t w14 = (uint32_t)(t & 0xFFFFFFFF);
            if (t >= (1L<<32)) break;

            /* W[16] = sig1(w14) + W16_base */
            uint32_t w16_n = sig1(w14) + W16_base;
            uint32_t w16_f = sig1(w14) + W16_base_f;

            /* Run rounds 14, 15, 16 from precomputed state */
            State sn2 = sn, sf2 = sf;

            /* Round 14: uses W[14] = w14 */
            sha_round_fn(&sn2, w14, K[14]);
            sha_round_fn(&sf2, w14, K[14]);

            /* Round 15: uses W[15] */
            sha_round_fn(&sn2, W[15], K[15]);
            sha_round_fn(&sf2, W[15], K[15]);

            /* Round 16: uses W[16] */
            sha_round_fn(&sn2, w16_n, K[16]);
            sha_round_fn(&sf2, w16_f, K[16]);

            /* De17 = e_n ^ e_f at position [4] */
            uint32_t de17 = sn2.s[4] ^ sf2.s[4];

            int h = popcount32(de17);
            if (h < best_hw) {
                best_hw = h;
                best_de17 = de17;
                best_msg = msg;
                best_w14 = w14;
                printf("  msg=%d w14=0x%08x De17=0x%08x HW=%d\n",
                       msg, w14, de17, h);
                fflush(stdout);
            }
            if (h == 0) {
                found = 1;
                break;
            }
        }

        if ((msg+1) % 100 == 0) {
            printf("  ... %d/%d messages done, best HW=%d\n", msg+1, M, best_hw);
            fflush(stdout);
        }
    }

    printf("\n");
    if (found) {
        printf("*** De17 = 0 FOUND! ***\n");
        printf("  Message #%d, W14 = 0x%08x\n", best_msg, best_w14);
        printf("  De17 = 0x%08x\n", best_de17);
    } else {
        printf("Best: msg=%d, W14=0x%08x, De17=0x%08x, HW=%d\n",
               best_msg, best_w14, best_de17, best_hw);
    }

    return 0;
}
