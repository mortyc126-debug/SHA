/*
 * step10_dual_barrier.c — Search for De17=De18=0 simultaneously.
 *
 * Strategy: sweep (W[14], W[15]) as 2 free words (64 bits).
 * Target: 64-bit condition (De17=0 AND De18=0).
 * Expected: ~2^64 / 2^64 = 1 solution per message.
 * With M messages: birthday gives ~sqrt(2^64) = 2^32 per message.
 * But direct: need ~2^64 trials total → too much.
 *
 * Better: use ALL 4 free words (W[12..15] = 128 bits).
 * Target: 64-bit condition. Expected: 2^64 solutions per message.
 * Random search: 1 hit per 2^64/2^128 = instant? No.
 * P(hit) = 2^{-64} per trial. Need 2^64 trials. Still too much.
 *
 * ACTUAL approach: fix W[12,13], sweep (W[14],W[15]).
 * 64-bit search space, 64-bit target → ~1 solution per combination.
 * Birthday on halves: store De17(W14) for 2^20 W14 values.
 * For each W15, compute De17 and De18. Look up.
 * This doesn't directly work because De17 depends on BOTH W14 and W15.
 *
 * REAL approach: 2-phase search
 * Phase 1: For 2^24 random (W14,W15) pairs, store De17 in hash table.
 *          Key=De17, Value=(W14,W15).
 * Phase 2: For 2^24 random (W14,W15) pairs, compute De17.
 *          If De17 matches any Phase 1 entry AND De18=0 → WIN.
 *          ... this doesn't work either.
 *
 * SIMPLEST: just sweep (W14, W15) with W14 in [0, 2^20), W15 in [0, 2^20).
 * Total: 2^40 trials. 64-bit target. P(hit) = 2^40/2^64 = 2^-24. Not enough.
 *
 * With M=2^24 messages: P = 1. Total work: 2^24 × 2^40 = 2^64. Too much.
 *
 * BETTER APPROACH: partial target.
 * Find De17=0 first (costs 2^32).
 * Given De17=0, how constrained is De18?
 * If De18 is "free" given De17=0 → need another 2^32 search.
 * But W[14] is already fixed by De17=0 → only W[15] free for De18.
 * W[15] → De18: 32→32 map. Need De18=0. Costs 2^32.
 * But W[15] also affects De17! So fixing De17=0 AND varying W[15]
 * might break De17=0.
 *
 * SOLUTION: Use W[14] for De17, W[15] for De18, check if they decouple.
 * If De17 ≈ f(W14) + g(W15) (approximately additive):
 *   Fix W14 for De17=0, then vary W15 for De18=0.
 * If not additive: need joint search.
 *
 * Let's just test: sweep W14 low range, for each check multiple W15.
 * Find (W14, W15) with De17=0 AND min HW(De18).
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define MASK 0xFFFFFFFFU

static inline uint32_t rotr(uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }
static inline uint32_t Sig0(uint32_t x) { return rotr(x,2) ^ rotr(x,13) ^ rotr(x,22); }
static inline uint32_t Sig1(uint32_t x) { return rotr(x,6) ^ rotr(x,11) ^ rotr(x,25); }
static inline uint32_t sig0f(uint32_t x) { return rotr(x,7) ^ rotr(x,18) ^ (x>>3); }
static inline uint32_t sig1f(uint32_t x) { return rotr(x,17) ^ rotr(x,19) ^ (x>>10); }
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

typedef struct { uint32_t s[8]; } State;

static inline void sha_round_fn(State *st, uint32_t W_r, uint32_t K_r) {
    uint32_t a=st->s[0], b=st->s[1], c=st->s[2], d=st->s[3];
    uint32_t e=st->s[4], f=st->s[5], g=st->s[6], h=st->s[7];
    uint32_t T1 = h + Sig1(e) + Ch(e,f,g) + K_r + W_r;
    uint32_t T2 = Sig0(a) + Maj(a,b,c);
    st->s[0]=T1+T2; st->s[1]=a; st->s[2]=b; st->s[3]=c;
    st->s[4]=d+T1;  st->s[5]=e; st->s[6]=f; st->s[7]=g;
}

static uint32_t xrand(void) {
    return ((uint32_t)rand() << 16) ^ (uint32_t)rand() ^ ((uint32_t)rand() << 28);
}

int main() {
    srand(time(NULL));

    printf("=== DUAL BARRIER SEARCH: De17=0 AND De18=0 ===\n\n");

    /* Strategy: for each random message, precompute state at round 13.
     * Then sweep (W14, W15) to find De17=De18=0.
     *
     * Optimization: For each W15, first sweep W14 to find De17=0.
     * If De17=0, check De18. This costs 2^32 per W15.
     *
     * Better: for each message, sweep W14 in [0, 2^22).
     * For each De17=0 hit, record (W14, De18).
     * Expected De17=0 hits per message: 2^22/2^32 = 2^-10 ≈ 0.001
     * Need ~2^10 messages to get one De17=0 hit.
     * Then check De18: P(De18=0 | De17=0) = 1/2^32.
     * Total: need 2^10 × 2^32 = 2^42 → too much.
     *
     * ACTUALLY: let's think differently.
     * Sweep W14 in full 32-bit range per message.
     * For the ~1 hit where De17=0, record W14 and De18.
     * Then across M messages, collect M values of De18.
     * If any De18=0 → we win.
     * Need M = 2^32 messages. Each with 2^32 W14 sweep. Total: 2^64. Bad.
     *
     * REAL SOLUTION: use BOTH W14 and W15 jointly.
     * (W14, W15) = 64 bits. Target (De17, De18) = 64 bits.
     * Random search: 2^64 per message. With M messages: 2^64/M.
     * For M=1: sweep 2^64 → too slow in C (~years).
     *
     * PRACTICAL: sweep moderate range, find NEAR misses.
     * Use W14 in [0, 2^20), W15 in [0, 2^20). Total: 2^40 per message.
     * Best De17+De18: expect HW ≈ 64 - 40 = 24 bits remaining.
     * With 100 messages: 100 × 2^40 = 2^46.6 → still big.
     *
     * LET'S JUST DO: 10 messages × W14 in [0, 2^20) × W15 in [0, 2^12)
     * = 10 × 2^32 = ~4×10^10. This is feasible (~minutes).
     */

    int best_hw17 = 32, best_hw18 = 32, best_joint = 64;
    uint32_t best_w14 = 0, best_w15 = 0;
    int best_msg = -1;
    int found = 0;

    int M = 20;
    uint32_t w14_range = 1U << 21;  /* 2^21 W14 values */
    uint32_t w15_range = 1U << 11;  /* 2^11 W15 values */

    printf("Search: M=%d msgs, W14 range=2^21, W15 range=2^11\n", M);
    printf("Total per msg: %u evals\n", w14_range * w15_range);
    printf("Total: %.2e evals\n\n", (double)M * w14_range * w15_range);

    for (int msg = 0; msg < M && !found; msg++) {
        uint32_t W[64];
        for (int i = 0; i < 16; i++) W[i] = xrand();

        /* Precompute state through round 13 */
        State sn, sf;
        memcpy(sn.s, IV, 32);
        memcpy(sf.s, IV, 32);

        uint32_t Wf0 = W[0] ^ DW_BIT;

        /* Full schedule for rounds 0..13 */
        for (int i = 16; i < 64; i++)
            W[i] = sig1f(W[i-2]) + W[i-7] + sig0f(W[i-15]) + W[i-16];

        for (int r = 0; r < 13; r++) {
            sha_round_fn(&sn, W[r], K[r]);
            sha_round_fn(&sf, (r==0 ? Wf0 : W[r]), K[r]);
        }

        /* Round 13 uses W[13] */
        sha_round_fn(&sn, W[13], K[13]);
        sha_round_fn(&sf, W[13], K[13]);

        /* Now sn, sf are states after round 13 (14 rounds done: 0..13) */

        /* Schedule constants that don't depend on W[14,15]:
         * W[16] = sig1(W[14]) + W[9] + sig0(W[1]) + W[0]
         * W[17] = sig1(W[15]) + W[10] + sig0(W[2]) + W[1]
         * W[18] = sig1(W[16]) + W[11] + sig0(W[3]) + W[2]
         *       = sig1(sig1(W[14]) + base16) + W[11] + sig0(W[3]) + W[2]
         */
        uint32_t base16_n = W[9] + sig0f(W[1]) + W[0];
        uint32_t base16_f = W[9] + sig0f(W[1]) + Wf0;
        uint32_t base17_rest = W[10] + sig0f(W[2]) + W[1];
        uint32_t base18_rest = W[11] + sig0f(W[3]) + W[2];

        for (uint32_t w14 = 0; w14 < w14_range; w14++) {
            uint32_t w16_n = sig1f(w14) + base16_n;
            uint32_t w16_f = sig1f(w14) + base16_f;

            /* Run round 14 with w14 */
            State s14n = sn, s14f = sf;
            sha_round_fn(&s14n, w14, K[14]);
            sha_round_fn(&s14f, w14, K[14]);

            for (uint32_t w15 = 0; w15 < w15_range; w15++) {
                uint32_t w17_n = sig1f(w15) + base17_rest;
                uint32_t w17_f = sig1f(w15) + base17_rest; /* W[1] same */

                uint32_t w18_n = sig1f(w16_n) + base18_rest;
                uint32_t w18_f = sig1f(w16_f) + base18_rest;

                /* Round 15: uses W[15] = w15 */
                State s15n = s14n, s15f = s14f;
                sha_round_fn(&s15n, w15, K[15]);
                sha_round_fn(&s15f, w15, K[15]);

                /* Round 16: uses W[16] */
                State s16n = s15n, s16f = s15f;
                sha_round_fn(&s16n, w16_n, K[16]);
                sha_round_fn(&s16f, w16_f, K[16]);

                uint32_t de17 = s16n.s[4] ^ s16f.s[4];

                /* Round 17: uses W[17] */
                State s17n = s16n, s17f = s16f;
                sha_round_fn(&s17n, w17_n, K[17]);
                sha_round_fn(&s17f, w17_f, K[17]);

                uint32_t de18 = s17n.s[4] ^ s17f.s[4];

                int hw17 = __builtin_popcount(de17);
                int hw18 = __builtin_popcount(de18);
                int joint = hw17 + hw18;

                if (joint < best_joint) {
                    best_joint = joint;
                    best_hw17 = hw17;
                    best_hw18 = hw18;
                    best_w14 = w14;
                    best_w15 = w15;
                    best_msg = msg;
                    printf("  msg=%d w14=0x%08x w15=0x%08x "
                           "De17=0x%08x(HW=%d) De18=0x%08x(HW=%d) sum=%d\n",
                           msg, w14, w15, de17, hw17, de18, hw18, joint);
                    fflush(stdout);
                    if (joint == 0) { found = 1; break; }
                }
            }
            if (found) break;
        }
        printf("  msg %d/%d done, best joint HW = %d\n", msg+1, M, best_joint);
        fflush(stdout);
    }

    printf("\n");
    if (found) {
        printf("*** De17=De18=0 FOUND! ***\n");
    } else {
        printf("Best: msg=%d w14=0x%x w15=0x%x joint_HW=%d (De17 HW=%d, De18 HW=%d)\n",
               best_msg, best_w14, best_w15, best_joint, best_hw17, best_hw18);
    }

    return 0;
}
