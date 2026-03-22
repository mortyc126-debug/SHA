/*
 * de17_verify.c — Find De17=0 AND output full message + verify.
 * Compile: gcc -O3 -o de17_verify de17_verify.c
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

static void schedule(uint32_t *W) {
    for (int i = 16; i < 64; i++)
        W[i] = sig1f(W[i-2]) + W[i-7] + sig0f(W[i-15]) + W[i-16];
}

static uint32_t xrand(void) {
    return ((uint32_t)rand() << 16) ^ (uint32_t)rand() ^ ((uint32_t)rand() << 28);
}

int main() {
    srand(time(NULL) ^ getpid());

    int found = 0;
    int best_hw = 32;
    uint32_t best_W[16], best_w14 = 0;

    printf("Searching for De17=0...\n");
    fflush(stdout);

    for (int msg = 0; msg < 100000 && !found; msg++) {
        uint32_t W[64], Wf[64];
        for (int i = 0; i < 16; i++) W[i] = xrand();

        /* Precompute state through round 13 */
        State sn, sf;
        memcpy(sn.s, IV, 32);
        memcpy(sf.s, IV, 32);

        /* Schedule for rounds 0..13 doesn't depend on W[14] */
        memcpy(Wf, W, 64*4);
        Wf[0] ^= DW_BIT;
        schedule(W); schedule(Wf);

        for (int r = 0; r < 14; r++) {
            sha_round_fn(&sn, W[r], K[r]);
            sha_round_fn(&sf, Wf[r], K[r]);
        }

        /* Sweep W[14] */
        uint32_t w9 = W[9], w1 = W[1], w0_n = W[0], w0_f = Wf[0];
        uint32_t w15 = W[15];
        uint32_t w16_base_n = w9 + sig0f(w1) + w0_n;
        uint32_t w16_base_f = w9 + sig0f(w1) + w0_f;

        for (uint32_t w14 = 0; w14 < (1U<<22); w14++) {
            uint32_t w16_n = sig1f(w14) + w16_base_n;
            uint32_t w16_f = sig1f(w14) + w16_base_f;

            State s2n = sn, s2f = sf;
            sha_round_fn(&s2n, w14, K[14]);
            sha_round_fn(&s2f, w14, K[14]);
            sha_round_fn(&s2n, w15, K[15]);
            sha_round_fn(&s2f, w15, K[15]);
            sha_round_fn(&s2n, w16_n, K[16]);
            sha_round_fn(&s2f, w16_f, K[16]);

            uint32_t de17 = s2n.s[4] ^ s2f.s[4];
            int h = __builtin_popcount(de17);

            if (h < best_hw) {
                best_hw = h;
                memcpy(best_W, W, 16*4);
                best_W[14] = w14;
                /* Undo schedule on best_W to get original W16 */
                best_w14 = w14;
                if (h <= 2)
                    printf("  msg=%d w14=0x%08x De17=0x%08x HW=%d\n", msg, w14, de17, h);
            }
            if (h == 0) {
                found = 1;
                /* Store the actual W[14] value */
                W[14] = w14;
                /* Re-schedule */
                for (int i = 16; i < 64; i++)
                    W[i] = sig1f(W[i-2]) + W[i-7] + sig0f(W[i-15]) + W[i-16];

                printf("\n*** De17=0 FOUND ***\n\n");
                printf("Message W[0..15]:\n");
                for (int i = 0; i < 16; i++)
                    printf("  W[%2d] = 0x%08x\n", i, W[i]);

                printf("\nDifferential: W'[0] = W[0] ^ 0x%08x = 0x%08x\n",
                       DW_BIT, W[0] ^ DW_BIT);

                /* Full verification: 64 rounds */
                memcpy(Wf, W, 64*4);
                Wf[0] ^= DW_BIT;
                for (int i = 16; i < 64; i++)
                    Wf[i] = sig1f(Wf[i-2]) + Wf[i-7] + sig0f(Wf[i-15]) + Wf[i-16];

                State vn, vf;
                memcpy(vn.s, IV, 32);
                memcpy(vf.s, IV, 32);

                printf("\nDifferential profile (all 64 rounds):\n");
                printf("  r |     De(XOR)  HW |     Da(XOR)  HW | total_HW\n");
                printf("----+------------------+------------------+----------\n");

                for (int r = 0; r < 64; r++) {
                    sha_round_fn(&vn, W[r], K[r]);
                    sha_round_fn(&vf, Wf[r], K[r]);

                    uint32_t de = vn.s[4] ^ vf.s[4];
                    uint32_t da = vn.s[0] ^ vf.s[0];
                    int total = 0;
                    for (int i = 0; i < 8; i++)
                        total += __builtin_popcount(vn.s[i] ^ vf.s[i]);

                    char marker[20] = "";
                    if (__builtin_popcount(de) == 0) strcpy(marker, " * De=0");
                    if (total == 0) strcpy(marker, " *** COLLISION");

                    printf(" %2d | 0x%08x  %2d | 0x%08x  %2d | %3d%s\n",
                           r+1, de, __builtin_popcount(de),
                           da, __builtin_popcount(da), total, marker);
                }

                /* Final hash */
                printf("\nFinal hash difference:\n");
                int hash_diff_total = 0;
                for (int i = 0; i < 8; i++) {
                    uint32_t hn = IV[i] + vn.s[i];
                    uint32_t hf = IV[i] + vf.s[i];
                    uint32_t d = hn ^ hf;
                    hash_diff_total += __builtin_popcount(d);
                    printf("  H[%d] = 0x%08x vs 0x%08x  diff=0x%08x HW=%d\n",
                           i, hn, hf, d, __builtin_popcount(d));
                }
                printf("Total hash diff: %d bits\n", hash_diff_total);

                break;
            }
        }

        if ((msg+1) % 500 == 0)
            printf("  ... %d messages, best HW=%d\n", msg+1, best_hw);
    }

    if (!found) {
        printf("Not found. Best HW=%d, W14=0x%08x\n", best_hw, best_w14);
    }

    return 0;
}
