/*
 * Step 17f: Nikolic-Biryukov style reduced-round collision search.
 *
 * Try multiple DW patterns and search for reduced-round collisions
 * by brute force. For N-round reduced SHA-256, just compute both
 * hashes and compare.
 *
 * DW patterns to try:
 * Pattern A: DW[i]=+1, rest=0  (simplest)
 * Pattern B: DW[i]=+1, DW[i+1]=-1 (NB-inspired)
 * Pattern C: DW[i]=+1, DW[i+1]=-1, DW[i+8]=-1
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
static inline uint32_t Sigma0(uint32_t x) { return rotr(x,2)^rotr(x,13)^rotr(x,22); }
static inline uint32_t Sigma1(uint32_t x) { return rotr(x,6)^rotr(x,11)^rotr(x,25); }
static inline uint32_t sigma0(uint32_t x) { return rotr(x,7)^rotr(x,18)^(x>>3); }
static inline uint32_t sigma1(uint32_t x) { return rotr(x,17)^rotr(x,19)^(x>>10); }
static inline uint32_t Ch(uint32_t e, uint32_t f, uint32_t g) { return (e&f)^(~e&g); }
static inline uint32_t Maj(uint32_t a,uint32_t b,uint32_t c) { return (a&b)^(a&c)^(b&c); }
static inline int hw(uint32_t x) { return __builtin_popcount(x); }

static const uint32_t K[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,
    0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,
    0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,
    0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,
    0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,
    0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,
    0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,
    0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,
    0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};
static const uint32_t IV[8] = {
    0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
    0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19
};

static uint64_t rng = 0xACE1ACE2ACE3ULL;
static inline uint32_t rand32(void) {
    rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
    return (uint32_t)(rng);
}

static inline void sha_round(uint32_t s[8], uint32_t w, uint32_t k) {
    uint32_t T1 = s[7]+Sigma1(s[4])+Ch(s[4],s[5],s[6])+k+w;
    uint32_t T2 = Sigma0(s[0])+Maj(s[0],s[1],s[2]);
    s[7]=s[6]; s[6]=s[5]; s[5]=s[4]; s[4]=s[3]+T1;
    s[3]=s[2]; s[2]=s[1]; s[1]=s[0]; s[0]=T1+T2;
}

static void expand(const uint32_t w16[16], uint32_t w64[64]) {
    for(int i=0;i<16;i++) w64[i]=w16[i];
    for(int t=16;t<64;t++)
        w64[t]=sigma1(w64[t-2])+w64[t-7]+sigma0(w64[t-15])+w64[t-16];
}

/* Compute n-round reduced hash (with feedforward) */
static void hash_n(const uint32_t w16[16], int n, uint32_t out[8]) {
    uint32_t w64[64], s[8];
    expand(w16, w64);
    for(int i=0;i<8;i++) s[i]=IV[i];
    for(int t=0;t<n;t++) sha_round(s, w64[t], K[t]);
    for(int i=0;i<8;i++) out[i]=IV[i]+s[i];
}

int main(void) {
    rng = (uint64_t)time(NULL) ^ 0xFEEDFACEULL;

    printf("Step 17f: Reduced-round SHA-256 collision search\n\n");

    /* Try each DW pattern and multiple start positions */
    struct {
        const char *name;
        int dw_pos[5];   /* positions with nonzero DW */
        int32_t dw_val[5]; /* DW values (signed) */
        int n_dw;        /* number of nonzero DWs */
    } patterns[] = {
        {"DW[i]=+1", {0}, {1}, 1},
        {"DW[i]=+1, DW[i+1]=-1", {0,1}, {1,-1}, 2},
        {"DW[i]=+1, DW[i+4]=-1", {0,4}, {1,-1}, 2},
        {"DW[i]=+1, DW[i+1]=-1, DW[i+8]=-1", {0,1,8}, {1,-1,-1}, 3},
        {"DW[i]=+1, DW[i+1]=-1, DW[i+4]=+1, DW[i+8]=-1", {0,1,4,8}, {1,-1,1,-1}, 4},
    };
    int n_patterns = 5;

    uint64_t total_trials = 50000000ULL;  /* 50M per config */

    for (int pi = 0; pi < n_patterns; pi++) {
        for (int start = 0; start <= 7; start++) {
            /* Check pattern fits in 16 words */
            int max_pos = start;
            for(int k=0;k<patterns[pi].n_dw;k++) {
                int p = start + patterns[pi].dw_pos[k];
                if(p > max_pos) max_pos = p;
            }
            if(max_pos >= 16) continue;

            int n_rounds = max_pos + 13; /* local collision + 4 rounds */
            if(n_rounds > 24) n_rounds = 24;

            int best_hw = 999;
            uint64_t collisions = 0;

            for(uint64_t trial = 0; trial < total_trials; trial++) {
                uint32_t W[16], W2[16];
                for(int i=0;i<16;i++) { W[i]=rand32(); W2[i]=W[i]; }

                /* Apply DW pattern */
                for(int k=0;k<patterns[pi].n_dw;k++) {
                    int p = start + patterns[pi].dw_pos[k];
                    W2[p] = W[p] + (uint32_t)patterns[pi].dw_val[k];
                }

                uint32_t h1[8], h2[8];
                hash_n(W, n_rounds, h1);
                hash_n(W2, n_rounds, h2);

                int thw = 0;
                for(int i=0;i<8;i++) thw += hw(h1[i]^h2[i]);

                if(thw < best_hw) {
                    best_hw = thw;
                    if(thw == 0) {
                        collisions++;
                        if(collisions <= 3) {
                            printf("*** %d-ROUND COLLISION! Pattern=%s start=%d trial=%lu\n",
                                   n_rounds, patterns[pi].name, start, (unsigned long)trial);
                            printf("  W  = {");
                            for(int i=0;i<16;i++) printf("0x%08x%s",W[i],i<15?",":"}\n");
                            printf("  W' = {");
                            for(int i=0;i<16;i++) printf("0x%08x%s",W2[i],i<15?",":"}\n");
                            printf("  H  = {");
                            for(int i=0;i<8;i++) printf("0x%08x%s",h1[i],i<7?",":"}\n");
                        }
                    }
                }

                if(trial % 10000000 == 0 && trial > 0) {
                    printf("  [%s s=%d %dr] %luM: best=%d coll=%lu\n",
                           patterns[pi].name, start, n_rounds,
                           (unsigned long)(trial/1000000), best_hw,
                           (unsigned long)collisions);
                }
            }

            printf("  RESULT: %s start=%d %d-rounds: best_hw=%d collisions=%lu\n",
                   patterns[pi].name, start, n_rounds, best_hw,
                   (unsigned long)collisions);
        }
    }

    return 0;
}
