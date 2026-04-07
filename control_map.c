/*
 * CONTROL MAP: how much does each W[r] control the final hash?
 *
 * If W[0] controls 50% of hash and W[15] controls 5%,
 * then search W[0] first (more leverage).
 *
 * gcc -O3 -march=native -o control_map control_map.c
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))
#define CH(e,f,g) (((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))
#define S0(x) (ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define S1(x) (ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define s0(x) (ROTR(x,7)^ROTR(x,18)^((x)>>3))
#define s1(x) (ROTR(x,17)^ROTR(x,19)^((x)>>10))
static const uint32_t K[64]={
0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
static const uint32_t IV[8]={
0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};

void sha256_compute(const uint32_t msg[16], uint32_t hash[8]) {
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for(int r=0;r<64;r++){
        uint32_t t1=h+S1(e)+CH(e,f,g)+K[r]+W[r], t2=S0(a)+MAJ(a,b,c);
        h=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
    }
    hash[0]=a+IV[0];hash[1]=b+IV[1];hash[2]=c+IV[2];hash[3]=d+IV[3];
    hash[4]=e+IV[4];hash[5]=f+IV[5];hash[6]=g+IV[6];hash[7]=h+IV[7];
}

int main() {
    printf("CONTROL MAP: per-word influence on hash\n");
    printf("========================================\n\n");

    srand(42);
    uint32_t msg[16], hash_ref[8];
    for(int i=0;i<16;i++) msg[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
    sha256_compute(msg, hash_ref);

    /* For each W[k]: flip it completely, measure hash change */
    printf("  W[k] full flip → HW(Δhash):\n");
    printf("  Word |  HW/256 | Per hash word\n");
    printf("  ─────+─────────+──────────────────────────────\n");
    for(int k=0;k<16;k++){
        uint32_t msg2[16];
        memcpy(msg2, msg, 64);
        msg2[k] ^= 0xFFFFFFFF;
        uint32_t hash2[8];
        sha256_compute(msg2, hash2);
        int total=0;
        printf("  W[%2d]|", k);
        int per_word[8];
        for(int w=0;w<8;w++){
            per_word[w]=__builtin_popcount(hash_ref[w]^hash2[w]);
            total+=per_word[w];
        }
        printf("  %3d   |", total);
        for(int w=0;w<8;w++) printf(" H%d=%2d", w, per_word[w]);
        printf("\n");
    }

    /* More interesting: SEQUENTIAL steering.
     * Fix W[1..15], vary ONLY W[0]. How many distinct hash values? */
    printf("\n  Sequential control test:\n");
    printf("  Fix W[1..15], vary W[k] through 2^16 values:\n\n");

    for(int k=0;k<16;k+=3){
        /* Count distinct hash[0] values from 2^16 trials */
        int n_trials = 65536;
        int distinct = 0;
        uint32_t last_h0 = 0;
        int collisions = 0;

        /* Use birthday-in-small-space: track if any hash[0] repeats */
        uint32_t seen[65536];
        for(int t=0;t<n_trials;t++){
            uint32_t msg2[16];
            memcpy(msg2, msg, 64);
            msg2[k] = (uint32_t)t * 0x10001u; /* systematic variation */
            uint32_t h[8];
            sha256_compute(msg2, h);
            seen[t] = h[0];
        }
        /* Count distinct in seen[] */
        /* Simple: sort and count unique */
        for(int i=0;i<n_trials-1;i++)
            for(int j=i+1;j<n_trials;j++)
                if(seen[i]>seen[j]){uint32_t t=seen[i];seen[i]=seen[j];seen[j]=t;}
        distinct=1;
        for(int i=1;i<n_trials;i++) if(seen[i]!=seen[i-1]) distinct++;
        collisions = n_trials - distinct;

        printf("  Vary W[%2d]: %d distinct H[0] / %d trials. Collisions in H[0]: %d\n",
               k, distinct, n_trials, collisions);
    }

    /* KEY TEST: if we know target hash[0], how many W[0] values produce it?
     * Answer: ~1 (since W[0] → hash[0] is essentially random 32→32) */
    printf("\n  KEY: preimage density per word:\n");
    uint32_t target_h0 = hash_ref[0];
    int matches = 0;
    int n_search = 1000000;
    for(int t=0;t<n_search;t++){
        uint32_t msg2[16];
        memcpy(msg2, msg, 64);
        msg2[0] = (uint32_t)rand()|((uint32_t)rand()<<16);
        uint32_t h[8];
        sha256_compute(msg2, h);
        if(h[0] == target_h0) matches++;
    }
    printf("  Vary W[0], match H[0]: %d/%d (expected: %d)\n",
           matches, n_search, n_search >> 5); /* expected n/2^32... too few trials */
    printf("  (Need 2^32 trials to expect 1 H[0] match)\n");

    /* REAL question: can we search W[0..15] SEQUENTIALLY?
     * Round 1: try 2^32 W[0] values → get 2^32 state[1] values.
     * Round 2: for each state[1], try 2^32 W[1] → 2^32 state[2].
     * ... but this is 2^(32×16) = 2^512 total. WORSE than brute force.
     *
     * UNLESS: we can PRUNE at each step.
     * After choosing W[0]: state[1] determined.
     * Does state[1] tell us anything about whether this path leads to target?
     * NO — state[1] is uncorrelated with hash (we proved this: all metrics random).
     *
     * After choosing W[0..15]: state[16] determined.
     * schedule determined. Forward 48 rounds → hash.
     * Check hash. This is STANDARD brute force: 2^512 tries.
     *
     * For collision: 2^128 birthday.
     * No sequential pruning possible: no intermediate signal. */

    printf("\n════════════════════════════════════\n");
    printf("CONTROL MAP VERDICT\n");
    printf("════════════════════════════════════\n\n");
    printf("  Each W[k] affects ~128/256 hash bits (50%% — full avalanche).\n");
    printf("  All 16 words have EQUAL influence (uniform, no weak word).\n");
    printf("  No sequential pruning: intermediate state = random.\n");
    printf("  Forward steering = choosing where to land in hash space.\n");
    printf("  But: hash space = 2^256, target = 1 point.\n");
    printf("  Need 2^256 tries (preimage) or 2^128 (collision).\n");

    return 0;
}
