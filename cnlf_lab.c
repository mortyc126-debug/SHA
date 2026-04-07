/*
 * CARRY×NLF LABORATORY: study with known W (controlled conditions).
 *
 * With W known: all locks open. We can trace everything.
 * Goal: build ATLAS of self-cancellation patterns, find universals.
 *
 * Experiments:
 * 1. Erasure map: for fixed W, WHERE exactly are bits erased?
 * 2. Quiet point search: find ΔW with minimum Q(W,ΔW)
 * 3. Self-cancellation pattern: is erasure pattern UNIVERSAL across W?
 * 4. Deadpool chain: full recovery trace with known W
 * 5. Newton convergence in quiet points
 *
 * gcc -O3 -march=native -o cnlf_lab cnlf_lab.c -lm
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

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

void sha256_full(const uint32_t msg[16], uint32_t st[65][8], uint32_t W[64]) {
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
    for(int i=0;i<8;i++) st[0][i]=IV[i];
    for(int r=0;r<64;r++){
        uint32_t a=st[r][0],b=st[r][1],c=st[r][2],d=st[r][3];
        uint32_t e=st[r][4],f=st[r][5],g=st[r][6],h=st[r][7];
        uint32_t t1=h+S1(e)+CH(e,f,g)+K[r]+W[r], t2=S0(a)+MAJ(a,b,c);
        st[r+1][0]=t1+t2; st[r+1][1]=a; st[r+1][2]=b; st[r+1][3]=c;
        st[r+1][4]=d+t1;  st[r+1][5]=e; st[r+1][6]=f; st[r+1][7]=g;
    }
}

int hw32(uint32_t x) { return __builtin_popcount(x); }
int hw256(const uint32_t a[8], const uint32_t b[8]) {
    int h=0; for(int i=0;i<8;i++) h+=hw32(a[i]^b[i]); return h;
}

int main() {
    printf("CARRY×NLF LABORATORY (known W)\n");
    printf("===============================\n\n");

    srand(42);
    uint32_t msg[16], st[65][8], W[64];
    for(int i=0;i<16;i++) msg[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
    sha256_full(msg, st, W);

    /* ═══ EXP 1: ERASURE MAP ═══ */
    printf("EXP 1: ERASURE MAP (per round, per bit)\n");
    printf("────────────────────────────────────────\n");

    /* Per round: which bits of e are erased by Ch?
     * Erased = bit position k where f[k]≠g[k]. */
    int total_erased_per_bit[32] = {0};
    int total_erased = 0;

    for(int r=0; r<64; r++) {
        uint32_t f = st[r][5], g = st[r][6];
        uint32_t erase = f ^ g;
        int n = hw32(erase);
        total_erased += n;
        for(int k=0; k<32; k++)
            if((erase>>k)&1) total_erased_per_bit[k]++;
    }

    printf("  Total erasures: %d across 64 rounds (avg %.1f/round)\n",
           total_erased, total_erased/64.0);
    printf("  Per bit position (across 64 rounds):\n    ");
    for(int k=0; k<32; k++) printf("%d", total_erased_per_bit[k]/6);
    printf("\n    (digits = erasure count / 6, scale 0-10)\n");

    /* Is erasure pattern UNIFORM across bit positions? */
    double mean_per_bit = total_erased / 32.0;
    double var = 0;
    for(int k=0; k<32; k++) {
        double d = total_erased_per_bit[k] - mean_per_bit;
        var += d*d;
    }
    var /= 32;
    printf("  Mean per bit: %.1f, std: %.1f → %s\n",
           mean_per_bit, sqrt(var),
           sqrt(var) < mean_per_bit * 0.15 ? "UNIFORM" : "NON-UNIFORM ★");

    /* ═══ EXP 2: QUIET POINT SEARCH ═══ */
    printf("\nEXP 2: QUIET POINT SEARCH (min Q = min HW(N))\n");
    printf("────────────────────────────────────────\n");

    /* Q(W, ΔW) = HW of nonlinear correction = HW(SHA(W⊕ΔW) ⊕ SHA(W) ⊕ J·ΔW)
     * Simplified: Q ≈ HW(SHA(W⊕ΔW) ⊕ SHA(W)) since J·ΔW ≈ SHA diff at first order.
     * Actually Q = HW of the NONLINEAR PART. But for measurement:
     * just use hash diff as proxy. */

    uint32_t hash1[8];
    for(int i=0;i<8;i++) hash1[i] = st[64][i]+IV[i];

    int best_q = 256;
    uint32_t best_dw[16] = {0};

    /* Random search */
    for(int trial=0; trial<500000; trial++) {
        uint32_t dw[16] = {0};
        /* Random sparse ΔW: 1-3 bit flips */
        int n_flips = 1 + (rand()%3);
        for(int f=0; f<n_flips; f++) {
            int w = rand()%16, b = rand()%32;
            dw[w] ^= (1u<<b);
        }

        uint32_t msg2[16], st2[65][8], W2[64];
        for(int i=0;i<16;i++) msg2[i] = msg[i]^dw[i];
        sha256_full(msg2, st2, W2);

        uint32_t hash2[8];
        for(int i=0;i<8;i++) hash2[i] = st2[64][i]+IV[i];

        int q = hw256(hash1, hash2);
        if(q < best_q && q > 0) {
            best_q = q;
            memcpy(best_dw, dw, 64);
        }
    }

    printf("  Random search (500K, 1-3 bit flips): min Q = %d/256\n", best_q);

    /* Hill climbing from best */
    for(int iter=0; iter<100000; iter++) {
        uint32_t dw2[16];
        memcpy(dw2, best_dw, 64);
        int w = rand()%16, b = rand()%32;
        dw2[w] ^= (1u<<b);

        uint32_t msg2[16], st2[65][8], W2[64];
        for(int i=0;i<16;i++) msg2[i] = msg[i]^dw2[i];
        sha256_full(msg2, st2, W2);
        uint32_t hash2[8];
        for(int i=0;i<8;i++) hash2[i] = st2[64][i]+IV[i];
        int q = hw256(hash1, hash2);
        if(q < best_q && q > 0) {
            best_q = q;
            memcpy(best_dw, dw2, 64);
        }
    }
    printf("  After hill climbing (100K): min Q = %d/256\n", best_q);

    /* Analyze best quiet point */
    int dw_hw = 0;
    for(int i=0;i<16;i++) dw_hw += hw32(best_dw[i]);
    printf("  Best ΔW: HW = %d msg bits flipped\n", dw_hw);

    /* ═══ EXP 3: UNIVERSALITY — same erasure pattern across W? ═══ */
    printf("\nEXP 3: UNIVERSALITY (erasure pattern across messages)\n");
    printf("────────────────────────────────────────\n");

    /* For 10 different messages: measure erasure count per round.
     * If same pattern → universal. If different → message-specific. */
    double erasure_per_round[64] = {0};
    int n_msgs = 100;
    for(int m=0; m<n_msgs; m++) {
        uint32_t msg_t[16], st_t[65][8], W_t[64];
        for(int i=0;i<16;i++) msg_t[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
        sha256_full(msg_t, st_t, W_t);
        for(int r=0; r<64; r++) {
            uint32_t erase = st_t[r][5] ^ st_t[r][6]; /* f⊕g */
            erasure_per_round[r] += hw32(erase);
        }
    }

    printf("  Avg erasure per round (100 messages):\n  ");
    int universal = 1;
    for(int r=0; r<64; r++) {
        erasure_per_round[r] /= n_msgs;
        if(r<8 || r>56 || r%16==0)
            printf("r%d=%.1f ", r, erasure_per_round[r]);
        if(fabs(erasure_per_round[r] - 16.0) > 3.0) universal = 0;
    }
    printf("\n  All rounds ≈ 16? %s\n", universal ? "YES — UNIVERSAL ★" : "NO — varies");

    /* ═══ EXP 4: DEADPOOL FULL CHAIN with known W ═══ */
    printf("\nEXP 4: DEADPOOL CHAIN (full recovery with known W)\n");
    printf("────────────────────────────────────────\n");

    /* With W known: can we recover ALL erased bits?
     * For each round r: compute T1 = a[r+1]-T2[r].
     * h[r] = T1 - S1(e[r]) - Ch(e[r],f[r],g[r]) - K[r] - W[r].
     * ALL known → h[r] recovered → erased bits restored. */

    int recovered_total = 0;
    for(int r=0; r<64; r++) {
        uint32_t a_rp1 = st[r+1][0];
        uint32_t a_r = st[r][0], b_r = st[r][1], c_r = st[r][2];
        uint32_t e_r = st[r][4], f_r = st[r][5], g_r = st[r][6];

        uint32_t T2 = S0(a_r) + MAJ(a_r, b_r, c_r);
        uint32_t T1 = a_rp1 - T2;
        uint32_t h_recovered = T1 - S1(e_r) - CH(e_r, f_r, g_r) - K[r] - W[r];
        uint32_t h_actual = st[r][7];

        if(h_recovered == h_actual) recovered_total++;
    }
    printf("  Recovered h[r] correctly: %d/64 rounds\n", recovered_total);
    printf("  → %s\n", recovered_total==64 ? "★ 100%% PERFECT RECOVERY with known W!" : "ERRORS");

    /* ═══ EXP 5: NEWTON IN QUIET POINT ═══ */
    printf("\nEXP 5: NEWTON CONVERGENCE at quiet point\n");
    printf("────────────────────────────────────────\n");

    /* At the quiet point ΔW found above:
     * SHA(W) and SHA(W⊕ΔW) differ by Q bits.
     * Can Newton iteration: ΔW_{k+1} = update(ΔW_k) converge to Q=0?
     *
     * Method: ΔW → compute SHA diff → adjust ΔW to reduce diff.
     * Per-word Newton: for each word i, adjust ΔW[i] to cancel
     * the contribution of word i to hash diff. */

    uint32_t dw_iter[16];
    memcpy(dw_iter, best_dw, 64);

    printf("  Starting from quiet point Q=%d\n", best_q);
    for(int iter=0; iter<20; iter++) {
        uint32_t msg2[16], st2[65][8], W2[64], hash2[8];
        for(int i=0;i<16;i++) msg2[i] = msg[i]^dw_iter[i];
        sha256_full(msg2, st2, W2);
        for(int i=0;i<8;i++) hash2[i] = st2[64][i]+IV[i];
        int q = hw256(hash1, hash2);

        if(q == 0) {
            printf("  ★★★ COLLISION FOUND at iteration %d!\n", iter);
            break;
        }

        /* Per-word adjustment: for first differing hash word,
         * XOR the hash diff into the corresponding ΔW word */
        for(int w=0; w<8; w++) {
            uint32_t hdiff = hash1[w] ^ hash2[w];
            if(hdiff) {
                /* Adjust a message word that primarily affects this hash word */
                /* Hash word w comes from state[64][w].
                 * state[64][0] = a[64] → affected by all W.
                 * Simplest: adjust ΔW[w%16] by hash diff. */
                dw_iter[w % 16] ^= hdiff;
                break;
            }
        }

        if(iter < 10 || iter == 19)
            printf("  Iter %2d: Q=%d\n", iter, q);
    }

    /* Better Newton: XOR-schedule based */
    printf("\n  XOR-Newton (schedule-aware):\n");
    memcpy(dw_iter, best_dw, 64);
    for(int iter=0; iter<50; iter++) {
        uint32_t msg2[16], st2[65][8], W2[64], hash2[8];
        for(int i=0;i<16;i++) msg2[i] = msg[i]^dw_iter[i];
        sha256_full(msg2, st2, W2);
        for(int i=0;i<8;i++) hash2[i] = st2[64][i]+IV[i];
        int q = hw256(hash1, hash2);

        if(q == 0) {
            printf("  ★★★ COLLISION at iteration %d!\n", iter);
            break;
        }

        /* Compute which ΔW word to flip using XOR-schedule sensitivity */
        /* For now: flip the bit in ΔW that reduces Q most (greedy) */
        int best_flip_w = -1, best_flip_b = -1, best_flip_q = q;
        for(int test=0; test<64; test++) { /* sample 64 random flips */
            int fw = rand()%16, fb = rand()%32;
            uint32_t dw_test[16];
            memcpy(dw_test, dw_iter, 64);
            dw_test[fw] ^= (1u<<fb);
            uint32_t m3[16], s3[65][8], w3[64], h3[8];
            for(int i=0;i<16;i++) m3[i]=msg[i]^dw_test[i];
            sha256_full(m3, s3, w3);
            for(int i=0;i<8;i++) h3[i]=s3[64][i]+IV[i];
            int q3 = hw256(hash1, h3);
            if(q3 < best_flip_q) {
                best_flip_q = q3; best_flip_w = fw; best_flip_b = fb;
            }
        }

        if(best_flip_w >= 0) {
            dw_iter[best_flip_w] ^= (1u<<best_flip_b);
            q = best_flip_q;
        }

        if(iter<10 || iter%10==0 || q < 80)
            printf("  Iter %2d: Q=%d%s\n", iter, q, q==0?" ★★★ COLLISION!":"");
        if(q==0) break;
    }

    printf("\n═══════════════════════════════════════\n");
    printf("LABORATORY SUMMARY\n");
    printf("═══════════════════════════════════════\n");

    return 0;
}
