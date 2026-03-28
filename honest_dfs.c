#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/*
 * HONEST TEST: DFS preimage for MULTI-BIT target on SHA-256.
 *
 * Previous: 1-bit target → n nodes (trivial: balanced function).
 * Now: k-bit target → how many nodes?
 *
 * If k=1:  nodes ≈ n (confirmed)
 * If k=8:  nodes ≈ ? (need 8 specific output bits to match)
 * If k=32: nodes ≈ ? (full word match = real preimage)
 *
 * For COLLISION: need H(M1) = H(M2) for M1≠M2.
 * Start with partial: H[4](M1) = H[4](M2) (32-bit match).
 */

#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))
#define SHR(x,n)  ((x)>>(n))
#define CH(e,f,g) (((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))
#define SIG0(x) (ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define SIG1(x) (ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define sig0(x) (ROTR(x,7)^ROTR(x,18)^SHR(x,3))
#define sig1(x) (ROTR(x,17)^ROTR(x,19)^SHR(x,10))

static const uint32_t H0[8] = {
    0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
    0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19
};
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

static inline void sha256_compute(uint32_t w0, uint32_t out[8]) {
    uint32_t W[64];
    W[0] = w0;
    for (int i = 1; i < 16; i++) W[i] = 0;
    for (int i = 16; i < 64; i++)
        W[i] = sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16];
    uint32_t a=H0[0],b=H0[1],c=H0[2],d=H0[3],
             e=H0[4],f=H0[5],g=H0[6],h=H0[7];
    for (int r = 0; r < 64; r++) {
        uint32_t T1 = h + SIG1(e) + CH(e,f,g) + K[r] + W[r];
        uint32_t T2 = SIG0(a) + MAJ(a,b,c);
        h=g; g=f; f=e; e=d+T1; d=c; c=b; b=a; a=T1+T2;
    }
    out[0]=a+H0[0]; out[1]=b+H0[1]; out[2]=c+H0[2]; out[3]=d+H0[3];
    out[4]=e+H0[4]; out[5]=f+H0[5]; out[6]=g+H0[6]; out[7]=h+H0[7];
}

/* Check if k lowest bits of out[target_word] match target_val */
static inline int check_match(uint32_t w0, int tw, uint32_t target_val, uint32_t mask) {
    uint32_t out[8];
    sha256_compute(w0, out);
    return (out[tw] & mask) == target_val;
}

#define PROP_DEPTH 8

static uint64_t g_nodes;
static uint64_t g_max_nodes;

int dfs_multibit(int n, int depth, uint32_t partial, int tw,
                 uint32_t target_val, uint32_t mask) {
    g_nodes++;
    if (g_nodes > g_max_nodes) return -1;

    if (depth == n) {
        return check_match(partial, tw, target_val, mask);
    }

    int n_free = n - depth - 1;

    for (int v = 0; v <= 1; v++) {
        uint32_t new_partial = partial | ((uint32_t)v << depth);

        if (n_free < PROP_DEPTH) {
            uint32_t n_rem = 1U << n_free;
            int first_match = -1;
            int all_same = 1;

            for (uint32_t rem = 0; rem < n_rem; rem++) {
                uint32_t full = new_partial | (rem << (depth + 1));
                int m = check_match(full, tw, target_val, mask);
                if (first_match == -1) {
                    first_match = m;
                } else if (m != first_match) {
                    all_same = 0;
                    break;
                }
            }

            if (all_same) {
                if (first_match == 1) return 1;
                continue; /* prune */
            }
        }

        int result = dfs_multibit(n, depth+1, new_partial, tw, target_val, mask);
        if (result == 1) return 1;
        if (result == -1) return -1;
    }
    return 0;
}

int main() {
    int n = 32; /* full word W[0] */
    int tw = 4; /* target word H[4] */
    int n_trials = 5;

    printf("HONEST MULTI-BIT DFS TEST on SHA-256\n");
    printf("n=%d bits of W[0], target word H[%d]\n\n", n, tw);

    srand(42);

    /* Test for k = 1, 2, 4, 8, 12, 16, 20, 24, 28, 32 target bits */
    int k_values[] = {1, 2, 4, 8, 12, 16, 20, 24, 28, 32};
    int n_k = 10;

    printf("  k bits | mask       |  mean nodes |     brute |      eps | speedup      | time\n");
    printf("  -------+------------+-------------+-----------+----------+--------------+------\n");

    for (int ki = 0; ki < n_k; ki++) {
        int k = k_values[ki];
        uint32_t mask = (k == 32) ? 0xFFFFFFFF : ((1U << k) - 1);

        uint64_t total_nodes = 0;
        int timeouts = 0;
        clock_t t0 = clock();

        g_max_nodes = (k <= 16) ? 100000000ULL : 10000000ULL;
        if (k >= 28) g_max_nodes = 1000000ULL;

        for (int trial = 0; trial < n_trials; trial++) {
            uint32_t target_input = (uint32_t)rand();
            uint32_t out[8];
            sha256_compute(target_input, out);
            uint32_t target_val = out[tw] & mask;

            g_nodes = 0;
            int result = dfs_multibit(n, 0, 0, tw, target_val, mask);

            if (result == -1) {
                timeouts++;
                total_nodes += g_max_nodes;
            } else {
                total_nodes += g_nodes;
            }
        }

        double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
        double mean = (double)total_nodes / n_trials;
        double eps = (mean > 1) ? (1.0 - log2(mean) / n) : 1.0;
        uint64_t brute_equiv = (uint64_t)1 << (k < 32 ? k : 31);
        /* For k-bit target: brute = 2^n / 2^(32-k) = 2^(n-32+k) preimages expected */
        /* But DFS searches 2^n space for any match */
        double speedup = (double)(1ULL << n) / mean;

        printf("  %5d  | 0x%08x | %11.0f | 2^%-7d | %8.4f | %12.1f | %.2fs",
               k, mask, mean, n, eps, speedup, elapsed);
        if (timeouts > 0) printf(" (%d TO)", timeouts);
        if (eps > 0.5) printf(" ★");
        printf("\n");
        fflush(stdout);
    }

    printf("\nINTERPRETATION:\n");
    printf("  k=1:  finding W[0] with H[4][b0]=target → trivial (balanced)\n");
    printf("  k=32: finding W[0] with H[4]=target     → preimage of 1 word\n");
    printf("  If nodes grow as 2^k → DFS = brute force (eps=0)\n");
    printf("  If nodes grow as n×k → DFS beats brute force (eps>0)\n");

    return 0;
}
