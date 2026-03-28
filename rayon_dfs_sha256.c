#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

/*
 * RAYON-style DFS with constant propagation on SHA-256
 * Tests: does ε > 0 hold at n=32 (full word)?
 *
 * DFS: fix bits of W[0] one by one (LSB first).
 * At each level: try 0 and 1.
 * Propagation: when few bits remain free (≤PROP_DEPTH),
 * exhaustively check if output bit is determined.
 * If determined and matches target → found.
 * If determined and doesn't match → prune.
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

/* SHA-256 output bit for W[0]=w0, W[1..15]=0 */
static inline int sha256_bit(uint32_t w0, int word, int bit) {
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
    uint32_t out[8] = {a+H0[0],b+H0[1],c+H0[2],d+H0[3],
                       e+H0[4],f+H0[5],g+H0[6],h+H0[7]};
    return (out[word] >> bit) & 1;
}

#define PROP_DEPTH 8  /* exhaustive check when ≤8 bits free */

static uint64_t g_nodes;
static uint64_t g_max_nodes;

/* DFS with propagation */
int dfs(int n, int depth, uint32_t partial, int target_val, int tw, int tb) {
    g_nodes++;
    if (g_nodes > g_max_nodes) return -1; /* timeout */

    if (depth == n) {
        return sha256_bit(partial, tw, tb) == target_val;
    }

    int n_free = n - depth - 1;

    for (int v = 0; v <= 1; v++) {
        uint32_t new_partial = partial | ((uint32_t)v << depth);

        /* Propagation: if few bits free, check exhaustively */
        if (n_free < PROP_DEPTH) {
            int first_val = -1;
            int determined = 1;
            uint32_t n_rem = 1U << n_free;

            for (uint32_t rem = 0; rem < n_rem; rem++) {
                uint32_t full = new_partial | (rem << (depth + 1));
                int oval = sha256_bit(full, tw, tb);
                if (first_val == -1) {
                    first_val = oval;
                } else if (oval != first_val) {
                    determined = 0;
                    break;
                }
            }

            if (determined) {
                if (first_val == target_val) return 1; /* found */
                continue; /* prune: determined but wrong */
            }
        }

        int result = dfs(n, depth + 1, new_partial, target_val, tw, tb);
        if (result == 1) return 1;
        if (result == -1) return -1;
    }
    return 0;
}

int main(int argc, char *argv[]) {
    int n = 24; /* default: 24 bits of W[0] */
    int n_trials = 10;
    int tw = 4, tb = 0; /* target: H[4][b0] */

    if (argc > 1) n = atoi(argv[1]);
    if (argc > 2) n_trials = atoi(argv[2]);

    uint64_t brute = 1ULL << n;
    g_max_nodes = brute < 100000000ULL ? brute : 100000000ULL;

    printf("RAYON DFS on SHA-256 — n=%d bits of W[0]\n", n);
    printf("Target: H[%d][b%d], %d trials, brute=2^%d=%lu\n", tw, tb, n_trials, n, brute);
    printf("Propagation depth: %d bits\n\n", PROP_DEPTH);

    srand(42);
    uint64_t total_nodes = 0;
    int timeouts = 0;

    for (int trial = 0; trial < n_trials; trial++) {
        uint32_t target_input = rand() & ((1U << n) - 1);
        int target_val = sha256_bit(target_input, tw, tb);

        g_nodes = 0;
        clock_t t0 = clock();
        int result = dfs(n, 0, 0, target_val, tw, tb);
        double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;

        if (result == -1) {
            timeouts++;
            printf("  Trial %2d: TIMEOUT (>%lu nodes) %.1fs\n", trial, g_max_nodes, elapsed);
            total_nodes += g_max_nodes;
        } else {
            total_nodes += g_nodes;
            printf("  Trial %2d: %8lu nodes (%.1fs) target=%d\n", trial, g_nodes, elapsed, target_val);
        }
    }

    double mean_nodes = (double)total_nodes / n_trials;
    double eps = (mean_nodes > 1) ? (1.0 - log2(mean_nodes) / n) : 1.0;

    printf("\nRESULTS:\n");
    printf("  Mean nodes: %.0f\n", mean_nodes);
    printf("  Brute force: %lu\n", brute);
    printf("  epsilon: %.4f\n", eps);
    printf("  Speedup: %.1fx\n", brute / mean_nodes);
    printf("  Timeouts: %d/%d\n", timeouts, n_trials);

    if (timeouts == 0 && eps > 0.1) {
        printf("\n  *** epsilon > 0 at n=%d! DFS beats brute force! ***\n", n);
    } else if (timeouts > n_trials / 2) {
        printf("\n  epsilon ~ 0 at n=%d. DFS = brute force.\n", n);
    }

    return 0;
}
