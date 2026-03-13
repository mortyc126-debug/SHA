/*
 * p_adic_tower.cu
 *
 * GPU experiment: measure the p-adic tower height of SHA-256.
 *
 * Based on П-47..П-53 from methodology_v15.md.
 *
 * Background (T_SHA_PADIC_TOWER):
 *   Define F: (Z/2^32)^16 -> (Z/2^32)^15 as the differential system:
 *     F(DW)[k] = De_{k+2}(DW)  for k=1..15  (additive differences De3..De17)
 *
 *   Sol_k = {DW in (Z/2^k)^16 : F(DW) ≡ 0 mod 2^k}
 *
 *   height_2(SHA-256) = sup{k : Sol_k ≠ ∅}
 *
 * Known results (from methodology):
 *   P(Sol_1 ≠ ∅) ≈ 9%  (vs 63% for random function)
 *   Sol_2 found via greedy cascade from Sol_1: P ≈ 75%
 *   height_2 ≥ 11  (П-53 result — T_CASCADE_UNIQUENESS)
 *
 * This experiment:
 *   1. For N random seeds, check if DW is in Sol_1 (P≈9%)
 *   2. For Sol_1 seeds, greedily lift to Sol_2, Sol_3, ...
 *   3. Report distribution of maximum height achieved
 *   4. Test hypothesis T_INFINITE_TOWER: height_2 = ∞ ?
 *
 * The greedy cascade (mod 2^k):
 *   For each coordinate i=0..15: try both DW_i and DW_i + 2^{k-1}
 *   Choose whichever gives F(DW) ≡ 0 mod 2^k for component i.
 *   (Greedy = one pass, not backtracking.)
 *
 * Build:
 *   nvcc -O3 -arch=sm_80 -o p_adic_tower p_adic_tower.cu
 *
 * Usage:
 *   ./p_adic_tower [num_seeds] [max_height] [out_file]
 *   ./p_adic_tower 1000000 20 tower_results.txt
 */

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ============================================================
 * SHA-256 CONSTANTS (same as birthday_search_17.cu)
 * ============================================================ */
static const uint32_t HOST_K[64] = {
    0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,
    0x3956c25bu,0x59f111f1u,0x923f82a4u,0xab1c5ed5u,
    0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,
    0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,
    0xe49b69c1u,0xefbe4786u,0x0fc19dc6u,0x240ca1ccu,
    0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
    0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,
    0xc6e00bf3u,0xd5a79147u,0x06ca6351u,0x14292967u,
    0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,
    0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,
    0xa2bfe8a1u,0xa81a664bu,0xc24b8b70u,0xc76c51a3u,
    0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
    0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,
    0x391c0cb3u,0x4ed8aa4au,0x5b9cca4fu,0x682e6ff3u,
    0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,
    0x90beffaau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u
};

#define IV_A 0x6a09e667u
#define IV_B 0xbb67ae85u
#define IV_C 0x3c6ef372u
#define IV_D 0xa54ff53au
#define IV_E 0x510e527fu
#define IV_F 0x9b05688cu
#define IV_G 0x1f83d9abu
#define IV_H 0x5be0cd19u

#define ROTR32(x,n)  (((x)>>(n))|((x)<<(32-(n))))
#define SIG0(x)  (ROTR32(x,2)^ROTR32(x,13)^ROTR32(x,22))
#define SIG1(x)  (ROTR32(x,6)^ROTR32(x,11)^ROTR32(x,25))
#define sig0(x)  (ROTR32(x,7)^ROTR32(x,18)^((x)>>3))
#define sig1(x)  (ROTR32(x,17)^ROTR32(x,19)^((x)>>10))
#define CH(e,f,g)  (((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))

#define SHA_ROUND(a,b,c,d,e,f,g,h,w,k) do { \
    uint32_t _T1=(h)+SIG1(e)+CH(e,f,g)+(k)+(w); \
    uint32_t _T2=SIG0(a)+MAJ(a,b,c); \
    (h)=(g);(g)=(f);(f)=(e);(e)=(d)+_T1; \
    (d)=(c);(c)=(b);(b)=(a);(a)=_T1+_T2; \
} while(0)

__constant__ uint32_t d_K[64];

/* ============================================================
 * DIFFERENTIAL EVALUATION
 *
 * Computes De3..De17 for given message pair (W, W+DW).
 * W[0..15] is fixed (all zeros except W[0] from seed).
 * DW[0..15] is the difference vector.
 *
 * Returns: bit mask indicating which De_r ≡ 0 mod 2^bits
 *          bit r-3 (r=3..17) = 1 if De_r ≡ 0 mod 2^bits, else 0
 * ============================================================ */
__device__ uint32_t eval_differential_mod(
        uint32_t W0_seed,    /* W[0] for base message */
        uint32_t DW[16],     /* difference vector */
        int      bits)       /* check mod 2^bits */
{
    /* Build W[0..15] = {W0_seed, 0, 0, ..., 0} */
    uint32_t mask = (bits >= 32) ? 0xFFFFFFFFu : ((1u << bits) - 1u);

    /* SHA state for W and W' */
    uint32_t a1=IV_A,b1=IV_B,c1=IV_C,d1=IV_D,e1=IV_E,f1=IV_F,g1=IV_G,h1=IV_H;
    uint32_t a2=IV_A,b2=IV_B,c2=IV_C,d2=IV_D,e2=IV_E,f2=IV_F,g2=IV_G,h2=IV_H;

    /* Round 1: W[0]=W0_seed, W'[0]=W0_seed+DW[0] */
    SHA_ROUND(a1,b1,c1,d1,e1,f1,g1,h1, W0_seed,         d_K[0]);
    SHA_ROUND(a2,b2,c2,d2,e2,f2,g2,h2, W0_seed+DW[0],   d_K[0]);

    /* Round 2: W[1]=W'[1]=0 +DW[1] */
    SHA_ROUND(a1,b1,c1,d1,e1,f1,g1,h1, 0u,       d_K[1]);
    SHA_ROUND(a2,b2,c2,d2,e2,f2,g2,h2, DW[1],    d_K[1]);

    uint32_t result = 0;

    #pragma unroll
    for (int r = 3; r <= 17; r++) {
        /* W[r-1] = 0 for r=3..16; W[16] from schedule for r=17 */
        uint32_t w1_r, w2_r;
        if (r <= 16) {
            w1_r = 0u;
            w2_r = DW[r-1];
        } else {
            /* W[16] = sig1(W[14])+W[9]+sig0(W[1])+W[0] = sig0(0)+W0_seed = W0_seed */
            /* W'[16] = sig1(DW[14])+DW[9]+sig0(DW[1])+(W0_seed+DW[0]) */
            w1_r = W0_seed;  /* sig0(0)=0, sig1(0)=0 */
            w2_r = (W0_seed + DW[0]) +
                   sig0(DW[1]) + DW[9] + sig1(DW[14]);
        }

        uint32_t T1_1 = h1 + SIG1(e1) + CH(e1,f1,g1) + d_K[r-1] + w1_r;
        uint32_t T2_1 = SIG0(a1) + MAJ(a1,b1,c1);
        h1=g1;g1=f1;f1=e1;e1=d1+T1_1;d1=c1;c1=b1;b1=a1;a1=T1_1+T2_1;

        uint32_t T1_2 = h2 + SIG1(e2) + CH(e2,f2,g2) + d_K[r-1] + w2_r;
        uint32_t T2_2 = SIG0(a2) + MAJ(a2,b2,c2);
        h2=g2;g2=f2;f2=e2;e2=d2+T1_2;d2=c2;c2=b2;b2=a2;a2=T1_2+T2_2;

        /* De_r = e2 - e1. Check mod 2^bits */
        uint32_t De_r = e2 - e1;
        if ((De_r & mask) == 0) result |= (1u << (r-3));
    }
    return result;
}

/* ============================================================
 * KERNEL: Test p-adic height for each seed
 *
 * Each thread:
 *   1. Uses seed to derive initial DW (1-bit per coordinate from seed bits)
 *   2. Checks Sol_1: F(DW) ≡ 0 mod 2 for all 15 components
 *   3. If in Sol_1: greedily lifts to Sol_2, Sol_3, ..., Sol_{max_k}
 *   4. Records maximum height achieved
 * ============================================================ */
#define MAX_HEIGHT 32

__global__ void p_adic_tower_kernel(
        uint64_t  n_seeds,
        uint64_t  seed_base,
        uint32_t  max_k,
        uint32_t* d_height_counts,   /* histogram: [0..max_k] */
        uint32_t* d_sol1_count,
        uint64_t* d_best_seed)       /* seed achieving max height */
{
    uint64_t idx = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_seeds) return;

    /* Derive seed */
    uint64_t seed = seed_base + idx;

    /* --- Generate initial DW from seed (1-bit per 16 coords) ---
     * DW[i] = low bit of (seed >> i)  for i=0..15
     * This gives DW in {0,1}^16 — the Sol_1 search space.
     */
    uint32_t DW[16];
    uint32_t W0_seed = (uint32_t)(seed >> 32) ^ (uint32_t)seed;
    for (int i = 0; i < 16; i++)
        DW[i] = (uint32_t)((seed >> i) & 1u);

    /* --- Check Sol_1: all 15 components of F(DW) even --- */
    uint32_t bits_set = eval_differential_mod(W0_seed, DW, 1);
    /* bits_set: bit k-3 = 1 if De_k ≡ 0 mod 2 */
    /* We need De3..De17 all 0 mod 2 → bits 0..14 all set */
    if ((bits_set & 0x7FFFu) != 0x7FFFu) {
        /* Not in Sol_1 — no height */
        atomicAdd(&d_height_counts[0], 1u);
        return;
    }

    /* In Sol_1! */
    atomicAdd(d_sol1_count, 1u);

    /* --- Greedy lifting from Sol_1 to Sol_{max_k} ---
     * At each step k (lifting to Sol_{k+1}):
     *   For each coordinate i=0..15:
     *     Try DW[i] and DW[i] + 2^k
     *     Keep whichever gives all 15 components of F(DW) ≡ 0 mod 2^{k+1}
     */
    int height = 1; /* already in Sol_1 */
    for (uint32_t k = 1; k <= max_k && k < MAX_HEIGHT; k++) {
        uint32_t DW_try[16];
        uint32_t bit_k = 1u << (k - 1);  /* 2^{k-1}: the bit to try flipping */

        /* Make a copy */
        for (int i = 0; i < 16; i++) DW_try[i] = DW[i];

        /* Greedy: for each coordinate, choose the value that satisfies mod 2^{k+1} */
        int lifted = 1;
        for (int i = 0; i < 16 && lifted; i++) {
            /* Try DW[i] as-is */
            uint32_t r0 = eval_differential_mod(W0_seed, DW_try, k+1);
            if ((r0 & 0x7FFFu) == 0x7FFFu) continue; /* coordinate already OK */

            /* Try DW[i] + bit_k */
            DW_try[i] ^= bit_k;
            uint32_t r1 = eval_differential_mod(W0_seed, DW_try, k+1);
            if ((r1 & 0x7FFFu) == 0x7FFFu) continue; /* flipped version works */

            /* Neither works: greedy fails at this coordinate */
            DW_try[i] ^= bit_k; /* restore */
            lifted = 0;
        }

        if (!lifted) break; /* can't lift to Sol_{k+1} */

        /* Successfully lifted */
        height = k + 1;
        for (int i = 0; i < 16; i++) DW[i] = DW_try[i];

        /* Track best seed */
        if (height >= (int)max_k) {
            atomicMax((unsigned long long*)d_best_seed, (unsigned long long)seed);
        }
    }

    atomicAdd(&d_height_counts[height], 1u);
}

/* ============================================================
 * MAIN
 * ============================================================ */
int main(int argc, char* argv[]) {
    uint64_t n_seeds   = 1000000ULL;  /* default: 1M seeds */
    uint32_t max_k     = 20;           /* default: test up to Sol_20 */
    const char* outfile = "tower_results.txt";

    if (argc > 1) n_seeds  = (uint64_t)strtoull(argv[1], NULL, 10);
    if (argc > 2) max_k    = (uint32_t)atoi(argv[2]);
    if (argc > 3) outfile  = argv[3];

    printf("=== SHA-256 p-Adic Tower Height Test ===\n");
    printf("Seeds     = %llu\n", (unsigned long long)n_seeds);
    printf("Max k     = %u\n", max_k);
    printf("Output    = %s\n\n", outfile);

    /* GPU setup */
    cudaMemcpyToSymbol(d_K, HOST_K, sizeof(HOST_K));

    uint32_t* d_height_counts = NULL;
    uint32_t* d_sol1_count    = NULL;
    uint64_t* d_best_seed     = NULL;
    cudaMalloc(&d_height_counts, (MAX_HEIGHT+1) * sizeof(uint32_t));
    cudaMalloc(&d_sol1_count, sizeof(uint32_t));
    cudaMalloc(&d_best_seed, sizeof(uint64_t));
    cudaMemset(d_height_counts, 0, (MAX_HEIGHT+1) * sizeof(uint32_t));
    cudaMemset(d_sol1_count, 0, sizeof(uint32_t));
    cudaMemset(d_best_seed, 0, sizeof(uint64_t));

    const int    BLOCK_SIZE = 256;
    const uint64_t CHUNK    = (1ULL << 22);  /* 4M seeds per launch */
    uint64_t seed_base      = (uint64_t)time(NULL) << 20;

    clock_t t0 = clock();

    for (uint64_t start = 0; start < n_seeds; start += CHUNK) {
        uint64_t n = (start + CHUNK <= n_seeds) ? CHUNK : (n_seeds - start);
        uint64_t blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
        p_adic_tower_kernel<<<(uint32_t)blocks, BLOCK_SIZE>>>(
                n, seed_base + start, max_k,
                d_height_counts, d_sol1_count, d_best_seed);
        cudaDeviceSynchronize();
        printf("  Progress: %.1f%%\r", 100.0*start/n_seeds);
        fflush(stdout);
    }
    printf("\n");

    double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;

    /* Copy results */
    uint32_t h_counts[MAX_HEIGHT+1];
    uint32_t h_sol1;
    uint64_t h_best;
    cudaMemcpy(h_counts, d_height_counts, (MAX_HEIGHT+1)*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_sol1, d_sol1_count, sizeof(uint32_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_best, d_best_seed, sizeof(uint64_t), cudaMemcpyDeviceToHost);

    /* Report */
    printf("\n=== RESULTS (%llu seeds, %.2f s) ===\n",
           (unsigned long long)n_seeds, elapsed);
    printf("Seeds in Sol_1: %u / %llu  (%.2f%%)\n",
           h_sol1, (unsigned long long)n_seeds,
           100.0*h_sol1/n_seeds);
    printf("  (Expected ~9%% for SHA-256, ~63%% for random function)\n\n");

    printf("Height distribution:\n");
    printf("  k | count   | P(height=k)\n");
    printf("  --+---------+------------\n");
    for (int k = 0; k <= (int)max_k && k <= MAX_HEIGHT; k++) {
        if (h_counts[k] > 0) {
            printf("  %2d| %7u | %.4f%%\n", k, h_counts[k],
                   100.0*h_counts[k]/n_seeds);
        }
    }

    /* Check for infinite tower */
    int max_observed = 0;
    for (int k = MAX_HEIGHT; k >= 0; k--) {
        if (h_counts[k] > 0) { max_observed = k; break; }
    }
    printf("\nMax observed height: %d\n", max_observed);
    if (max_observed >= (int)max_k) {
        printf("HYPOTHESIS: height_2 >= %u (reached test limit — try higher max_k)\n", max_k);
        printf("Best seed achieving max height: 0x%016llx\n", (unsigned long long)h_best);
    }

    /* Save to file */
    FILE* fp = fopen(outfile, "w");
    if (fp) {
        fprintf(fp, "# SHA-256 p-adic tower heights\n");
        fprintf(fp, "# Seeds: %llu, Time: %.2f s\n",
                (unsigned long long)n_seeds, elapsed);
        fprintf(fp, "# Sol_1 rate: %.4f%% (expected ~9%%)\n",
                100.0*h_sol1/n_seeds);
        fprintf(fp, "# height,count,probability\n");
        for (int k = 0; k <= (int)max_k && k <= MAX_HEIGHT; k++) {
            fprintf(fp, "%d,%u,%.6f\n", k, h_counts[k],
                    (double)h_counts[k]/n_seeds);
        }
        fprintf(fp, "# max_height_seed=0x%016llx\n", (unsigned long long)h_best);
        fclose(fp);
    }
    printf("\nResults saved to: %s\n", outfile);

    cudaFree(d_height_counts);
    cudaFree(d_sol1_count);
    cudaFree(d_best_seed);
    return 0;
}
