/*
 * birthday_search_17.cu
 *
 * GPU birthday search: find SHA-256 message pairs (W, W') with
 *   De3 = De4 = ... = De17 = 0  (15 consecutive zero differentials)
 *
 * Based on T_CASCADE_17 (П-13) and T_BIRTHDAY_COST17 (П-32)
 * from methodology_v15.md.
 *
 * Algorithm (T_CASCADE_17):
 *   Fix W[0]=W0, W[1]=w1 (iterated), W[2..15]=0.
 *   W'[0] = W[0] + DW0  (additive difference, default DW0=1)
 *   For rounds r=3..16: adaptively set W'[r-1] to cancel De_r.
 *   Check f17 = Da13 + DW16 == 0  (T_DE17_DECOMPOSITION).
 *
 * Expected performance:
 *   A100 (80GB): ~0.04 sec for one pair (E[cost] = 2^32)
 *   RTX 3090   : ~0.3  sec for one pair
 *   RTX 4090   : ~0.15 sec for one pair
 *
 * Known verification pairs (from methodology):
 *   П-15: W0=0xe82222c7, w1=0x516cfb41, DW0=1  -> De3..De17=0
 *   П-16: W0=0xd4254551, w1=0x679ea4de, DW0=1  -> De3..De17=0
 *
 * Build:
 *   nvcc -O3 -arch=sm_80 -o birthday_search_17 birthday_search_17.cu
 *
 * Usage:
 *   ./birthday_search_17 [W0_hex] [DW0_hex] [num_W0_values] [out_file]
 *   ./birthday_search_17                          # defaults: W0=random, DW0=1, n=1, pairs.txt
 *   ./birthday_search_17 e82222c7 1 1 pairs.txt  # reproduces П-15
 */

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ============================================================
 * SHA-256 CONSTANTS
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

/* SHA-256 IV */
#define IV_A 0x6a09e667u
#define IV_B 0xbb67ae85u
#define IV_C 0x3c6ef372u
#define IV_D 0xa54ff53au
#define IV_E 0x510e527fu
#define IV_F 0x9b05688cu
#define IV_G 0x1f83d9abu
#define IV_H 0x5be0cd19u

/* ============================================================
 * SHA-256 MACROS
 * ============================================================ */
#define ROTR32(x,n)  (((x)>>(n))|((x)<<(32-(n))))
#define SIG0(x)      (ROTR32(x,2)  ^ ROTR32(x,13) ^ ROTR32(x,22))
#define SIG1(x)      (ROTR32(x,6)  ^ ROTR32(x,11) ^ ROTR32(x,25))
#define sig0(x)      (ROTR32(x,7)  ^ ROTR32(x,18) ^ ((x)>>3))
#define sig1(x)      (ROTR32(x,17) ^ ROTR32(x,19) ^ ((x)>>10))
#define CH(e,f,g)    (((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c)   (((a)&(b))^((a)&(c))^((b)&(c)))

/* One SHA-256 round, updating (a..h) using word w and constant k */
#define SHA_ROUND(a,b,c,d,e,f,g,h,w,k) do { \
    uint32_t _T1 = (h) + SIG1(e) + CH(e,f,g) + (k) + (w); \
    uint32_t _T2 = SIG0(a) + MAJ(a,b,c); \
    (h)=(g); (g)=(f); (f)=(e); (e)=(d)+_T1; \
    (d)=(c); (c)=(b); (b)=(a); (a)=_T1+_T2; \
} while(0)

/* Constant memory for K (faster access on GPU) */
__constant__ uint32_t d_K[64];

/* ============================================================
 * RESULT STRUCTURE
 * ============================================================ */
#define MAX_RESULTS  500000

typedef struct {
    uint32_t W0;
    uint32_t w1;    /* W[1][0] — the iterated parameter */
    uint32_t DW0;   /* additive difference at W[0] */
    uint32_t Da13;  /* a-register diff at round 13: a'13 - a13 */
    uint32_t DW16;  /* W'[16] - W[16] (from schedule) */
    /* Verification: Da13 + DW16 == 0 */
} Result;

/* ============================================================
 * GPU KERNEL: Birthday search for f17 = Da13 + DW16 = 0
 *
 * Differential convention: De_r = e_r(W') - e_r(W)  (W' minus W)
 *                           DW_k = W'[k] - W[k]
 *
 * Cascade insight (T_CASCADE_17):
 *   Setting DW[r-1] = -De_r_nat makes De_r = 0 deterministically.
 *   After cascade: De3..De16 = 0 for any (W0, w1).
 *   Then f17 = Da13 + DW16 is a pseudo-random 32-bit function of w1.
 *   Expect f17 = 0 with probability 2^{-32}.
 *
 * Key formula (T_DEk_DECOMPOSITION):
 *   De_r = De_r_nat + DW[r-1]   where De_r_nat = what De_r would be
 *   with DW[r-1]=0.  Setting DW[r-1] = -De_r_nat  =>  De_r = 0.
 *
 * Cascade math:
 *   After De3=0, subsequent rounds give:
 *     T1_2 = T1_1 + d1 - d2   (ensures e2_new = e1_new = De_r = 0)
 *   where d1, d2 are d-registers just before round r for W and W'.
 * ============================================================ */
__global__ void birthday_kernel(
        uint32_t W0,
        uint32_t DW0,
        uint64_t start,
        uint64_t n,
        Result*  d_results,
        uint32_t* d_count)
{
    uint64_t idx = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;

    uint32_t w1 = (uint32_t)(start + idx);

    /* --- Message W[0..15] ---
     *   W[0] = W0   (fixed)
     *   W[1] = w1   (iterated)
     *   W[2..15] = 0
     */

    /* --- SHA state for W (1) and W' (2) --- */
    uint32_t a1=IV_A, b1=IV_B, c1=IV_C, d1=IV_D;
    uint32_t e1=IV_E, f1=IV_F, g1=IV_G, h1=IV_H;
    uint32_t a2=IV_A, b2=IV_B, c2=IV_C, d2=IV_D;
    uint32_t e2=IV_E, f2=IV_F, g2=IV_G, h2=IV_H;

    /* ---- Round 1: uses W[0] and W'[0] = W[0]+DW0 ---- */
    SHA_ROUND(a1,b1,c1,d1,e1,f1,g1,h1, W0,       d_K[0]);
    SHA_ROUND(a2,b2,c2,d2,e2,f2,g2,h2, W0+DW0,   d_K[0]);

    /* ---- Round 2: uses W[1]=w1 (same for both, DW[1]=0) ---- */
    SHA_ROUND(a1,b1,c1,d1,e1,f1,g1,h1, w1, d_K[1]);
    SHA_ROUND(a2,b2,c2,d2,e2,f2,g2,h2, w1, d_K[1]);

    /* ---- Rounds 3..16: adaptive cascade ----
     *
     * For round r (r=3..16), W[r-1]=0.
     *
     * Natural differential:
     *   De_r_nat = (d2 + T1_2t) - (d1 + T1_1)
     *   where T1_1  = h1+SIG1(e1)+CH(e1,f1,g1)+K[r-1]  (W[r-1]=0)
     *         T1_2t = h2+SIG1(e2)+CH(e2,f2,g2)+K[r-1]  (W[r-1]=0, tentative)
     *
     * Cascade: set DW[r-1] = -De_r_nat => De_r = 0.
     *
     * Key simplification: after setting DW[r-1]:
     *   T1_2 = T1_2t + DW[r-1] = T1_2t - De_r_nat = T1_1 + d1 - d2
     * So for state updates we don't need T1_2t.
     * We only need T1_2t for r=10 and r=15 to compute Wp[9] and Wp[14]
     * (used in DW16 = W'[16] - W[16]).
     */

    uint32_t Da13 = 0;
    uint32_t Wp9  = 0;   /* W'[9]  set by cascade at round r=10 */
    uint32_t Wp14 = 0;   /* W'[14] set by cascade at round r=15 */

    /* Unroll the loop for GPU efficiency */
    #pragma unroll
    for (int r = 3; r <= 16; r++) {
        /* T1 for W (W[r-1]=0) */
        uint32_t k  = d_K[r-1];
        uint32_t T1_1 = h1 + SIG1(e1) + CH(e1,f1,g1) + k;
        /* T1 for W': T1_2 = T1_1 + d1 - d2  (cascade formula) */
        uint32_t T1_2 = T1_1 + d1 - d2;

        /* Compute Wp[r-1] only when needed for DW16 */
        if (r == 10 || r == 15) {
            uint32_t T1_2t = h2 + SIG1(e2) + CH(e2,f2,g2) + k;
            uint32_t De_nat = (d2 + T1_2t) - (d1 + T1_1);
            uint32_t wp_rm1 = (uint32_t)(-(int32_t)De_nat); /* = -De_nat mod 2^32 */
            if (r == 10) Wp9  = wp_rm1;
            if (r == 15) Wp14 = wp_rm1;
        }

        /* Update W state */
        {
            uint32_t T2 = SIG0(a1) + MAJ(a1,b1,c1);
            h1=g1; g1=f1; f1=e1; e1=d1+T1_1; d1=c1; c1=b1; b1=a1; a1=T1_1+T2;
        }
        /* Update W' state */
        {
            uint32_t T2 = SIG0(a2) + MAJ(a2,b2,c2);
            h2=g2; g2=f2; f2=e2; e2=d2+T1_2; d2=c2; c2=b2; b2=a2; a2=T1_2+T2;
        }

        /* Save Da13 after round 13 */
        if (r == 13) Da13 = a2 - a1;
    }

    /* ---- Compute DW16 ----
     *
     * W[16]  = sig1(W[14]=0)  + W[9]=0  + sig0(W[1]=w1) + W[0]=W0
     *        = sig0(w1) + W0
     * W'[16] = sig1(Wp14)     + Wp9     + sig0(w1)       + (W0+DW0)
     *
     * DW16 = W'[16] - W[16] = sig1(Wp14) + Wp9 + DW0
     */
    uint32_t DW16 = sig1(Wp14) + Wp9 + DW0;

    /* ---- Check f17 = Da13 + DW16 == 0 ---- */
    if (Da13 + DW16 == 0) {
        uint32_t pos = atomicAdd(d_count, 1u);
        if (pos < MAX_RESULTS) {
            d_results[pos].W0   = W0;
            d_results[pos].w1   = w1;
            d_results[pos].DW0  = DW0;
            d_results[pos].Da13 = Da13;
            d_results[pos].DW16 = DW16;
        }
    }
}

/* ============================================================
 * CPU VERIFIER: check one pair using full SHA-256 computation
 * Returns 1 if De3..De17=0, 0 otherwise.
 * ============================================================ */
static int verify_pair_cpu(uint32_t W0, uint32_t w1, uint32_t DW0) {
    /* Build W[0..15] */
    uint32_t W[16]  = {0};
    uint32_t Wp[16] = {0};
    W[0]  = W0;    W[1]  = w1;
    Wp[0] = W0+DW0; Wp[1] = w1;

    /* SHA-256 state */
    uint32_t a1=IV_A,b1=IV_B,c1=IV_C,d1=IV_D,e1=IV_E,f1=IV_F,g1=IV_G,h1=IV_H;
    uint32_t a2=IV_A,b2=IV_B,c2=IV_C,d2=IV_D,e2=IV_E,f2=IV_F,g2=IV_G,h2=IV_H;

    SHA_ROUND(a1,b1,c1,d1,e1,f1,g1,h1, W[0],  HOST_K[0]);
    SHA_ROUND(a2,b2,c2,d2,e2,f2,g2,h2, Wp[0], HOST_K[0]);
    SHA_ROUND(a1,b1,c1,d1,e1,f1,g1,h1, W[1],  HOST_K[1]);
    SHA_ROUND(a2,b2,c2,d2,e2,f2,g2,h2, Wp[1], HOST_K[1]);

    uint32_t saved_a1_13=0, saved_a2_13=0;
    uint32_t Wp9=0, Wp14=0;

    for (int r = 3; r <= 16; r++) {
        uint32_t k   = HOST_K[r-1];
        uint32_t T1_1 = h1 + SIG1(e1) + CH(e1,f1,g1) + k;
        uint32_t T1_2t = h2 + SIG1(e2) + CH(e2,f2,g2) + k;
        uint32_t De_nat = (d2 + T1_2t) - (d1 + T1_1);
        uint32_t wp_rm1 = (uint32_t)(-(int32_t)De_nat);
        Wp[r-1] = wp_rm1;
        if (r==10) Wp9  = wp_rm1;
        if (r==15) Wp14 = wp_rm1;

        uint32_t T1_2 = T1_2t + wp_rm1; /* = T1_1 + d1 - d2 */
        {
            uint32_t T2 = SIG0(a1)+MAJ(a1,b1,c1);
            h1=g1;g1=f1;f1=e1;e1=d1+T1_1;d1=c1;c1=b1;b1=a1;a1=T1_1+T2;
        }
        {
            uint32_t T2 = SIG0(a2)+MAJ(a2,b2,c2);
            h2=g2;g2=f2;f2=e2;e2=d2+T1_2;d2=c2;c2=b2;b2=a2;a2=T1_2+T2;
        }
        /* Check De_r == 0 */
        if (e1 != e2) {
            printf("  [CPU verify] FAIL: De%d = 0x%08x\n", r, e2-e1);
            return 0;
        }
        if (r==13) { saved_a1_13=a1; saved_a2_13=a2; }
    }

    uint32_t Da13 = saved_a2_13 - saved_a1_13;
    uint32_t DW16 = sig1(Wp14) + Wp9 + DW0;
    int ok = (Da13 + DW16 == 0);
    if (!ok)
        printf("  [CPU verify] FAIL: f17 = Da13+DW16 = 0x%08x\n", Da13+DW16);
    return ok;
}

/* ============================================================
 * MAIN
 * ============================================================ */
int main(int argc, char* argv[]) {
    /* --- Parse arguments --- */
    uint32_t W0_start  = 0;            /* first W0 to test */
    uint32_t DW0       = 1;            /* default: ΔW[0] = 1 */
    int      num_W0    = 1;            /* number of W0 values to scan */
    const char* outfile = "pairs.txt"; /* output file */

    if (argc > 1) W0_start = (uint32_t)strtoul(argv[1], NULL, 16);
    if (argc > 2) DW0      = (uint32_t)strtoul(argv[2], NULL, 16);
    if (argc > 3) num_W0   = atoi(argv[3]);
    if (argc > 4) outfile  = argv[4];

    if (W0_start == 0) {
        /* Random W0 by default */
        srand((unsigned)time(NULL));
        W0_start = ((uint32_t)rand() << 16) ^ (uint32_t)rand();
    }

    printf("=== SHA-256 GPU Birthday Search (De3..De17=0) ===\n");
    printf("DW0    = 0x%08x\n", DW0);
    printf("W0     = 0x%08x .. 0x%08x  (%d values)\n",
           W0_start, W0_start + num_W0 - 1, num_W0);
    printf("Output = %s\n\n", outfile);

    /* --- Self-test: verify known pair П-15 --- */
    printf("[Self-test] Verifying known pair П-15: ");
    printf("W0=0xe82222c7, w1=0x516cfb41, DW0=1 ... ");
    fflush(stdout);
    if (verify_pair_cpu(0xe82222c7u, 0x516cfb41u, 1)) {
        printf("OK\n");
    } else {
        printf("FAIL - check implementation!\n");
        return 1;
    }

    /* --- GPU setup --- */
    cudaError_t err;
    err = cudaMemcpyToSymbol(d_K, HOST_K, sizeof(HOST_K));
    if (err != cudaSuccess) {
        fprintf(stderr, "cudaMemcpyToSymbol failed: %s\n", cudaGetErrorString(err));
        return 1;
    }

    Result*  d_results = NULL;
    uint32_t* d_count  = NULL;
    cudaMalloc(&d_results, MAX_RESULTS * sizeof(Result));
    cudaMalloc(&d_count, sizeof(uint32_t));

    Result*  h_results = (Result*)malloc(MAX_RESULTS * sizeof(Result));

    /* Output file */
    FILE* fp = fopen(outfile, "w");
    if (!fp) { perror("fopen"); return 1; }
    fprintf(fp, "# W0,w1,DW0,Da13,DW16\n");

    const int    BLOCK_SIZE  = 256;
    const uint64_t CHUNK     = (1ULL << 26);  /* 64M candidates per launch */
    const uint64_t FULL_SCAN = (1ULL << 32);  /* full W1[0] space */

    uint64_t total_pairs = 0;
    clock_t  t_start_all = clock();

    /* --- Scan num_W0 values of W0 --- */
    for (int wi = 0; wi < num_W0; wi++) {
        uint32_t W0 = W0_start + (uint32_t)wi;
        uint32_t h_count = 0;
        cudaMemset(d_count, 0, sizeof(uint32_t));

        clock_t t0 = clock();
        uint64_t pairs_this_W0 = 0;

        printf("[W0=0x%08x] Scanning 2^32 candidates ...\n", W0);
        fflush(stdout);

        for (uint64_t start = 0; start < FULL_SCAN; start += CHUNK) {
            uint64_t n = CHUNK;
            if (start + n > FULL_SCAN) n = FULL_SCAN - start;

            uint64_t blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
            birthday_kernel<<<(uint32_t)blocks, BLOCK_SIZE>>>(
                    W0, DW0, start, n, d_results, d_count);
            cudaDeviceSynchronize();

            /* Check for new results */
            cudaMemcpy(&h_count, d_count, sizeof(uint32_t), cudaMemcpyDeviceToHost);
            if (h_count > pairs_this_W0) {
                uint32_t new_count = h_count < MAX_RESULTS ? h_count : MAX_RESULTS;
                cudaMemcpy(h_results, d_results,
                           new_count * sizeof(Result), cudaMemcpyDeviceToHost);
                for (uint64_t i = pairs_this_W0; i < new_count; i++) {
                    fprintf(fp, "0x%08x,0x%08x,0x%08x,0x%08x,0x%08x\n",
                            h_results[i].W0, h_results[i].w1,
                            h_results[i].DW0, h_results[i].Da13,
                            h_results[i].DW16);
                    fflush(fp);
                    printf("  FOUND pair #%llu: w1=0x%08x  Da13=0x%08x  DW16=0x%08x\n",
                           (unsigned long long)(total_pairs + i + 1),
                           h_results[i].w1, h_results[i].Da13, h_results[i].DW16);
                }
                pairs_this_W0 = new_count;
            }

            /* Progress every 16 chunks (1G candidates) */
            if (((start / CHUNK) & 15) == 0) {
                double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
                double rate = start > 0 ? start / elapsed / 1e9 : 0;
                printf("  %.1f%%  %.2f G/s  %.1f s elapsed\r",
                       100.0*(double)start/FULL_SCAN, rate, elapsed);
                fflush(stdout);
            }
        }

        double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
        printf("\n  W0=0x%08x done: %llu pairs in %.2f s\n",
               W0, (unsigned long long)pairs_this_W0, elapsed);
        total_pairs += pairs_this_W0;

        /* Reset for next W0 */
        cudaMemset(d_count, 0, sizeof(uint32_t));
    }

    double total_elapsed = (double)(clock()-t_start_all)/CLOCKS_PER_SEC;
    printf("\n=== DONE: %llu total pairs found in %.2f s ===\n",
           (unsigned long long)total_pairs, total_elapsed);
    printf("Results saved to: %s\n", outfile);

    fclose(fp);
    free(h_results);
    cudaFree(d_results);
    cudaFree(d_count);
    return 0;
}
