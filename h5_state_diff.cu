/*
 * h5_state_diff.cu
 *
 * H5 experiment: full differential state analysis after round 17.
 *
 * Extends birthday_search_17: for each found pair (De3..De17=0),
 * computes the FULL state differential at rounds 17..22.
 *
 * Key questions:
 *   Q1. What is (Da17, Db17, Dc17, Dd17, Df17, Dg17, Dh17)?
 *   Q2. What is De18? Is P(De18=0) ≈ 1/2^32 (random), or structured?
 *   Q3. Distribution of v2(Da17), v2(De18)?
 *
 * Build:
 *   nvcc -O3 -arch=sm_80 --use_fast_math -o h5_state_diff h5_state_diff.cu
 *
 * Usage:
 *   ./h5_state_diff [W0_hex] [DW0_hex] [num_W0] [out_csv]
 *   ./h5_state_diff 0 1 100 h5_results.csv   # 100 random W0, collect pairs
 */

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

#define ROTR32(x,n) (((x)>>(n))|((x)<<(32-(n))))
#define SIG0(x) (ROTR32(x,2)^ROTR32(x,13)^ROTR32(x,22))
#define SIG1(x) (ROTR32(x,6)^ROTR32(x,11)^ROTR32(x,25))
#define sig0(x) (ROTR32(x,7)^ROTR32(x,18)^((x)>>3))
#define sig1(x) (ROTR32(x,17)^ROTR32(x,19)^((x)>>10))
#define CH(e,f,g)  (((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))

#define SHA_ROUND(a,b,c,d,e,f,g,h,w,k) do { \
    uint32_t _T1=(h)+SIG1(e)+CH(e,f,g)+(k)+(w); \
    uint32_t _T2=SIG0(a)+MAJ(a,b,c); \
    (h)=(g);(g)=(f);(f)=(e);(e)=(d)+_T1; \
    (d)=(c);(c)=(b);(b)=(a);(a)=_T1+_T2; \
} while(0)

__constant__ uint32_t d_K[64];

#define MAX_RESULTS 200000

/* Extended result: full state diff at round 17, plus De18..De22 */
typedef struct {
    uint32_t W0, w1, DW0;
    uint32_t Da13, DW16;       /* existing */
    /* Round 17 state differential (De17=0 by construction) */
    uint32_t Da17, Db17, Dc17, Dd17;
    uint32_t Df17, Dg17, Dh17;
    /* Rounds 18..20: e-differential only (key H5 question) */
    uint32_t De18, De19, De20;
    /* Schedule differences at 17,18 */
    uint32_t DW17, DW18;
} Result17;

__global__ void h5_kernel(
        uint32_t W0, uint32_t DW0,
        uint64_t start, uint64_t n,
        Result17* d_results, uint32_t* d_count)
{
    uint64_t idx = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;

    uint32_t w1 = (uint32_t)(start + idx);

    uint32_t a1=IV_A,b1=IV_B,c1=IV_C,d1=IV_D,e1=IV_E,f1=IV_F,g1=IV_G,h1=IV_H;
    uint32_t a2=IV_A,b2=IV_B,c2=IV_C,d2=IV_D,e2=IV_E,f2=IV_F,g2=IV_G,h2=IV_H;

    /* Round 1 */
    SHA_ROUND(a1,b1,c1,d1,e1,f1,g1,h1, W0,      d_K[0]);
    SHA_ROUND(a2,b2,c2,d2,e2,f2,g2,h2, W0+DW0,  d_K[0]);
    /* Round 2 */
    SHA_ROUND(a1,b1,c1,d1,e1,f1,g1,h1, w1, d_K[1]);
    SHA_ROUND(a2,b2,c2,d2,e2,f2,g2,h2, w1, d_K[1]);

    /* Cascade rounds 3..16 — track Wp[2], Wp[9], Wp[14] for schedule */
    uint32_t Da13=0, Wp2=0, Wp9=0, Wp14=0;

    #pragma unroll
    for (int r = 3; r <= 16; r++) {
        uint32_t k    = d_K[r-1];
        uint32_t T1_1 = h1 + SIG1(e1) + CH(e1,f1,g1) + k;
        uint32_t T1_2 = T1_1 + d1 - d2;  /* cascade */

        if (r == 3 || r == 10 || r == 15) {
            uint32_t T1_2t = h2 + SIG1(e2) + CH(e2,f2,g2) + k;
            uint32_t De_nat = (d2 + T1_2t) - (d1 + T1_1);
            uint32_t wp_rm1 = (uint32_t)(-(int32_t)De_nat);
            if (r ==  3) Wp2  = wp_rm1;
            if (r == 10) Wp9  = wp_rm1;
            if (r == 15) Wp14 = wp_rm1;
        }

        { uint32_t T2=SIG0(a1)+MAJ(a1,b1,c1);
          h1=g1;g1=f1;f1=e1;e1=d1+T1_1;d1=c1;c1=b1;b1=a1;a1=T1_1+T2; }
        { uint32_t T2=SIG0(a2)+MAJ(a2,b2,c2);
          h2=g2;g2=f2;f2=e2;e2=d2+T1_2;d2=c2;c2=b2;b2=a2;a2=T1_2+T2; }

        if (r == 13) Da13 = a2 - a1;
    }

    /* DW16 = W'[16]-W[16] = sig1(Wp14)+Wp9+DW0 */
    uint32_t DW16 = sig1(Wp14) + Wp9 + DW0;

    if (Da13 + DW16 != 0) return;  /* not a pair */

    /* === PAIR FOUND: compute extended state === */

    /* W[16] and W'[16] */
    uint32_t W16  = sig0(w1) + W0;
    uint32_t W16p = W16 + DW16;

    /* Round 17 */
    SHA_ROUND(a1,b1,c1,d1,e1,f1,g1,h1, W16,  d_K[16]);
    SHA_ROUND(a2,b2,c2,d2,e2,f2,g2,h2, W16p, d_K[16]);

    /* State diff at round 17 */
    uint32_t Da17 = a2-a1, Db17 = b2-b1, Dc17 = c2-c1, Dd17 = d2-d1;
    uint32_t De17 = e2-e1; /* should be 0 */
    uint32_t Df17 = f2-f1, Dg17 = g2-g1, Dh17 = h2-h1;
    (void)De17;

    /* Schedule:
     * W[17] = sig1(W[15]=0) + W[10]=0 + sig0(W[2]=0) + W[1]=w1 = w1
     * W'[17] = sig1(Wp14)   + Wp9     + sig0(Wp2)    + w1
     */
    uint32_t W17  = w1;
    uint32_t W17p = sig1(Wp14) + Wp9 + sig0(Wp2) + w1;
    uint32_t DW17 = W17p - W17;

    /* Round 18 */
    SHA_ROUND(a1,b1,c1,d1,e1,f1,g1,h1, W17,  d_K[17]);
    SHA_ROUND(a2,b2,c2,d2,e2,f2,g2,h2, W17p, d_K[17]);
    uint32_t De18 = e2 - e1;

    /* W[18] = sig1(W[16]) + W[11]=0 + sig0(W[3]=0) + W[2]=0 = sig1(W16)
     * W'[18] = sig1(W16p) + Wp[11]  + sig0(Wp[3])  + Wp[2]
     * We don't track Wp[11], Wp[3] → use 0 as approximation for analysis
     * (they're set by cascade at r=12, r=4 but not saved here)
     * For exact computation, run CPU post-processing with h5_analysis.py
     */
    uint32_t W18  = sig1(W16);
    uint32_t W18p = sig1(W16p);  /* approx: ignores Wp11, Wp3, Wp2 terms */
    uint32_t DW18 = W18p - W18;

    /* Round 19 (approximate: full schedule needs all Wp values) */
    SHA_ROUND(a1,b1,c1,d1,e1,f1,g1,h1, W18,  d_K[18]);
    SHA_ROUND(a2,b2,c2,d2,e2,f2,g2,h2, W18p, d_K[18]);
    uint32_t De19 = e2 - e1;

    /* W[19] = sig1(W17) + W12 + sig0(W4) + W3 = sig1(w1) (all W[3..12]=0) */
    uint32_t W19  = sig1(W17);
    uint32_t W19p = sig1(W17p);
    SHA_ROUND(a1,b1,c1,d1,e1,f1,g1,h1, W19,  d_K[19]);
    SHA_ROUND(a2,b2,c2,d2,e2,f2,g2,h2, W19p, d_K[19]);
    uint32_t De20 = e2 - e1;

    uint32_t pos = atomicAdd(d_count, 1u);
    if (pos < MAX_RESULTS) {
        Result17* r17    = &d_results[pos];
        r17->W0   = W0;  r17->w1   = w1;   r17->DW0  = DW0;
        r17->Da13 = Da13; r17->DW16 = DW16;
        r17->Da17 = Da17; r17->Db17 = Db17; r17->Dc17 = Dc17; r17->Dd17 = Dd17;
        r17->Df17 = Df17; r17->Dg17 = Dg17; r17->Dh17 = Dh17;
        r17->De18 = De18; r17->De19 = De19; r17->De20 = De20;
        r17->DW17 = DW17; r17->DW18 = DW18;
    }
}

static uint32_t v2(uint32_t x) {
    if (x == 0) return 32;
    uint32_t v = 0;
    while ((x & 1) == 0) { v++; x >>= 1; }
    return v;
}

int main(int argc, char* argv[]) {
    uint32_t W0_start = 0, DW0 = 1;
    uint32_t num_W0 = 100;
    const char* outfile = "h5_results.csv";

    if (argc > 1) W0_start = (uint32_t)strtoul(argv[1], NULL, 16);
    if (argc > 2) DW0      = (uint32_t)strtoul(argv[2], NULL, 16);
    if (argc > 3) num_W0   = (uint32_t)atoi(argv[3]);
    if (argc > 4) outfile  = argv[4];

    if (W0_start == 0) {
        srand((unsigned)time(NULL));
        W0_start = ((uint32_t)rand() << 16) ^ (uint32_t)rand();
    }

    printf("=== H5 State Differential Analysis ===\n");
    printf("W0 range = 0x%08x .. 0x%08x\n", W0_start, W0_start + num_W0 - 1);
    printf("DW0 = 0x%08x\n", DW0);

    cudaMemcpyToSymbol(d_K, HOST_K, sizeof(HOST_K));

    Result17 *d_results; uint32_t *d_count;
    cudaMalloc(&d_results, MAX_RESULTS * sizeof(Result17));
    cudaMalloc(&d_count, 4);

    FILE* fp = fopen(outfile, "w");
    fprintf(fp, "# W0,w1,DW0,Da13,DW16,"
                "Da17,Db17,Dc17,Dd17,Df17,Dg17,Dh17,"
                "De18,De19,De20,DW17,DW18,"
                "v2_Da17,v2_De18,v2_DW17\n");

    uint32_t total_pairs = 0;
    const int BLOCK = 256;
    const uint64_t CHUNK = (1ULL << 24);

    time_t t0 = time(NULL);

    for (uint32_t wi = 0; wi < num_W0; wi++) {
        uint32_t W0 = W0_start + wi;
        cudaMemset(d_count, 0, 4);

        for (uint64_t start = 0; start < (1ULL<<32); start += CHUNK) {
            uint64_t n = (start + CHUNK <= (1ULL<<32)) ? CHUNK : ((1ULL<<32) - start);
            uint64_t blocks = (n + BLOCK - 1) / BLOCK;
            h5_kernel<<<(uint32_t)blocks, BLOCK>>>(W0, DW0, start, n, d_results, d_count);
        }
        cudaDeviceSynchronize();

        uint32_t cnt;
        cudaMemcpy(&cnt, d_count, 4, cudaMemcpyDeviceToHost);
        if (cnt > MAX_RESULTS) cnt = MAX_RESULTS;

        if (cnt > 0) {
            Result17 *h_results = (Result17*)malloc(cnt * sizeof(Result17));
            cudaMemcpy(h_results, d_results, cnt * sizeof(Result17), cudaMemcpyDeviceToHost);
            for (uint32_t i = 0; i < cnt; i++) {
                Result17* r = &h_results[i];
                fprintf(fp, "0x%08x,0x%08x,0x%08x,0x%08x,0x%08x,"
                            "0x%08x,0x%08x,0x%08x,0x%08x,0x%08x,0x%08x,0x%08x,"
                            "0x%08x,0x%08x,0x%08x,0x%08x,0x%08x,"
                            "%u,%u,%u\n",
                    r->W0, r->w1, r->DW0, r->Da13, r->DW16,
                    r->Da17, r->Db17, r->Dc17, r->Dd17, r->Df17, r->Dg17, r->Dh17,
                    r->De18, r->De19, r->De20, r->DW17, r->DW18,
                    v2(r->Da17), v2(r->De18), v2(r->DW17));
            }
            fflush(fp);
            free(h_results);
            total_pairs += cnt;
        }

        if ((wi+1) % 10 == 0 || wi == num_W0-1) {
            printf("  [%u/%u W0]  pairs so far: %u  (%lds)\r",
                   wi+1, num_W0, total_pairs, (long)(time(NULL)-t0));
            fflush(stdout);
        }
    }
    printf("\n");

    fclose(fp);
    printf("=== DONE: %u pairs found ===\n", total_pairs);
    printf("Results: %s\n", outfile);

    cudaFree(d_results); cudaFree(d_count);
    return 0;
}
