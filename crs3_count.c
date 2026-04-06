/*
 * CRS-3: КОЛИЧЕСТВО vs ДЛИНА циклов.
 *
 * g32: {0,1}^32 → {0,1}^32. Input space = 2^32.
 * Random function 2^32→2^32:
 *   Expected #cycles ≈ 0.5 * ln(2^32) ≈ 11
 *   Expected avg λ ≈ sqrt(π * 2^32 / 8) ≈ 23K
 *   Expected max λ ≈ 2^16 = 65K
 *
 * SHA-256 g32: 3 cycles. Is 3 << 11 significant?
 *
 * Test: many random functions 2^N → 2^N for small N,
 * count cycles, compare with SHA-256-like construction.
 *
 * Compile: gcc -O3 -march=native -o crs3_count crs3_count.c
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))
#define CH(e,f,g) (((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))
#define SIG0(x) (ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define SIG1(x) (ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define sig0(x) (ROTR(x,7)^ROTR(x,18)^((x)>>3))
#define sig1(x) (ROTR(x,17)^ROTR(x,19)^((x)>>10))

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

/*
 * SMALL N exact cycle counting.
 * For N-bit space (2^N elements), build full function table,
 * follow every element to its cycle, count distinct cycles.
 */

/* SHA-256-like mini hash for N-bit */
static inline uint32_t mini_sha(uint32_t x, int N) {
    uint32_t mask = (1u << N) - 1;
    uint32_t W[64]; memset(W,0,sizeof(W)); W[0] = x & mask;
    for(int i=16;i<64;i++) W[i]=(sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16])&mask;
    uint32_t a=IV[0]&mask,b=IV[1]&mask,c=IV[2]&mask,d=IV[3]&mask;
    uint32_t e=IV[4]&mask,f=IV[5]&mask,g=IV[6]&mask,h=IV[7]&mask;
    for(int r=0;r<64;r++){
        uint32_t t1=(h+SIG1(e)+CH(e,f,g)+K[r]+W[r])&mask;
        uint32_t t2=(SIG0(a)+MAJ(a,b,c))&mask;
        h=g;g=f;f=e;e=(d+t1)&mask;d=c;c=b;b=a;a=(t1+t2)&mask;
    }
    return (a+(IV[0]&mask))&mask;
}

/* Random function for N-bit (precomputed table) */
static uint32_t *rand_table = NULL;

static inline uint32_t rand_func(uint32_t x, int N) {
    return rand_table[x];
}

/* Count cycles by visiting every node */
int count_all_cycles(uint32_t (*f)(uint32_t, int), int N) {
    uint32_t size = 1u << N;
    uint8_t *visited = calloc(size, 1);  /* 0=unvisited, 1=in-progress, 2=done */
    int n_cycles = 0;

    for (uint32_t start = 0; start < size; start++) {
        if (visited[start] == 2) continue;

        /* Follow chain from start */
        uint32_t x = start;
        while (visited[x] == 0) {
            visited[x] = 1; /* mark in-progress */
            x = f(x, N);
        }

        if (visited[x] == 1) {
            /* Found a new cycle! x is on the cycle. */
            n_cycles++;
            /* Mark entire cycle as done */
            uint32_t p = x;
            do {
                visited[p] = 2;
                p = f(p, N);
            } while (p != x);
        }

        /* Mark tail (non-cycle) as done */
        x = start;
        while (visited[x] == 1) {
            visited[x] = 2;
            x = f(x, N);
        }
    }

    free(visited);
    return n_cycles;
}

int main() {
    printf("CRS-3: CYCLE COUNT — SHA-256 vs RANDOM\n");
    printf("=========================================\n\n");

    /* Small N where we can enumerate ALL 2^N elements */
    printf("EXACT cycle count (full enumeration):\n");
    printf("%-4s | %-8s %-8s | %-20s\n", "N", "SHA-mini", "Random", "Expected ≈ 0.5*ln(2^N)");
    printf("-----+-------------------+---------------------\n");

    for (int N = 8; N <= 24; N += 2) {
        uint32_t size = 1u << N;

        /* SHA-mini */
        clock_t t0 = clock();
        int sha_cycles = count_all_cycles(mini_sha, N);
        double sha_time = (double)(clock()-t0)/CLOCKS_PER_SEC;

        /* Random function */
        rand_table = malloc(size * sizeof(uint32_t));
        srand(N * 12345 + 67);
        for (uint32_t i = 0; i < size; i++)
            rand_table[i] = ((uint32_t)rand() ^ ((uint32_t)rand()<<16)) & (size-1);

        t0 = clock();
        int rand_cycles = count_all_cycles(rand_func, N);
        double rand_time = (double)(clock()-t0)/CLOCKS_PER_SEC;
        free(rand_table); rand_table = NULL;

        double expected = 0.5 * N * log(2);

        printf("%-4d | %-8d %-8d | expected ≈ %.1f",
               N, sha_cycles, rand_cycles, expected);

        double ratio = (double)sha_cycles / rand_cycles;
        if (ratio < 0.5)
            printf("  ★ SHA %.1fx FEWER", 1.0/ratio);
        else if (ratio > 2.0)
            printf("  SHA %.1fx MORE", ratio);
        printf("  (%.1fs+%.1fs)\n", sha_time, rand_time);

        if (N >= 22 && sha_time > 30) break; /* too slow */
    }

    /* Multiple random functions for statistics at N=16 */
    printf("\nStatistics at N=16 (100 random functions):\n");
    {
        int N = 16;
        uint32_t size = 1u << N;

        int sha_c = count_all_cycles(mini_sha, N);

        int rand_counts[100];
        double sum = 0, sum2 = 0;
        rand_table = malloc(size * sizeof(uint32_t));
        for (int trial = 0; trial < 100; trial++) {
            srand(trial * 9999 + 42);
            for (uint32_t i = 0; i < size; i++)
                rand_table[i] = ((uint32_t)rand() ^ ((uint32_t)rand()<<16)) & (size-1);
            rand_counts[trial] = count_all_cycles(rand_func, N);
            sum += rand_counts[trial];
            sum2 += rand_counts[trial] * rand_counts[trial];
        }
        free(rand_table); rand_table = NULL;

        double mean = sum / 100;
        double std = sqrt(sum2/100 - mean*mean);

        /* Sort for percentiles */
        for (int i=0;i<99;i++) for(int j=i+1;j<100;j++)
            if(rand_counts[i]>rand_counts[j]){int t=rand_counts[i];rand_counts[i]=rand_counts[j];rand_counts[j]=t;}

        int count_below_sha = 0;
        for (int i=0;i<100;i++) if(rand_counts[i] <= sha_c) count_below_sha++;

        printf("  SHA-mini N=16: %d cycles\n", sha_c);
        printf("  Random N=16:   mean=%.1f, std=%.1f, min=%d, max=%d\n",
               mean, std, rand_counts[0], rand_counts[99]);
        printf("  P(random ≤ SHA): %d/100 = %.0f%%\n", count_below_sha, count_below_sha*1.0);
        printf("  Z-score: %.2f\n", (sha_c - mean) / std);

        if (count_below_sha <= 5)
            printf("  ★★★ STATISTICALLY SIGNIFICANT: SHA has anomalously FEW cycles!\n");
        else if (count_below_sha <= 15)
            printf("  ★★ Unusual but not extreme.\n");
        else
            printf("  Within normal range.\n");
    }

    return 0;
}
