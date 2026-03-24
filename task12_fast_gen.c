/*
 * ЗАДАНИЕ 12: Fast amortized pair generator
 * Uses amortized approach: sweep W0 once, test against N backgrounds.
 * Outputs Wn/Wf pairs in hex format for Python analysis.
 *
 * Compile: gcc -O3 -march=native -o task12_fast_gen task12_fast_gen.c -lpthread
 * Usage:   ./task12_fast_gen [num_sets] [num_threads] [max_seconds]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#define MASK 0xFFFFFFFFU

static const uint32_t K[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
};
static const uint32_t H0[8] = {
    0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
    0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19
};

static inline uint32_t rotr(uint32_t x, int n) { return (x>>n)|(x<<(32-n)); }
static inline uint32_t sig0(uint32_t x) { return rotr(x,7)^rotr(x,18)^(x>>3); }
static inline uint32_t sig1(uint32_t x) { return rotr(x,17)^rotr(x,19)^(x>>10); }
static inline uint32_t Sig0(uint32_t x) { return rotr(x,2)^rotr(x,13)^rotr(x,22); }
static inline uint32_t Sig1(uint32_t x) { return rotr(x,6)^rotr(x,11)^rotr(x,25); }
static inline uint32_t Ch(uint32_t e,uint32_t f,uint32_t g) { return (e&f)^((~e)&g); }
static inline uint32_t Maj(uint32_t a,uint32_t b,uint32_t c) { return (a&b)^(a&c)^(b&c); }

/* Wang chain: compute δe[17]. Fill Wf. */
static uint32_t wang_de17(const uint32_t Wn[16], uint32_t Wf[16]) {
    memcpy(Wf, Wn, 64);
    Wf[0] = Wn[0] ^ 0x8000;
    uint32_t sn[8], sf[8];
    memcpy(sn, H0, 32); memcpy(sf, H0, 32);

    uint32_t T1n = sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[0]+Wn[0];
    uint32_t T2n = Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
    uint32_t T1f = sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[0]+Wf[0];
    uint32_t T2f = Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);
    uint32_t nn[8]={T1n+T2n,sn[0],sn[1],sn[2],sn[3]+T1n,sn[4],sn[5],sn[6]};
    uint32_t nf[8]={T1f+T2f,sf[0],sf[1],sf[2],sf[3]+T1f,sf[4],sf[5],sf[6]};
    memcpy(sn,nn,32); memcpy(sf,nf,32);

    for(int r=1;r<16;r++) {
        uint32_t dd=sf[3]-sn[3], dh=sf[7]-sn[7];
        uint32_t dS=Sig1(sf[4])-Sig1(sn[4]);
        uint32_t dC=Ch(sf[4],sf[5],sf[6])-Ch(sn[4],sn[5],sn[6]);
        Wf[r] = Wn[r] + (uint32_t)(-(int32_t)(dd+dh+dS+dC));
        T1n = sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[r]+Wn[r];
        T2n = Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
        T1f = sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[r]+Wf[r];
        T2f = Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);
        uint32_t tn[8]={T1n+T2n,sn[0],sn[1],sn[2],sn[3]+T1n,sn[4],sn[5],sn[6]};
        uint32_t tf[8]={T1f+T2f,sf[0],sf[1],sf[2],sf[3]+T1f,sf[4],sf[5],sf[6]};
        memcpy(sn,tn,32); memcpy(sf,tf,32);
    }

    uint32_t Wn16 = sig1(Wn[14])+Wn[9]+sig0(Wn[1])+Wn[0];
    uint32_t Wf16 = sig1(Wf[14])+Wf[9]+sig0(Wf[1])+Wf[0];
    T1n = sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[16]+Wn16;
    T1f = sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[16]+Wf16;
    return (sf[3]+T1f)-(sn[3]+T1n);
}

/* Thread state */
typedef struct {
    int thread_id, num_threads;
    int num_sets;
    uint32_t (*backgrounds)[16]; /* backgrounds[i][0..15] */
    int *found;                  /* found[i] = 1 if hit */
    uint32_t (*result_Wn)[16];
    uint32_t (*result_Wf)[16];
    volatile int *total_found;
    volatile int *stop_flag;
    pthread_mutex_t *mutex;
    uint64_t evals;
    int time_limit;
} targ_t;

static void *worker(void *arg) {
    targ_t *ta = (targ_t *)arg;
    ta->evals = 0;
    uint32_t Wn[16], Wf[16];
    time_t t0 = time(NULL);

    uint64_t start = (uint64_t)ta->thread_id * (0x100000000ULL / ta->num_threads);
    uint64_t end = (ta->thread_id == ta->num_threads - 1) ?
                   0x100000000ULL :
                   (uint64_t)(ta->thread_id + 1) * (0x100000000ULL / ta->num_threads);

    for (uint64_t w0 = start; w0 < end; w0++) {
        if (*ta->stop_flag || *ta->total_found >= ta->num_sets) break;
        if (ta->time_limit > 0 && (w0 & 0xFFFFF) == 0) {
            if (time(NULL) - t0 >= ta->time_limit) break;
        }

        for (int s = 0; s < ta->num_sets; s++) {
            if (ta->found[s]) continue;

            for (int i = 1; i < 16; i++) Wn[i] = ta->backgrounds[s][i];
            Wn[0] = (uint32_t)w0;

            uint32_t de17 = wang_de17(Wn, Wf);
            ta->evals++;

            if (de17 == 0) {
                pthread_mutex_lock(ta->mutex);
                if (!ta->found[s]) {
                    ta->found[s] = 1;
                    memcpy(ta->result_Wn[s], Wn, 64);
                    memcpy(ta->result_Wf[s], Wf, 64);
                    (*ta->total_found)++;
                    fprintf(stderr, "\r  Found %d/%d pairs (%.1fM evals)   ",
                            *ta->total_found, ta->num_sets, ta->evals/1e6);
                }
                pthread_mutex_unlock(ta->mutex);
            }
        }
    }
    return NULL;
}

int main(int argc, char **argv) {
    int num_sets = 200;
    int num_threads = 4;
    int time_limit = 600; /* seconds */

    if (argc > 1) num_sets = atoi(argv[1]);
    if (argc > 2) num_threads = atoi(argv[2]);
    if (argc > 3) time_limit = atoi(argv[3]);

    fprintf(stderr, "Task12 fast generator: %d sets, %d threads, %ds limit\n",
            num_sets, num_threads, time_limit);

    /* Allocate */
    uint32_t (*backgrounds)[16] = malloc(num_sets * 16 * sizeof(uint32_t));
    int *found = calloc(num_sets, sizeof(int));
    uint32_t (*result_Wn)[16] = malloc(num_sets * 16 * sizeof(uint32_t));
    uint32_t (*result_Wf)[16] = malloc(num_sets * 16 * sizeof(uint32_t));
    volatile int total_found = 0;
    volatile int stop_flag = 0;
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

    /* Generate random backgrounds */
    FILE *rng = fopen("/dev/urandom", "rb");
    for (int s = 0; s < num_sets; s++) {
        fread(backgrounds[s], 4, 16, rng);
    }
    fclose(rng);

    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);

    pthread_t threads[64];
    targ_t targs[64];
    for (int t = 0; t < num_threads; t++) {
        targs[t].thread_id = t;
        targs[t].num_threads = num_threads;
        targs[t].num_sets = num_sets;
        targs[t].backgrounds = backgrounds;
        targs[t].found = found;
        targs[t].result_Wn = result_Wn;
        targs[t].result_Wf = result_Wf;
        targs[t].total_found = &total_found;
        targs[t].stop_flag = &stop_flag;
        targs[t].mutex = &mutex;
        targs[t].evals = 0;
        targs[t].time_limit = time_limit;
        pthread_create(&threads[t], NULL, worker, &targs[t]);
    }

    for (int t = 0; t < num_threads; t++)
        pthread_join(threads[t], NULL);

    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double elapsed = (ts1.tv_sec-ts0.tv_sec)+(ts1.tv_nsec-ts0.tv_nsec)*1e-9;

    uint64_t total_evals = 0;
    for (int t = 0; t < num_threads; t++) total_evals += targs[t].evals;

    fprintf(stderr, "\nDone: %d pairs in %.1fs (%.2e evals, %.1fM/s)\n",
            (int)total_found, elapsed, (double)total_evals,
            total_evals/elapsed/1e6);

    /* Output pairs */
    for (int s = 0; s < num_sets; s++) {
        if (!found[s]) continue;
        printf("PAIR\n");
        printf("Wn =");
        for (int w = 0; w < 16; w++) printf(" %08x", result_Wn[s][w]);
        printf("\n");
        printf("Wf =");
        for (int w = 0; w < 16; w++) printf(" %08x", result_Wf[s][w]);
        printf("\n");
    }

    free(backgrounds); free(found);
    free(result_Wn); free(result_Wf);
    return 0;
}
