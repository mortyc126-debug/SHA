/*
 * ЗАДАНИЕ 12: Probabilistic pair generator
 * Strategy: N backgrounds × M random W0 samples each.
 * With N*M total trials and P(hit)=1/2^32 per trial,
 * expected hits = N*M / 2^32.
 *
 * For 10000 hits: need N*M ≈ 10000 * 2^32 ≈ 4.3e13.
 * At 20M/s: ~2.15e6 seconds. Too slow.
 *
 * BETTER: use large N (many backgrounds), moderate M.
 * Each background acts as independent trial.
 * With N=4e9 backgrounds and M=1 sample each: 4e9 trials, ~1 hit.
 * With N=1e6 and M=4300: N*M=4.3e9, ~1 hit.
 *
 * BEST: for each background, compute δe[17] for ONE W0 value.
 * P(hit)=2^{-32}. Need 2^{32}*target backgrounds.
 * At 20M/s: 2^{32}*10000/20M = forever.
 *
 * ACTUAL BEST: Amortize by keeping MANY backgrounds active.
 * Use hash table: for each random (bg_id, W0), compute δe[17].
 * Store in hash map: δe[17] → (bg_id, W0).
 * Look for δe[17] = 0.
 *
 * This is effectively just random sampling. Speed = 20M evals/s.
 * Expected time per hit: 2^32 / 20M ≈ 215 seconds.
 * For 50 hits: ~10750s ≈ 3 hours.
 *
 * OPTIMIZATION: use SIMD-like approach with 4 independent streams.
 * Each thread: random bg, random W0, compute wang_de17, check.
 * No amortization overhead. Pure random sampling.
 *
 * gcc -O3 -march=native -o task12_sample_gen task12_sample_gen.c -lpthread
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
static inline uint32_t Sig0(uint32_t x) { return rotr(x,2)^rotr(x,13)^rotr(x,22); }
static inline uint32_t Sig1(uint32_t x) { return rotr(x,6)^rotr(x,11)^rotr(x,25); }
static inline uint32_t sig0(uint32_t x) { return rotr(x,7)^rotr(x,18)^(x>>3); }
static inline uint32_t sig1(uint32_t x) { return rotr(x,17)^rotr(x,19)^(x>>10); }
static inline uint32_t Ch(uint32_t e,uint32_t f,uint32_t g) { return (e&f)^((~e)&g); }
static inline uint32_t Maj(uint32_t a,uint32_t b,uint32_t c) { return (a&b)^(a&c)^(b&c); }

/* Fast PRNG: xoshiro128** */
typedef struct { uint32_t s[4]; } rng_t;
static inline uint32_t rotl(uint32_t x, int k) { return (x<<k)|(x>>(32-k)); }
static inline uint32_t rng_next(rng_t *r) {
    uint32_t res = rotl(r->s[1]*5,7)*9;
    uint32_t t = r->s[1]<<9;
    r->s[2]^=r->s[0]; r->s[3]^=r->s[1]; r->s[1]^=r->s[2]; r->s[0]^=r->s[3];
    r->s[2]^=t; r->s[3]=rotl(r->s[3],11);
    return res;
}
static void rng_seed(rng_t *r, uint64_t s) {
    r->s[0]=(uint32_t)s; r->s[1]=(uint32_t)(s>>32);
    r->s[2]=r->s[0]^0xDEADBEEF; r->s[3]=r->s[1]^0xCAFEBABE;
    for(int i=0;i<20;i++) rng_next(r);
}

static uint32_t wang_de17(const uint32_t Wn[16], uint32_t Wf[16]) {
    for(int i=0;i<16;i++) Wf[i]=Wn[i];
    Wf[0] = Wn[0] ^ 0x8000;
    uint32_t sn[8], sf[8];
    for(int i=0;i<8;i++){sn[i]=H0[i];sf[i]=H0[i];}

    uint32_t T1n=sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[0]+Wn[0];
    uint32_t T2n=Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
    uint32_t T1f=sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[0]+Wf[0];
    uint32_t T2f=Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);
    {uint32_t t[8]={T1n+T2n,sn[0],sn[1],sn[2],sn[3]+T1n,sn[4],sn[5],sn[6]};memcpy(sn,t,32);}
    {uint32_t t[8]={T1f+T2f,sf[0],sf[1],sf[2],sf[3]+T1f,sf[4],sf[5],sf[6]};memcpy(sf,t,32);}

    for(int r=1;r<16;r++){
        uint32_t dd=sf[3]-sn[3],dh=sf[7]-sn[7];
        uint32_t dS=Sig1(sf[4])-Sig1(sn[4]),dC=Ch(sf[4],sf[5],sf[6])-Ch(sn[4],sn[5],sn[6]);
        Wf[r]=Wn[r]+(uint32_t)(-(int32_t)(dd+dh+dS+dC));
        T1n=sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[r]+Wn[r];
        T2n=Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
        T1f=sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[r]+Wf[r];
        T2f=Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);
        {uint32_t t[8]={T1n+T2n,sn[0],sn[1],sn[2],sn[3]+T1n,sn[4],sn[5],sn[6]};memcpy(sn,t,32);}
        {uint32_t t[8]={T1f+T2f,sf[0],sf[1],sf[2],sf[3]+T1f,sf[4],sf[5],sf[6]};memcpy(sf,t,32);}
    }
    uint32_t Wn16=sig1(Wn[14])+Wn[9]+sig0(Wn[1])+Wn[0];
    uint32_t Wf16=sig1(Wf[14])+Wf[9]+sig0(Wf[1])+Wf[0];
    T1n=sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[16]+Wn16;
    T1f=sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[16]+Wf16;
    return (sf[3]+T1f)-(sn[3]+T1n);
}

typedef struct {
    int thread_id;
    uint64_t seed;
    int time_limit;
    /* Results */
    uint32_t (*out_Wn)[16];
    uint32_t (*out_Wf)[16];
    int found;
    int max_results;
    uint64_t evals;
    volatile int *global_found;
    volatile int *target;
    pthread_mutex_t *out_mutex;
} targ_t;

static void *worker(void *arg) {
    targ_t *ta = (targ_t *)arg;
    rng_t rng;
    rng_seed(&rng, ta->seed);
    ta->found = 0;
    ta->evals = 0;
    time_t t0 = time(NULL);

    uint32_t Wn[16], Wf[16];

    while(*ta->global_found < *ta->target) {
        if((ta->evals & 0xFFFFF) == 0 && time(NULL) - t0 >= ta->time_limit) break;

        /* Generate random background + random W0 */
        for(int i=0;i<16;i++) Wn[i] = rng_next(&rng);

        uint32_t de17 = wang_de17(Wn, Wf);
        ta->evals++;

        if(de17 == 0) {
            pthread_mutex_lock(ta->out_mutex);
            if(ta->found < ta->max_results) {
                memcpy(ta->out_Wn[ta->found], Wn, 64);
                memcpy(ta->out_Wf[ta->found], Wf, 64);
                ta->found++;
                (*ta->global_found)++;
                /* Output immediately with flush */
                printf("PAIR\n");
                printf("Wn =");
                for(int w=0;w<16;w++) printf(" %08x", Wn[w]);
                printf("\n");
                printf("Wf =");
                for(int w=0;w<16;w++) printf(" %08x", Wf[w]);
                printf("\n");
                fflush(stdout);
                fprintf(stderr, "\r  Found %d pairs (%luM evals, %.0fs)   ",
                        *ta->global_found, ta->evals/1000000,
                        difftime(time(NULL), t0));
            }
            pthread_mutex_unlock(ta->out_mutex);
        }
    }
    return NULL;
}

int main(int argc, char **argv) {
    int target = 10000;
    int num_threads = 4;
    int time_limit = 600;
    if(argc>1) target = atoi(argv[1]);
    if(argc>2) num_threads = atoi(argv[2]);
    if(argc>3) time_limit = atoi(argv[3]);

    fprintf(stderr, "Task12 random sampler: target=%d, threads=%d, limit=%ds\n",
            target, num_threads, time_limit);
    fprintf(stderr, "Expected rate: ~1 pair per 2^32/rate seconds\n");

    volatile int global_found = 0;
    volatile int vol_target = target;
    pthread_mutex_t out_mutex = PTHREAD_MUTEX_INITIALIZER;

    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);

    int max_per_thread = target/num_threads + 100;
    pthread_t threads[64];
    targ_t targs[64];

    for(int t=0; t<num_threads; t++) {
        targs[t].thread_id = t;
        targs[t].seed = time(NULL) ^ ((uint64_t)(t+1) << 32) ^ 0xABCDEF0123456789ULL;
        targs[t].time_limit = time_limit;
        targs[t].out_Wn = malloc(max_per_thread * 16 * sizeof(uint32_t));
        targs[t].out_Wf = malloc(max_per_thread * 16 * sizeof(uint32_t));
        targs[t].found = 0;
        targs[t].max_results = max_per_thread;
        targs[t].evals = 0;
        targs[t].global_found = &global_found;
        targs[t].target = &vol_target;
        targs[t].out_mutex = &out_mutex;
        pthread_create(&threads[t], NULL, worker, &targs[t]);
    }

    for(int t=0; t<num_threads; t++)
        pthread_join(threads[t], NULL);

    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double elapsed = (ts1.tv_sec-ts0.tv_sec)+(ts1.tv_nsec-ts0.tv_nsec)*1e-9;

    uint64_t total_evals = 0;
    int total_pairs = 0;
    for(int t=0; t<num_threads; t++) {
        total_evals += targs[t].evals;
        for(int i=0; i<targs[t].found; i++) {
            printf("PAIR\n");
            printf("Wn =");
            for(int w=0;w<16;w++) printf(" %08x", targs[t].out_Wn[i][w]);
            printf("\n");
            printf("Wf =");
            for(int w=0;w<16;w++) printf(" %08x", targs[t].out_Wf[i][w]);
            printf("\n");
            total_pairs++;
        }
        free(targs[t].out_Wn);
        free(targs[t].out_Wf);
    }

    fprintf(stderr, "\nDone: %d pairs in %.1fs (%.2e evals, %.1fM/s, %.0f evals/pair)\n",
            total_pairs, elapsed, (double)total_evals,
            total_evals/elapsed/1e6,
            total_pairs > 0 ? (double)total_evals/total_pairs : 0);
    return 0;
}
