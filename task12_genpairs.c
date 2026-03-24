/*
 * ЗАДАНИЕ 12: Mass δe[17]=0 pair generator
 * Outputs Wn/Wf pairs to stdout for Python analysis.
 *
 * Compile: gcc -O3 -march=native -o task12_genpairs task12_genpairs.c -lpthread
 * Usage:   ./task12_genpairs [num_sets] [num_threads]
 *          Default: 5000 sets, 4 threads
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

/* PRNG: xoshiro128** */
typedef struct { uint32_t s[4]; } rng_t;
static inline uint32_t rotl(uint32_t x, int k) { return (x<<k)|(x>>(32-k)); }
static uint32_t rng_next(rng_t *r) {
    uint32_t res = rotl(r->s[1]*5,7)*9;
    uint32_t t = r->s[1]<<9;
    r->s[2] ^= r->s[0]; r->s[3] ^= r->s[1];
    r->s[1] ^= r->s[2]; r->s[0] ^= r->s[3];
    r->s[2] ^= t; r->s[3] = rotl(r->s[3],11);
    return res;
}
static void rng_seed(rng_t *r, uint64_t s) {
    r->s[0]=(uint32_t)s; r->s[1]=(uint32_t)(s>>32);
    r->s[2]=r->s[0]^0x12345678; r->s[3]=r->s[1]^0x9abcdef0;
    for(int i=0;i<20;i++) rng_next(r);
}

/* Wang chain: returns δe[17]. Also fills Wf[0..15]. */
static uint32_t wang_de17(const uint32_t Wn[16], uint32_t Wf[16]) {
    memcpy(Wf, Wn, 64);
    Wf[0] = Wn[0] ^ 0x8000;
    uint32_t sn[8], sf[8];
    memcpy(sn, H0, 32); memcpy(sf, H0, 32);

    /* Round 0 */
    uint32_t T1n = sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[0]+Wn[0];
    uint32_t T2n = Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
    uint32_t T1f = sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[0]+Wf[0];
    uint32_t T2f = Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);
    uint32_t nn[8]={T1n+T2n,sn[0],sn[1],sn[2],sn[3]+T1n,sn[4],sn[5],sn[6]};
    uint32_t nf[8]={T1f+T2f,sf[0],sf[1],sf[2],sf[3]+T1f,sf[4],sf[5],sf[6]};
    memcpy(sn,nn,32); memcpy(sf,nf,32);

    /* Rounds 1..15 */
    for(int r=1;r<16;r++) {
        uint32_t dd=sf[3]-sn[3], dh=sf[7]-sn[7];
        uint32_t dS=Sig1(sf[4])-Sig1(sn[4]);
        uint32_t dC=Ch(sf[4],sf[5],sf[6])-Ch(sn[4],sn[5],sn[6]);
        uint32_t dWr = -(dd+dh+dS+dC);
        Wf[r] = Wn[r]+dWr;
        T1n = sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[r]+Wn[r];
        T2n = Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
        T1f = sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[r]+Wf[r];
        T2f = Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);
        uint32_t tn[8]={T1n+T2n,sn[0],sn[1],sn[2],sn[3]+T1n,sn[4],sn[5],sn[6]};
        uint32_t tf[8]={T1f+T2f,sf[0],sf[1],sf[2],sf[3]+T1f,sf[4],sf[5],sf[6]};
        memcpy(sn,tn,32); memcpy(sf,tf,32);
    }

    /* Round 16: compute δe[17] */
    uint32_t Wn16 = sig1(Wn[14])+Wn[9]+sig0(Wn[1])+Wn[0];
    uint32_t Wf16 = sig1(Wf[14])+Wf[9]+sig0(Wf[1])+Wf[0];
    T1n = sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[16]+Wn16;
    T1f = sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[16]+Wf16;
    return (sf[3]+T1f)-(sn[3]+T1n);
}

/* Thread data */
typedef struct {
    int thread_id;
    int num_threads;
    int num_sets;
    uint64_t seed;
    /* Results buffer */
    uint32_t (*results)[32]; /* pairs of (Wn[16], Wf[16]) */
    int result_count;
    int max_results;
    uint64_t total_evals;
} thread_data_t;

static void *search_thread(void *arg) {
    thread_data_t *td = (thread_data_t *)arg;
    rng_t rng;
    rng_seed(&rng, td->seed);
    td->result_count = 0;
    td->total_evals = 0;

    for(int s = td->thread_id; s < td->num_sets; s += td->num_threads) {
        uint32_t Wn[16], Wf[16];
        for(int i=1;i<16;i++) Wn[i] = rng_next(&rng);

        /* Sweep W[0] through full 2^32 range */
        for(uint64_t w0 = 0; w0 < (uint64_t)0x100000000ULL; w0++) {
            Wn[0] = (uint32_t)w0;
            uint32_t de17 = wang_de17(Wn, Wf);
            td->total_evals++;
            if(de17 == 0 && td->result_count < td->max_results) {
                memcpy(td->results[td->result_count], Wn, 64);
                memcpy(td->results[td->result_count]+16, Wf, 64);
                td->result_count++;
                break; /* next set */
            }
        }
    }
    return NULL;
}

int main(int argc, char **argv) {
    int num_sets = 5000;
    int num_threads = 4;
    if(argc > 1) num_sets = atoi(argv[1]);
    if(argc > 2) num_threads = atoi(argv[2]);

    int max_per_thread = (num_sets / num_threads) + 2;

    fprintf(stderr, "Task12 pair generator: %d sets, %d threads\n", num_sets, num_threads);

    pthread_t threads[num_threads];
    thread_data_t tdata[num_threads];

    for(int t=0; t<num_threads; t++) {
        tdata[t].thread_id = t;
        tdata[t].num_threads = num_threads;
        tdata[t].num_sets = num_sets;
        tdata[t].seed = time(NULL) ^ ((uint64_t)t << 32) ^ 0xDEADBEEF12ULL;
        tdata[t].max_results = max_per_thread;
        tdata[t].results = malloc(max_per_thread * 32 * sizeof(uint32_t));
        tdata[t].result_count = 0;
        tdata[t].total_evals = 0;
        pthread_create(&threads[t], NULL, search_thread, &tdata[t]);
    }

    for(int t=0; t<num_threads; t++) {
        pthread_join(threads[t], NULL);
    }

    /* Output all pairs */
    int total_pairs = 0;
    uint64_t total_evals = 0;
    for(int t=0; t<num_threads; t++) {
        total_evals += tdata[t].total_evals;
        for(int i=0; i<tdata[t].result_count; i++) {
            printf("PAIR\n");
            printf("Wn =");
            for(int w=0;w<16;w++) printf(" %08x", tdata[t].results[i][w]);
            printf("\n");
            printf("Wf =");
            for(int w=0;w<16;w++) printf(" %08x", tdata[t].results[i][16+w]);
            printf("\n");
            total_pairs++;
        }
        free(tdata[t].results);
    }

    fprintf(stderr, "Found %d pairs from %d sets (%.2e evaluations)\n",
            total_pairs, num_sets, (double)total_evals);
    return 0;
}
