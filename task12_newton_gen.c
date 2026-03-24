/*
 * ЗАДАНИЕ 12: Newton/Hensel pair generator
 * Uses iterative refinement to solve δe[17](W0) = 0 for each background.
 * Much faster than brute force sweep.
 *
 * gcc -O3 -march=native -o task12_newton_gen task12_newton_gen.c -lpthread
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#define MASK 0xFFFFFFFFU

static const uint32_t K_const[64] = {
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

/* Compute δe[17] for given Wn[0..15] with Wf[0] = Wn[0]^0x8000 */
static uint32_t wang_de17(const uint32_t Wn[16], uint32_t Wf[16]) {
    for(int i=0;i<16;i++) Wf[i]=Wn[i];
    Wf[0] = Wn[0] ^ 0x8000;
    uint32_t sn[8], sf[8];
    for(int i=0;i<8;i++){sn[i]=H0[i];sf[i]=H0[i];}

    uint32_t T1n=sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K_const[0]+Wn[0];
    uint32_t T2n=Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
    uint32_t T1f=sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K_const[0]+Wf[0];
    uint32_t T2f=Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);
    {uint32_t t[8]={T1n+T2n,sn[0],sn[1],sn[2],sn[3]+T1n,sn[4],sn[5],sn[6]};memcpy(sn,t,32);}
    {uint32_t t[8]={T1f+T2f,sf[0],sf[1],sf[2],sf[3]+T1f,sf[4],sf[5],sf[6]};memcpy(sf,t,32);}

    for(int r=1;r<16;r++){
        uint32_t dd=sf[3]-sn[3],dh=sf[7]-sn[7];
        uint32_t dS=Sig1(sf[4])-Sig1(sn[4]),dC=Ch(sf[4],sf[5],sf[6])-Ch(sn[4],sn[5],sn[6]);
        Wf[r]=Wn[r]+(uint32_t)(-(int32_t)(dd+dh+dS+dC));
        T1n=sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K_const[r]+Wn[r];
        T2n=Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
        T1f=sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K_const[r]+Wf[r];
        T2f=Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);
        {uint32_t t[8]={T1n+T2n,sn[0],sn[1],sn[2],sn[3]+T1n,sn[4],sn[5],sn[6]};memcpy(sn,t,32);}
        {uint32_t t[8]={T1f+T2f,sf[0],sf[1],sf[2],sf[3]+T1f,sf[4],sf[5],sf[6]};memcpy(sf,t,32);}
    }
    uint32_t Wn16=sig1(Wn[14])+Wn[9]+sig0(Wn[1])+Wn[0];
    uint32_t Wf16=sig1(Wf[14])+Wf[9]+sig0(Wf[1])+Wf[0];
    T1n=sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K_const[16]+Wn16;
    T1f=sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K_const[16]+Wf16;
    return (sf[3]+T1f)-(sn[3]+T1n);
}

/* Modular inverse mod 2^32 using extended Euclidean / Hensel lifting */
static uint32_t mod_inv(uint32_t a) {
    /* a must be odd. Returns a^{-1} mod 2^32. */
    uint32_t x = a; /* a*a ≡ 1 mod 2 if a is odd, start with x=a */
    /* Newton: x = x*(2 - a*x) mod 2^k, doubling k each step */
    for(int i = 0; i < 5; i++) {
        x = x * (2 - a * x);
    }
    return x;
}

/* Solve f(W0) = δe[17](bg, W0) = 0 using Hensel lifting bit-by-bit.
 * Returns 1 if found, 0 if failed. Sets *out_w0. */
static int hensel_solve(const uint32_t bg[16], uint32_t *out_w0, uint64_t *evals) {
    uint32_t Wn[16], Wf[16];
    for(int i=1;i<16;i++) Wn[i]=bg[i];

    /* Try Hensel lifting: solve f(x) = 0 mod 2^k for k=1,2,...,32
     * f(x) mod 2 depends only on x mod 2.
     * Start with x=0 and x=1, pick whichever gives f(x)≡0 mod 2.
     * Then lift: x_{k+1} = x_k or x_k + 2^k, pick the one with f≡0 mod 2^{k+1}. */
    uint32_t x = 0;
    for(int k = 0; k < 32; k++) {
        /* Test x */
        Wn[0] = x;
        uint32_t f0 = wang_de17(Wn, Wf);
        (*evals)++;

        /* Check if f(x) ≡ 0 mod 2^{k+1} */
        uint32_t mask = (k < 31) ? ((1U << (k+1)) - 1) : 0xFFFFFFFFU;
        if((f0 & mask) == 0) {
            /* x is good, continue */
        } else {
            /* Try x + 2^k */
            x ^= (1U << k);
            Wn[0] = x;
            uint32_t f1 = wang_de17(Wn, Wf);
            (*evals)++;

            if((f1 & mask) != 0) {
                /* Neither works — function is not 2-adically smooth at this point.
                 * This can happen due to XOR/rotation operations.
                 * Fallback: try both branches and see if either leads to a solution. */
                /* We'll just continue with whichever has more trailing zeros */
                int tz0 = __builtin_ctz(f0 | (1U<<31));
                int tz1 = __builtin_ctz(f1 | (1U<<31));
                if(tz0 > tz1) {
                    x ^= (1U << k); /* revert */
                }
                /* Continue anyway — may or may not converge */
            }
        }
    }

    /* Final check */
    Wn[0] = x;
    uint32_t final_f = wang_de17(Wn, Wf);
    (*evals)++;
    if(final_f == 0) {
        *out_w0 = x;
        return 1;
    }
    return 0;
}

/* Multi-start Hensel with random perturbations */
static int multi_hensel_solve(const uint32_t bg[16], uint32_t *out_w0, uint64_t *evals, int max_starts) {
    uint32_t Wn[16], Wf[16];
    for(int i=1;i<16;i++) Wn[i]=bg[i];

    /* Strategy 1: pure Hensel from bit 0 */
    if(hensel_solve(bg, out_w0, evals)) return 1;

    /* Strategy 2: try random starting points with Newton correction */
    /* f(x) = 0: Newton step x' = x - f(x) * (f'(x))^{-1} mod 2^32
     * f'(x) ≈ f(x+1) - f(x) */
    for(int start = 0; start < max_starts; start++) {
        /* Random start */
        uint32_t x;
        {
            /* Simple LCG-like random */
            static uint32_t seed = 12345;
            seed = seed * 1103515245 + 12345 + start * 7919;
            x = seed ^ (uint32_t)start;
        }

        for(int iter = 0; iter < 100; iter++) {
            Wn[0] = x;
            uint32_t fx = wang_de17(Wn, Wf);
            (*evals)++;
            if(fx == 0) { *out_w0 = x; return 1; }

            Wn[0] = x + 1;
            uint32_t fx1 = wang_de17(Wn, Wf);
            (*evals)++;
            uint32_t deriv = fx1 - fx;

            if(deriv & 1) {
                /* Derivative is odd = invertible mod 2^32 */
                uint32_t inv = mod_inv(deriv);
                x = x - fx * inv;
            } else {
                /* Not invertible, perturb and try again */
                x += 0x10001 * (iter + 1);
            }
        }
    }
    return 0;
}

/* Thread data */
typedef struct {
    int thread_id, num_threads;
    int target;
    uint32_t (*backgrounds)[16];
    int bg_count;
    /* Output buffer */
    uint32_t (*out_Wn)[16];
    uint32_t (*out_Wf)[16];
    int found_count;
    int max_output;
    uint64_t total_evals;
    volatile int *global_found;
    pthread_mutex_t *mutex;
} targ_t;

static void *worker(void *arg) {
    targ_t *ta = (targ_t *)arg;
    ta->found_count = 0;
    ta->total_evals = 0;

    for(int i = ta->thread_id; i < ta->bg_count && *ta->global_found < ta->target; i += ta->num_threads) {
        uint32_t w0;
        uint64_t evals = 0;
        int ok = multi_hensel_solve(ta->backgrounds[i], &w0, &evals, 200);
        ta->total_evals += evals;

        if(ok && ta->found_count < ta->max_output) {
            uint32_t Wn[16], Wf[16];
            for(int j=1;j<16;j++) Wn[j] = ta->backgrounds[i][j];
            Wn[0] = w0;
            wang_de17(Wn, Wf);

            pthread_mutex_lock(ta->mutex);
            memcpy(ta->out_Wn[ta->found_count], Wn, 64);
            memcpy(ta->out_Wf[ta->found_count], Wf, 64);
            ta->found_count++;
            (*ta->global_found)++;
            if(*ta->global_found % 100 == 0 || *ta->global_found <= 20)
                fprintf(stderr, "\r  Found %d pairs (%.2eM evals)   ",
                        *ta->global_found, ta->total_evals/1e6);
            pthread_mutex_unlock(ta->mutex);
        }
    }
    return NULL;
}

int main(int argc, char **argv) {
    int target = 10000;
    int num_threads = 4;
    int bg_mult = 5; /* try bg_mult × target backgrounds */

    if(argc > 1) target = atoi(argv[1]);
    if(argc > 2) num_threads = atoi(argv[2]);
    if(argc > 3) bg_mult = atoi(argv[3]);

    int bg_count = target * bg_mult;

    fprintf(stderr, "Task12 Newton generator: target=%d, threads=%d, backgrounds=%d\n",
            target, num_threads, bg_count);

    /* Generate random backgrounds */
    uint32_t (*backgrounds)[16] = malloc(bg_count * 16 * sizeof(uint32_t));
    FILE *rng = fopen("/dev/urandom", "rb");
    for(int i = 0; i < bg_count; i++)
        fread(backgrounds[i], 4, 16, rng);
    fclose(rng);

    volatile int global_found = 0;
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);

    int max_per_thread = target; /* generous buffer */
    pthread_t threads[64];
    targ_t targs[64];

    for(int t = 0; t < num_threads; t++) {
        targs[t].thread_id = t;
        targs[t].num_threads = num_threads;
        targs[t].target = target;
        targs[t].backgrounds = backgrounds;
        targs[t].bg_count = bg_count;
        targs[t].out_Wn = malloc(max_per_thread * 16 * sizeof(uint32_t));
        targs[t].out_Wf = malloc(max_per_thread * 16 * sizeof(uint32_t));
        targs[t].found_count = 0;
        targs[t].max_output = max_per_thread;
        targs[t].total_evals = 0;
        targs[t].global_found = &global_found;
        targs[t].mutex = &mutex;
        pthread_create(&threads[t], NULL, worker, &targs[t]);
    }

    for(int t = 0; t < num_threads; t++)
        pthread_join(threads[t], NULL);

    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double elapsed = (ts1.tv_sec-ts0.tv_sec)+(ts1.tv_nsec-ts0.tv_nsec)*1e-9;

    uint64_t total_evals = 0;
    int total_pairs = 0;
    for(int t = 0; t < num_threads; t++) {
        total_evals += targs[t].total_evals;
        for(int i = 0; i < targs[t].found_count; i++) {
            printf("PAIR\n");
            printf("Wn =");
            for(int w = 0; w < 16; w++) printf(" %08x", targs[t].out_Wn[i][w]);
            printf("\n");
            printf("Wf =");
            for(int w = 0; w < 16; w++) printf(" %08x", targs[t].out_Wf[i][w]);
            printf("\n");
            total_pairs++;
        }
        free(targs[t].out_Wn);
        free(targs[t].out_Wf);
    }

    fprintf(stderr, "\nDone: %d pairs in %.1fs (%.2e evals, %.0f evals/pair)\n",
            total_pairs, elapsed, (double)total_evals,
            total_pairs > 0 ? (double)total_evals/total_pairs : 0);

    free(backgrounds);
    return 0;
}
