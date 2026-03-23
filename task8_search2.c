/*
 * ЗАДАНИЕ 8: 17-round SHA-256 Wang collision search (v2)
 * Strategy: sweep Wn[0] with multiple Wn[1..15] sets per thread.
 * Each thread handles different Wn[0] ranges, all share the same
 * set of N_SETS random Wn[1..15] backgrounds.
 *
 * Compile: gcc -O3 -march=native -o task8_search2 task8_search2.c -lpthread
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#define MASK 0xFFFFFFFFU
#define N_SETS 10  /* number of Wn[1..15] background sets */

static const uint32_t K[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,
    0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,
    0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,
    0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,
    0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,
    0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,
    0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,
    0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,
    0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
};

static const uint32_t H0[8] = {
    0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
    0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19,
};

static inline uint32_t rotr(uint32_t x, int n) {
    return (x >> n) | (x << (32 - n));
}
static inline uint32_t sig0(uint32_t x) { return rotr(x,7)^rotr(x,18)^(x>>3); }
static inline uint32_t sig1(uint32_t x) { return rotr(x,17)^rotr(x,19)^(x>>10); }
static inline uint32_t Sig0(uint32_t x) { return rotr(x,2)^rotr(x,13)^rotr(x,22); }
static inline uint32_t Sig1(uint32_t x) { return rotr(x,6)^rotr(x,11)^rotr(x,25); }
static inline uint32_t Ch(uint32_t e,uint32_t f,uint32_t g) { return (e&f)^((~e)&g); }
static inline uint32_t Maj(uint32_t a,uint32_t b,uint32_t c) { return (a&b)^(a&c)^(b&c); }
static inline int popcount32(uint32_t x) { return __builtin_popcount(x); }

/*
 * Wang chain check for δe[17]=0.
 * Returns δe[17]. Writes Wf[0..15] if Wf_out != NULL.
 */
static uint32_t wang_de17(const uint32_t Wn[16], uint32_t Wf_out[16]) {
    uint32_t Wf[16];
    for (int i = 0; i < 16; i++) Wf[i] = Wn[i];
    Wf[0] = Wn[0] ^ 0x8000U;

    uint32_t sn[8], sf[8];
    for (int i = 0; i < 8; i++) { sn[i] = H0[i]; sf[i] = H0[i]; }

    /* Round 0 */
    uint32_t T1n = sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[0]+Wn[0];
    uint32_t T2n = Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
    uint32_t T1f = sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[0]+Wf[0];
    uint32_t T2f = Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);
    {
        uint32_t a=T1n+T2n,b=sn[0],c=sn[1],d=sn[2],e=sn[3]+T1n,f2=sn[4],g=sn[5],h=sn[6];
        sn[0]=a;sn[1]=b;sn[2]=c;sn[3]=d;sn[4]=e;sn[5]=f2;sn[6]=g;sn[7]=h;
    }
    {
        uint32_t a=T1f+T2f,b=sf[0],c=sf[1],d=sf[2],e=sf[3]+T1f,f2=sf[4],g=sf[5],h=sf[6];
        sf[0]=a;sf[1]=b;sf[2]=c;sf[3]=d;sf[4]=e;sf[5]=f2;sf[6]=g;sf[7]=h;
    }

    /* Rounds 1..15 */
    for (int r = 1; r < 16; r++) {
        uint32_t dd=sf[3]-sn[3], dh=sf[7]-sn[7];
        uint32_t dSig1=Sig1(sf[4])-Sig1(sn[4]);
        uint32_t dCh=Ch(sf[4],sf[5],sf[6])-Ch(sn[4],sn[5],sn[6]);
        Wf[r] = Wn[r] + (uint32_t)(-(int32_t)(dd+dh+dSig1+dCh));

        T1n = sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[r]+Wn[r];
        T2n = Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
        T1f = sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[r]+Wf[r];
        T2f = Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);

        uint32_t a,b,c,d,e,f2,g,h;
        a=T1n+T2n;b=sn[0];c=sn[1];d=sn[2];e=sn[3]+T1n;f2=sn[4];g=sn[5];h=sn[6];
        sn[0]=a;sn[1]=b;sn[2]=c;sn[3]=d;sn[4]=e;sn[5]=f2;sn[6]=g;sn[7]=h;
        a=T1f+T2f;b=sf[0];c=sf[1];d=sf[2];e=sf[3]+T1f;f2=sf[4];g=sf[5];h=sf[6];
        sf[0]=a;sf[1]=b;sf[2]=c;sf[3]=d;sf[4]=e;sf[5]=f2;sf[6]=g;sf[7]=h;
    }

    /* Round 16: compute δe[17] */
    uint32_t Wn16=sig1(Wn[14])+Wn[9]+sig0(Wn[1])+Wn[0];
    uint32_t Wf16=sig1(Wf[14])+Wf[9]+sig0(Wf[1])+Wf[0];
    T1n = sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[16]+Wn16;
    T1f = sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[16]+Wf16;
    uint32_t en17=sn[3]+T1n, ef17=sf[3]+T1f;

    if (Wf_out) {
        for (int i = 0; i < 16; i++) Wf_out[i] = Wf[i];
    }
    return ef17 - en17;
}

/* Full hash for verification */
void sha256_hash(const uint32_t M[16], uint32_t hash_out[8]) {
    uint32_t W[64];
    for (int i=0;i<16;i++) W[i]=M[i];
    for (int i=16;i<64;i++) W[i]=sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16];
    uint32_t a=H0[0],b=H0[1],c=H0[2],d=H0[3],e=H0[4],f=H0[5],g=H0[6],h=H0[7];
    for (int r=0;r<64;r++){
        uint32_t T1=h+Sig1(e)+Ch(e,f,g)+K[r]+W[r];
        uint32_t T2=Sig0(a)+Maj(a,b,c);
        h=g;g=f;f=e;e=d+T1;d=c;c=b;b=a;a=T1+T2;
    }
    hash_out[0]=a+H0[0];hash_out[1]=b+H0[1];hash_out[2]=c+H0[2];hash_out[3]=d+H0[3];
    hash_out[4]=e+H0[4];hash_out[5]=f+H0[5];hash_out[6]=g+H0[6];hash_out[7]=h+H0[7];
}

/* Verify and print a found pair */
void verify_pair(const uint32_t Wn[16], const uint32_t Wf[16], int pair_id) {
    /* Full round-by-round verification */
    uint32_t Wn_full[64], Wf_full[64];
    for(int i=0;i<16;i++){Wn_full[i]=Wn[i];Wf_full[i]=Wf[i];}
    for(int i=16;i<64;i++){
        Wn_full[i]=sig1(Wn_full[i-2])+Wn_full[i-7]+sig0(Wn_full[i-15])+Wn_full[i-16];
        Wf_full[i]=sig1(Wf_full[i-2])+Wf_full[i-7]+sig0(Wf_full[i-15])+Wf_full[i-16];
    }
    uint32_t sn[8],sf[8];
    for(int i=0;i<8;i++){sn[i]=H0[i];sf[i]=H0[i];}
    int de_zero=0;
    printf("  Round-by-round δe:\n");
    for(int r=0;r<64;r++){
        uint32_t T1n=sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[r]+Wn_full[r];
        uint32_t T2n=Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
        uint32_t T1f=sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[r]+Wf_full[r];
        uint32_t T2f=Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);
        uint32_t nsn[8]={T1n+T2n,sn[0],sn[1],sn[2],sn[3]+T1n,sn[4],sn[5],sn[6]};
        uint32_t nsf[8]={T1f+T2f,sf[0],sf[1],sf[2],sf[3]+T1f,sf[4],sf[5],sf[6]};
        for(int i=0;i<8;i++){sn[i]=nsn[i];sf[i]=nsf[i];}
        uint32_t de=sf[4]-sn[4];
        if(r>=1 && r<=16 && de==0) de_zero++;
        if(r<=20 || r==63)
            printf("    δe[%2d]=0x%08x HW=%d%s\n",r+1,de,popcount32(de),de==0?" <<<":"");
    }
    printf("  δe[2..17] zeros: %d/16\n", de_zero);

    /* Hash difference */
    uint32_t hn[8],hf[8];
    sha256_hash(Wn,hn); sha256_hash(Wf,hf);
    int hw=0;
    printf("  hash_n = ");
    for(int i=0;i<8;i++) printf("%08x",hn[i]);
    printf("\n  hash_f = ");
    for(int i=0;i<8;i++) printf("%08x",hf[i]);
    printf("\n  δhash  = ");
    for(int i=0;i<8;i++){uint32_t d=hn[i]^hf[i];hw+=popcount32(d);printf("%08x",d);}
    printf("\n  HW(δhash) = %d/256\n", hw);

    /* Print messages */
    printf("  Wn = ");
    for(int i=0;i<16;i++) printf("%08x ",Wn[i]);
    printf("\n  Wf = ");
    for(int i=0;i<16;i++) printf("%08x ",Wf[i]);
    printf("\n");
}

/* Shared data */
typedef struct {
    uint32_t sets[N_SETS][16];  /* Wn[0..15] templates (Wn[0] unused) */
    int found[N_SETS];
    uint32_t found_Wn[N_SETS][16];
    uint32_t found_Wf[N_SETS][16];
    uint32_t found_w0[N_SETS];
    pthread_mutex_t mutex;
    int total_found;
} shared_t;

typedef struct {
    shared_t *shared;
    uint64_t start, end;
    uint64_t iters;
} targ_t;

void *worker(void *arg) {
    targ_t *ta = (targ_t *)arg;
    shared_t *sh = ta->shared;
    ta->iters = 0;

    uint32_t Wn[16], Wf[16];

    for (uint64_t w0 = ta->start; w0 < ta->end; w0++) {
        if (sh->total_found >= N_SETS) break;

        for (int s = 0; s < N_SETS; s++) {
            if (sh->found[s]) continue;

            for (int i = 1; i < 16; i++) Wn[i] = sh->sets[s][i];
            Wn[0] = (uint32_t)w0;
            uint32_t de17 = wang_de17(Wn, Wf);
            ta->iters++;

            if (de17 == 0) {
                pthread_mutex_lock(&sh->mutex);
                if (!sh->found[s]) {
                    sh->found[s] = 1;
                    sh->found_w0[s] = (uint32_t)w0;
                    for (int i=0;i<16;i++){sh->found_Wn[s][i]=Wn[i];sh->found_Wf[s][i]=Wf[i];}
                    sh->total_found++;
                    fprintf(stderr, "  [HIT] Set %d: w0=0x%08x (total %d/%d)\n",
                            s, (uint32_t)w0, sh->total_found, N_SETS);
                }
                pthread_mutex_unlock(&sh->mutex);
            }
        }
    }
    return NULL;
}

int main(int argc, char **argv) {
    int num_threads = 4;
    if (argc > 1) num_threads = atoi(argv[1]);

    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    printf("=== TASK 8: 17-round Wang collision search (v2) ===\n");
    printf("Searching %d sets × 2^32 candidates, %d threads\n", N_SETS, num_threads);
    printf("Expected hits: ~%d (Poisson, some sets may miss)\n\n", N_SETS);

    shared_t shared;
    memset(&shared, 0, sizeof(shared));
    pthread_mutex_init(&shared.mutex, NULL);

    /* Generate N_SETS random backgrounds */
    FILE *rng = fopen("/dev/urandom", "rb");
    for (int s = 0; s < N_SETS; s++) {
        fread(shared.sets[s], 4, 16, rng);
    }
    fclose(rng);

    uint64_t total = (uint64_t)1 << 32;
    uint64_t chunk = total / num_threads;

    pthread_t threads[64];
    targ_t targs[64];

    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);

    for (int t = 0; t < num_threads; t++) {
        targs[t].shared = &shared;
        targs[t].start = (uint64_t)t * chunk;
        targs[t].end = (t == num_threads-1) ? total : (uint64_t)(t+1)*chunk;
        targs[t].iters = 0;
        pthread_create(&threads[t], NULL, worker, &targs[t]);
    }

    for (int t = 0; t < num_threads; t++)
        pthread_join(threads[t], NULL);

    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double elapsed = (ts1.tv_sec-ts0.tv_sec)+(ts1.tv_nsec-ts0.tv_nsec)*1e-9;

    uint64_t total_iters = 0;
    for (int t = 0; t < num_threads; t++) total_iters += targs[t].iters;
    double rate = total_iters / elapsed / 1e6;

    printf("\n=== RESULTS ===\n");
    printf("Total iterations: %lu (%.1fs, %.1f M/s)\n", total_iters, elapsed, rate);
    printf("Pairs found: %d/%d\n\n", shared.total_found, N_SETS);

    for (int s = 0; s < N_SETS; s++) {
        if (shared.found[s]) {
            printf("--- PAIR %d (w0=0x%08x) ---\n", s, shared.found_w0[s]);
            verify_pair(shared.found_Wn[s], shared.found_Wf[s], s);
            printf("\n");
        }
    }

    /* Summary */
    if (shared.total_found > 0) {
        printf("=== SUMMARY ===\n");
        int min_hw = 999, max_hw = 0;
        for (int s = 0; s < N_SETS; s++) {
            if (!shared.found[s]) continue;
            uint32_t hn[8],hf[8];
            sha256_hash(shared.found_Wn[s],hn);
            sha256_hash(shared.found_Wf[s],hf);
            int hw=0;
            for(int i=0;i<8;i++) hw+=popcount32(hn[i]^hf[i]);
            printf("  Set %d: HW(δhash) = %d\n", s, hw);
            if(hw<min_hw) min_hw=hw;
            if(hw>max_hw) max_hw=hw;
        }
        printf("  Min HW = %d, Max HW = %d\n", min_hw, max_hw);
        printf("  Random expectation: 128\n");
    }

    pthread_mutex_destroy(&shared.mutex);
    return 0;
}
