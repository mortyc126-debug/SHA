/*
 * ЗАДАНИЕ 8: 17-round SHA-256 Wang collision search
 * Brute-force Wn[0] to find δe[17]=0 with Wang chain.
 *
 * Compile: gcc -O3 -march=native -o task8_search task8_search.c -lpthread
 * Usage:   ./task8_search [num_pairs] [num_threads]
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#define MASK 0xFFFFFFFFU

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

/* Full 64-round SHA-256 compression (state only, no IV add) */
void sha256_compress(const uint32_t M[16], const uint32_t iv[8],
                     uint32_t state_out[8]) {
    uint32_t W[64];
    for (int i = 0; i < 16; i++) W[i] = M[i];
    for (int i = 16; i < 64; i++)
        W[i] = sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16];

    uint32_t a=iv[0],b=iv[1],c=iv[2],d=iv[3];
    uint32_t e=iv[4],f=iv[5],g=iv[6],h=iv[7];
    for (int r = 0; r < 64; r++) {
        uint32_t T1 = h + Sig1(e) + Ch(e,f,g) + K[r] + W[r];
        uint32_t T2 = Sig0(a) + Maj(a,b,c);
        h=g; g=f; f=e; e=d+T1; d=c; c=b; b=a; a=T1+T2;
    }
    state_out[0]=a; state_out[1]=b; state_out[2]=c; state_out[3]=d;
    state_out[4]=e; state_out[5]=f; state_out[6]=g; state_out[7]=h;
}

/* SHA-256 hash = compress + IV add */
void sha256_hash(const uint32_t M[16], uint32_t hash_out[8]) {
    sha256_compress(M, H0, hash_out);
    for (int i = 0; i < 8; i++) hash_out[i] += H0[i];
}

/*
 * Wang chain: given IV and Wn[0..15], compute Wf[0..15] such that δe[2..16]=0.
 * Returns δe[17] (the first uncontrolled difference).
 * Also returns state_n[17] and state_f[17] if pointers provided.
 */
uint32_t wang_chain_de17(const uint32_t iv[8],
                         const uint32_t Wn[16],
                         uint32_t Wf_out[16]) {
    /* Copy Wn to Wf, flip bit 15 of W[0] */
    for (int i = 0; i < 16; i++) Wf_out[i] = Wn[i];
    Wf_out[0] = Wn[0] ^ 0x8000U;

    /* State arrays: [a,b,c,d,e,f,g,h] */
    uint32_t sn[8], sf[8];
    for (int i = 0; i < 8; i++) { sn[i] = iv[i]; sf[i] = iv[i]; }

    /* Round 0 */
    {
        uint32_t T1n = sn[7] + Sig1(sn[4]) + Ch(sn[4],sn[5],sn[6]) + K[0] + Wn[0];
        uint32_t T2n = Sig0(sn[0]) + Maj(sn[0],sn[1],sn[2]);
        uint32_t T1f = sf[7] + Sig1(sf[4]) + Ch(sf[4],sf[5],sf[6]) + K[0] + Wf_out[0];
        uint32_t T2f = Sig0(sf[0]) + Maj(sf[0],sf[1],sf[2]);
        uint32_t new_sn[8] = {T1n+T2n, sn[0], sn[1], sn[2], sn[3]+T1n, sn[4], sn[5], sn[6]};
        uint32_t new_sf[8] = {T1f+T2f, sf[0], sf[1], sf[2], sf[3]+T1f, sf[4], sf[5], sf[6]};
        for (int i = 0; i < 8; i++) { sn[i] = new_sn[i]; sf[i] = new_sf[i]; }
    }

    /* Rounds 1..15: adaptive δW[r] to ensure δe[r+1]=0 */
    for (int r = 1; r < 16; r++) {
        uint32_t dd = sf[3] - sn[3];
        uint32_t dh = sf[7] - sn[7];
        uint32_t dSig1 = Sig1(sf[4]) - Sig1(sn[4]);
        uint32_t dCh = Ch(sf[4],sf[5],sf[6]) - Ch(sn[4],sn[5],sn[6]);
        uint32_t dWr = (uint32_t)(-(int32_t)(dd + dh + dSig1 + dCh));
        Wf_out[r] = Wn[r] + dWr;

        uint32_t T1n = sn[7] + Sig1(sn[4]) + Ch(sn[4],sn[5],sn[6]) + K[r] + Wn[r];
        uint32_t T2n = Sig0(sn[0]) + Maj(sn[0],sn[1],sn[2]);
        uint32_t T1f = sf[7] + Sig1(sf[4]) + Ch(sf[4],sf[5],sf[6]) + K[r] + Wf_out[r];
        uint32_t T2f = Sig0(sf[0]) + Maj(sf[0],sf[1],sf[2]);
        uint32_t new_sn[8] = {T1n+T2n, sn[0], sn[1], sn[2], sn[3]+T1n, sn[4], sn[5], sn[6]};
        uint32_t new_sf[8] = {T1f+T2f, sf[0], sf[1], sf[2], sf[3]+T1f, sf[4], sf[5], sf[6]};
        for (int i = 0; i < 8; i++) { sn[i] = new_sn[i]; sf[i] = new_sf[i]; }
    }

    /* Now compute round 16 to get e[17] */
    /* Need W[16] from schedule */
    uint32_t Wn16 = sig1(Wn[14]) + Wn[9] + sig0(Wn[1]) + Wn[0];
    uint32_t Wf16 = sig1(Wf_out[14]) + Wf_out[9] + sig0(Wf_out[1]) + Wf_out[0];

    uint32_t T1n = sn[7] + Sig1(sn[4]) + Ch(sn[4],sn[5],sn[6]) + K[16] + Wn16;
    uint32_t T1f = sf[7] + Sig1(sf[4]) + Ch(sf[4],sf[5],sf[6]) + K[16] + Wf16;

    uint32_t en17 = sn[3] + T1n;
    uint32_t ef17 = sf[3] + T1f;

    return ef17 - en17;  /* δe[17] */
}

/* Result structure */
typedef struct {
    uint32_t Wn[16];
    uint32_t Wf[16];
    uint32_t w0_found;
    uint64_t iters;
    int found;
} search_result_t;

/* Thread argument */
typedef struct {
    uint32_t Wn_base[16];  /* Wn[1..15] fixed, Wn[0] varies */
    uint64_t start;
    uint64_t end;  /* exclusive */
    search_result_t result;
    volatile int *global_found;
} thread_arg_t;

void *search_thread(void *arg) {
    thread_arg_t *ta = (thread_arg_t *)arg;
    uint32_t Wn[16], Wf[16];
    for (int i = 0; i < 16; i++) Wn[i] = ta->Wn_base[i];

    ta->result.found = 0;
    ta->result.iters = 0;

    for (uint64_t w0 = ta->start; w0 < ta->end; w0++) {
        if (*(ta->global_found)) break;

        Wn[0] = (uint32_t)w0;
        uint32_t de17 = wang_chain_de17(H0, Wn, Wf);
        ta->result.iters++;

        if (de17 == 0) {
            ta->result.found = 1;
            ta->result.w0_found = (uint32_t)w0;
            for (int i = 0; i < 16; i++) {
                ta->result.Wn[i] = Wn[i];
                ta->result.Wf[i] = Wf[i];
            }
            *(ta->global_found) = 1;
            break;
        }
    }
    return NULL;
}

/* Generate random uint32_t */
uint32_t rand32(void) {
    uint32_t r;
    FILE *f = fopen("/dev/urandom", "rb");
    fread(&r, 4, 1, f);
    fclose(f);
    return r;
}

void rand_words(uint32_t *buf, int n) {
    FILE *f = fopen("/dev/urandom", "rb");
    fread(buf, 4, n, f);
    fclose(f);
}

int main(int argc, char **argv) {
    int num_pairs = 10;
    int num_threads = 4;
    if (argc > 1) num_pairs = atoi(argv[1]);
    if (argc > 2) num_threads = atoi(argv[2]);

    printf("=== TASK 8: 17-round Wang collision search ===\n");
    printf("Target: %d pairs, %d threads\n\n", num_pairs, num_threads);

    for (int pair = 0; pair < num_pairs; pair++) {
        printf("--- Pair %d/%d ---\n", pair+1, num_pairs);

        /* Random Wn[1..15] */
        uint32_t Wn_base[16];
        rand_words(Wn_base, 16);
        /* Wn[0] will be varied by threads */

        volatile int global_found = 0;
        pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
        thread_arg_t *args = malloc(num_threads * sizeof(thread_arg_t));

        /* Split 2^32 range among threads */
        uint64_t total = (uint64_t)1 << 32;
        uint64_t chunk = total / num_threads;

        struct timespec ts_start, ts_end;
        clock_gettime(CLOCK_MONOTONIC, &ts_start);

        for (int t = 0; t < num_threads; t++) {
            for (int i = 0; i < 16; i++) args[t].Wn_base[i] = Wn_base[i];
            args[t].start = (uint64_t)t * chunk;
            args[t].end = (t == num_threads-1) ? total : (uint64_t)(t+1) * chunk;
            args[t].global_found = &global_found;
            pthread_create(&threads[t], NULL, search_thread, &args[t]);
        }

        search_result_t best;
        best.found = 0;
        uint64_t total_iters = 0;

        for (int t = 0; t < num_threads; t++) {
            pthread_join(threads[t], NULL);
            total_iters += args[t].result.iters;
            if (args[t].result.found && !best.found) {
                best = args[t].result;
            }
        }

        clock_gettime(CLOCK_MONOTONIC, &ts_end);
        double elapsed = (ts_end.tv_sec - ts_start.tv_sec)
                       + (ts_end.tv_nsec - ts_start.tv_nsec) * 1e-9;
        double rate = total_iters / elapsed / 1e6;

        if (best.found) {
            printf("FOUND! Wn[0] = 0x%08x after %lu iters (%.1fs, %.1f M/s)\n",
                   best.w0_found, total_iters, elapsed, rate);

            /* Verify: compute full hashes */
            uint32_t hash_n[8], hash_f[8];
            sha256_hash(best.Wn, hash_n);
            sha256_hash(best.Wf, hash_f);

            /* Verify δe[2..17]=0 */
            uint32_t Wn_full[64], Wf_full[64];
            for (int i = 0; i < 16; i++) { Wn_full[i] = best.Wn[i]; Wf_full[i] = best.Wf[i]; }
            for (int i = 16; i < 64; i++) {
                Wn_full[i] = sig1(Wn_full[i-2]) + Wn_full[i-7] + sig0(Wn_full[i-15]) + Wn_full[i-16];
                Wf_full[i] = sig1(Wf_full[i-2]) + Wf_full[i-7] + sig0(Wf_full[i-15]) + Wf_full[i-16];
            }
            uint32_t sn[8], sf[8];
            for (int i = 0; i < 8; i++) { sn[i] = H0[i]; sf[i] = H0[i]; }
            int de_zero_count = 0;
            for (int r = 0; r < 64; r++) {
                uint32_t T1n = sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[r]+Wn_full[r];
                uint32_t T2n = Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
                uint32_t T1f = sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[r]+Wf_full[r];
                uint32_t T2f = Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);
                uint32_t nsn[8] = {T1n+T2n,sn[0],sn[1],sn[2],sn[3]+T1n,sn[4],sn[5],sn[6]};
                uint32_t nsf[8] = {T1f+T2f,sf[0],sf[1],sf[2],sf[3]+T1f,sf[4],sf[5],sf[6]};
                for (int i=0;i<8;i++){sn[i]=nsn[i];sf[i]=nsf[i];}
                uint32_t de = sf[4] - sn[4];
                if (r >= 1 && r <= 16 && de == 0) de_zero_count++;
                if (r >= 16 && r <= 20) {
                    printf("  δe[%d] = 0x%08x (HW=%d)\n", r+1, de, popcount32(de));
                }
            }
            printf("  δe[2..17] all zero: %s (%d/16)\n",
                   de_zero_count >= 16 ? "YES" : "NO", de_zero_count);

            /* Hash difference */
            int hw = 0;
            printf("  δhash = ");
            for (int i = 0; i < 8; i++) {
                uint32_t d = hash_n[i] ^ hash_f[i];
                hw += popcount32(d);
                printf("%08x ", d);
            }
            printf("\n  HW(δhash) = %d/256\n", hw);

            /* Print Wn */
            printf("  Wn = ");
            for (int i = 0; i < 16; i++) printf("%08x ", best.Wn[i]);
            printf("\n  Wf = ");
            for (int i = 0; i < 16; i++) printf("%08x ", best.Wf[i]);
            printf("\n\n");
        } else {
            printf("NOT FOUND in %lu iters (%.1fs)\n\n", total_iters, elapsed);
        }

        free(threads);
        free(args);
    }

    return 0;
}
