/*
 * ЗАДАНИЕ 11, Experiment A: SA maximizing carry_agreement(18..40)
 * High-performance C implementation.
 *
 * Outputs best pair data + 1000 random pairs for correlation analysis.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MASK 0xFFFFFFFFU
#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))

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

static inline uint32_t sig0(uint32_t x) { return ROTR(x,7)^ROTR(x,18)^(x>>3); }
static inline uint32_t sig1(uint32_t x) { return ROTR(x,17)^ROTR(x,19)^(x>>10); }
static inline uint32_t Sig0(uint32_t x) { return ROTR(x,2)^ROTR(x,13)^ROTR(x,22); }
static inline uint32_t Sig1(uint32_t x) { return ROTR(x,6)^ROTR(x,11)^ROTR(x,25); }
static inline uint32_t Ch(uint32_t e,uint32_t f,uint32_t g) { return (e&f)^((~e)&g); }
static inline uint32_t Maj(uint32_t a,uint32_t b,uint32_t c) { return (a&b)^(a&c)^(b&c); }
static inline int popcount32(uint32_t x) { return __builtin_popcount(x); }

/* Fast PRNG (xoshiro128**) */
static uint32_t rng_s[4];
static inline uint32_t rotl(uint32_t x, int k) { return (x<<k)|(x>>(32-k)); }
static uint32_t rng_next(void) {
    uint32_t r = rotl(rng_s[1]*5,7)*9;
    uint32_t t = rng_s[1]<<9;
    rng_s[2] ^= rng_s[0]; rng_s[3] ^= rng_s[1];
    rng_s[1] ^= rng_s[2]; rng_s[0] ^= rng_s[3];
    rng_s[2] ^= t; rng_s[3] = rotl(rng_s[3],11);
    return r;
}
static void rng_seed(uint64_t s) {
    rng_s[0] = (uint32_t)(s); rng_s[1] = (uint32_t)(s>>32);
    rng_s[2] = rng_s[0]^0x12345678; rng_s[3] = rng_s[1]^0x9abcdef0;
    for(int i=0;i<20;i++) rng_next();
}

static void schedule(const uint32_t M[16], uint32_t W[64]) {
    memcpy(W, M, 16*sizeof(uint32_t));
    for(int i=16;i<64;i++)
        W[i] = sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16];
}

/* Compute carry bits for x + y */
static uint32_t compute_carries(uint32_t x, uint32_t y) {
    uint32_t carries = 0;
    uint32_t c = 0;
    for(int i=0;i<32;i++) {
        uint32_t xi = (x>>i)&1, yi = (y>>i)&1;
        c = (xi&yi)|(xi&c)|(yi&c);
        if(c) carries |= (1u<<i);
    }
    return carries;
}

/* Wang chain: given Wn[0..15], produce Wf[0..15] with δe[r]=0 for r=2..17 */
static void wang_chain(const uint32_t Wn[16], uint32_t Wf[16]) {
    memcpy(Wf, Wn, 16*sizeof(uint32_t));
    Wf[0] = Wn[0] ^ 0x8000;

    uint32_t sn[8], sf[8];
    memcpy(sn, H0, 32); memcpy(sf, H0, 32);

    /* Round 0 */
    uint32_t T1n = sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[0]+Wn[0];
    uint32_t T2n = Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
    uint32_t T1f = sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[0]+Wf[0];
    uint32_t T2f = Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);
    uint32_t new_sn[8] = {T1n+T2n,sn[0],sn[1],sn[2],sn[3]+T1n,sn[4],sn[5],sn[6]};
    uint32_t new_sf[8] = {T1f+T2f,sf[0],sf[1],sf[2],sf[3]+T1f,sf[4],sf[5],sf[6]};
    memcpy(sn,new_sn,32); memcpy(sf,new_sf,32);

    /* Rounds 1..15: adapt Wf[r] to keep δe[r+1]=0 */
    for(int r=1;r<16;r++) {
        uint32_t dd = sf[3]-sn[3];
        uint32_t dh = sf[7]-sn[7];
        uint32_t dS = Sig1(sf[4])-Sig1(sn[4]);
        uint32_t dC = Ch(sf[4],sf[5],sf[6])-Ch(sn[4],sn[5],sn[6]);
        uint32_t dWr = -(dd+dh+dS+dC);
        Wf[r] = Wn[r]+dWr;

        T1n = sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[r]+Wn[r];
        T2n = Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
        T1f = sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[r]+Wf[r];
        T2f = Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);
        uint32_t ns[8] = {T1n+T2n,sn[0],sn[1],sn[2],sn[3]+T1n,sn[4],sn[5],sn[6]};
        uint32_t nf[8] = {T1f+T2f,sf[0],sf[1],sf[2],sf[3]+T1f,sf[4],sf[5],sf[6]};
        memcpy(sn,ns,32); memcpy(sf,nf,32);
    }
}

typedef struct {
    uint32_t states_n[65][8];
    uint32_t states_f[65][8];
    uint32_t carries_n[64][5]; /* 5 sub-additions per round */
    uint32_t carries_f[64][5];
    uint32_t T1_n[64], T1_f[64], T2_n[64], T2_f[64];
    uint32_t Ch_n[64], Ch_f[64], Maj_n[64], Maj_f[64];
    uint32_t Sig1_n[64], Sig1_f[64], Sig0_n[64], Sig0_f[64];
    uint32_t Wn_sched[64], Wf_sched[64];
    uint32_t hash_n[8], hash_f[8];
} FullData;

static void full_sha256_pair(const uint32_t Wn[16], const uint32_t Wf[16], FullData *d) {
    schedule(Wn, d->Wn_sched);
    schedule(Wf, d->Wf_sched);

    /* Init states */
    memcpy(d->states_n[0], H0, 32);
    memcpy(d->states_f[0], H0, 32);

    for(int which=0;which<2;which++) {
        uint32_t s[8];
        memcpy(s, H0, 32);
        uint32_t *W = (which==0) ? d->Wn_sched : d->Wf_sched;
        uint32_t (*states)[8] = (which==0) ? d->states_n : d->states_f;
        uint32_t (*carries)[5] = (which==0) ? d->carries_n : d->carries_f;
        uint32_t *T1arr = (which==0) ? d->T1_n : d->T1_f;
        uint32_t *T2arr = (which==0) ? d->T2_n : d->T2_f;
        uint32_t *Charr = (which==0) ? d->Ch_n : d->Ch_f;
        uint32_t *Marr = (which==0) ? d->Maj_n : d->Maj_f;
        uint32_t *S1arr = (which==0) ? d->Sig1_n : d->Sig1_f;
        uint32_t *S0arr = (which==0) ? d->Sig0_n : d->Sig0_f;
        uint32_t *hash = (which==0) ? d->hash_n : d->hash_f;

        memcpy(states[0], H0, 32);
        for(int r=0;r<64;r++) {
            uint32_t a=s[0],b=s[1],c=s[2],dd=s[3],e=s[4],f=s[5],g=s[6],h=s[7];
            uint32_t s1e = Sig1(e); S1arr[r] = s1e;
            uint32_t ch = Ch(e,f,g); Charr[r] = ch;
            uint32_t s0a = Sig0(a); S0arr[r] = s0a;
            uint32_t maj = Maj(a,b,c); Marr[r] = maj;

            uint32_t sum1 = h + s1e;
            uint32_t sum2 = sum1 + ch;
            uint32_t sum3 = sum2 + K[r];
            uint32_t sum4 = sum3 + W[r];
            uint32_t sum5 = s0a + maj;

            carries[r][0] = compute_carries(h, s1e);
            carries[r][1] = compute_carries(sum1, ch);
            carries[r][2] = compute_carries(sum2, K[r]);
            carries[r][3] = compute_carries(sum3, W[r]);
            carries[r][4] = compute_carries(s0a, maj);

            T1arr[r] = sum4;
            T2arr[r] = sum5;

            s[0] = sum4+sum5; s[1]=a; s[2]=b; s[3]=c;
            s[4] = dd+sum4; s[5]=e; s[6]=f; s[7]=g;
            memcpy(states[r+1], s, 32);
        }
        for(int i=0;i<8;i++) hash[i] = s[i] + H0[i];
    }
}

/* Compute carry_agreement for one round */
static double carry_agreement_round(const FullData *d, int r) {
    int agree = 0;
    for(int si=0;si<5;si++) {
        uint32_t dc = d->carries_n[r][si] ^ d->carries_f[r][si];
        agree += 32 - popcount32(dc);
    }
    return agree / 160.0;
}

/* Average carry_agreement over rounds r1..r2 (inclusive) */
static double avg_carry_agreement(const FullData *d, int r1, int r2) {
    double sum = 0;
    for(int r=r1;r<=r2;r++) sum += carry_agreement_round(d, r);
    return sum / (r2-r1+1);
}

/* HW of hash difference */
static int hw_dhash(const FullData *d) {
    int hw = 0;
    for(int i=0;i<8;i++) hw += popcount32(d->hash_n[i] ^ d->hash_f[i]);
    return hw;
}

/* Linear residual for round r */
static int linear_residual(const FullData *d, int r) {
    uint32_t dh = d->states_n[r][7] ^ d->states_f[r][7];
    uint32_t dSig1 = d->Sig1_n[r] ^ d->Sig1_f[r];
    uint32_t dCh = d->Ch_n[r] ^ d->Ch_f[r];
    uint32_t dW = d->Wn_sched[r] ^ d->Wf_sched[r];
    uint32_t dT1 = d->T1_n[r] ^ d->T1_f[r];
    uint32_t lin = dh ^ dSig1 ^ dCh ^ dW;
    return popcount32(dT1 ^ lin);
}

/* SA optimization */
static FullData best_data;
static double global_best = 0;
static uint32_t global_best_Wn[16], global_best_Wf[16];

static void sa_restart(int restart_id) {
    uint32_t Wn[16], Wf[16];
    FullData d;

    /* Random initial Wn */
    for(int i=0;i<16;i++) Wn[i] = rng_next();
    wang_chain(Wn, Wf);
    full_sha256_pair(Wn, Wf, &d);
    double score = avg_carry_agreement(&d, 18, 40);
    double best = score;
    uint32_t bestWn[16];
    memcpy(bestWn, Wn, 64);
    FullData bestD;
    memcpy(&bestD, &d, sizeof(d));

    double T = 0.05;
    double T_min = 0.0001;
    double alpha = pow(T_min/T, 1.0/200000.0);

    for(int step=0;step<200000;step++) {
        /* Mutate: flip 1 random bit in Wn[0..15] */
        uint32_t saved[16];
        memcpy(saved, Wn, 64);
        int word = rng_next() % 16;
        int bit = rng_next() % 32;
        Wn[word] ^= (1u << bit);

        wang_chain(Wn, Wf);
        full_sha256_pair(Wn, Wf, &d);
        double new_score = avg_carry_agreement(&d, 18, 40);

        double delta = new_score - score;
        if(delta > 0 || (double)(rng_next()&0xFFFFFF)/0xFFFFFF < exp(delta/T)) {
            score = new_score;
            if(score > best) {
                best = score;
                memcpy(bestWn, Wn, 64);
                memcpy(&bestD, &d, sizeof(d));
            }
        } else {
            memcpy(Wn, saved, 64);
        }
        T *= alpha;
    }

    /* Check global best */
    if(best > global_best) {
        global_best = best;
        memcpy(global_best_Wn, bestWn, 64);
        wang_chain(bestWn, global_best_Wf);
        memcpy(&best_data, &bestD, sizeof(bestD));
    }
    fprintf(stderr, "  Restart %2d: best=%.6f (global=%.6f)\n", restart_id, best, global_best);
}

int main(void) {
    rng_seed(time(NULL) ^ (uint64_t)42);

    fprintf(stderr, "=== Experiment A: SA for carry_agreement(18..40) ===\n");
    fprintf(stderr, "Running 50 restarts x 200K steps...\n");

    /* SA phase */
    global_best = 0;
    for(int r=0;r<50;r++) sa_restart(r);

    /* Output best pair data */
    printf("=== BEST PAIR ===\n");
    printf("carry_agreement(18..40) = %.6f\n", global_best);
    printf("Wn =");
    for(int i=0;i<16;i++) printf(" %08x", global_best_Wn[i]);
    printf("\n");
    printf("Wf =");
    for(int i=0;i<16;i++) printf(" %08x", global_best_Wf[i]);
    printf("\n");
    printf("HW(dhash) = %d\n", hw_dhash(&best_data));
    printf("\n");

    /* Full profiles for best pair */
    printf("=== CARRY AGREEMENT PROFILE (best pair) ===\n");
    for(int r=0;r<64;r++) {
        printf("r=%2d carry_agree=%.4f", r, carry_agreement_round(&best_data, r));
        printf(" lin_resid=%2d", linear_residual(&best_data, r));
        /* HW(dstate) */
        int hw_ds = 0;
        for(int i=0;i<8;i++) hw_ds += popcount32(best_data.states_n[r][i] ^ best_data.states_f[r][i]);
        printf(" HW(dstate)=%3d", hw_ds);
        /* HW(de) */
        printf(" HW(de)=%2d", popcount32(best_data.states_n[r][4] ^ best_data.states_f[r][4]));
        printf("\n");
    }
    /* r=64 state */
    {
        int hw_ds = 0;
        for(int i=0;i<8;i++) hw_ds += popcount32(best_data.states_n[64][i] ^ best_data.states_f[64][i]);
        printf("r=64 HW(dstate)=%3d HW(de)=%2d\n", hw_ds,
               popcount32(best_data.states_n[64][4] ^ best_data.states_f[64][4]));
    }

    /* Per-word HW(dhash) */
    printf("\nPer-word HW(dhash): ");
    for(int i=0;i<8;i++) printf("%d ", popcount32(best_data.hash_n[i] ^ best_data.hash_f[i]));
    printf("\n");

    /* Correlation analysis: 1000 random Wang pairs */
    printf("\n=== CORRELATION: carry_agreement(18..40) vs HW(dhash) ===\n");
    printf("N=1000 random Wang pairs\n");

    double sum_ca=0, sum_hw=0, sum_ca2=0, sum_hw2=0, sum_cahw=0;
    int N = 1000;
    FullData cd;

    /* Also collect for scatter plot output */
    printf("--- SCATTER DATA ---\n");
    for(int i=0;i<N;i++) {
        uint32_t Wn[16], Wf[16];
        for(int j=0;j<16;j++) Wn[j] = rng_next();
        wang_chain(Wn, Wf);
        full_sha256_pair(Wn, Wf, &cd);
        double ca = avg_carry_agreement(&cd, 18, 40);
        int hw_h = hw_dhash(&cd);

        if(i < 1000) printf("%.6f %d\n", ca, hw_h); /* all 1000 */

        sum_ca += ca; sum_hw += hw_h;
        sum_ca2 += ca*ca; sum_hw2 += (double)hw_h*hw_h;
        sum_cahw += ca*hw_h;
    }
    printf("--- END SCATTER ---\n");

    double mean_ca = sum_ca/N, mean_hw = sum_hw/N;
    double var_ca = sum_ca2/N - mean_ca*mean_ca;
    double var_hw = sum_hw2/N - mean_hw*mean_hw;
    double cov = sum_cahw/N - mean_ca*mean_hw;
    double pearson = cov / sqrt(var_ca * var_hw);

    printf("\nBaseline mean carry_agreement(18..40) = %.4f\n", mean_ca);
    printf("Baseline mean HW(dhash) = %.2f\n", mean_hw);
    printf("Pearson r(carry_agree, HW_dhash) = %.4f\n", pearson);
    printf("SA best = %.6f (%.1fx above baseline)\n", global_best, global_best/mean_ca);

    /* Histogram of carry_agreement */
    printf("\n=== CARRY_AGREEMENT HISTOGRAM (baseline) ===\n");
    int hist[20] = {0};
    /* recount from scratch */
    for(int i=0;i<10000;i++) {
        uint32_t Wn[16], Wf[16];
        for(int j=0;j<16;j++) Wn[j] = rng_next();
        wang_chain(Wn, Wf);
        full_sha256_pair(Wn, Wf, &cd);
        double ca = avg_carry_agreement(&cd, 18, 40);
        int bin = (int)(ca * 20);
        if(bin >= 20) bin = 19;
        hist[bin]++;
    }
    for(int i=0;i<20;i++)
        printf("[%.3f-%.3f): %d\n", i/20.0, (i+1)/20.0, hist[i]);

    return 0;
}
