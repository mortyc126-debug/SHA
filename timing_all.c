/*
 * TIMING ALL PATHS: micro-measurement of correct vs wrong.
 *
 * For each "tупик" that requires birthday: measure CPU cycles
 * for correct answer vs wrong answer. Even 1 cycle difference
 * over 2^128 candidates = meaningful.
 *
 * Uses rdtsc for cycle-accurate timing.
 *
 * gcc -O3 -march=native -o timing_all timing_all.c
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))
#define CH(e,f,g) (((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))
#define S0(x) (ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define S1(x) (ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define s0(x) (ROTR(x,7)^ROTR(x,18)^((x)>>3))
#define s1(x) (ROTR(x,17)^ROTR(x,19)^((x)>>10))
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

static inline uint64_t rdtsc(void) {
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ((uint64_t)lo) | (((uint64_t)hi) << 32);
}

void sha256_compute(const uint32_t msg[16], uint32_t hash[8]) {
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for(int r=0;r<64;r++){
        uint32_t t1=h+S1(e)+CH(e,f,g)+K[r]+W[r], t2=S0(a)+MAJ(a,b,c);
        h=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
    }
    hash[0]=a+IV[0];hash[1]=b+IV[1];hash[2]=c+IV[2];hash[3]=d+IV[3];
    hash[4]=e+IV[4];hash[5]=f+IV[5];hash[6]=g+IV[6];hash[7]=h+IV[7];
}

/* Measure cycles for N iterations, return average */
double measure_cycles(void (*func)(const uint32_t*, uint32_t*),
                      const uint32_t *input, uint32_t *output, int N) {
    /* Warmup */
    for(int i=0;i<100;i++) func(input, output);

    uint64_t total = 0;
    for(int i=0;i<N;i++) {
        uint64_t start = rdtsc();
        func(input, output);
        uint64_t end = rdtsc();
        total += (end - start);
    }
    return (double)total / N;
}

int main() {
    printf("TIMING ALL PATHS: cycle-accurate measurement\n");
    printf("=============================================\n\n");

    srand(42);
    uint32_t msg_correct[16], msg_wrong[16], hash_target[8], hash_out[8];
    for(int i=0;i<16;i++) msg_correct[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
    sha256_compute(msg_correct, hash_target);

    int N = 100000;

    printf("N = %d iterations per measurement\n\n", N);

    /* ═══ TEST 1: Full SHA-256 — correct vs wrong msg ═══ */
    printf("TEST 1: Full SHA-256 forward\n");
    printf("────────────────────────────\n");
    {
        double cyc_correct = measure_cycles(sha256_compute, msg_correct, hash_out, N);

        /* Multiple wrong messages */
        double cyc_wrong_total = 0;
        int n_wrong = 10;
        for(int w=0; w<n_wrong; w++) {
            for(int i=0;i<16;i++) msg_wrong[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
            cyc_wrong_total += measure_cycles(sha256_compute, msg_wrong, hash_out, N);
        }
        double cyc_wrong = cyc_wrong_total / n_wrong;
        double diff = cyc_correct - cyc_wrong;
        double pct = diff / cyc_wrong * 100;

        printf("  Correct msg:  %.1f cycles\n", cyc_correct);
        printf("  Wrong msgs:   %.1f cycles (avg of %d)\n", cyc_wrong, n_wrong);
        printf("  Difference:   %+.1f cycles (%+.3f%%)\n", diff, pct);
        printf("  → %s\n", (diff > 5 || diff < -5) ? "★ TIMING SIGNAL!" : "No signal (within noise)");
    }

    /* ═══ TEST 2: Near-miss vs far-miss ═══ */
    printf("\nTEST 2: Near-miss (1 bit diff) vs far-miss (random)\n");
    printf("────────────────────────────────────────────────────\n");
    {
        /* Near miss: correct msg with 1 bit flipped */
        uint32_t msg_near[16];
        memcpy(msg_near, msg_correct, 64);
        msg_near[0] ^= 1;

        double cyc_near = measure_cycles(sha256_compute, msg_near, hash_out, N);

        /* Far miss */
        for(int i=0;i<16;i++) msg_wrong[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
        double cyc_far = measure_cycles(sha256_compute, msg_wrong, hash_out, N);

        double cyc_correct = measure_cycles(sha256_compute, msg_correct, hash_out, N);

        printf("  Correct:    %.1f cycles\n", cyc_correct);
        printf("  Near-miss:  %.1f cycles (1 bit diff)\n", cyc_near);
        printf("  Far-miss:   %.1f cycles (random)\n", cyc_far);
        printf("  Near-correct: %+.1f cycles\n", cyc_near - cyc_correct);
        printf("  Far-correct:  %+.1f cycles\n", cyc_far - cyc_correct);
    }

    /* ═══ TEST 3: Schedule computation ═══ */
    printf("\nTEST 3: Schedule only (no rounds)\n");
    printf("──────────────────────────────────\n");
    {
        uint32_t W1[64], W2[64];
        volatile uint32_t sink;

        uint64_t total_correct = 0, total_wrong = 0;
        for(int iter=0;iter<N;iter++){
            for(int i=0;i<16;i++) W1[i]=msg_correct[i];
            uint64_t t0=rdtsc();
            for(int i=16;i<64;i++) W1[i]=s1(W1[i-2])+W1[i-7]+s0(W1[i-15])+W1[i-16];
            sink=W1[63];
            total_correct += rdtsc()-t0;

            for(int i=0;i<16;i++) W2[i]=msg_wrong[i];
            t0=rdtsc();
            for(int i=16;i<64;i++) W2[i]=s1(W2[i-2])+W2[i-7]+s0(W2[i-15])+W2[i-16];
            sink=W2[63];
            total_wrong += rdtsc()-t0;
        }
        printf("  Correct msg schedule:  %.1f cycles\n", (double)total_correct/N);
        printf("  Wrong msg schedule:    %.1f cycles\n", (double)total_wrong/N);
        printf("  Diff: %+.1f\n", (double)(total_correct-total_wrong)/N);
    }

    /* ═══ TEST 4: Backward chain ═══ */
    printf("\nTEST 4: Backward chain (створочне)\n");
    printf("───────────────────────────────────\n");
    {
        volatile uint32_t sink;
        uint64_t total_correct = 0, total_wrong = 0;

        for(int iter=0;iter<N;iter++){
            /* Correct: backward from correct hash */
            uint32_t a[65];
            a[64]=hash_target[0]-IV[0]; a[63]=hash_target[1]-IV[1];
            a[62]=hash_target[2]-IV[2]; a[61]=hash_target[3]-IV[3];
            uint32_t e64=hash_target[4]-IV[4], e63=hash_target[5]-IV[5];
            uint32_t e62=hash_target[6]-IV[6], e61=hash_target[7]-IV[7];

            uint64_t t0=rdtsc();
            for(int r=63;r>=60;r--){
                uint32_t T2=S0(a[r])+MAJ(a[r],a[r-1],a[r-2]);
                uint32_t T1=a[r+1]-T2;
                uint32_t e_rp1;
                if(r==63)e_rp1=e64; else if(r==62)e_rp1=e63;
                else if(r==61)e_rp1=e62; else e_rp1=e61;
                a[r-3]=e_rp1-T1;
            }
            sink=a[57];
            total_correct += rdtsc()-t0;

            /* Wrong: backward from random hash */
            uint32_t fake_hash[8];
            for(int i=0;i<8;i++) fake_hash[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
            a[64]=fake_hash[0]-IV[0]; a[63]=fake_hash[1]-IV[1];
            a[62]=fake_hash[2]-IV[2]; a[61]=fake_hash[3]-IV[3];
            e64=fake_hash[4]-IV[4]; e63=fake_hash[5]-IV[5];
            e62=fake_hash[6]-IV[6]; e61=fake_hash[7]-IV[7];

            t0=rdtsc();
            for(int r=63;r>=60;r--){
                uint32_t T2=S0(a[r])+MAJ(a[r],a[r-1],a[r-2]);
                uint32_t T1=a[r+1]-T2;
                uint32_t e_rp1;
                if(r==63)e_rp1=e64; else if(r==62)e_rp1=e63;
                else if(r==61)e_rp1=e62; else e_rp1=e61;
                a[r-3]=e_rp1-T1;
            }
            sink=a[57];
            total_wrong += rdtsc()-t0;
        }
        printf("  Correct backward: %.1f cycles\n", (double)total_correct/N);
        printf("  Wrong backward:   %.1f cycles\n", (double)total_wrong/N);
        printf("  Diff: %+.1f\n", (double)(total_correct-total_wrong)/N);
    }

    /* ═══ TEST 5: Hash comparison (the ONLY difference) ═══ */
    printf("\nTEST 5: Hash comparison (match vs mismatch)\n");
    printf("────────────────────────────────────────────\n");
    {
        uint32_t h1[8], h2[8], h3[8];
        memcpy(h1, hash_target, 32);
        memcpy(h2, hash_target, 32); /* match */
        for(int i=0;i<8;i++) h3[i]=(uint32_t)rand()|((uint32_t)rand()<<16); /* mismatch */

        volatile int sink;
        uint64_t total_match = 0, total_mismatch = 0;

        for(int iter=0;iter<N;iter++){
            uint64_t t0=rdtsc();
            int match = (h1[0]==h2[0] && h1[1]==h2[1] && h1[2]==h2[2] && h1[3]==h2[3] &&
                         h1[4]==h2[4] && h1[5]==h2[5] && h1[6]==h2[6] && h1[7]==h2[7]);
            sink=match;
            total_match += rdtsc()-t0;

            t0=rdtsc();
            int nomatch = (h1[0]==h3[0] && h1[1]==h3[1] && h1[2]==h3[2] && h1[3]==h3[3] &&
                           h1[4]==h3[4] && h1[5]==h3[5] && h1[6]==h3[6] && h1[7]==h3[7]);
            sink=nomatch;
            total_mismatch += rdtsc()-t0;
        }
        printf("  Hash match:    %.1f cycles\n", (double)total_match/N);
        printf("  Hash mismatch: %.1f cycles\n", (double)total_mismatch/N);
        printf("  Diff: %+.1f cycles\n", (double)(total_match-total_mismatch)/N);
        if((double)(total_match-total_mismatch)/N > 2.0)
            printf("  ★ MATCH IS SLOWER (checks all 8 words)!\n");
        else if((double)(total_match-total_mismatch)/N < -2.0)
            printf("  ★ MISMATCH IS SLOWER!\n");
        else
            printf("  → No significant difference.\n");
    }

    /* ═══ TEST 6: Carry chain length as signal ═══ */
    printf("\nTEST 6: Carry chain length (correct vs wrong)\n");
    printf("──────────────────────────────────────────────\n");
    {
        /* Correct msg: measure total carry propagation length */
        uint32_t W[64], st[65][8];
        for(int i=0;i<16;i++) W[i]=msg_correct[i];
        for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
        for(int i=0;i<8;i++) st[0][i]=IV[i];

        long long carry_correct = 0, carry_wrong = 0;

        /* For correct msg: count carry bits across all rounds */
        for(int r=0;r<64;r++){
            uint32_t a=st[r][0],b=st[r][1],c=st[r][2],d=st[r][3];
            uint32_t e=st[r][4],f=st[r][5],g=st[r][6],h=st[r][7];
            uint32_t sig1e=S1(e), che=CH(e,f,g);
            /* T1 = h + sig1e + che + K[r] + W[r] (4 additions) */
            uint32_t s1_val = h + sig1e;
            carry_correct += __builtin_popcount(s1_val ^ (h ^ sig1e));
            uint32_t s2_val = s1_val + che;
            carry_correct += __builtin_popcount(s2_val ^ (s1_val ^ che));
            uint32_t s3_val = s2_val + K[r];
            carry_correct += __builtin_popcount(s3_val ^ (s2_val ^ K[r]));
            uint32_t t1 = s3_val + W[r];
            carry_correct += __builtin_popcount(t1 ^ (s3_val ^ W[r]));

            uint32_t sig0a=S0(a), maja=MAJ(a,b,c);
            uint32_t t2 = sig0a + maja;
            carry_correct += __builtin_popcount(t2 ^ (sig0a ^ maja));

            uint32_t a_new = t1 + t2;
            carry_correct += __builtin_popcount(a_new ^ (t1 ^ t2));
            uint32_t e_new = d + t1;
            carry_correct += __builtin_popcount(e_new ^ (d ^ t1));

            st[r+1][0]=a_new;st[r+1][1]=a;st[r+1][2]=b;st[r+1][3]=c;
            st[r+1][4]=e_new;st[r+1][5]=e;st[r+1][6]=f;st[r+1][7]=g;
        }

        /* For wrong msg */
        for(int i=0;i<16;i++) W[i]=msg_wrong[i];
        for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
        for(int i=0;i<8;i++) st[0][i]=IV[i];
        for(int r=0;r<64;r++){
            uint32_t a=st[r][0],b=st[r][1],c=st[r][2],d=st[r][3];
            uint32_t e=st[r][4],f=st[r][5],g=st[r][6],h=st[r][7];
            uint32_t sig1e=S1(e), che=CH(e,f,g);
            uint32_t s1_val=h+sig1e; carry_wrong+=__builtin_popcount(s1_val^(h^sig1e));
            uint32_t s2_val=s1_val+che; carry_wrong+=__builtin_popcount(s2_val^(s1_val^che));
            uint32_t s3_val=s2_val+K[r]; carry_wrong+=__builtin_popcount(s3_val^(s2_val^K[r]));
            uint32_t t1=s3_val+W[r]; carry_wrong+=__builtin_popcount(t1^(s3_val^W[r]));
            uint32_t sig0a=S0(a), maja=MAJ(a,b,c);
            uint32_t t2=sig0a+maja; carry_wrong+=__builtin_popcount(t2^(sig0a^maja));
            uint32_t a_new=t1+t2; carry_wrong+=__builtin_popcount(a_new^(t1^t2));
            uint32_t e_new=d+t1; carry_wrong+=__builtin_popcount(e_new^(d^t1));
            st[r+1][0]=a_new;st[r+1][1]=a;st[r+1][2]=b;st[r+1][3]=c;
            st[r+1][4]=e_new;st[r+1][5]=e;st[r+1][6]=f;st[r+1][7]=g;
        }

        printf("  Correct msg total carry bits: %lld\n", carry_correct);
        printf("  Wrong msg total carry bits:   %lld\n", carry_wrong);
        printf("  Diff: %+lld (%.2f%%)\n", carry_correct-carry_wrong,
               100.0*(carry_correct-carry_wrong)/carry_correct);

        /* Stats over many messages */
        long long carry_stats[100];
        for(int m=0;m<100;m++){
            uint32_t msg_t[16];
            for(int i=0;i<16;i++) msg_t[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
            for(int i=0;i<16;i++) W[i]=msg_t[i];
            for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
            for(int i=0;i<8;i++) st[0][i]=IV[i];
            long long c=0;
            for(int r=0;r<64;r++){
                uint32_t a=st[r][0],b=st[r][1],cc=st[r][2],d=st[r][3];
                uint32_t e=st[r][4],f=st[r][5],g=st[r][6],h=st[r][7];
                uint32_t t1=h+S1(e)+CH(e,f,g)+K[r]+W[r], t2=S0(a)+MAJ(a,b,cc);
                c+=__builtin_popcount((h+S1(e))^(h^S1(e)));
                c+=__builtin_popcount(t2^(S0(a)^MAJ(a,b,cc)));
                st[r+1][0]=t1+t2;st[r+1][1]=a;st[r+1][2]=b;st[r+1][3]=cc;
                st[r+1][4]=d+t1;st[r+1][5]=e;st[r+1][6]=f;st[r+1][7]=g;
            }
            carry_stats[m]=c;
        }
        double mean=0,var=0;
        for(int i=0;i<100;i++) mean+=carry_stats[i]; mean/=100;
        for(int i=0;i<100;i++) {double d=carry_stats[i]-mean;var+=d*d;} var/=100;
        printf("  Carry stats (100 msgs): mean=%.0f, std=%.0f\n", mean, sqrt(var));
        printf("  Correct msg carry: %lld (Z=%.2f)\n", carry_correct,
               (carry_correct-mean)/sqrt(var));
    }

    printf("\n════════════════════════════════════\n");
    printf("SUMMARY\n");
    printf("════════════════════════════════════\n");

    return 0;
}
