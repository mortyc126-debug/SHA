/*
 * AMPLIFIED SIGNAL SPACE: custom projection that amplifies ch×carry.
 *
 * Strategy: instead of looking at hash (random) or full Φ (weak),
 * construct a CUSTOM projection P(state) designed to:
 *   1. AMPLIFY the ch×carry interaction signal
 *   2. SUPPRESS the random noise
 *   3. Be COMPUTABLE from forward state
 *
 * Approach: at each round, extract ONLY the bits where
 * Ch self-cancellation AND carry interact — and accumulate
 * across rounds. This is our "tuned antenna".
 *
 * gcc -O3 -march=native -o amp_space amplified_space.c -lm
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

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

/* Custom signal projections — different "antennas" */

typedef struct {
    int score;  /* higher = more "correct-like" */
} Signal;

/* Antenna 1: ch×carry at erased positions, accumulated */
Signal antenna_ch_carry(const uint32_t msg[16]) {
    Signal sig = {0};
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for(int r=0;r<64;r++){
        uint32_t erase = f ^ g;  /* Ch erased positions */
        uint32_t s1e=S1(e), che=CH(e,f,g);
        uint32_t v1 = h + s1e;
        uint32_t carry1 = v1 ^ (h ^ s1e);
        /* Carry AT erased positions */
        sig.score += __builtin_popcount(carry1 & erase);
        /* Also: carry OF Ch addition at erased positions */
        uint32_t v2 = v1 + che;
        uint32_t carry2 = v2 ^ (v1 ^ che);
        sig.score += __builtin_popcount(carry2 & erase);

        uint32_t t1=v2+K[r]+W[r]; /* simplified: skip intermediate carries */
        uint32_t mja=MAJ(a,b,c), s0a=S0(a);
        uint32_t t2=s0a+mja;
        uint32_t anew=t1+t2, enew=d+t1;
        h=g;g=f;f=e;e=enew;d=c;c=b;b=a;a=anew;
    }
    return sig;
}

/* Antenna 2: carry chain LENGTH at erased positions */
Signal antenna_chain_length(const uint32_t msg[16]) {
    Signal sig = {0};
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for(int r=0;r<64;r++){
        uint32_t erase = f ^ g;
        uint32_t s1e=S1(e), che=CH(e,f,g);
        uint32_t v1 = h + s1e;
        uint32_t carry = v1 ^ (h ^ s1e);
        /* Longest consecutive carry=1 AT erased positions */
        uint32_t masked = carry & erase;
        int chain=0, maxchain=0;
        for(int k=0;k<32;k++){
            if((masked>>k)&1){chain++;if(chain>maxchain)maxchain=chain;}
            else chain=0;
        }
        sig.score += maxchain;

        uint32_t t1=h+s1e+che+K[r]+W[r];
        uint32_t t2=S0(a)+MAJ(a,b,c);
        uint32_t anew=t1+t2, enew=d+t1;
        h=g;g=f;f=e;e=enew;d=c;c=b;b=a;a=anew;
    }
    return sig;
}

/* Antenna 3: T1 XOR T2 at erased positions (NLF interaction) */
Signal antenna_t1_xor_t2(const uint32_t msg[16]) {
    Signal sig = {0};
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for(int r=0;r<64;r++){
        uint32_t erase = f ^ g;
        uint32_t t1=h+S1(e)+CH(e,f,g)+K[r]+W[r];
        uint32_t t2=S0(a)+MAJ(a,b,c);
        sig.score += __builtin_popcount((t1 ^ t2) & erase);
        uint32_t anew=t1+t2, enew=d+t1;
        h=g;g=f;f=e;e=enew;d=c;c=b;b=a;a=anew;
    }
    return sig;
}

/* Antenna 4: Deadpool survivor count (how many erased bits survive 3 rounds) */
Signal antenna_deadpool(const uint32_t msg[16]) {
    Signal sig = {0};
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    uint32_t prev_erase[4] = {0,0,0,0};
    for(int r=0;r<64;r++){
        uint32_t erase = f ^ g;
        /* At round r: bits erased at r-3 are now at h position.
         * h = g[r-1] = f[r-2] = e[r-3].
         * Count: h bits that MATCH the erasure pattern from r-3. */
        if(r >= 3) {
            sig.score += __builtin_popcount(h & prev_erase[r%4]);
        }
        prev_erase[r%4] = erase;

        uint32_t t1=h+S1(e)+CH(e,f,g)+K[r]+W[r];
        uint32_t t2=S0(a)+MAJ(a,b,c);
        uint32_t anew=t1+t2, enew=d+t1;
        h=g;g=f;f=e;e=enew;d=c;c=b;b=a;a=anew;
    }
    return sig;
}

/* Antenna 5: COMBINED (weighted sum of all signals) */
Signal antenna_combined(const uint32_t msg[16]) {
    Signal s1 = antenna_ch_carry(msg);
    Signal s2 = antenna_chain_length(msg);
    Signal s3 = antenna_t1_xor_t2(msg);
    Signal s4 = antenna_deadpool(msg);
    return (Signal){s1.score * 3 + s2.score * 5 + s3.score * 2 + s4.score * 4};
}

typedef Signal (*antenna_fn)(const uint32_t msg[16]);

int main() {
    printf("AMPLIFIED SIGNAL SPACE: custom antennas\n");
    printf("========================================\n\n");

    srand(42);
    uint32_t msg_correct[16];
    for(int i=0;i<16;i++) msg_correct[i]=(uint32_t)rand()|((uint32_t)rand()<<16);

    antenna_fn antennas[] = {antenna_ch_carry, antenna_chain_length,
                             antenna_t1_xor_t2, antenna_deadpool, antenna_combined};
    const char *names[] = {"ch×carry","chain_len","t1⊕t2@erase","deadpool","COMBINED"};
    int n_ant = 5;

    int N = 20000;

    printf("%-14s | %8s | %8s | %6s | %6s | Signal?\n",
           "Antenna","Correct","Mean(rnd)","Std","Z");
    printf("───────────────+──────────+──────────+────────+────────+─────────\n");

    for(int ai=0;ai<n_ant;ai++){
        int correct_val = antennas[ai](msg_correct).score;
        double sum=0,sum2=0;
        for(int t=0;t<N;t++){
            uint32_t msg_w[16];
            for(int i=0;i<16;i++) msg_w[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
            int v = antennas[ai](msg_w).score;
            sum+=v; sum2+=(double)v*v;
        }
        double mean=sum/N, std=sqrt(sum2/N-mean*mean);
        double z=std>0?(correct_val-mean)/std:0;
        printf("%-14s | %8d | %8.1f | %6.1f | %+6.2f | %s\n",
               names[ai],correct_val,mean,std,z,
               fabs(z)>3?"★★★ YES!":fabs(z)>2?"★★ maybe":"no");
    }

    /* Test with PARTIAL correct words */
    printf("\nPartial correct words (K=0..16), COMBINED antenna:\n");
    printf("──────────────────────────────────────────────────\n");
    int correct_combined = antenna_combined(msg_correct).score;
    printf("  Full correct: %d\n\n", correct_combined);

    for(int K=0;K<=16;K+=2){
        double sum=0,sum2=0;
        int n=5000;
        for(int t=0;t<n;t++){
            uint32_t msg_t[16];
            for(int i=0;i<K;i++) msg_t[i]=msg_correct[i];
            for(int i=K;i<16;i++) msg_t[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
            int v=antenna_combined(msg_t).score;
            sum+=v;sum2+=(double)v*v;
        }
        double mean=sum/n, std=sqrt(sum2/n-mean*mean);
        double z_vs_random = std>0?(mean - antenna_combined(msg_correct).score*0+
                                    (correct_combined-mean))/std:0;
        /* Actually: compare mean(K correct) with mean(0 correct) */
        printf("  K=%2d: mean=%.1f, std=%.1f",K,mean,std);
        if(K>0){
            /* Does having K correct words change the antenna reading? */
            printf(", shift=%+.1f", mean-(sum/n)); /* self-reference bug, fix: */
        }
        printf("\n");
    }

    /* Better: measure shift from K=0 baseline */
    printf("\n  Signal vs K=0 baseline:\n");
    double baseline_mean=0, baseline_std=0;
    {
        double s=0,s2=0; int n=10000;
        for(int t=0;t<n;t++){
            uint32_t msg_t[16];
            for(int i=0;i<16;i++) msg_t[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
            int v=antenna_combined(msg_t).score;
            s+=v;s2+=(double)v*v;
        }
        baseline_mean=s/n; baseline_std=sqrt(s2/n-baseline_mean*baseline_mean);
    }
    printf("  K=0 baseline: mean=%.1f, std=%.1f\n\n", baseline_mean, baseline_std);

    for(int K=0;K<=16;K++){
        double s=0; int n=5000;
        for(int t=0;t<n;t++){
            uint32_t msg_t[16];
            for(int i=0;i<K;i++) msg_t[i]=msg_correct[i];
            for(int i=K;i<16;i++) msg_t[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
            s+=antenna_combined(msg_t).score;
        }
        double mean=s/n;
        double z=(mean-baseline_mean)/baseline_std;
        printf("  K=%2d: mean=%7.1f, Z vs baseline=%+6.2f %s\n",
               K,mean,z,fabs(z)>3?"★★★":fabs(z)>2?"★★":fabs(z)>1?"★":"");
    }

    return 0;
}
