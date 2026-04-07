/*
 * HASH-DERIVED ANTENNA: tuned to target hash.
 *
 * We know from hash: state[64], a[57..64], erasure[63].
 * Use this to BUILD antenna specific to this target.
 *
 * Key: erasure pattern at round 63 = e[62]⊕e[61] = KNOWN.
 * Carry pattern at round 63 = function of state[63] = partially known.
 *
 * For candidate msg: compute forward → extract features →
 * compare with hash-derived EXPECTED features.
 *
 * gcc -O3 -march=native -o hash_ant hash_antenna.c -lm
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

/* Forward SHA-256 returning state[64] */
void sha256_with_state64(const uint32_t msg[16], uint32_t state64[8], uint32_t hash[8]) {
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for(int r=0;r<64;r++){
        uint32_t t1=h+S1(e)+CH(e,f,g)+K[r]+W[r], t2=S0(a)+MAJ(a,b,c);
        h=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
    }
    state64[0]=a;state64[1]=b;state64[2]=c;state64[3]=d;
    state64[4]=e;state64[5]=f;state64[6]=g;state64[7]=h;
    for(int i=0;i<8;i++) hash[i]=state64[i]+IV[i];
}

/* Hash-derived antenna features */
typedef struct {
    /* Known from hash backward */
    uint32_t a[65]; /* a[57..64] */
    uint32_t erase_63; /* erasure pattern at round 63: e[62]⊕e[61] */
    uint32_t erase_62_mask; /* partial: known bits of e[62]⊕e[60] */
    uint32_t T2_63, T1_63; /* known from backward */
    uint32_t T2_62, T2_61, T2_60;
    /* Expected Ch and Maj outputs at round 63 */
    uint32_t ch_63, maj_63;
    /* Known carry pattern at round 63 */
    uint32_t carry_63_expected; /* carry of a[64] = T1+T2 */
} HashAntenna;

void build_hash_antenna(const uint32_t hash[8], HashAntenna *ant) {
    /* state[64] = hash - IV */
    uint32_t s64[8];
    for(int i=0;i<8;i++) s64[i] = hash[i]-IV[i];

    /* a[61..64] directly */
    ant->a[64]=s64[0]; ant->a[63]=s64[1]; ant->a[62]=s64[2]; ant->a[61]=s64[3];
    /* e[61..64] */
    uint32_t e64=s64[4],e63=s64[5],e62=s64[6],e61=s64[7];

    /* створочне: a[57..60] */
    for(int r=63;r>=60;r--){
        uint32_t T2=S0(ant->a[r])+MAJ(ant->a[r],ant->a[r-1],ant->a[r-2]);
        uint32_t T1=ant->a[r+1]-T2;
        uint32_t e_rp1;
        if(r==63)e_rp1=e64; else if(r==62)e_rp1=e63;
        else if(r==61)e_rp1=e62; else e_rp1=e61;
        ant->a[r-3]=e_rp1-T1;
    }

    /* Erasure at round 63: e[62]⊕e[61] */
    ant->erase_63 = e62 ^ e61;

    /* T2[63], T1[63] */
    ant->T2_63 = S0(ant->a[63])+MAJ(ant->a[63],ant->a[62],ant->a[61]);
    ant->T1_63 = ant->a[64]-ant->T2_63;

    /* Ch and Maj at round 63: uses e[63],e[62],e[61] for Ch; a[63],a[62],a[61] for Maj */
    ant->ch_63 = CH(e63,e62,e61);
    ant->maj_63 = MAJ(ant->a[63],ant->a[62],ant->a[61]);

    /* Carry at round 63: carry of T1+T2 → a[64] */
    ant->carry_63_expected = (ant->T1_63 + ant->T2_63) ^ (ant->T1_63 ^ ant->T2_63);

    /* T2 for rounds 60..62 */
    ant->T2_62 = S0(ant->a[62])+MAJ(ant->a[62],ant->a[61],ant->a[60]);
    ant->T2_61 = S0(ant->a[61])+MAJ(ant->a[61],ant->a[60],ant->a[59]);
    ant->T2_60 = S0(ant->a[60])+MAJ(ant->a[60],ant->a[59],ant->a[58]);
}

/* Score candidate msg against hash-derived antenna */
int score_candidate(const uint32_t msg[16], const HashAntenna *ant) {
    uint32_t state64[8], hash[8];
    sha256_with_state64(msg, state64, hash);

    int score = 0;

    /* 1. Match a[57..64] from candidate vs target */
    /* a[64] = state64[0], a[63] = state64[1], etc. */
    for(int i=0;i<8;i++){
        /* Count matching BITS between candidate state[64] and target */
        score += 32 - __builtin_popcount(state64[i] ^ ((uint32_t*)&ant->a[64-i])[0]);
        /* Wait: ant->a is indexed differently. Use proper values. */
    }
    /* Fix: compare a[57..64] properly */
    score = 0;
    score += 32 - __builtin_popcount(state64[0] ^ ant->a[64]); /* a[64] */
    score += 32 - __builtin_popcount(state64[1] ^ ant->a[63]); /* a[63] */
    score += 32 - __builtin_popcount(state64[2] ^ ant->a[62]);
    score += 32 - __builtin_popcount(state64[3] ^ ant->a[61]);

    /* 2. Erasure pattern match at round 63 */
    /* From candidate: e[62]⊕e[61] at round 63 = state64[6]⊕state64[7] */
    uint32_t cand_erase_63 = state64[6] ^ state64[7];
    score += 32 - __builtin_popcount(cand_erase_63 ^ ant->erase_63);

    /* 3. Carry pattern at final addition */
    uint32_t cand_T2 = S0(state64[1])+MAJ(state64[1],state64[2],state64[3]);
    uint32_t cand_T1 = state64[0]-cand_T2;
    uint32_t cand_carry = (cand_T1+cand_T2) ^ (cand_T1^cand_T2);
    score += 32 - __builtin_popcount(cand_carry ^ ant->carry_63_expected);

    /* 4. Ch output match */
    uint32_t cand_ch = CH(state64[5],state64[6],state64[7]);
    score += 32 - __builtin_popcount(cand_ch ^ ant->ch_63);

    return score;
}

int main() {
    printf("HASH-DERIVED ANTENNA\n");
    printf("====================\n\n");

    srand(42);
    uint32_t msg_correct[16];
    for(int i=0;i<16;i++) msg_correct[i]=(uint32_t)rand()|((uint32_t)rand()<<16);

    uint32_t state64_correct[8], target_hash[8];
    sha256_with_state64(msg_correct, state64_correct, target_hash);

    /* Build antenna from hash */
    HashAntenna ant;
    build_hash_antenna(target_hash, &ant);

    printf("Target hash: ");
    for(int i=0;i<8;i++) printf("%08x ", target_hash[i]);
    printf("\n");
    printf("Erasure[63] = %08x (HW=%d)\n", ant.erase_63,
           __builtin_popcount(ant.erase_63));
    printf("Known a[57..64]: ");
    for(int r=57;r<=64;r++) printf("%08x ", ant.a[r]);
    printf("\n\n");

    /* Score correct msg */
    int score_correct = score_candidate(msg_correct, &ant);
    printf("Correct msg score: %d\n", score_correct);

    /* Score random msgs */
    int N = 50000;
    double sum=0,sum2=0;
    int max_wrong = 0;
    for(int t=0;t<N;t++){
        uint32_t msg_w[16];
        for(int i=0;i<16;i++) msg_w[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
        int s = score_candidate(msg_w, &ant);
        sum+=s; sum2+=(double)s*s;
        if(s>max_wrong) max_wrong=s;
    }
    double mean=sum/N, std=sqrt(sum2/N-mean*mean);
    double z = std>0?(score_correct-mean)/std:0;

    printf("Random msg:  mean=%.1f, std=%.1f, max=%d\n", mean, std, max_wrong);
    printf("Z-score: %+.2f %s\n\n", z,
           fabs(z)>5?"★★★★★ EXTREME":fabs(z)>3?"★★★ STRONG":fabs(z)>2?"★★":"no");

    /* Test UNIVERSALITY: does this work for other correct messages? */
    printf("UNIVERSALITY: antenna from hash → score other correct msgs\n");
    printf("───────────────────────────────────────────────────────────\n");

    /* For 10 different messages: build antenna from THEIR hash,
     * then score correct vs random. */
    printf("  msg# | correct_score | random_mean | Z\n");
    printf("  ─────+───────────────+─────────────+──────\n");

    for(int m=0;m<10;m++){
        uint32_t msg_m[16], st64[8], h_m[8];
        srand(m*7777+13);
        for(int i=0;i<16;i++) msg_m[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
        sha256_with_state64(msg_m, st64, h_m);

        HashAntenna ant_m;
        build_hash_antenna(h_m, &ant_m);

        int sc = score_candidate(msg_m, &ant_m);
        double s=0,s2=0;
        for(int t=0;t<5000;t++){
            uint32_t msg_w[16];
            for(int i=0;i<16;i++) msg_w[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
            int v = score_candidate(msg_w, &ant_m);
            s+=v;s2+=(double)v*v;
        }
        double mn=s/5000, sd=sqrt(s2/5000-mn*mn);
        double zz=sd>0?(sc-mn)/sd:0;
        printf("  %4d | %13d | %11.1f | %+6.2f %s\n",
               m, sc, mn, zz,
               fabs(zz)>3?"★★★":fabs(zz)>2?"★★":"");
    }

    /* Partial correct words with hash antenna */
    printf("\nPartial correct words (hash-derived antenna):\n");
    printf("──────────────────────────────────────────────\n");
    for(int K=0;K<=16;K+=2){
        double s=0; int n=5000;
        for(int t=0;t<n;t++){
            uint32_t msg_t[16];
            for(int i=0;i<K;i++) msg_t[i]=msg_correct[i];
            for(int i=K;i<16;i++) msg_t[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
            s+=score_candidate(msg_t, &ant);
        }
        double mn=s/n;
        double zz=std>0?(mn-mean)/std:0;
        printf("  K=%2d: mean_score=%.1f, Z=%+.2f %s\n",K,mn,zz,
               fabs(zz)>2?"★★":fabs(zz)>1?"★":"");
    }

    return 0;
}
