/*
 * CARRY × NLF ALGEBRA: trace the exact interaction.
 *
 * We know: without carry×NLF = 128 deterministic hash bits.
 *          with carry×NLF = 0 deterministic bits.
 *
 * Questions:
 * 1. WHERE do the 128 bits die? Which round, which register?
 * 2. Is the death GRADUAL (1 bit per round) or SUDDEN?
 * 3. Can we ISOLATE the carry×NLF term at each round?
 * 4. What's the ALGEBRAIC STRUCTURE of carry(Ch) and carry(Maj)?
 *
 * gcc -O3 -march=native -o cnlf carry_nlf_algebra.c
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))
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

/* Full SHA-256 with per-round state */
void sha256_full(const uint32_t msg[16], uint32_t st[65][8]) {
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=(ROTR(W[i-2],17)^ROTR(W[i-2],19)^(W[i-2]>>10))+W[i-7]+(ROTR(W[i-15],7)^ROTR(W[i-15],18)^(W[i-15]>>3))+W[i-16];
    for(int i=0;i<8;i++) st[0][i]=IV[i];
    for(int r=0;r<64;r++){
        uint32_t a=st[r][0],b=st[r][1],c=st[r][2],d=st[r][3];
        uint32_t e=st[r][4],f=st[r][5],g=st[r][6],h=st[r][7];
        uint32_t t1=h+(ROTR(e,6)^ROTR(e,11)^ROTR(e,25))+((e&f)^(~e&g))+K[r]+W[r];
        uint32_t t2=(ROTR(a,2)^ROTR(a,13)^ROTR(a,22))+((a&b)^(a&c)^(b&c));
        st[r+1][7]=g;st[r+1][6]=f;st[r+1][5]=e;st[r+1][4]=d+t1;
        st[r+1][3]=c;st[r+1][2]=b;st[r+1][1]=a;st[r+1][0]=t1+t2;
    }
}

/* SHA-256 with NLF replaced by XOR3 (carry preserved) */
void sha256_no_nlf(const uint32_t msg[16], uint32_t st[65][8]) {
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=(ROTR(W[i-2],17)^ROTR(W[i-2],19)^(W[i-2]>>10))+W[i-7]+(ROTR(W[i-15],7)^ROTR(W[i-15],18)^(W[i-15]>>3))+W[i-16];
    for(int i=0;i<8;i++) st[0][i]=IV[i];
    for(int r=0;r<64;r++){
        uint32_t a=st[r][0],b=st[r][1],c=st[r][2],d=st[r][3];
        uint32_t e=st[r][4],f=st[r][5],g=st[r][6],h=st[r][7];
        uint32_t t1=h+(ROTR(e,6)^ROTR(e,11)^ROTR(e,25))+(e^f^g)+K[r]+W[r]; /* Ch→XOR3 */
        uint32_t t2=(ROTR(a,2)^ROTR(a,13)^ROTR(a,22))+(a^b^c);             /* Maj→XOR3 */
        st[r+1][7]=g;st[r+1][6]=f;st[r+1][5]=e;st[r+1][4]=d+t1;
        st[r+1][3]=c;st[r+1][2]=b;st[r+1][1]=a;st[r+1][0]=t1+t2;
    }
}

int main() {
    printf("CARRY × NLF ALGEBRA: where exactly do 128 bits die?\n");
    printf("====================================================\n\n");

    uint32_t msg[16], msg2[16];
    srand(42);
    for(int i=0;i<16;i++) msg[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
    memcpy(msg2, msg, 64);
    msg2[0] ^= 1;

    /* Test 1: per-round deterministic bit count.
     * For MANY messages, flip W[0][b0], track which state bits ALWAYS flip.
     * Do this for FULL SHA-256 and NO-NLF variant. */
    int N = 2000;
    int det_full[65] = {0};  /* deterministic bits per round, full SHA */
    int det_nonlf[65] = {0}; /* deterministic bits per round, no-NLF */

    /* Count per (round, register, bit): how many times this bit flips */
    /* Too much memory for all — just count per round total */
    int flip_count_full[65][256];
    int flip_count_nonlf[65][256];
    memset(flip_count_full, 0, sizeof(flip_count_full));
    memset(flip_count_nonlf, 0, sizeof(flip_count_nonlf));

    for(int trial=0; trial<N; trial++) {
        uint32_t m[16], m2[16], s1[65][8], s2[65][8];
        for(int i=0;i<16;i++) m[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
        memcpy(m2, m, 64);
        m2[0] ^= 1;

        /* Full SHA */
        sha256_full(m, s1);
        sha256_full(m2, s2);
        for(int r=0;r<=64;r++)
            for(int w=0;w<8;w++)
                for(int b=0;b<32;b++)
                    if(((s1[r][w]^s2[r][w])>>b)&1)
                        flip_count_full[r][w*32+b]++;

        /* No-NLF */
        sha256_no_nlf(m, s1);
        sha256_no_nlf(m2, s2);
        for(int r=0;r<=64;r++)
            for(int w=0;w<8;w++)
                for(int b=0;b<32;b++)
                    if(((s1[r][w]^s2[r][w])>>b)&1)
                        flip_count_nonlf[r][w*32+b]++;
    }

    printf("DETERMINISTIC BITS PER ROUND (N=%d messages)\n", N);
    printf("  = bits that flip in 100%% of messages when W[0][b0] flips\n\n");
    printf("%-6s | %-10s %-10s | %-10s %-10s\n",
           "Round", "Full:det", "Full:>90%", "NoNLF:det", "NoNLF:>90%");
    printf("-------+------------------------+------------------------\n");

    for(int r=0;r<=64;r++) {
        int df=0, hf=0, dn=0, hn=0;
        for(int i=0;i<256;i++) {
            if(flip_count_full[r][i]==N) df++;
            if(flip_count_full[r][i]>N*9/10) hf++;
            if(flip_count_nonlf[r][i]==N) dn++;
            if(flip_count_nonlf[r][i]>N*9/10) hn++;
        }
        if(r<=8 || r==16 || r==32 || r==48 || r==64 || df>0 || dn != det_nonlf[r>0?r-1:0]) {
            printf("  r=%2d | %6d     %6d     | %6d     %6d\n", r, df, hf, dn, hn);
        }
        det_full[r]=df; det_nonlf[r]=dn;
    }

    /* Test 2: WHICH bits are deterministic in No-NLF?
     * At round 64: 128 det bits. Which registers? Which bit positions? */
    printf("\nDETERMINISTIC BITS IN NO-NLF at round 64:\n");
    printf("  Register  | Det bits | Bit positions\n");
    printf("  ----------+----------+---------------------\n");
    const char* rnames[] = {"a","b","c","d","e","f","g","h"};
    for(int w=0;w<8;w++) {
        int count=0;
        printf("  %-10s| ", rnames[w]);
        for(int b=0;b<32;b++)
            if(flip_count_nonlf[64][w*32+b]==N) count++;
        printf("%3d      | ", count);
        for(int b=0;b<32;b++)
            if(flip_count_nonlf[64][w*32+b]==N) printf("%d,",b);
        printf("\n");
    }

    /* Test 3: NLF CROSS-TERM per round.
     * Define: cross_term[r] = (full SHA diff) XOR (no-NLF diff).
     * This isolates the CONTRIBUTION of Ch/Maj to the carry chain.
     * If cross_term = 0: NLF had no carry effect.
     * If cross_term != 0: NLF-carry interaction at this round. */
    printf("\nNLF CROSS-TERM: (Full diff) XOR (No-NLF diff) per round\n");
    printf("  = bits where carry×NLF interaction changes the outcome\n\n");

    sha256_full(msg, (uint32_t(*)[8])flip_count_full); /* reuse buffer */
    /* Actually recompute properly */
    {
        uint32_t sf1[65][8], sf2[65][8], sn1[65][8], sn2[65][8];
        sha256_full(msg, sf1); sha256_full(msg2, sf2);
        sha256_no_nlf(msg, sn1); sha256_no_nlf(msg2, sn2);

        printf("  Round | Full_diff | NoNLF_diff | Cross_term (XOR) | NLF effect\n");
        printf("  ------+-----------+------------+------------------+-----------\n");
        for(int r=0;r<=64;r++) {
            int df=0, dn=0, ct=0;
            for(int w=0;w<8;w++) {
                uint32_t full_d = sf1[r][w]^sf2[r][w];
                uint32_t nonlf_d = sn1[r][w]^sn2[r][w];
                uint32_t cross = full_d ^ nonlf_d;
                df += __builtin_popcount(full_d);
                dn += __builtin_popcount(nonlf_d);
                ct += __builtin_popcount(cross);
            }
            if(r<=8 || r==16 || r==32 || r==64)
                printf("  %5d | %6d    | %7d    | %9d        | %s\n",
                       r, df, dn, ct,
                       ct==0?"none":ct<10?"small":ct<50?"medium":"LARGE");
        }
    }

    return 0;
}
