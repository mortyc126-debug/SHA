/*
 * SIGNAL KILLER: какая из 3 функций глушит сигнал?
 *
 * 3 механизма SHA-256:
 *   1. CARRY (mod 2^32 addition vs XOR)
 *   2. ROTATION (Σ₀, Σ₁, σ₀, σ₁)
 *   3. NONLINEAR FUNCTIONS (Ch, Maj)
 *
 * Отключаем по одному, измеряем: насколько быстро сигнал
 * достигает 50% (saturation) и можно ли его ОБРАТИТЬ.
 *
 * gcc -O3 -march=native -o sigkill signal_killer.c
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

/* Flags for disabling components */
#define F_FULL     0
#define F_NO_CARRY 1   /* + → ^ */
#define F_NO_ROT   2   /* Σ₀=Σ₁=σ₀=σ₁=identity */
#define F_NO_NLF   4   /* Ch=XOR3, Maj=XOR3 */
#define F_NO_SCHED 8   /* W[16..63]=0 (no schedule expansion) */

void sha256_variant(const uint32_t msg[16], uint32_t states[65][8], int flags) {
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=msg[i];

    if(flags & F_NO_SCHED) {
        for(int i=16;i<64;i++) W[i]=0;
    } else if(flags & F_NO_CARRY) {
        for(int i=16;i<64;i++) {
            uint32_t s1 = (flags&F_NO_ROT) ? W[i-2] : (ROTR(W[i-2],17)^ROTR(W[i-2],19)^(W[i-2]>>10));
            uint32_t s0 = (flags&F_NO_ROT) ? W[i-15] : (ROTR(W[i-15],7)^ROTR(W[i-15],18)^(W[i-15]>>3));
            W[i] = s1 ^ W[i-7] ^ s0 ^ W[i-16];
        }
    } else {
        for(int i=16;i<64;i++) {
            uint32_t s1 = (flags&F_NO_ROT) ? W[i-2] : (ROTR(W[i-2],17)^ROTR(W[i-2],19)^(W[i-2]>>10));
            uint32_t s0 = (flags&F_NO_ROT) ? W[i-15] : (ROTR(W[i-15],7)^ROTR(W[i-15],18)^(W[i-15]>>3));
            W[i] = s1 + W[i-7] + s0 + W[i-16];
        }
    }

    for(int i=0;i<8;i++) states[0][i]=IV[i];

    for(int r=0;r<64;r++){
        uint32_t a=states[r][0],b=states[r][1],c=states[r][2],d=states[r][3];
        uint32_t e=states[r][4],f=states[r][5],g=states[r][6],h=states[r][7];

        uint32_t S1 = (flags&F_NO_ROT) ? e : (ROTR(e,6)^ROTR(e,11)^ROTR(e,25));
        uint32_t S0 = (flags&F_NO_ROT) ? a : (ROTR(a,2)^ROTR(a,13)^ROTR(a,22));
        uint32_t ch = (flags&F_NO_NLF) ? (e^f^g) : ((e&f)^(~e&g));
        uint32_t mj = (flags&F_NO_NLF) ? (a^b^c) : ((a&b)^(a&c)^(b&c));

        uint32_t t1, t2;
        if(flags & F_NO_CARRY) {
            t1 = h^S1^ch^K[r]^W[r];
            t2 = S0^mj;
            states[r+1][0]=t1^t2;
            states[r+1][4]=d^t1;
        } else {
            t1 = h+S1+ch+K[r]+W[r];
            t2 = S0+mj;
            states[r+1][0]=t1+t2;
            states[r+1][4]=d+t1;
        }
        states[r+1][1]=a; states[r+1][2]=b; states[r+1][3]=c;
        states[r+1][5]=e; states[r+1][6]=f; states[r+1][7]=g;
    }
}

int hw256(const uint32_t a[8], const uint32_t b[8]) {
    int h=0; for(int i=0;i<8;i++) h+=__builtin_popcount(a[i]^b[i]); return h;
}

void measure_signal(const char *name, int flags) {
    uint32_t msg[16], msg2[16], st1[65][8], st2[65][8];
    srand(42);
    for(int i=0;i<16;i++) msg[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
    memcpy(msg2, msg, 64);
    msg2[0] ^= 1; /* flip W[0] bit 0 */

    sha256_variant(msg, st1, flags);
    sha256_variant(msg2, st2, flags);

    printf("%-28s |", name);
    int sat_round = -1;
    for(int r=0;r<=64;r++){
        int d=hw256(st1[r],st2[r]);
        if(r<=1 || r==2 || r==4 || r==6 || r==8 || r==16 || r==32 || r==64) {
            printf(" r%d=%d", r, d);
        }
        if(d >= 120 && sat_round < 0) sat_round = r;
    }
    int final_d = hw256(st1[64], st2[64]);
    printf(" | sat@%d final=%d\n", sat_round, final_d);
}

int main() {
    printf("SIGNAL KILLER: which component kills the signal?\n");
    printf("=================================================\n\n");
    printf("Flip W[0] bit 0, track signal through 64 rounds.\n");
    printf("sat = first round where diff >= 120/256.\n\n");

    printf("%-28s | Signal trace                            | Summary\n", "Config");
    printf("------------------------------+------------------------------------------+---------\n");

    measure_signal("FULL SHA-256",               F_FULL);
    measure_signal("No carry (XOR only)",         F_NO_CARRY);
    measure_signal("No rotation (Sig=id)",        F_NO_ROT);
    measure_signal("No Ch/Maj (use XOR3)",        F_NO_NLF);
    measure_signal("No carry + no rotation",      F_NO_CARRY|F_NO_ROT);
    measure_signal("No carry + no NLF",           F_NO_CARRY|F_NO_NLF);
    measure_signal("No rotation + no NLF",        F_NO_ROT|F_NO_NLF);
    measure_signal("ONLY shift register",         F_NO_CARRY|F_NO_ROT|F_NO_NLF);
    measure_signal("No schedule expansion",       F_NO_SCHED);

    /* Now: SIGNAL REVERSIBILITY test.
     * For each config: can we RECOVER the input bit from the output?
     * Method: flip W[0] bit 0, look at hash diff.
     * Count: how many hash bits are DETERMINISTICALLY linked to input bit.
     * = how many hash bits flip EVERY TIME we flip W[0][b0] (across many messages) */
    printf("\n\nSIGNAL REVERSIBILITY: which hash bits ALWAYS flip with W[0][b0]?\n");
    printf("=================================================================\n\n");

    int configs[] = {F_FULL, F_NO_CARRY, F_NO_ROT, F_NO_NLF,
                     F_NO_CARRY|F_NO_ROT, F_NO_CARRY|F_NO_NLF,
                     F_NO_CARRY|F_NO_ROT|F_NO_NLF};
    const char *cnames[] = {"FULL","No carry","No rot","No NLF",
                            "No carry+rot","No carry+NLF","ONLY shift"};
    int nc = 7;

    for(int ci=0; ci<nc; ci++) {
        int always_flip[256] = {0};
        int N = 1000;
        srand(ci * 999 + 42);

        for(int trial=0; trial<N; trial++) {
            uint32_t m[16], m2[16], s1[65][8], s2[65][8];
            for(int i=0;i<16;i++) m[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
            memcpy(m2, m, 64);
            m2[0] ^= 1;

            sha256_variant(m, s1, configs[ci]);
            sha256_variant(m2, s2, configs[ci]);

            for(int w=0;w<8;w++)
                for(int b=0;b<32;b++)
                    if(((s1[64][w]^s2[64][w])>>b)&1)
                        always_flip[w*32+b]++;
        }

        int deterministic = 0;
        int high_corr = 0; /* >90% flip rate */
        int any_corr = 0;  /* >60% flip rate */
        for(int i=0;i<256;i++) {
            if(always_flip[i] == N) deterministic++;
            if(always_flip[i] > N*9/10) high_corr++;
            if(always_flip[i] > N*6/10) any_corr++;
        }

        printf("  %-20s: deterministic=%d, >90%%=%d, >60%%=%d\n",
               cnames[ci], deterministic, high_corr, any_corr);
    }

    printf("\n  'deterministic' = hash bits that ALWAYS flip when W[0][b0] flips.\n");
    printf("  These bits carry a RECOVERABLE signal from input to output.\n");
    printf("  If deterministic > 0: signal SURVIVES all 64 rounds!\n");

    return 0;
}
