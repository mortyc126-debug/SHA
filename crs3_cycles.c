/*
 * CRS-3 Cycle Anomaly: root cause analysis.
 * Counts cycles of g32(x) = SHA-256(x||0..0)[word] with variants.
 *
 * Compile: gcc -O3 -march=native -o cycles crs3_cycles.c
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))
#define CH(e,f,g) (((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))
#define SIG0(x) (ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define SIG1(x) (ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define sig0(x) (ROTR(x,7)^ROTR(x,18)^((x)>>3))
#define sig1(x) (ROTR(x,17)^ROTR(x,19)^((x)>>10))

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

/* g32: SHA-256(x||pad)[out_word], input at position in_pos */
static inline uint32_t g32(uint32_t x, int in_pos, int out_word,
                            int rounds, const uint32_t pad[15]) {
    uint32_t W[64];
    for (int i=0;i<16;i++) W[i] = (i==in_pos) ? x : pad[i < 1 ? 0 : i-1];
    /* Fix: build msg properly */
    uint32_t msg[16];
    for (int i=0;i<16;i++) msg[i]=0;
    msg[in_pos] = x;
    for (int i=0;i<15;i++) { int j=(in_pos+1+i)%16; msg[j]=pad[i]; }
    for (int i=0;i<16;i++) W[i]=msg[i];
    for (int i=16;i<64;i++) W[i]=sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16];

    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for (int r=0;r<rounds;r++){
        uint32_t t1=h+SIG1(e)+CH(e,f,g)+K[r]+W[r];
        uint32_t t2=SIG0(a)+MAJ(a,b,c);
        h=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
    }
    uint32_t hash[8]={a+IV[0],b+IV[1],c+IV[2],d+IV[3],e+IV[4],f+IV[5],g+IV[6],h+IV[7]};
    return hash[out_word];
}

/* Simple g32 with zeros padding */
static inline uint32_t g32_simple(uint32_t x) {
    uint32_t W[64]={0}; W[0]=x;
    for (int i=16;i<64;i++) W[i]=sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for (int r=0;r<64;r++){
        uint32_t t1=h+SIG1(e)+CH(e,f,g)+K[r]+W[r];
        uint32_t t2=SIG0(a)+MAJ(a,b,c);
        h=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
    }
    return a+IV[0];
}

/* g32 with XOR instead of + (no carry) */
static inline uint32_t g32_xor(uint32_t x) {
    uint32_t W[64]={0}; W[0]=x;
    for (int i=16;i<64;i++) W[i]=sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for (int r=0;r<64;r++){
        uint32_t t1=h^SIG1(e)^CH(e,f,g)^K[r]^W[r];
        uint32_t t2=SIG0(a)^MAJ(a,b,c);
        h=g;g=f;f=e;e=d^t1;d=c;c=b;b=a;a=t1^t2;
    }
    return a^IV[0];
}

/* g32 reduced rounds */
static inline uint32_t g32_rounds(uint32_t x, int R) {
    uint32_t W[64]={0}; W[0]=x;
    for (int i=16;i<64;i++) W[i]=sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for (int r=0;r<R;r++){
        uint32_t t1=h+SIG1(e)+CH(e,f,g)+K[r]+W[r];
        uint32_t t2=SIG0(a)+MAJ(a,b,c);
        h=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
    }
    return a+IV[0];
}

/* g32 output word selection */
static inline uint32_t g32_word(uint32_t x, int w) {
    uint32_t W[64]={0}; W[0]=x;
    for (int i=16;i<64;i++) W[i]=sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for (int r=0;r<64;r++){
        uint32_t t1=h+SIG1(e)+CH(e,f,g)+K[r]+W[r];
        uint32_t t2=SIG0(a)+MAJ(a,b,c);
        h=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
    }
    uint32_t hash[8]={a+IV[0],b+IV[1],c+IV[2],d+IV[3],e+IV[4],f+IV[5],g+IV[6],h+IV[7]};
    return hash[w];
}

/* g32 input at position pos */
static inline uint32_t g32_pos(uint32_t x, int pos) {
    uint32_t W[64]={0}; W[pos]=x;
    for (int i=16;i<64;i++) W[i]=sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for (int r=0;r<64;r++){
        uint32_t t1=h+SIG1(e)+CH(e,f,g)+K[r]+W[r];
        uint32_t t2=SIG0(a)+MAJ(a,b,c);
        h=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
    }
    return a+IV[0];
}

typedef struct { uint32_t lam; uint32_t canon; } Cycle;

int find_cycles(uint32_t (*f)(uint32_t), int n_starts, uint32_t max_steps,
                Cycle *out, int max_cycles) {
    int n_cycles = 0;
    srand(42);
    for (int s=0; s<n_starts && n_cycles<max_cycles; s++) {
        uint32_t x0 = (uint32_t)rand() | ((uint32_t)rand()<<16);
        /* Brent */
        uint32_t power=1, lam=1;
        uint32_t tort=x0, hare=f(x0);
        while (tort != hare && lam < max_steps) {
            if (power==lam) { tort=hare; power*=2; lam=0; }
            hare=f(hare); lam++;
        }
        if (lam >= max_steps) continue;
        /* Find canonical */
        uint32_t p=tort, canon=p;
        for (uint32_t i=0;i<lam;i++) { p=f(p); if(p<canon) canon=p; }
        /* Check if new */
        int found=0;
        for (int c=0;c<n_cycles;c++)
            if (out[c].canon==canon) { found=1; break; }
        if (!found) {
            out[n_cycles].lam=lam;
            out[n_cycles].canon=canon;
            n_cycles++;
        }
    }
    return n_cycles;
}

/* Wrappers for function pointers */
static uint32_t wrap_simple(uint32_t x) { return g32_simple(x); }
static uint32_t wrap_xor(uint32_t x) { return g32_xor(x); }
static uint32_t wrap_r8(uint32_t x) { return g32_rounds(x,8); }
static uint32_t wrap_r16(uint32_t x) { return g32_rounds(x,16); }
static uint32_t wrap_r32(uint32_t x) { return g32_rounds(x,32); }
static uint32_t wrap_r48(uint32_t x) { return g32_rounds(x,48); }
static uint32_t wrap_w1(uint32_t x) { return g32_word(x,1); }
static uint32_t wrap_w4(uint32_t x) { return g32_word(x,4); }
static uint32_t wrap_w7(uint32_t x) { return g32_word(x,7); }
static uint32_t wrap_pos1(uint32_t x) { return g32_pos(x,1); }
static uint32_t wrap_pos7(uint32_t x) { return g32_pos(x,7); }
static uint32_t wrap_pos15(uint32_t x) { return g32_pos(x,15); }

int main() {
    printf("CRS-3 CYCLE ANOMALY: ROOT CAUSE (C implementation)\n");
    printf("====================================================\n\n");

    Cycle cycles[100];
    int n;
    int starts=30;
    uint32_t max_s=2000000;

    typedef struct { const char* name; uint32_t (*func)(uint32_t); } Test;
    Test tests[] = {
        {"Standard g32 (R=64)", wrap_simple},
        {"No carry (XOR, R=64)", wrap_xor},
        {"R=8 rounds", wrap_r8},
        {"R=16 rounds", wrap_r16},
        {"R=32 rounds", wrap_r32},
        {"R=48 rounds", wrap_r48},
        {"Output H[1]", wrap_w1},
        {"Output H[4]", wrap_w4},
        {"Output H[7]", wrap_w7},
        {"Input at W[1]", wrap_pos1},
        {"Input at W[7]", wrap_pos7},
        {"Input at W[15]", wrap_pos15},
    };
    int n_tests = sizeof(tests)/sizeof(tests[0]);

    printf("%-25s | Cycles | Lengths\n", "Config");
    printf("------------------------------------------------------\n");

    for (int t=0; t<n_tests; t++) {
        clock_t start = clock();
        n = find_cycles(tests[t].func, starts, max_s, cycles, 100);
        double elapsed = (double)(clock()-start)/CLOCKS_PER_SEC;

        printf("%-25s | %4d   | ", tests[t].name, n);
        for (int i=0; i<n && i<5; i++) printf("%u ", cycles[i].lam);
        printf(" (%.1fs)\n", elapsed);
    }

    return 0;
}
