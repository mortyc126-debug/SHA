#include <stdio.h>
#include <stdint.h>
#include <time.h>

#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))
#define SHR(x,n)  ((x)>>(n))
#define CH(e,f,g) (((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))
#define SIG0(x) (ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define SIG1(x) (ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define sig0(x) (ROTR(x,7)^ROTR(x,18)^SHR(x,3))
#define sig1(x) (ROTR(x,17)^ROTR(x,19)^SHR(x,10))

static const uint32_t H0[8] = {
    0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
    0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19
};
static const uint32_t K[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

/* Compute specific output bit of SHA-256(W0||0..0) for 64 rounds */
static inline int sha256_bit(uint32_t w0, int word, int bit) {
    uint32_t W[64];
    W[0] = w0;
    for (int i = 1; i < 16; i++) W[i] = 0;
    for (int i = 16; i < 64; i++)
        W[i] = sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16];

    uint32_t a=H0[0],b=H0[1],c=H0[2],d=H0[3],
             e=H0[4],f=H0[5],g=H0[6],h=H0[7];
    for (int r = 0; r < 64; r++) {
        uint32_t T1 = h + SIG1(e) + CH(e,f,g) + K[r] + W[r];
        uint32_t T2 = SIG0(a) + MAJ(a,b,c);
        h=g; g=f; f=e; e=d+T1; d=c; c=b; b=a; a=T1+T2;
    }
    uint32_t out[8] = {a+H0[0],b+H0[1],c+H0[2],d+H0[3],
                       e+H0[4],f+H0[5],g+H0[6],h+H0[7]};
    return (out[word] >> bit) & 1;
}

/* Compute top coefficient: XOR of f(x) over all 2^n inputs */
int main(int argc, char *argv[]) {
    int n = 24;
    if (argc > 1) n = atoi(argv[1]);

    /* Targets: (word, bit) */
    int targets[][2] = {{7,23},{4,0},{7,0},{0,15}};
    int n_targets = 4;

    printf("DEGREE n=%d — C implementation\n", n);
    printf("2^%d = %lu evaluations per bit\n\n", n, 1UL<<n);

    uint64_t N = 1UL << n;

    for (int t = 0; t < n_targets; t++) {
        int w = targets[t][0], bi = targets[t][1];
        clock_t t0 = clock();

        int xor_all = 0;
        for (uint64_t x = 0; x < N; x++) {
            xor_all ^= sha256_bit((uint32_t)x, w, bi);
        }

        double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
        const char *status = (xor_all == 0) ? "DEFICIT *" : "FULL";
        printf("  H[%d][b%2d]: top=%d  %s  (%.1fs)\n", w, bi, xor_all, status, elapsed);
        fflush(stdout);
    }

    return 0;
}
