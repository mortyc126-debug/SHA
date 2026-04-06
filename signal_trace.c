/*
 * SIGNAL TRACING: track ONE input bit through all 64 rounds.
 *
 * Method: flip bit (w,b) of message. At EACH round, record
 * which state bits changed. This is the "signal wavefront".
 *
 * The wavefront shows HOW the signal propagates:
 * - Which rounds it's in
 * - How many bits it affects
 * - WHERE in the state it lives (a,b,c,d,e,f,g,h)
 * - Whether it GROWS or SHRINKS
 *
 * Then: from HASH backward — can we trace which input bits
 * contributed to each hash bit?
 *
 * Compile: gcc -O3 -march=native -o sigtrace signal_trace.c
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

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

void sha256_trace(const uint32_t msg[16], uint32_t states[65][8]) {
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
    for(int i=0;i<8;i++) states[0][i]=IV[i];
    for(int r=0;r<64;r++){
        uint32_t a=states[r][0],b=states[r][1],c=states[r][2],d=states[r][3];
        uint32_t e=states[r][4],f=states[r][5],g=states[r][6],h=states[r][7];
        uint32_t t1=h+S1(e)+CH(e,f,g)+K[r]+W[r], t2=S0(a)+MAJ(a,b,c);
        states[r+1][0]=t1+t2; states[r+1][1]=a; states[r+1][2]=b; states[r+1][3]=c;
        states[r+1][4]=d+t1;  states[r+1][5]=e; states[r+1][6]=f; states[r+1][7]=g;
    }
}

int main() {
    printf("SIGNAL TRACING: one bit through 64 rounds\n");
    printf("==========================================\n\n");

    uint32_t msg[16], msg2[16];
    uint32_t st1[65][8], st2[65][8];

    /* Fixed message */
    srand(42);
    for(int i=0;i<16;i++) msg[i]=(uint32_t)rand()|((uint32_t)rand()<<16);

    sha256_trace(msg, st1);

    /* Test several input bits */
    int test_bits[][2] = {{0,0},{0,15},{0,31},{7,0},{15,0},{15,31}};
    int n_tests = 6;

    for(int t=0; t<n_tests; t++) {
        int mw = test_bits[t][0], mb = test_bits[t][1];
        memcpy(msg2, msg, 64);
        msg2[mw] ^= (1u << mb);
        sha256_trace(msg2, st2);

        printf("INPUT: W[%d] bit %d\n", mw, mb);
        printf("  Round |  a   b   c   d   e   f   g   h | Total | Signal path\n");
        printf("  ------+--------------------------------+-------+------------\n");

        for(int r=0; r<=64; r++) {
            int per_reg[8], total=0;
            for(int i=0;i<8;i++) {
                per_reg[i] = __builtin_popcount(st1[r][i]^st2[r][i]);
                total += per_reg[i];
            }

            /* Show interesting rounds */
            if(r <= 5 || r >= 60 || r == mw || r == mw+1 || r == mw+4 ||
               r == 16 || r == 32 || total <= 5 || total >= 125) {
                printf("  %5d | %2d  %2d  %2d  %2d  %2d  %2d  %2d  %2d | %5d | ",
                       r, per_reg[0],per_reg[1],per_reg[2],per_reg[3],
                       per_reg[4],per_reg[5],per_reg[6],per_reg[7], total);

                /* Describe signal path */
                if(r == 0) printf("IV (no diff)");
                else if(r == mw+1 && mw < 16) {
                    printf("W[%d] enters → a,e affected", mw);
                } else if(total <= 3) {
                    /* Which regs? */
                    for(int i=0;i<8;i++) if(per_reg[i]>0) {
                        const char* names[]={"a","b","c","d","e","f","g","h"};
                        printf("%s(%d) ",names[i],per_reg[i]);
                    }
                } else if(total >= 120) printf("SATURATED");
                else if(per_reg[0]==0 && per_reg[4]==0) printf("only shifted regs");
                printf("\n");
            }
        }

        /* Hash diff */
        int hash_diff = 0;
        uint32_t h1[8], h2[8];
        for(int i=0;i<8;i++) {
            h1[i]=st1[64][i]+IV[i];
            h2[i]=st2[64][i]+IV[i];
            hash_diff += __builtin_popcount(h1[i]^h2[i]);
        }
        printf("  Hash diff: %d/256 bits\n\n", hash_diff);
    }

    /* REVERSE TRACE: from hash bit, which input bits affect it? */
    printf("══════════════════════════════════════════\n");
    printf("REVERSE SIGNAL: hash bit → input bits\n");
    printf("══════════════════════════════════════════\n\n");

    /* For hash bit H[0][b0]: flip each of 512 msg bits,
     * check if H[0][b0] flips. Count and identify. */
    sha256_trace(msg, st1);
    uint32_t ref_h0 = st1[64][0] + IV[0];

    for(int hw=0; hw<8; hw+=3) {
        for(int hb=0; hb<32; hb+=15) {
            uint32_t ref_hash[8];
            for(int i=0;i<8;i++) ref_hash[i]=st1[64][i]+IV[i];
            int ref_bit = (ref_hash[hw] >> hb) & 1;

            int sources[16] = {0}; /* per-word count */
            int total = 0;

            for(int mw=0;mw<16;mw++) for(int mb=0;mb<32;mb++) {
                memcpy(msg2, msg, 64);
                msg2[mw] ^= (1u << mb);
                sha256_trace(msg2, st2);
                uint32_t new_hash = st2[64][hw] + IV[hw];
                if(((new_hash >> hb) & 1) != ref_bit) {
                    sources[mw]++;
                    total++;
                }
            }

            printf("  H[%d][b%2d] depends on %d/512 msg bits: ", hw, hb, total);
            for(int mw=0;mw<16;mw++) if(sources[mw]>0)
                printf("W%d:%d ", mw, sources[mw]);
            printf("\n");
        }
    }

    /* KEY EXPERIMENT: signal SURVIVAL through rounds.
     * For input W[0][b0]: at which round does it first
     * affect EVERY state bit? (= full diffusion round) */
    printf("\n══════════════════════════════════════════\n");
    printf("SIGNAL DIFFUSION SPEED\n");
    printf("══════════════════════════════════════════\n\n");

    memcpy(msg2, msg, 64);
    msg2[0] ^= 1; /* flip W[0] bit 0 */
    sha256_trace(msg, st1);
    sha256_trace(msg2, st2);

    int first_full = -1;
    int first_half = -1;
    printf("  W[0][b0] signal propagation:\n");
    for(int r=0;r<=64;r++) {
        int total = 0;
        for(int i=0;i<8;i++) total+=__builtin_popcount(st1[r][i]^st2[r][i]);
        if(total >= 128 && first_half < 0) first_half = r;
        if(total >= 240 && first_full < 0) first_full = r;
        if(r<=8 || r==16 || r==first_half || r==first_full || r==64)
            printf("    r=%2d: %3d/256 bits affected%s%s\n", r, total,
                   r==first_half?" ← HALF":"", r==first_full?" ← FULL":"");
    }
    printf("  First 50%% affected: round %d\n", first_half);
    printf("  First 94%% affected: round %d\n", first_full);

    /* SIGNAL DECAY from known side:
     * From hash backward, state[64] is known.
     * If we flip msg bit: state[64] changes.
     * How many of those changes are in the 128-bit GAP (a[53..56])
     * vs the KNOWN zone (a[57..64])? */
    printf("\n══════════════════════════════════════════\n");
    printf("SIGNAL PARTITION: known zone vs gap zone\n");
    printf("══════════════════════════════════════════\n\n");

    printf("  When flipping W[k] bit 0, how does δstate[64] distribute?\n");
    printf("  Known zone: H[1..3,5..7] (192 bits, from shift of a[57..63],e[57..63])\n");
    printf("  Active zone: H[0,4] (64 bits, from round 63 output a[64],e[64])\n\n");

    for(int mw=0; mw<16; mw+=3) {
        memcpy(msg2, msg, 64);
        msg2[mw] ^= 1;
        sha256_trace(msg2, st2);

        int known_diff = 0, active_diff = 0;
        for(int i=0;i<8;i++) {
            int d = __builtin_popcount(st1[64][i]^st2[64][i]);
            if(i==0 || i==4) active_diff += d;
            else known_diff += d;
        }
        printf("  W[%2d]: known_zone=%3d/192, active_zone=%2d/64, total=%3d/256\n",
               mw, known_diff, active_diff, known_diff+active_diff);
    }

    printf("\n  If active_zone ≈ 16 (= 64*0.25) consistently:\n");
    printf("  → signal distributes UNIFORMLY across all hash words.\n");
    printf("  → No concentration in gap/known zone.\n");
    printf("  → SHA-256 mixes perfectly even at state[64] level.\n");

    return 0;
}
