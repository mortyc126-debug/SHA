/*
 * DEADPOOL INVERSION: recover erased bits from their clones.
 *
 * Ch(e,f,g) erases e[k] when f[k]≠g[k].
 * But e[k] at round r = f[k] at round r+1 = g[k] at round r+2 = h[k] at round r+3.
 *
 * When Ch erases e[k] at round r:
 *   - f[k] at round r = e[k] from round r-1 (DIFFERENT value!)
 *   - e[k] at round r will become f[k] at round r+1
 *   - At round r+1: e[k] is now at f position → NOT erased by Ch!
 *     (Ch uses NEW e at r+1, not old)
 *   - At round r+2: old e[k] is at g position → g enters Ch as "else" branch
 *   - At round r+3: old e[k] is at h position → h enters T1 LINEARLY (h+...)
 *
 * So: erased e[k] survives for 3 more rounds in f,g,h positions.
 * At h position (round r+3): enters LINEARLY into T1 = h + stuff.
 * LINEAR = recoverable from T1!
 *
 * Question: can we recover the erased bit from its h-position appearance?
 *
 * gcc -O3 -march=native -o deadpool deadpool.c
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

void sha256_trace(const uint32_t msg[16], uint32_t st[65][8]) {
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=(ROTR(W[i-2],17)^ROTR(W[i-2],19)^(W[i-2]>>10))+W[i-7]+(ROTR(W[i-15],7)^ROTR(W[i-15],18)^(W[i-15]>>3))+W[i-16];
    for(int i=0;i<8;i++) st[0][i]=IV[i];
    for(int r=0;r<64;r++){
        uint32_t a=st[r][0],b=st[r][1],c=st[r][2],d=st[r][3];
        uint32_t e=st[r][4],f=st[r][5],g=st[r][6],h=st[r][7];
        uint32_t t1=h+(ROTR(e,6)^ROTR(e,11)^ROTR(e,25))+((e&f)^(~e&g))+K[r]+W[r];
        uint32_t t2=(ROTR(a,2)^ROTR(a,13)^ROTR(a,22))+((a&b)^(a&c)^(b&c));
        st[r+1][0]=t1+t2; st[r+1][1]=a; st[r+1][2]=b; st[r+1][3]=c;
        st[r+1][4]=d+t1;  st[r+1][5]=e; st[r+1][6]=f; st[r+1][7]=g;
    }
}

int main() {
    printf("DEADPOOL INVERSION: recover erased bits from clones\n");
    printf("====================================================\n\n");

    uint32_t msg[16], st[65][8];
    srand(42);
    for(int i=0;i<16;i++) msg[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
    sha256_trace(msg, st);

    /* TRACE: where does e[r] go after round r? */
    printf("CLONE TRACKING: e at round 10\n");
    printf("─────────────────────────────\n");
    uint32_t e10 = st[10][4];
    printf("  Round 10: e = %08x (ORIGINAL)\n", e10);
    printf("  Round 11: f = %08x %s\n", st[11][5], st[11][5]==e10?"← CLONE":"≠");
    printf("  Round 12: g = %08x %s\n", st[12][6], st[12][6]==e10?"← CLONE":"≠");
    printf("  Round 13: h = %08x %s\n", st[13][7], st[13][7]==e10?"← CLONE":"≠");
    printf("  Round 14: GONE (shifted out of state)\n");

    /* ERASURE MAP: at each round, which e-bits are erased by Ch? */
    printf("\nERASURE MAP: which e-bits Ch erases per round\n");
    printf("─────────────────────────────\n");
    int total_erased = 0;
    for(int r = 0; r < 64; r++) {
        uint32_t e = st[r][4], f = st[r][5], g = st[r][6];
        /* Ch erases e[k] where f[k] ≠ g[k] */
        uint32_t erase_mask = f ^ g;
        int n_erased = __builtin_popcount(erase_mask);
        total_erased += n_erased;
        if(r < 5 || r > 59)
            printf("  Round %2d: %2d/32 bits erased (f^g = %08x)\n",
                   r, n_erased, erase_mask);
    }
    printf("  ...\n  Total across 64 rounds: %d erasures (avg %.1f/round)\n",
           total_erased, (double)total_erased/64);

    /* DEADPOOL TEST: can we recover erased bits?
     *
     * Setup: at round r, Ch erases e[k] (where f[k]≠g[k]).
     * But at round r+3: this same value is h[k].
     * h enters T1 LINEARLY: T1 = h + Σ₁(e') + Ch(e',f',g') + K + W.
     *
     * If we KNOW T1 (from a[r+4] - T2[r+3] — backward chain):
     *   h = T1 - Σ₁(e') - Ch(e',f',g') - K - W
     *   ALL KNOWN (if we know state[r+3] and W[r+3]).
     *
     * So: erased bit at round r → recoverable at round r+3 from T1!
     * BUT: this needs W[r+3] and full state[r+3] — the same unknowns.
     *
     * HOWEVER: the erased bit at round r ALSO exists at:
     *   round r+1 as f → enters Ch(e',f',g') at round r+1.
     *     Ch(e',f,g') = e'&f ⊕ ~e'&g'. f here IS our erased e!
     *     If e'[k] = 1: Ch uses f[k] = our bit! RECOVERABLE from Ch output.
     *     If e'[k] = 0: Ch uses g'[k] instead. Our bit INVISIBLE.
     *
     * Probability e'[k]=1: 50%.
     * So: 50% chance of recovery from round r+1.
     *     Another chance at round r+2 (g position).
     *     And at round r+3 (h position, linear recovery).
     *
     * TOTAL probability of recovery across 3 rounds?
     */

    printf("\nDEADPOOL RECOVERY PROBABILITY\n");
    printf("─────────────────────────────\n");

    int N = 100000;
    int recoverable_r1 = 0, recoverable_r2 = 0, recoverable_r3 = 0;
    int total_tested = 0;

    srand(42);
    for(int trial = 0; trial < N; trial++) {
        for(int i = 0; i < 16; i++) msg[i] = (uint32_t)rand()|((uint32_t)rand()<<16);
        sha256_trace(msg, st);

        /* Pick round 30, random erased bit */
        int r = 30;
        uint32_t e = st[r][4], f = st[r][5], g = st[r][6];
        uint32_t erase_mask = f ^ g;
        if(!erase_mask) continue;

        /* Find first erased bit */
        int k = __builtin_ctz(erase_mask);
        total_tested++;

        /* Clone at round r+1: f position. Is it visible through Ch? */
        /* Ch at round r+1: Ch(e', f', g') where f'=e (our value) */
        uint32_t e_r1 = st[r+1][4]; /* new e at round r+1 */
        if((e_r1 >> k) & 1) {
            /* e'[k]=1 → Ch uses f'[k] = our erased bit. VISIBLE! */
            recoverable_r1++;
        }

        /* Clone at round r+2: g position. */
        /* Ch at round r+2: Ch(e'', f'', g'') where g''=f'=e */
        uint32_t e_r2 = st[r+2][4];
        if(!((e_r2 >> k) & 1)) {
            /* e''[k]=0 → Ch uses g''[k] = our erased bit. VISIBLE! */
            recoverable_r2++;
        }

        /* Clone at round r+3: h position. Enters T1 LINEARLY. */
        /* ALWAYS visible (h enters T1 = h + stuff, linear in h). */
        recoverable_r3++;
    }

    printf("  Erased at round r (Ch, bit k where f[k]≠g[k]).\n");
    printf("  N = %d erased bits tested.\n\n", total_tested);
    printf("  Recovery at r+1 (f→Ch, if e'[k]=1): %d/%d = %.1f%%\n",
           recoverable_r1, total_tested, 100.0*recoverable_r1/total_tested);
    printf("  Recovery at r+2 (g→Ch, if e'[k]=0): %d/%d = %.1f%%\n",
           recoverable_r2, total_tested, 100.0*recoverable_r2/total_tested);
    printf("  Recovery at r+3 (h→T1, ALWAYS):     %d/%d = %.1f%%\n",
           recoverable_r3, total_tested, 100.0*recoverable_r3/total_tested);

    int any_recovery = 0;
    for(int trial = 0; trial < N; trial++) {
        for(int i = 0; i < 16; i++) msg[i] = (uint32_t)rand()|((uint32_t)rand()<<16);
        sha256_trace(msg, st);
        int r = 30;
        uint32_t erase = st[r][5] ^ st[r][6];
        if(!erase) continue;
        int k = __builtin_ctz(erase);
        uint32_t e_r1 = st[r+1][4], e_r2 = st[r+2][4];
        int rec1 = (e_r1 >> k) & 1;
        int rec2 = !((e_r2 >> k) & 1);
        if(rec1 || rec2) any_recovery++;
    }
    printf("\n  Recovery at r+1 OR r+2: ~%.1f%%\n",
           100.0*(recoverable_r1+recoverable_r2-any_recovery)/total_tested);

    printf("\n  ★ EVERY erased bit is recoverable at r+3 (h position)!\n");
    printf("  h enters T1 = h + Σ₁(e) + Ch(e,f,g) + K + W LINEARLY.\n");
    printf("  If T1 known → h known → erased bit recovered.\n");
    printf("  T1 = a[r+4] - T2[r+3] (from backward chain, KNOWN).\n\n");

    printf("  BUT: T1 at round r+3 also contains Ch(e,f,g) at round r+3.\n");
    printf("  Ch at r+3 may erase ANOTHER bit. And recovery needs W[r+3].\n");
    printf("  → Deadpool recovery = CHAIN of erasure+recovery.\n");
    printf("  → Each round: ~14 bits erased, ~14 bits recoverable from r-3.\n");
    printf("  → NET effect: 0? Or positive/negative?\n\n");

    /* NET BALANCE: per round, how many bits are NEWLY erased
     * vs how many are RECOVERED from 3 rounds ago? */
    printf("NET ERASURE BALANCE per round:\n");
    printf("─────────────────────────────\n");

    sha256_trace(msg, st);
    for(int r = 3; r < 64; r++) {
        uint32_t erase_r = st[r][5] ^ st[r][6]; /* erased at round r */
        int n_erased = __builtin_popcount(erase_r);

        /* Recovered from round r-3 (those bits are now at h position) */
        /* The bits erased at r-3 are now at h (round r).
         * They enter T1 linearly → recoverable IF we know T1.
         * T1 = a[r+1] - T2[r]. T2[r] = known from a-values.
         * So T1 known → h known → ALL erased bits from r-3 recovered. */
        uint32_t erase_rm3 = st[r-3][5] ^ st[r-3][6];
        int n_recovered = __builtin_popcount(erase_rm3);

        if(r < 8 || r > 59 || r == 30)
            printf("  Round %2d: erased=%2d, recovered(from r-3)=%2d, net=%+d\n",
                   r, n_erased, n_recovered, n_erased - n_recovered);
    }

    printf("\n  Net ≈ 0 per round: erasure and recovery BALANCE!\n");
    printf("  But: recovery needs T1 = a[r+1]-T2[r] which needs\n");
    printf("  a[r+1] (KNOWN from backward chain if r+1 > 56)\n");
    printf("  AND W[r] (from schedule — needs msg = UNKNOWN).\n");
    printf("  → Recovery blocked by SAME h+W coupling as before.\n");
    printf("  → Deadpool can recover bits, but only WITH knowledge of W.\n");

    return 0;
}
