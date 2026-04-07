/*
 * TIMING FILTER: does wrong a[56] fail FASTER than right one?
 *
 * For each candidate a[56]:
 *   1. Compute h[63] = e[60] = a[60] + a[56] - T2[59]
 *   2. W[63] = C[63] - h[63]
 *   3. Schedule backward: does W[63] lead to consistent W[0..15]?
 *   4. If inconsistent: HOW QUICKLY does contradiction appear?
 *
 * If wrong candidates fail at step 2-3 (O(1)):
 *   total cost = 2^32 × O(1) per candidate = 2^32.
 *
 * gcc -O3 -march=native -o timing timing_filter.c
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

void sha256_full(const uint32_t msg[16], uint32_t st[65][8], uint32_t W[64]) {
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
    for(int i=0;i<8;i++) st[0][i]=IV[i];
    for(int r=0;r<64;r++){
        uint32_t a=st[r][0],b=st[r][1],c=st[r][2],d=st[r][3];
        uint32_t e=st[r][4],f=st[r][5],g=st[r][6],h=st[r][7];
        uint32_t t1=h+S1(e)+CH(e,f,g)+K[r]+W[r], t2=S0(a)+MAJ(a,b,c);
        st[r+1][0]=t1+t2; st[r+1][1]=a; st[r+1][2]=b; st[r+1][3]=c;
        st[r+1][4]=d+t1;  st[r+1][5]=e; st[r+1][6]=f; st[r+1][7]=g;
    }
}

int main() {
    printf("TIMING FILTER: do wrong candidates fail faster?\n");
    printf("================================================\n\n");

    srand(42);
    uint32_t msg[16], st[65][8], W[64];
    for(int i=0;i<16;i++) msg[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
    sha256_full(msg, st, W);

    uint32_t hash[8];
    for(int i=0;i<8;i++) hash[i]=st[64][i]+IV[i];

    /* Known from backward chain */
    uint32_t a_known[65];
    a_known[64]=st[64][0]; a_known[63]=st[64][1];
    a_known[62]=st[64][2]; a_known[61]=st[64][3];
    /* створочне: a[57..60] */
    for(int r=63;r>=60;r--) {
        uint32_t T2=S0(a_known[r])+MAJ(a_known[r],a_known[r-1],a_known[r-2]);
        uint32_t T1=a_known[r+1]-T2;
        uint32_t e_rp1;
        if(r==63) e_rp1=st[64][4]; else if(r==62) e_rp1=st[64][5];
        else if(r==61) e_rp1=st[64][6]; else e_rp1=st[64][7];
        a_known[r-3]=e_rp1-T1;
    }

    uint32_t actual_a56 = st[56][0];
    printf("Actual a[56] = %08x\n\n", actual_a56);

    /* h[63]+W[63] = C (known from state[64]) */
    uint32_t T2_63 = S0(a_known[63])+MAJ(a_known[63],a_known[62],a_known[61]);
    uint32_t T1_63 = a_known[64]-T2_63;
    uint32_t C_63 = T1_63 - S1(st[64][5]) - CH(st[64][5],st[64][6],st[64][7]) - K[63];
    /* C_63 = h[63] + W[63] */

    uint32_t T2_59 = S0(a_known[59])+MAJ(a_known[59],a_known[58],a_known[57]);

    /* For each candidate a[56]:
     * h[63] = a[60] + a56_guess - T2[59]
     * W[63] = C_63 - h[63]
     *
     * CHECK: does W[63] match the schedule?
     * W[63] = s1(W[61]) + W[56] + s0(W[48]) + W[47]
     * We DON'T know W[61],W[56],W[48],W[47] directly.
     *
     * BUT: we can check if our W[63] is CONSISTENT
     * by extending backward and forward and looking for
     * SELF-CONSISTENCY in the message.
     *
     * Simple check: compute msg from schedule backward,
     * then forward SHA-256, compare hash.
     * Cost: O(64) per candidate. Not O(1).
     *
     * ALTERNATIVE: partial consistency check.
     * W[63] from our chain = specific value.
     * W[63] from schedule = s1(W[61])+W[56]+s0(W[48])+W[47].
     * If we also guess a[55] → W[62], etc.
     * After 16 guesses: full msg recoverable → check hash.
     *
     * Can we check EARLIER? */

    /* TEST 1: for candidate a[56], how many backward chain steps
     * are CONSISTENT with the hash?
     *
     * Approach: extend the chain.
     * a[56] guess → h[63] → W[63]
     * Continue: a[55] guess → h[62] → W[62]
     * ...until W[48..63] → schedule backward → W[0..15] → forward → hash check.
     *
     * BUT: each step needs ANOTHER guess (a[55], a[54]...).
     * We can't verify a[56] alone.
     *
     * HOWEVER: for the CORRECT a[56], the chain is INTERNALLY consistent
     * at every step. For WRONG a[56]:
     * W[63] is wrong → schedule backward gives wrong msg →
     * forward SHA gives wrong hash. But we can't check until we have full msg.
     *
     * Is there ANY intermediate check that filters without full msg? */

    /* TEST 2: STRUCTURAL consistency of W[63].
     * W[63] must be a valid schedule output.
     * Are there values of W[63] that are IMPOSSIBLE?
     * No — any 32-bit value can be W[63] for some msg.
     * → No structural filter on W[63] alone. */

    /* TEST 3: CARRY consistency.
     * When computing h[63] from a[56]:
     * h[63] = a[60] + a[56] - T2[59]
     * This uses mod 2^32 arithmetic. The carry in this addition
     * depends on a[60] and a[56].
     * But the carry is ALWAYS consistent (mod arithmetic is always valid).
     * → No carry inconsistency filter. */

    /* TEST 4: MULTI-ROUND consistency.
     * For candidate a[56]:
     * h[63] determined → W[63] = C-h[63]
     * Now: a[57] = KNOWN. T2[56] = S0(a[56]) + MAJ(a[56],a[55],a[54]).
     * a[55],a[54] unknown → T2[56] unknown → can't check a[57] consistency.
     * → Need a[55] too. Can't filter with a[56] alone. */

    /* CONCLUSION: there is NO cheap filter for a[56] alone.
     * Every 32-bit value of a[56] produces a valid h[63] and W[63].
     * No structural, carry, or multi-round consistency check exists
     * that discriminates correct from wrong without more unknowns.
     *
     * The information about "correct vs wrong" is spread across
     * ALL 128 unknown bits (a[53..56]).
     * Checking requires ALL of them → cost 2^128. */

    printf("TEST: filtering power of a[56] alone\n");
    printf("─────────────────────────────────────\n");

    /* For 1000 random a[56] candidates: compute W[63], check
     * if it "looks different" from correct W[63]. */
    int n_test = 1000;
    uint32_t correct_W63 = W[63];

    printf("  Correct W[63] = %08x\n", correct_W63);
    printf("  Checking %d random a[56] candidates...\n", n_test);

    int hw_sum = 0;
    for(int i=0; i<n_test; i++) {
        uint32_t a56g = (uint32_t)rand()|((uint32_t)rand()<<16);
        uint32_t h63 = a_known[60] + a56g - T2_59;
        uint32_t w63 = C_63 - h63;
        int hw = __builtin_popcount(w63 ^ correct_W63);
        hw_sum += hw;
    }
    double avg_hw = (double)hw_sum / n_test;

    printf("  Avg HW(W63_candidate ⊕ W63_correct) = %.1f (expected 16 for random)\n", avg_hw);
    printf("  → %s\n", avg_hw > 15 && avg_hw < 17 ?
           "RANDOM — no filtering possible from W[63] alone" :
           "STRUCTURED — potential filter!");

    /* TEST 5: FORWARD PROPAGATION SPEED.
     * For correct msg: forward SHA computes in 64 steps, all consistent.
     * For wrong msg: also 64 steps, also "consistent" (SHA is a function).
     * There is NO "inconsistency" — SHA always produces SOME output.
     * The question is only: does output = target hash?
     * That check requires FULL computation. No shortcut. */

    printf("\nTEST: does correct msg propagate DIFFERENTLY?\n");
    printf("─────────────────────────────────────────────\n");

    /* Compare: intermediate state HW(δstate) between correct and wrong msg.
     * If correct msg gives "smoother" propagation → timing signal. */

    /* Correct msg */
    uint32_t st_correct[65][8], W_correct[64];
    sha256_full(msg, st_correct, W_correct);

    /* Wrong msg (random) */
    uint32_t msg_wrong[16], st_wrong[65][8], W_wrong[64];
    for(int i=0;i<16;i++) msg_wrong[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
    sha256_full(msg_wrong, st_wrong, W_wrong);

    printf("  Forward propagation profile (HW of state at each round):\n");
    printf("  Round | HW(correct) | HW(wrong) | Diff\n");
    printf("  ------+-------------+-----------+------\n");
    for(int r=0;r<=64;r+=8) {
        int hw_c=0, hw_w=0;
        for(int i=0;i<8;i++) {
            hw_c+=__builtin_popcount(st_correct[r][i]);
            hw_w+=__builtin_popcount(st_wrong[r][i]);
        }
        printf("  %5d | %11d | %9d | %+d\n", r, hw_c, hw_w, hw_c-hw_w);
    }
    printf("  → Both propagate at ~128 HW. No timing difference.\n");

    printf("\n════════════════════════════════════════\n");
    printf("TIMING FILTER VERDICT\n");
    printf("════════════════════════════════════════\n\n");
    printf("  SHA-256 is a TOTAL FUNCTION: every input produces valid output.\n");
    printf("  There is no 'inconsistency' or 'contradiction' for wrong inputs.\n");
    printf("  Wrong msg → wrong hash, but computation is equally fast.\n");
    printf("  No intermediate state distinguishes correct from wrong.\n\n");
    printf("  FILTERING requires comparing output with TARGET hash.\n");
    printf("  That comparison needs FULL computation (64 rounds).\n");
    printf("  Cost per candidate: O(64). Number of candidates: 2^128.\n");
    printf("  Total: O(64 × 2^128) = O(2^134). Standard brute force.\n\n");
    printf("  The ONLY 'timing' difference: correct msg matches hash.\n");
    printf("  That's a RESULT check, not a propagation speed difference.\n");

    return 0;
}
