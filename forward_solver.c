/*
 * FORWARD SOLVER: choose W[r] to steer state toward target.
 *
 * Rounds 0..15: W[r] = free (msg word). Can choose to steer state.
 * Rounds 16..63: W[r] = schedule(W[0..15]). Fixed once msg chosen.
 *
 * Strategy: at each round r (0..15), choose W[r] so that
 * a[r+1] or e[r+1] matches some target value.
 *
 * Key insight: T1 = h + S1(e) + Ch(e,f,g) + K[r] + W[r]
 * W[r] enters T1 LINEARLY (addition). So:
 * For ANY desired T1_target: W[r] = T1_target - h - S1(e) - Ch - K[r]
 * UNIQUE solution! W[r] is DETERMINED by desired T1.
 *
 * And: a[r+1] = T1 + T2. T2 = S0(a) + Maj(a,b,c) = known from state.
 * So: for desired a[r+1] = target: T1 = target - T2.
 * → W[r] = (target - T2) - h - S1(e) - Ch - K[r]
 *
 * We can set a[r+1] to ANY value by choosing W[r]!
 * BUT: e[r+1] = d + T1 = d + (target - T2) = also determined.
 * We can't independently control BOTH a and e.
 *
 * Question: with 16 free W-words, can we steer 256-bit state[16]
 * to match what backward chain needs?
 *
 * gcc -O3 -march=native -o fwd_solve forward_solver.c -lm
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

int main() {
    printf("FORWARD SOLVER: steer state with W[r]\n");
    printf("======================================\n\n");

    /* Generate target hash */
    srand(42);
    uint32_t msg_orig[16];
    for(int i=0;i<16;i++) msg_orig[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
    uint32_t target_hash[8];
    sha256_compute(msg_orig, target_hash);

    printf("Target hash: ");
    for(int i=0;i<8;i++) printf("%08x ", target_hash[i]);
    printf("\n\n");

    /* TEST 1: Can we choose W[r] to set a[r+1] to ANY target? */
    printf("TEST 1: Set a[1] to target value via W[0]\n");
    printf("──────────────────────────────────────────\n");
    {
        uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3];
        uint32_t e=IV[4],f=IV[5],g=IV[6],h=IV[7];

        uint32_t target_a1 = 0xDEADBEEF;

        /* T2 = S0(a) + Maj(a,b,c) — known from state[0] */
        uint32_t t2 = S0(a) + MAJ(a,b,c);
        /* T1 = target_a1 - T2 */
        uint32_t t1_needed = target_a1 - t2;
        /* W[0] = T1 - h - S1(e) - Ch(e,f,g) - K[0] */
        uint32_t w0 = t1_needed - h - S1(e) - CH(e,f,g) - K[0];

        /* Verify */
        uint32_t t1_check = h + S1(e) + CH(e,f,g) + K[0] + w0;
        uint32_t a1_check = t1_check + t2;
        uint32_t e1_check = d + t1_check;

        printf("  Target a[1] = %08x\n", target_a1);
        printf("  Computed W[0] = %08x\n", w0);
        printf("  Resulting a[1] = %08x %s\n", a1_check,
               a1_check==target_a1 ? "✓ MATCH!" : "✗");
        printf("  Resulting e[1] = %08x (NOT controlled — determined by T1)\n", e1_check);
        printf("  → We control a[r+1] but NOT e[r+1] independently.\n");
    }

    /* TEST 2: Over 16 rounds, how much of state[16] can we control? */
    printf("\nTEST 2: Controllable dimensions over 16 rounds\n");
    printf("───────────────────────────────────────────────\n");
    {
        /* Each W[r] gives us control over a[r+1] (32 bits).
         * But e[r+1] = d + T1 is determined (not independently controlled).
         *
         * State[16] = (a16, a15, a14, a13, e16, e15, e14, e13).
         * a13..a16 = last 4 values of a-sequence. We set a[1]..a[16]
         * by choosing W[0]..W[15]. So a[13..16] = fully controlled!
         *
         * e13..e16: e[r] = d[r-1] + T1[r-1] = a[r-4] + T1[r-1].
         * T1[r-1] = a[r] - T2[r-1]. a[r] = our choice!
         * So T1[r-1] = chosen_a[r] - T2[r-1] = determined by our choices.
         * e[r] = a[r-4] + T1[r-1] = determined.
         *
         * Actually e[r+1] = d[r] + T1[r]. d[r] = a[r-3] (chosen!).
         * T1[r] = a[r+1] - T2[r] = chosen - T2 = determined.
         * So e[r+1] = a[r-3] + a[r+1] - T2[r].
         *
         * e[r+1] depends on a[r-3] and a[r+1] — BOTH chosen by us!
         * → e[r+1] IS controllable through our choice of a-values!
         *
         * BUT: T2[r] = S0(a[r]) + Maj(a[r], a[r-1], a[r-2]) — nonlinear in a!
         * We choose a[r+1] = target. This determines T1 = target - T2.
         * Then e[r+1] = d + T1 = a[r-3] + target - T2.
         * T2 depends on a[r], a[r-1], a[r-2] which we ALREADY chose.
         * So e[r+1] is a SPECIFIC function of our previous choices + current target.
         *
         * → e is NOT independently controllable. It's a CONSEQUENCE of a-choices.
         */

        /* Concrete: set a[1..16] to targets. What e[13..16] do we get? */
        uint32_t target_a[17]; /* a[0]=IV, a[1..16]=targets */
        target_a[0] = IV[0];
        for(int r=1;r<=16;r++) target_a[r] = 0x11111111 * r; /* arbitrary targets */

        uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3];
        uint32_t e_val=IV[4],f_val=IV[5],g_val=IV[6],h_val=IV[7];
        uint32_t W_chosen[16];

        for(int r=0;r<16;r++) {
            uint32_t t2 = S0(a) + MAJ(a,b,c);
            uint32_t t1 = target_a[r+1] - t2;
            W_chosen[r] = t1 - h_val - S1(e_val) - CH(e_val,f_val,g_val) - K[r];

            uint32_t e_new = d + t1;
            h_val=g_val;g_val=f_val;f_val=e_val;e_val=e_new;
            d=c;c=b;b=a;a=target_a[r+1];
        }

        printf("  Set a[1..16] to arbitrary targets:\n");
        printf("  a[13..16] = %08x %08x %08x %08x (CONTROLLED)\n",
               target_a[13], target_a[14], target_a[15], target_a[16]);
        printf("  e[13..16] = %08x %08x %08x %08x (DETERMINED by a-choices)\n",
               /* These are the actual e values from computation */
               f_val, g_val, h_val, e_val); /* current state shifted */

        /* Verify: compute forward with chosen W, check state[16] */
        uint32_t hash_check[8];
        sha256_compute(W_chosen, hash_check);
        /* state[16] from full computation */
        uint32_t W_full[64], st[65][8];
        for(int i=0;i<16;i++) W_full[i]=W_chosen[i];
        for(int i=16;i<64;i++) W_full[i]=s1(W_full[i-2])+W_full[i-7]+s0(W_full[i-15])+W_full[i-16];
        for(int i=0;i<8;i++) st[0][i]=IV[i];
        for(int r=0;r<64;r++){
            uint32_t aa=st[r][0],bb=st[r][1],cc=st[r][2],dd=st[r][3];
            uint32_t ee=st[r][4],ff=st[r][5],gg=st[r][6],hh=st[r][7];
            uint32_t t1=hh+S1(ee)+CH(ee,ff,gg)+K[r]+W_full[r];
            uint32_t t2=S0(aa)+MAJ(aa,bb,cc);
            st[r+1][0]=t1+t2;st[r+1][1]=aa;st[r+1][2]=bb;st[r+1][3]=cc;
            st[r+1][4]=dd+t1;st[r+1][5]=ee;st[r+1][6]=ff;st[r+1][7]=gg;
        }

        printf("\n  Actual state[16] = ");
        for(int i=0;i<8;i++) printf("%08x ", st[16][i]);
        printf("\n  a[13..16] match: %s\n",
               st[16][3]==target_a[13] && st[16][2]==target_a[14] &&
               st[16][1]==target_a[15] && st[16][0]==target_a[16] ? "✓ YES" : "✗ NO");

        /* What hash does this produce? */
        printf("\n  Hash from chosen W: ");
        for(int i=0;i<8;i++) printf("%08x ", hash_check[i]);
        printf("\n  Target hash:        ");
        for(int i=0;i<8;i++) printf("%08x ", target_hash[i]);
        int hw = 0;
        for(int i=0;i<8;i++) hw += __builtin_popcount(hash_check[i]^target_hash[i]);
        printf("\n  Hash diff: %d/256 bits\n", hw);
    }

    /* TEST 3: MEET IN THE MIDDLE concept */
    printf("\nTEST 3: Forward steering + backward target\n");
    printf("────────────────────────────────────────────\n");
    {
        /* Backward: from target hash, we know a[57..64].
         * Forward: we can set a[1..16] to anything by choosing W[0..15].
         * GAP: rounds 17..56. State at round 16 → forward with schedule → state[57].
         *
         * If state[16] is chosen → schedule determined → state[57] determined.
         * We need state[57] to match backward target.
         *
         * But: state[16] = 256 bits. We control a[13..16] = 128 bits (4 words).
         * e[13..16] = determined (not free). So we have 128 bits of control.
         *
         * Target at state[57]: a[57..60] = known from backward (128 bits).
         * Match: 128 control bits → 128 target bits. SQUARE SYSTEM!
         *
         * But: the mapping state[16] → state[57] (41 rounds with schedule) is
         * highly nonlinear. The 128 control bits affect ALL 256 bits of state[57].
         * And state[57] has 256 bits but we target only 128 (a[57..60]).
         *
         * Birthday: need 2^64 forward trials to match 128-bit target.
         * BUT: each trial is O(41 rounds) = O(41).
         * Total: 2^64 × O(41) ≈ 2^70.
         *
         * WAIT: is this sub-birthday for PREIMAGE?
         * Standard preimage: 2^256. Our approach: 2^70? That's wrong.
         *
         * The issue: we control only 128 bits (a[13..16]).
         * But msg = W[0..15] = 512 bits. 128 bits of state = covers 128 msg bits.
         * Other 384 msg bits = also chosen by us (they determine a[1..12]).
         * So we actually control ALL 512 msg bits.
         *
         * The question: can we choose ALL 512 msg bits to make state[57]
         * match backward target? If yes → preimage found.
         *
         * Degrees of freedom: 512 (msg bits).
         * Constraints: state[57..64] = 256 bits (from hash).
         * Free: 512 - 256 = 256. Standard preimage.
         */

        printf("  Forward control: 16 W-words = 512 msg bits.\n");
        printf("  These determine a[1..16] (512 bits → 512 a-values).\n");
        printf("  a[1..16] → schedule → forward 41 rounds → state[57].\n");
        printf("  Target: state[57..64] from backward = 256 bits.\n");
        printf("  512 DOF - 256 constraints = 256 free.\n");
        printf("  → STANDARD PREIMAGE (2^256). No improvement.\n\n");

        printf("  For COLLISION: two different W giving same hash.\n");
        printf("  Both go through same process. Birthday = 2^128.\n\n");

        printf("  The forward solver can SET state[16] to anything,\n");
        printf("  but state[16]→state[64] is FIXED by schedule.\n");
        printf("  No way to adjust after round 16.\n");
    }

    printf("\n════════════════════════════════════\n");
    printf("FORWARD SOLVER VERDICT\n");
    printf("════════════════════════════════════\n\n");
    printf("  W[r] enters T1 LINEARLY → can set a[r+1] to any target.\n");
    printf("  But e[r+1] = determined by a-choices (not independent).\n");
    printf("  16 free W-words → 512 msg bits → full control.\n");
    printf("  BUT: after round 16, schedule fixes remaining W.\n");
    printf("  state[16]→state[64] = deterministic, no further control.\n");
    printf("  Finding msg with correct hash = standard preimage/birthday.\n");

    return 0;
}
