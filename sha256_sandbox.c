/*
 * SHA-256 SANDBOX: instrumented solver with signal hunting.
 *
 * Fixed: target hash → state[57..64] (backward chain).
 * Variable: W[0..15] (msg).
 * Forward: W[0..15] → rounds 0..56 → state[57].
 * MATCH: state[56] output → does round 57 reach target state[57]?
 *
 * SIGNAL HUNT: at intermediate rounds, compare with what
 * the CORRECT msg produces. Look for ANY metric that
 * correlates with "number of correct W-words".
 *
 * Key idea: if W[0] is correct (but W[1..15] wrong),
 * does ANY intermediate state show a signal?
 *
 * gcc -O3 -march=native -o sandbox sha256_sandbox.c -lm
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

/* Full forward with state trace */
void sha256_trace(const uint32_t msg[16], uint32_t st[65][8], uint32_t W[64]) {
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

int hw256(const uint32_t a[8], const uint32_t b[8]) {
    int h=0; for(int i=0;i<8;i++) h+=__builtin_popcount(a[i]^b[i]); return h;
}

int main() {
    printf("SHA-256 SANDBOX: signal hunt with partial correct words\n");
    printf("=======================================================\n\n");

    srand(42);
    uint32_t msg_correct[16], st_correct[65][8], W_correct[64];
    for(int i=0;i<16;i++) msg_correct[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
    sha256_trace(msg_correct, st_correct, W_correct);

    uint32_t target_hash[8];
    for(int i=0;i<8;i++) target_hash[i]=st_correct[64][i]+IV[i];

    /* EXPERIMENT: fix K correct words, randomize the rest.
     * K=0: all random (baseline).
     * K=1: W[0] correct, W[1..15] random.
     * ...
     * K=16: all correct (= match).
     *
     * For each K: measure HW(state[r] vs correct state[r]) at various rounds.
     * SIGNAL = does having more correct words make intermediate state CLOSER? */

    printf("Signal hunt: K correct words (from W[0]), rest random.\n");
    printf("Measure: avg HW(δstate[r]) vs correct at rounds 1,4,8,16,32,57,64.\n\n");

    int N = 5000; /* trials per K */

    printf("K correct |");
    int checkpoints[] = {1, 2, 4, 8, 16, 32, 48, 57, 64};
    int n_cp = 9;
    for(int c=0;c<n_cp;c++) printf(" r=%-3d|", checkpoints[c]);
    printf(" hash\n");
    printf("──────────+");
    for(int c=0;c<n_cp;c++) printf("──────+");
    printf("──────\n");

    for(int K=0; K<=16; K++) {
        double avg_hw[65] = {0};
        double avg_hash_hw = 0;

        for(int trial=0; trial<N; trial++) {
            uint32_t msg_test[16], st_test[65][8], W_test[64];

            /* First K words correct, rest random */
            for(int i=0;i<K;i++) msg_test[i] = msg_correct[i];
            for(int i=K;i<16;i++) msg_test[i] = (uint32_t)rand()|((uint32_t)rand()<<16);

            sha256_trace(msg_test, st_test, W_test);

            for(int c=0;c<n_cp;c++) {
                int r = checkpoints[c];
                avg_hw[r] += hw256(st_test[r], st_correct[r]);
            }

            uint32_t h[8];
            for(int i=0;i<8;i++) h[i]=st_test[64][i]+IV[i];
            avg_hash_hw += hw256(h, target_hash);
        }

        printf("  K=%-2d    |", K);
        for(int c=0;c<n_cp;c++) {
            int r = checkpoints[c];
            double v = avg_hw[r] / N;
            /* Signal: deviation from 128 (random baseline) */
            if(v < 120)
                printf(" %5.1f★|", v);
            else if(v < 125)
                printf(" %5.1f◆|", v);
            else
                printf(" %5.1f |", v);
        }
        printf(" %5.1f", avg_hash_hw / N);
        if(avg_hash_hw / N < 120) printf("★");
        printf("\n");
    }

    printf("\n★ = significant deviation from 128 (random)\n");
    printf("◆ = slight deviation\n");

    /* DEEP DIVE: for K=1 (only W[0] correct), look at EVERY round */
    printf("\nDEEP: K=1 (only W[0] correct), δstate per round:\n");
    {
        double avg_hw[65] = {0};
        for(int trial=0; trial<N; trial++) {
            uint32_t msg_test[16], st_test[65][8], W_test[64];
            msg_test[0] = msg_correct[0];
            for(int i=1;i<16;i++) msg_test[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
            sha256_trace(msg_test, st_test, W_test);
            for(int r=0;r<=64;r++)
                avg_hw[r] += hw256(st_test[r], st_correct[r]);
        }

        printf("  r: ");
        for(int r=0;r<=64;r++) {
            double v = avg_hw[r]/N;
            if(r<=3 || r==8 || r==16 || r==32 || r==57 || r==64)
                printf("%d:%.0f ", r, v);
        }
        printf("\n");

        /* KEY: at round 1, state should be CLOSER because W[0] correct.
         * state[1] = f(IV, W[0]). Same W[0] → same state[1]!
         * δstate[1] = 0 when K≥1! */
        printf("  r=1: δstate = %.1f (should be 0 if W[0] fully determines state[1])\n",
               avg_hw[1]/N);
    }

    /* K=1: state[1] depends ONLY on W[0] and IV. If W[0] correct → state[1] = correct! */
    /* So for K=1: δstate[1] = 0 EXACTLY. */
    /* For K=2: δstate[2] = 0 EXACTLY. Etc. */
    /* This is TRIVIAL: correct W[0..K-1] → state[K] = correct. */
    /* The question: does state[K] being correct help with state[64]? */

    printf("\n\nKEY TEST: knowing state[K] is correct → does it help with hash?\n");
    printf("──────────────────────────────────────────────────────────\n");
    printf("  If W[0..K-1] correct → state[K] = correct state[K].\n");
    printf("  But state[K] → hash goes through 64-K more rounds.\n");
    printf("  With WRONG W[K..15] → schedule WRONG → hash WRONG.\n");
    printf("  Having correct state[K] doesn't help with wrong W[K..15].\n\n");

    printf("  UNLESS: we can CHECK state[K] against target without\n");
    printf("  computing all remaining rounds. That's the SIGNAL.\n\n");

    /* FINAL TEST: correlate HW(δstate[r]) with HW(δhash) */
    printf("CORRELATION: does closer intermediate state → closer hash?\n");
    printf("──────────────────────────────────────────────────────────\n");
    {
        /* For K=8 (half correct): many trials.
         * For each: measure δstate[16] and δhash. Correlate. */
        int n = 10000;
        double sum_x=0,sum_y=0,sum_xy=0,sum_x2=0,sum_y2=0;
        for(int trial=0;trial<n;trial++){
            uint32_t msg_test[16], st_test[65][8], W_test[64];
            for(int i=0;i<8;i++) msg_test[i]=msg_correct[i];
            for(int i=8;i<16;i++) msg_test[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
            sha256_trace(msg_test, st_test, W_test);

            double x = hw256(st_test[16], st_correct[16]); /* δstate[16] */
            uint32_t h[8];
            for(int i=0;i<8;i++) h[i]=st_test[64][i]+IV[i];
            double y = hw256(h, target_hash); /* δhash */

            sum_x+=x; sum_y+=y; sum_xy+=x*y; sum_x2+=x*x; sum_y2+=y*y;
        }
        double mx=sum_x/n, my=sum_y/n;
        double cov=sum_xy/n-mx*my;
        double sx=sqrt(sum_x2/n-mx*mx), sy=sqrt(sum_y2/n-my*my);
        double corr=cov/(sx*sy+1e-10);

        printf("  K=8: corr(δstate[16], δhash) = %.4f\n", corr);
        printf("  mean δstate[16] = %.1f, mean δhash = %.1f\n", mx, my);
        printf("  → %s\n", fabs(corr) > 0.05 ? "★ CORRELATION EXISTS!" :
               fabs(corr) > 0.02 ? "◆ weak signal" : "no correlation");

        /* Same for different checkpoints */
        for(int cp=1;cp<=57;cp+=7){
            sum_x=sum_y=sum_xy=sum_x2=sum_y2=0;
            for(int trial=0;trial<n;trial++){
                uint32_t msg_test[16], st_test[65][8], W_test[64];
                for(int i=0;i<8;i++) msg_test[i]=msg_correct[i];
                for(int i=8;i<16;i++) msg_test[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
                sha256_trace(msg_test, st_test, W_test);

                double x = hw256(st_test[cp], st_correct[cp]);
                uint32_t h[8];
                for(int i=0;i<8;i++) h[i]=st_test[64][i]+IV[i];
                double y = hw256(h, target_hash);
                sum_x+=x;sum_y+=y;sum_xy+=x*y;sum_x2+=x*x;sum_y2+=y*y;
            }
            mx=sum_x/n;my=sum_y/n;
            cov=sum_xy/n-mx*my;
            sx=sqrt(sum_x2/n-mx*mx);sy=sqrt(sum_y2/n-my*my);
            corr=cov/(sx*sy+1e-10);
            printf("  K=8: corr(δstate[%2d], δhash) = %+.4f %s\n",
                   cp, corr, fabs(corr)>0.05?"★":"");
        }
    }

    return 0;
}
