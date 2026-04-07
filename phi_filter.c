/*
 * Φ-FILTER: look for signal in Φ-space, not hash-space.
 *
 * Φ-coordinates per round:
 *   carry_hw    = HW of carry correction in T1 addition
 *   margin      = min distance of operands to carry threshold
 *   t1_t2_ratio = HW(T1) vs HW(T2)
 *   phi_corr    = HW of total carry correction (φ = Σ carry corrections)
 *
 * Question: for correct msg vs wrong msg, is there signal in Φ-space?
 * We proved: NO signal in hash. But Φ is 4.73× structured (methodology).
 *
 * gcc -O3 -march=native -o phi_filter phi_filter.c -lm
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

typedef struct {
    int carry_hw_total;     /* total HW of carry corrections */
    int t1_hw_total;        /* sum HW(T1) */
    int t2_hw_total;        /* sum HW(T2) */
    int margin_total;       /* sum of carry margins */
    int erasure_total;      /* sum HW(f^g) */
    int ch_carry_interact;  /* carry bits AT erased positions */
    int phi_correction;     /* HW(SHA output ^ XOR-SHA output) */
} PhiCoords;

void sha256_phi(const uint32_t msg[16], PhiCoords *phi, uint32_t hash[8]) {
    memset(phi, 0, sizeof(*phi));
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];

    /* Also compute XOR-schedule for phi_correction */
    uint32_t Wx[64];
    for(int i=0;i<16;i++) Wx[i]=msg[i];
    for(int i=16;i<64;i++) Wx[i]=s1(Wx[i-2])^Wx[i-7]^s0(Wx[i-15])^Wx[i-16];

    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    /* XOR path */
    uint32_t ax=IV[0],bx=IV[1],cx=IV[2],dx=IV[3],ex=IV[4],fx=IV[5],gx=IV[6],hx=IV[7];

    for(int r=0;r<64;r++){
        uint32_t s1e=S1(e), che=CH(e,f,g), s0a=S0(a), mja=MAJ(a,b,c);

        /* T1 additions: h + s1e + che + K[r] + W[r] */
        uint32_t v1=h+s1e;
        uint32_t c1=v1^(h^s1e); /* carry correction */
        uint32_t v2=v1+che;
        uint32_t c2=v2^(v1^che);
        uint32_t v3=v2+K[r];
        uint32_t c3=v3^(v2^K[r]);
        uint32_t t1=v3+W[r];
        uint32_t c4=t1^(v3^W[r]);
        uint32_t t2=s0a+mja;
        uint32_t c5=t2^(s0a^mja);
        uint32_t anew=t1+t2;
        uint32_t c6=anew^(t1^t2);
        uint32_t enew=d+t1;
        uint32_t c7=enew^(d^t1);

        phi->carry_hw_total += __builtin_popcount(c1)+__builtin_popcount(c2)+
            __builtin_popcount(c3)+__builtin_popcount(c4)+
            __builtin_popcount(c5)+__builtin_popcount(c6)+__builtin_popcount(c7);

        phi->t1_hw_total += __builtin_popcount(t1);
        phi->t2_hw_total += __builtin_popcount(t2);

        /* Margin: for each addition, how close were operands to carry threshold? */
        /* Simplified: count bits where both operands = 1 (guaranteed carry) */
        phi->margin_total += __builtin_popcount(h & s1e); /* generate bits in first add */

        /* Erasure */
        uint32_t erase = f ^ g;
        phi->erasure_total += __builtin_popcount(erase);

        /* Ch-carry interaction: carry bits AT erased positions */
        phi->ch_carry_interact += __builtin_popcount(c1 & erase);

        /* XOR path */
        uint32_t s1ex=S1(ex), chex=CH(ex,fx,gx), s0ax=S0(ax), mjax=MAJ(ax,bx,cx);
        uint32_t t1x=hx^s1ex^chex^K[r]^Wx[r];
        uint32_t t2x=s0ax^mjax;

        h=g;g=f;f=e;e=enew;d=c;c=b;b=a;a=anew;
        hx=gx;gx=fx;fx=ex;ex=dx^t1x;dx=cx;cx=bx;bx=ax;ax=t1x^t2x;
    }

    hash[0]=a+IV[0];hash[1]=b+IV[1];hash[2]=c+IV[2];hash[3]=d+IV[3];
    hash[4]=e+IV[4];hash[5]=f+IV[5];hash[6]=g+IV[6];hash[7]=h+IV[7];

    uint32_t hx_out[8]={ax^IV[0],bx^IV[1],cx^IV[2],dx^IV[3],
                         ex^IV[4],fx^IV[5],gx^IV[6],hx^IV[7]};
    phi->phi_correction = 0;
    for(int i=0;i<8;i++) phi->phi_correction += __builtin_popcount(hash[i]^hx_out[i]);
}

int main() {
    printf("Φ-FILTER: signal in Φ-space for correct vs wrong\n");
    printf("=================================================\n\n");

    srand(42);
    uint32_t msg_correct[16];
    for(int i=0;i<16;i++) msg_correct[i]=(uint32_t)rand()|((uint32_t)rand()<<16);

    PhiCoords phi_correct;
    uint32_t target_hash[8];
    sha256_phi(msg_correct, &phi_correct, target_hash);

    /* Collect Φ-coordinates for N wrong messages */
    int N = 50000;
    double sum[7]={0}, sum2[7]={0};
    int correct_vals[7];
    correct_vals[0]=phi_correct.carry_hw_total;
    correct_vals[1]=phi_correct.t1_hw_total;
    correct_vals[2]=phi_correct.t2_hw_total;
    correct_vals[3]=phi_correct.margin_total;
    correct_vals[4]=phi_correct.erasure_total;
    correct_vals[5]=phi_correct.ch_carry_interact;
    correct_vals[6]=phi_correct.phi_correction;

    for(int trial=0;trial<N;trial++){
        uint32_t msg_w[16], h[8];
        for(int i=0;i<16;i++) msg_w[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
        PhiCoords phi_w;
        sha256_phi(msg_w, &phi_w, h);

        int vals[7]={phi_w.carry_hw_total, phi_w.t1_hw_total, phi_w.t2_hw_total,
                     phi_w.margin_total, phi_w.erasure_total,
                     phi_w.ch_carry_interact, phi_w.phi_correction};
        for(int i=0;i<7;i++){sum[i]+=vals[i];sum2[i]+=(double)vals[i]*vals[i];}
    }

    const char *names[7]={"carry_hw","t1_hw","t2_hw","margin","erasure","ch×carry","phi_corr"};

    printf("%-12s | %10s | %10s | %8s | %6s | Signal?\n",
           "Φ-metric","Correct","Mean(wrong)","Std","Z");
    printf("─────────────+────────────+────────────+──────────+────────+────────\n");

    for(int i=0;i<7;i++){
        double mean=sum[i]/N;
        double std=sqrt(sum2[i]/N-mean*mean);
        double z=std>0?(correct_vals[i]-mean)/std:0;
        printf("%-12s | %10d | %10.1f | %8.1f | %+6.2f | %s\n",
               names[i],correct_vals[i],mean,std,z,
               fabs(z)>3?"★★★ YES":fabs(z)>2?"★★ maybe":"no");
    }

    /* Correlation between each Φ-metric and hash distance */
    printf("\nΦ-metric vs hash distance CORRELATION (K=8 correct words):\n");
    printf("──────────────────────────────────────────────────────────\n");

    for(int mi=0;mi<7;mi++){
        double sx=0,sy=0,sxy=0,sx2=0,sy2=0;
        int n=10000;
        for(int trial=0;trial<n;trial++){
            uint32_t msg_t[16],h_t[8];
            for(int i=0;i<8;i++) msg_t[i]=msg_correct[i]; /* 8 correct */
            for(int i=8;i<16;i++) msg_t[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
            PhiCoords phi_t;
            sha256_phi(msg_t, &phi_t, h_t);

            int vals[7]={phi_t.carry_hw_total, phi_t.t1_hw_total, phi_t.t2_hw_total,
                         phi_t.margin_total, phi_t.erasure_total,
                         phi_t.ch_carry_interact, phi_t.phi_correction};

            double x=vals[mi];
            int hd=0; for(int i=0;i<8;i++) hd+=__builtin_popcount(h_t[i]^target_hash[i]);
            double y=hd;

            sx+=x;sy+=y;sxy+=x*y;sx2+=x*x;sy2+=y*y;
        }
        double mx=sx/n,my=sy/n;
        double cov=sxy/n-mx*my;
        double stdx=sqrt(sx2/n-mx*mx),stdy=sqrt(sy2/n-my*my);
        double corr=cov/(stdx*stdy+1e-10);

        printf("  %-12s vs δhash: corr=%+.4f %s\n",
               names[mi],corr,fabs(corr)>0.05?"★ SIGNAL!":fabs(corr)>0.02?"◆ weak":"no");
    }

    return 0;
}
