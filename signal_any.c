/*
 * ANY SIGNAL: measure EVERYTHING about correct vs wrong msg.
 *
 * For correct msg and 1000 wrong msgs:
 * measure every possible internal metric.
 * If ANY metric distinguishes correct → signal exists.
 *
 * gcc -O3 -march=native -o signal_any signal_any.c -lm
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
    int carry_total;      /* total carry bits across all additions */
    int state_hw_sum;     /* sum of HW(state) across all rounds */
    int ch_hw_sum;        /* sum of HW(Ch output) */
    int maj_hw_sum;       /* sum of HW(Maj output) */
    int t1_hw_sum;        /* sum of HW(T1) */
    int t2_hw_sum;        /* sum of HW(T2) */
    int carry_max_chain;  /* longest consecutive carry=1 across all rounds */
    int state_hw_var;     /* variance proxy of state HW */
    int w_hw_sum;         /* sum of HW(W[r]) */
    int erasure_total;    /* total erased bits (f⊕g) */
    uint32_t hash[8];
} Metrics;

void sha256_metrics(const uint32_t msg[16], Metrics *m) {
    memset(m, 0, sizeof(*m));
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];

    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];

    for(int r=0;r<64;r++){
        /* Metrics */
        int state_hw = __builtin_popcount(a)+__builtin_popcount(b)+
                        __builtin_popcount(c)+__builtin_popcount(d)+
                        __builtin_popcount(e)+__builtin_popcount(f)+
                        __builtin_popcount(g)+__builtin_popcount(h);
        m->state_hw_sum += state_hw;
        m->state_hw_var += (state_hw - 128)*(state_hw - 128);

        uint32_t ch_val = CH(e,f,g);
        uint32_t maj_val = MAJ(a,b,c);
        m->ch_hw_sum += __builtin_popcount(ch_val);
        m->maj_hw_sum += __builtin_popcount(maj_val);

        /* Erasure */
        m->erasure_total += __builtin_popcount(f ^ g);

        /* W */
        m->w_hw_sum += __builtin_popcount(W[r]);

        /* Carry in T1 = h + S1(e) + ch + K[r] + W[r] */
        uint32_t s1e = S1(e);
        uint32_t v1 = h + s1e;
        m->carry_total += __builtin_popcount(v1 ^ (h ^ s1e));
        uint32_t v2 = v1 + ch_val;
        m->carry_total += __builtin_popcount(v2 ^ (v1 ^ ch_val));
        uint32_t v3 = v2 + K[r];
        m->carry_total += __builtin_popcount(v3 ^ (v2 ^ K[r]));
        uint32_t t1 = v3 + W[r];
        m->carry_total += __builtin_popcount(t1 ^ (v3 ^ W[r]));

        /* T2 */
        uint32_t s0a = S0(a);
        uint32_t t2 = s0a + maj_val;
        m->carry_total += __builtin_popcount(t2 ^ (s0a ^ maj_val));

        /* a_new, e_new */
        uint32_t a_new = t1 + t2;
        m->carry_total += __builtin_popcount(a_new ^ (t1 ^ t2));
        uint32_t e_new = d + t1;
        m->carry_total += __builtin_popcount(e_new ^ (d ^ t1));

        m->t1_hw_sum += __builtin_popcount(t1);
        m->t2_hw_sum += __builtin_popcount(t2);

        /* Carry max chain */
        uint32_t carry_bits = v1 ^ (h ^ s1e); /* carry of first addition */
        int chain = 0, max_chain = 0;
        for(int k=0;k<32;k++){
            if((carry_bits>>k)&1) { chain++; if(chain>max_chain)max_chain=chain; }
            else chain=0;
        }
        if(max_chain > m->carry_max_chain) m->carry_max_chain = max_chain;

        h=g;g=f;f=e;e=e_new;d=c;c=b;b=a;a=a_new;
    }

    m->hash[0]=a+IV[0];m->hash[1]=b+IV[1];m->hash[2]=c+IV[2];m->hash[3]=d+IV[3];
    m->hash[4]=e+IV[4];m->hash[5]=f+IV[5];m->hash[6]=g+IV[6];m->hash[7]=h+IV[7];
}

int main() {
    printf("ANY SIGNAL: correct vs wrong, ALL metrics\n");
    printf("==========================================\n\n");

    srand(42);
    uint32_t msg_correct[16];
    for(int i=0;i<16;i++) msg_correct[i]=(uint32_t)rand()|((uint32_t)rand()<<16);

    Metrics m_correct;
    sha256_metrics(msg_correct, &m_correct);

    /* Collect metrics for N wrong messages */
    int N = 10000;
    double sum_carry=0, sum2_carry=0;
    double sum_shw=0, sum2_shw=0;
    double sum_ch=0, sum2_ch=0;
    double sum_maj=0, sum2_maj=0;
    double sum_t1=0, sum2_t1=0;
    double sum_t2=0, sum2_t2=0;
    double sum_maxc=0, sum2_maxc=0;
    double sum_var=0, sum2_var=0;
    double sum_whw=0, sum2_whw=0;
    double sum_erase=0, sum2_erase=0;

    for(int trial=0; trial<N; trial++) {
        uint32_t msg_w[16];
        for(int i=0;i<16;i++) msg_w[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
        Metrics m;
        sha256_metrics(msg_w, &m);

        sum_carry+=m.carry_total; sum2_carry+=(double)m.carry_total*m.carry_total;
        sum_shw+=m.state_hw_sum; sum2_shw+=(double)m.state_hw_sum*m.state_hw_sum;
        sum_ch+=m.ch_hw_sum; sum2_ch+=(double)m.ch_hw_sum*m.ch_hw_sum;
        sum_maj+=m.maj_hw_sum; sum2_maj+=(double)m.maj_hw_sum*m.maj_hw_sum;
        sum_t1+=m.t1_hw_sum; sum2_t1+=(double)m.t1_hw_sum*m.t1_hw_sum;
        sum_t2+=m.t2_hw_sum; sum2_t2+=(double)m.t2_hw_sum*m.t2_hw_sum;
        sum_maxc+=m.carry_max_chain; sum2_maxc+=(double)m.carry_max_chain*m.carry_max_chain;
        sum_var+=m.state_hw_var; sum2_var+=(double)m.state_hw_var*m.state_hw_var;
        sum_whw+=m.w_hw_sum; sum2_whw+=(double)m.w_hw_sum*m.w_hw_sum;
        sum_erase+=m.erasure_total; sum2_erase+=(double)m.erasure_total*m.erasure_total;
    }

    printf("%-20s | %10s | %10s | %8s | %6s | Signal?\n",
           "Metric", "Correct", "Mean(wrong)", "Std", "Z-score");
    printf("─────────────────────+────────────+────────────+──────────+────────+────────\n");

    #define REPORT(name, val_c, sum_w, sum2_w) { \
        double mean = sum_w / N; \
        double std = sqrt(sum2_w/N - mean*mean); \
        double z = std > 0 ? (val_c - mean) / std : 0; \
        printf("%-20s | %10d | %10.1f | %8.1f | %+6.2f | %s\n", \
               name, val_c, mean, std, z, \
               fabs(z) > 3.0 ? "★★★ YES" : fabs(z) > 2.0 ? "★★ maybe" : "no"); \
    }

    REPORT("Carry total",     m_correct.carry_total,     sum_carry, sum2_carry);
    REPORT("State HW sum",    m_correct.state_hw_sum,    sum_shw,   sum2_shw);
    REPORT("Ch HW sum",       m_correct.ch_hw_sum,       sum_ch,    sum2_ch);
    REPORT("Maj HW sum",      m_correct.maj_hw_sum,      sum_maj,   sum2_maj);
    REPORT("T1 HW sum",       m_correct.t1_hw_sum,       sum_t1,    sum2_t1);
    REPORT("T2 HW sum",       m_correct.t2_hw_sum,       sum_t2,    sum2_t2);
    REPORT("Max carry chain",  m_correct.carry_max_chain, sum_maxc,  sum2_maxc);
    REPORT("State HW variance",m_correct.state_hw_var,    sum_var,   sum2_var);
    REPORT("W HW sum",        m_correct.w_hw_sum,        sum_whw,   sum2_whw);
    REPORT("Erasure total",   m_correct.erasure_total,    sum_erase, sum2_erase);

    printf("\n");
    printf("ANY metric with |Z| > 3 would indicate the correct message\n");
    printf("has a STATISTICALLY DISTINGUISHABLE internal computation.\n");

    return 0;
}
