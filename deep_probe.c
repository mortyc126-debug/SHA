/*
 * deep_probe.c — глубокий зонд для hw59<90
 * Задача: получить достаточно образцов hw59<90 для φ-анализа
 * Вопросы:
 *   Q-I: Откуда lag=44? c[8]↔c[52] самая сильная пара?
 *   Q-J: Сколько образцов нужно для надёжной статистики при hw59<90?
 *   Q-K: Есть ли монотонный переход phi(hw59<k) от k=100 до k=85?
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

/* ===== SHA-256 core ===== */
static const uint32_t K256[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,
    0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,
    0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,
    0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,
    0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,
    0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,
    0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,
    0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,
    0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};
static const uint32_t IV[8]={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
                               0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))
#define S0(x) (ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define S1(x) (ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define s0(x) (ROTR(x,7)^ROTR(x,18)^((x)>>3))
#define s1(x) (ROTR(x,17)^ROTR(x,19)^((x)>>10))
#define CH(e,f,g) (((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))

static int sha256_hw59(const uint32_t W0[16], const uint32_t W1[16]) {
    uint32_t W[2][64];
    for (int i=0;i<16;i++){W[0][i]=W0[i];W[1][i]=W1[i];}
    for (int i=16;i<64;i++){
        W[0][i]=s1(W[0][i-2])+W[0][i-7]+s0(W[0][i-15])+W[0][i-16];
        W[1][i]=s1(W[1][i-2])+W[1][i-7]+s0(W[1][i-15])+W[1][i-16];
    }
    uint32_t s[2][8];
    for(int j=0;j<2;j++) for(int k=0;k<8;k++) s[j][k]=IV[k];
    uint32_t a,b,c,d,e,f,g,h,T1,T2;
    for (int r=0;r<64;r++){
        for(int j=0;j<2;j++){
            a=s[j][0];b=s[j][1];c=s[j][2];d=s[j][3];
            e=s[j][4];f=s[j][5];g=s[j][6];h=s[j][7];
            T1=h+S1(e)+CH(e,f,g)+K256[r]+W[j][r];
            T2=S0(a)+MAJ(a,b,c);
            s[j][7]=g;s[j][6]=f;s[j][5]=e;s[j][4]=d+T1;
            s[j][3]=c;s[j][2]=b;s[j][1]=a;s[j][0]=T1+T2;
        }
    }
    int hw=0;
    for(int k=0;k<8;k++) hw+=__builtin_popcount(s[0][k]^s[1][k]);
    return hw;
}

/* carry + T1+T2 glass box */
static void sha256_carry(const uint32_t W0[16], uint8_t carry[64]) {
    uint32_t W[64];
    for (int i=0;i<16;i++) W[i]=W0[i];
    for (int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for (int r=0;r<64;r++){
        uint32_t T1=h+S1(e)+CH(e,f,g)+K256[r]+W[r];
        uint32_t T2=S0(a)+MAJ(a,b,c);
        uint64_t sum=(uint64_t)T1+T2;
        carry[r]=(sum>>32)&1;
        h=g;g=f;f=e;e=d+T1;d=c;c=b;b=a;a=(uint32_t)sum;
    }
}

/* ===== RNG ===== */
static __thread uint64_t rng_state;
static inline uint64_t xrand(void){
    rng_state^=rng_state<<13;rng_state^=rng_state>>7;rng_state^=rng_state<<17;
    return rng_state;
}

/* ===== Статистика φ ===== */
#define NPAIRS_FOCUS 20   /* фокусные пары для мониторинга */
typedef struct {
    long n;
    double sx[64], sy[64], sxy[64*64];
} PhiMat;

static void phimat_add(PhiMat *m, const uint8_t carry[64]) {
    m->n++;
    for (int i=0;i<64;i++){
        m->sx[i]+=carry[i];
        for(int j=i+1;j<64;j++)
            m->sxy[i*64+j]+=carry[i]*carry[j];
    }
}
static double phimat_phi(const PhiMat *m, int i, int j){
    if(m->n<30) return NAN;
    double N=m->n;
    double mi=m->sx[i]/N, mj=m->sx[j]/N;
    double vi=mi*(1-mi), vj=mj*(1-mj);
    if(vi<1e-9||vj<1e-9) return NAN;
    double mij=m->sxy[i*64+j]/N;
    double cov=mij-mi*mj;
    return cov/sqrt(vi*vj);
}

/* ===== SA для hw59<90 ===== */
#define DW0 0x80000000u

typedef struct {
    int tid; long target; uint64_t seed;
    PhiMat mat85, mat90, mat95;
    long cnt85, cnt90, cnt95;
} WArg;

static void *worker(void *varg){
    WArg *wa=(WArg*)varg;
    rng_state=wa->seed^(uint64_t)wa->tid*0xdeadbeef;

    uint32_t W[16], W1[16];
    for(int i=0;i<16;i++) W[i]=(uint32_t)xrand();
    W1[0]=W[0]^DW0;
    for(int i=1;i<16;i++) W1[i]=W[i];

    double T=50.0;
    uint32_t bestW[16]; int besthw=256;
    memcpy(bestW,W,64);

    long iters=0, max_iters=wa->target*20000LL;
    long cnt85=0,cnt90=0,cnt95=0;

    while(iters<max_iters && cnt85<wa->target){
        iters++;

        W1[0]=W[0]^DW0;
        int hw=sha256_hw59(W,W1);

        uint8_t carry[64];
        sha256_carry(W,carry);

        if(hw<95){
            phimat_add(&wa->mat95,carry);
            cnt95++;
        }
        if(hw<90){
            phimat_add(&wa->mat90,carry);
            cnt90++;
        }
        if(hw<85){
            phimat_add(&wa->mat85,carry);
            cnt85++;
        }

        if(hw<besthw){besthw=hw;memcpy(bestW,W,64);}

        /* мутация */
        int k=(int)(xrand()%16);
        int bit=(int)(xrand()%32);
        uint32_t saved=W[k];
        W[k]^=(1u<<bit);
        W1[0]=W[0]^DW0;
        int nhw=sha256_hw59(W,W1);

        double dE=(nhw-hw);
        double prob=(dE<0)?1.0:exp(-dE/T);
        if((double)(xrand()%1000000)/1000000.0 > prob) W[k]=saved;

        /* Адаптивная температура */
        T*=0.999997;
        if(T<0.5) T=0.5;

        /* Периодический рестарт к лучшему */
        if(iters%200000==0 && besthw>90){
            memcpy(W,bestW,64);
            T=30.0+(double)(xrand()%20);
        }
        /* Большой рестарт при зависании */
        if(iters%2000000==0 && besthw>95){
            for(int i=0;i<16;i++) W[i]=(uint32_t)xrand();
            besthw=256; T=50.0;
        }
    }

    wa->cnt85=cnt85; wa->cnt90=cnt90; wa->cnt95=cnt95;
    return NULL;
}

/* ===== Анализ пар ===== */
static void analyze_mat(const PhiMat *m, const char *label, long N, int topN){
    printf("\n=== %s (N=%ld) ===\n",label,N);
    if(N<30){printf("  Недостаточно образцов.\n");return;}

    /* Топ пар */
    typedef struct{int i,j;double phi;} Pair;
    Pair pairs[64*64/2];
    int np=0;
    for(int i=0;i<64;i++) for(int j=i+1;j<64;j++){
        double phi=phimat_phi(m,i,j);
        if(!isnan(phi)){pairs[np].i=i;pairs[np].j=j;pairs[np].phi=phi;np++;}
    }
    /* Sort by |phi| desc */
    for(int a=0;a<np-1;a++) for(int b=a+1;b<np;b++){
        if(fabs(pairs[b].phi)>fabs(pairs[a].phi)){Pair tmp=pairs[a];pairs[a]=pairs[b];pairs[b]=tmp;}
    }

    printf("  Топ-%d пар:\n",topN<np?topN:np);
    for(int k=0;k<topN&&k<np;k++){
        int i=pairs[k].i,j=pairs[k].j;
        printf("    c[%2d]↔c[%2d] lag=%2d: phi=%+.4f\n",
               i,j,j-i,pairs[k].phi);
    }

    /* Лаговая структура */
    printf("  Лаговая структура (avg|phi|, max|phi|):\n");
    for(int lag=1;lag<=20;lag++){
        double sum=0,mx=0; int cnt=0;
        for(int i=0;i+lag<64;i++){
            double phi=phimat_phi(m,i,i+lag);
            if(!isnan(phi)){sum+=fabs(phi);if(fabs(phi)>mx)mx=fabs(phi);cnt++;}
        }
        if(cnt>0) printf("    lag=%2d: avg|phi|=%.4f max|phi|=%.4f\n",
                         lag,sum/cnt,mx);
    }

    /* Зонный анализ */
    double E1=0,E2=0,E3=0; int n1=0,n2=0,n3=0;
    for(int i=0;i<64;i++) for(int j=i+1;j<64;j++){
        double phi=phimat_phi(m,i,j);
        if(isnan(phi)) continue;
        int zi=(i<18)?0:(i<59)?1:2;
        int zj=(j<18)?0:(j<59)?1:2;
        /* пара в зоне = обе координаты */
        if(zi==0&&zj==0){E1+=fabs(phi);n1++;}
        else if(zi==1&&zj==1){E2+=fabs(phi);n2++;}
        else if(zi==2&&zj==2){E3+=fabs(phi);n3++;}
    }
    printf("  Зонная энергия: Z1=%.4f(%d) Z2=%.4f(%d) Z3=%.4f(%d)\n",
           n1?E1/n1:0,n1, n2?E2/n2:0,n2, n3?E3/n3:0,n3);
}

int main(int argc, char **argv){
    int nthreads=4;
    long target_per_thread=500;
    if(argc>1) target_per_thread=atol(argv[1]);
    if(argc>2) nthreads=atoi(argv[2]);

    printf("deep_probe: target=%ld образцов hw59<85 на поток, потоков=%d\n",
           target_per_thread,nthreads);
    printf("\nВопросы:\n");
    printf("  Q-I: Откуда lag=44? Почему c[8]↔c[52] самая сильная пара?\n");
    printf("  Q-J: Достаточно ли N образцов для phi при hw59<90?\n");
    printf("  Q-K: Монотонный переход phi(hw59<k) от k=100 до k=85?\n\n");

    time_t t0=time(NULL);

    WArg *args=calloc(nthreads,sizeof(WArg));
    pthread_t *threads=malloc(nthreads*sizeof(pthread_t));

    for(int t=0;t<nthreads;t++){
        args[t].tid=t;
        args[t].target=target_per_thread;
        args[t].seed=(uint64_t)time(NULL)^(uint64_t)t*0x123456789abcdefLL;
        pthread_create(&threads[t],NULL,worker,&args[t]);
    }
    for(int t=0;t<nthreads;t++) pthread_join(threads[t],NULL);

    /* Объединяем статистику */
    PhiMat mat85,mat90,mat95;
    memset(&mat85,0,sizeof(mat85));
    memset(&mat90,0,sizeof(mat90));
    memset(&mat95,0,sizeof(mat95));
    long tot85=0,tot90=0,tot95=0;
    for(int t=0;t<nthreads;t++){
        tot85+=args[t].cnt85; tot90+=args[t].cnt90; tot95+=args[t].cnt95;
        mat85.n+=args[t].mat85.n; mat90.n+=args[t].mat90.n; mat95.n+=args[t].mat95.n;
        for(int i=0;i<64;i++){
            mat85.sx[i]+=args[t].mat85.sx[i];
            mat90.sx[i]+=args[t].mat90.sx[i];
            mat95.sx[i]+=args[t].mat95.sx[i];
            for(int j=i+1;j<64;j++){
                mat85.sxy[i*64+j]+=args[t].mat85.sxy[i*64+j];
                mat90.sxy[i*64+j]+=args[t].mat90.sxy[i*64+j];
                mat95.sxy[i*64+j]+=args[t].mat95.sxy[i*64+j];
            }
        }
    }

    printf("Время: %ld сек\n",(long)(time(NULL)-t0));
    printf("Образцы: hw59<85=%ld hw59<90=%ld hw59<95=%ld\n",tot85,tot90,tot95);

    /* Анализ */
    analyze_mat(&mat95,"hw59<95",tot95,20);
    analyze_mat(&mat90,"hw59<90",tot90,20);
    analyze_mat(&mat85,"hw59<85",tot85,15);

    /* Q-I: Специальный анализ c[8]↔c[52] */
    printf("\n=== Q-I: Анализ c[8]↔c[52] lag=44 ===\n");
    printf("Гипотеза: lag=44 = 32+12? или 7×6+2? Нет стандартной структуры.\n");
    printf("Проверим окрестность: корреляции c[8] со всеми c[j], j=40..63:\n");
    for(int j=40;j<64;j++){
        if(j<=8) continue;
        double phi95=phimat_phi(&mat95,8,j);
        double phi90=phimat_phi(&mat90,8,j);
        if(!isnan(phi95)){
            printf("  c[8]↔c[%2d] lag=%2d: phi95=%+.4f phi90=%s\n",
                   j,j-8,phi95,isnan(phi90)?"-":
                   (char[]){phi90>=0?'+':'-',
                   '0'+(int)(fabs(phi90)*10)%10,'.',
                   '0'+(int)(fabs(phi90)*1000)%10,0});
        }
    }

    free(args); free(threads);
    return 0;
}
