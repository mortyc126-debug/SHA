/*
 * anomaly_probe.c — анализ аномалии c[43]↔c[44] (положительная lag=1)
 * и получение hw59<85 образцов
 *
 * Q182: Почему c[43]↔c[44] lag=1 ПОЛОЖИТЕЛЬНАЯ?
 *       Все другие lag=1 пары — отрицательные.
 *
 * Q181: Измерить φ при hw59<85
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

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
    for (int r=0;r<64;r++){
        for(int j=0;j<2;j++){
            uint32_t a=s[j][0],b=s[j][1],c=s[j][2],d=s[j][3];
            uint32_t e=s[j][4],f=s[j][5],g=s[j][6],h=s[j][7];
            uint32_t T1=h+S1(e)+CH(e,f,g)+K256[r]+W[j][r];
            uint32_t T2=S0(a)+MAJ(a,b,c);
            s[j][7]=g;s[j][6]=f;s[j][5]=e;s[j][4]=d+T1;
            s[j][3]=c;s[j][2]=b;s[j][1]=a;s[j][0]=T1+T2;
        }
    }
    int hw=0;
    for(int k=0;k<8;k++) hw+=__builtin_popcount(s[0][k]^s[1][k]);
    return hw;
}

/* Glass box SHA с T1,T2 детальным выводом */
typedef struct {
    uint32_t T1[64], T2[64], A[64];
    uint8_t  carry[64];
    uint64_t sumT[64];  /* T1+T2 полная */
} GlassBox;

static void sha256_glass(const uint32_t W0[16], GlassBox *g) {
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=W0[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];

    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3];
    uint32_t e=IV[4],f=IV[5],gg=IV[6],h=IV[7];
    for(int r=0;r<64;r++){
        g->T1[r]=h+S1(e)+CH(e,f,gg)+K256[r]+W[r];
        g->T2[r]=S0(a)+MAJ(a,b,c);
        uint64_t sum=(uint64_t)g->T1[r]+g->T2[r];
        g->sumT[r]=sum;
        g->carry[r]=(sum>>32)&1;
        g->A[r]=(uint32_t)sum;
        h=gg;gg=f;f=e;e=d+g->T1[r];d=c;c=b;b=a;a=g->A[r];
    }
}

/* Ключевой анализ: почему carry[43]=1 → carry[44]=1 (положительная) */
static void analyze_lag1_signs(void) {
    printf("=== Q182: Анализ lag=1 пар — знаки φ ===\n");
    printf("Теория: carry[r]=1 → a[r+1] мало → T2[r+1] мало → carry[r+1]=0 (анти)\n");
    printf("НО: carry[43]=1 → carry[44]=1 (?)\n\n");

    /* Сгенерируем 100 образцов hw59<95 и изучим распределение */
    printf("Изучаем P(carry[44]=1 | carry[43]=1) vs P(carry[44]=1 | carry[43]=0)\n");
    printf("для случайных W:\n\n");

    /* Случайный W */
    srand(42);
    long cnt[4]={0}; /* [c43][c44] */
    for(int s=0;s<1000000;s++){
        uint32_t W[16];
        for(int i=0;i<16;i++) W[i]=(uint32_t)(rand()^(rand()<<16));
        GlassBox g;
        sha256_glass(W,&g);
        cnt[g.carry[43]*2+g.carry[44]]++;
    }
    long tot=cnt[0]+cnt[1]+cnt[2]+cnt[3];
    printf("Случайный W (1M образцов):\n");
    printf("  P(c44=1|c43=0) = %.4f (%ld/%ld)\n",
           (double)cnt[1]/(cnt[0]+cnt[1]),cnt[1],cnt[0]+cnt[1]);
    printf("  P(c44=1|c43=1) = %.4f (%ld/%ld)\n",
           (double)cnt[3]/(cnt[2]+cnt[3]),cnt[3],cnt[2]+cnt[3]);
    printf("  phi(c43,c44) theoretical = %.6f\n",
           ((double)cnt[3]/tot - (double)(cnt[2]+cnt[3])/tot * (double)(cnt[1]+cnt[3])/tot) /
           sqrt((double)(cnt[2]+cnt[3])/tot*(double)(cnt[0]+cnt[1])/tot*
                (double)(cnt[1]+cnt[3])/tot*(double)(cnt[0]+cnt[2])/tot));

    /* Проверим: при каком именно условии на состояние (e,a) в раунде 43 
     * возникает позитивная корреляция? */
    printf("\nАнализ T1[43]+T2[43] при carry[43]=1:\n");
    printf("Если carry[43]=1: a[44] = (T1[43]+T2[43]) mod 2^32\n");
    printf("                   т.е. a[44] = (большое число) mod 2^32 — СЛУЧАЙНЫЙ 32-bit\n");
    printf("                   T2[44] = S0(a[44])+MAJ(a[44],a[43],a[42])\n");
    printf("                   Зависит от конкретного a[44], а не от его 'малости'\n\n");
    
    printf("ВАЖНО: Причина анти-корреляции carry[r]→carry[r+1]:\n");
    printf("  a[r+1] = T1[r]+T2[r] mod 2^32\n");
    printf("  Если carry[r]=1: T1[r]+T2[r] ≥ 2^32\n");
    printf("  → a[r+1] = T1[r]+T2[r] - 2^32 (меньше, чем было бы без переноса)\n");
    printf("  → S0(a[r+1]) другой, но НЕ обязательно меньший!\n\n");
    printf("  S0 = ROTR(x,2)^ROTR(x,13)^ROTR(x,22) — это XOR вращений\n");
    printf("  Нет монотонной связи: большой x не → большой S0(x)\n\n");
    printf("  Механизм анти-корреляции более тонкий:\n");
    printf("  При W близком к решению (Wang SC), a[r+1] ограничен конкретными\n");
    printf("  битовыми условиями. Carry[r]=1 означает определённый паттерн битов,\n");
    printf("  который через Wang SC КОРРЕЛИРУЕТ с carry[r+1].\n");
    printf("  Знак зависит от конкретных SC в раунде r+1.\n");
}

/* RNG */
static uint64_t rng_s;
static inline uint64_t xr(void){rng_s^=rng_s<<13;rng_s^=rng_s>>7;rng_s^=rng_s<<17;return rng_s;}

typedef struct {
    long target85, cnt85, cnt90;
    uint64_t seed; int tid;
    /* Accumulate: all lag=1 pairs and global matrix */
    double sx[64], sy[64], sxy_lag1[64];  /* lag=1: i→i+1 */
    double sx_all[64], sxy_all[64*64];
    long N85, N90;
} Deep85Arg;

static void *deep85_worker(void *v){
    Deep85Arg *a=(Deep85Arg*)v;
    rng_s=a->seed^(uint64_t)a->tid*0xdeadbeef1234LL;
    
    uint32_t W[16], W1[16];
    for(int i=0;i<16;i++) W[i]=(uint32_t)xr();
    
    double T=60.0;
    uint32_t bestW[16]; int besthw=256;
    memcpy(bestW,W,64);
    
    long iters=0;
    while(a->cnt85<a->target85){
        iters++;
        W1[0]=W[0]^0x80000000u;
        for(int i=1;i<16;i++) W1[i]=W[i];
        int hw=sha256_hw59(W,W1);
        
        uint8_t carry[64];
        GlassBox g;
        sha256_glass(W,&g);
        for(int i=0;i<64;i++) carry[i]=g.carry[i];

        if(hw<90){
            a->N90++;
            for(int i=0;i<64;i++){
                a->sx[i]+=carry[i];
                if(i+1<64) a->sxy_lag1[i]+=carry[i]*carry[i+1];
            }
        }
        if(hw<85){
            a->cnt85++;
            a->N85++;
            for(int i=0;i<64;i++){
                a->sx_all[i]+=carry[i];
                for(int j=i+1;j<64;j++)
                    a->sxy_all[i*64+j]+=carry[i]*carry[j];
            }
        }
        
        if(hw<besthw){besthw=hw;memcpy(bestW,W,64);}
        
        /* SA */
        int k=(int)(xr()%16), bit=(int)(xr()%32);
        uint32_t saved=W[k];
        W[k]^=(1u<<bit);
        W1[0]=W[0]^0x80000000u;
        int nhw=sha256_hw59(W,W1);
        double dE=nhw-hw;
        if(dE>0 && (double)(xr()%1000000)/1000000.0>exp(-dE/T)) W[k]=saved;
        T*=0.9999975;
        if(T<0.3) T=0.3;
        if(iters%300000==0){memcpy(W,bestW,64);T=40.0+(double)(xr()%30);}
        if(iters%5000000==0&&besthw>92){for(int i=0;i<16;i++)W[i]=(uint32_t)xr();besthw=256;T=60.0;}
    }
    return NULL;
}

int main(int argc, char **argv){
    int nthreads=6;
    long target=100;  /* per thread */
    if(argc>1) target=atol(argv[1]);
    if(argc>2) nthreads=atoi(argv[2]);
    
    printf("anomaly_probe: target=%ld образцов hw59<85 на поток, потоков=%d\n\n",target,nthreads);
    
    analyze_lag1_signs();
    
    printf("\n=== Получение hw59<85 образцов ===\n");
    time_t t0=time(NULL);
    
    Deep85Arg *args=calloc(nthreads,sizeof(Deep85Arg));
    pthread_t *th=malloc(nthreads*sizeof(pthread_t));
    for(int t=0;t<nthreads;t++){
        args[t].tid=t;
        args[t].target85=target;
        args[t].seed=(uint64_t)time(NULL)^(uint64_t)t*0x987654321LL;
        pthread_create(&th[t],NULL,deep85_worker,&args[t]);
    }
    for(int t=0;t<nthreads;t++) pthread_join(th[t],NULL);
    
    /* Объединяем */
    long tot85=0,tot90=0;
    double sx[64]={0}, sxy_lag1[64]={0};
    double sx_all[64]={0}, sxy_all[64*64];
    memset(sxy_all,0,sizeof(sxy_all));
    
    for(int t=0;t<nthreads;t++){
        tot85+=args[t].N85; tot90+=args[t].N90;
        for(int i=0;i<64;i++){
            sx[i]+=args[t].sx[i];
            sxy_lag1[i]+=args[t].sxy_lag1[i];
            sx_all[i]+=args[t].sx_all[i];
            for(int j=i+1;j<64;j++)
                sxy_all[i*64+j]+=args[t].sxy_all[i*64+j];
        }
    }
    
    printf("Время: %ld сек\n",(long)(time(NULL)-t0));
    printf("Образцы: hw59<85=%ld hw59<90=%ld\n\n",tot85,tot90);
    
    /* Lag=1 анализ при hw59<90 */
    if(tot90>50){
        printf("=== lag=1 корреляции при hw59<90 (N=%ld) ===\n",tot90);
        double N=(double)tot90;
        for(int i=0;i<63;i++){
            double mi=sx[i]/N, mj=sx[i+1]/N;
            double vi=mi*(1-mi), vj=mj*(1-mj);
            if(vi<1e-9||vj<1e-9){printf("  c[%2d]↔c[%2d]: degenerate (mi=%.3f mj=%.3f)\n",i,i+1,mi,mj);continue;}
            double mij=sxy_lag1[i]/N;
            double phi=(mij-mi*mj)/sqrt(vi*vj);
            printf("  c[%2d]↔c[%2d]: phi=%+.4f  (E[c%d]=%.3f E[c%d]=%.3f)\n",
                   i,i+1,phi,i,mi,i+1,mj);
        }
    }
    
    /* Топ пар при hw59<85 */
    if(tot85>50){
        printf("\n=== Топ-15 пар при hw59<85 (N=%ld) ===\n",tot85);
        typedef struct{int i,j;double phi;}Pair;
        Pair best[15]; int nb=0;
        double Esum=0; int npairs=0;
        for(int i=0;i<64;i++) for(int j=i+1;j<64;j++){
            if(tot85<30) continue;
            double N=(double)tot85;
            double mi=sx_all[i]/N, mj=sx_all[j]/N;
            double vi=mi*(1-mi), vj=mj*(1-mj);
            if(vi<1e-9||vj<1e-9) continue;
            double mij=sxy_all[i*64+j]/N;
            double phi=(mij-mi*mj)/sqrt(vi*vj);
            Esum+=fabs(phi); npairs++;
            if(nb<15){best[nb].i=i;best[nb].j=j;best[nb].phi=phi;nb++;}
            else{
                int worst=0;
                for(int k=1;k<15;k++) if(fabs(best[k].phi)<fabs(best[worst].phi)) worst=k;
                if(fabs(phi)>fabs(best[worst].phi)){best[worst].i=i;best[worst].j=j;best[worst].phi=phi;}
            }
        }
        /* Sort */
        for(int a=0;a<nb-1;a++) for(int b=a+1;b<nb;b++){
            if(fabs(best[b].phi)>fabs(best[a].phi)){Pair tmp=best[a];best[a]=best[b];best[b]=tmp;}
        }
        for(int k=0;k<nb;k++)
            printf("  c[%2d]↔c[%2d] lag=%2d: phi=%+.4f\n",
                   best[k].i,best[k].j,best[k].j-best[k].i,best[k].phi);
        if(npairs>0)
            printf("\nЭнергия при hw59<85: E=%.4f (норм=%.4f)\n",Esum,(double)Esum/npairs);
    }
    
    free(args); free(th);
    return 0;
}
