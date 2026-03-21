/*
 * dsc_rank_sweep.c — DSC + точно такой же SA как rank_sweep
 *
 * Цель: получить 1500+ образцов при hw59<85 (как rank_sweep)
 * и измерить: сколько бит delta[64] зафиксированы → оставшаяся сложность
 *
 * Q-D2: нефиксированные биты delta[64] = оставшаяся сложность
 * Q-D5: профиль rank(delta-code) = новый вид линейного кода
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

static const uint32_t K256[64]={
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
static const uint32_t IV[8]={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
                               0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
#define ROTR(x,n)(((x)>>(n))|((x)<<(32-(n))))
#define S0(x)(ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define S1(x)(ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define s0(x)(ROTR(x,7)^ROTR(x,18)^((x)>>3))
#define s1(x)(ROTR(x,17)^ROTR(x,19)^((x)>>10))
#define CH(e,f,g)(((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c)(((a)&(b))^((a)&(c))^((b)&(c)))

/* Вычислить hw59 И delta[64] */
static int sha256_full(const uint32_t W0[16], uint32_t diff64[8]){
    uint32_t W[2][64];
    for(int i=0;i<16;i++){W[0][i]=W0[i];W[1][i]=W0[i];}
    W[1][0]^=0x80000000u;
    for(int i=16;i<64;i++){
        W[0][i]=s1(W[0][i-2])+W[0][i-7]+s0(W[0][i-15])+W[0][i-16];
        W[1][i]=s1(W[1][i-2])+W[1][i-7]+s0(W[1][i-15])+W[1][i-16];
    }
    uint32_t s[2][8];
    for(int j=0;j<2;j++) for(int k=0;k<8;k++) s[j][k]=IV[k];
    for(int r=0;r<64;r++) for(int j=0;j<2;j++){
        uint32_t a=s[j][0],b=s[j][1],c=s[j][2],d=s[j][3],
                 e=s[j][4],f=s[j][5],g=s[j][6],h=s[j][7];
        uint32_t T1=h+S1(e)+CH(e,f,g)+K256[r]+W[j][r], T2=S0(a)+MAJ(a,b,c);
        s[j][7]=g;s[j][6]=f;s[j][5]=e;s[j][4]=d+T1;
        s[j][3]=c;s[j][2]=b;s[j][1]=a;s[j][0]=T1+T2;
    }
    int hw=0;
    for(int k=0;k<8;k++){
        uint32_t d=s[0][k]^s[1][k];
        if(diff64) diff64[k]=d;
        hw+=__builtin_popcount(d);
    }
    return hw;
}

/* Быстрая версия только hw59 */
static int hw59_fast(const uint32_t W0[16]){
    return sha256_full(W0,NULL);
}

/* RNG */
static __thread uint64_t rng_s;
static inline uint64_t xr(){rng_s^=rng_s<<13;rng_s^=rng_s>>7;rng_s^=rng_s<<17;return rng_s;}

/* Аккумулятор DSC для одного уровня */
typedef struct {
    int thr;
    long target, n;
    long zero[8][32];   /* счётчик delta64[w][b]==0 */
    double sum1[8][32]; /* среднее delta64[w][b] */
    double hw_sum;
    pthread_mutex_t mu;
} DSCLevel;

static void dsc_add(DSCLevel *lv, const uint32_t diff64[8], int hw){
    pthread_mutex_lock(&lv->mu);
    if(lv->n<lv->target){
        lv->n++;
        lv->hw_sum+=hw;
        for(int w=0;w<8;w++){
            for(int b=0;b<32;b++){
                int bit=(diff64[w]>>b)&1;
                lv->sum1[w][b]+=bit;
                if(!bit) lv->zero[w][b]++;
            }
        }
    }
    pthread_mutex_unlock(&lv->mu);
}

#define NLEVELS 4
static int THRESHOLDS[NLEVELS]={95,90,85,80};

static DSCLevel levels[NLEVELS];

/* SA worker — ТОЧНАЯ КОПИЯ rank_sweep.c */
typedef struct { int tid; long max_iters; uint64_t seed; } TArg;

static void *sa_worker(void *v){
    TArg *a=(TArg*)v;
    rng_s=a->seed^(uint64_t)a->tid*0xdeadbeef12LL;
    uint32_t W[16]; for(int i=0;i<16;i++) W[i]=(uint32_t)xr();
    double T=60.0; uint32_t bW[16]; int bh=256; memcpy(bW,W,64);
    long iter=0;
    while(iter<a->max_iters){
        iter++;
        int hw=hw59_fast(W);

        /* Проверяем, нужно ли собрать образец */
        int need=0;
        for(int lv=0;lv<NLEVELS;lv++)
            if(hw<THRESHOLDS[lv] && levels[lv].n<levels[lv].target){need=1;break;}
        if(need){
            uint32_t diff64[8];
            sha256_full(W,diff64);  /* дополнительный вызов только когда нужен */
            for(int lv=0;lv<NLEVELS;lv++)
                if(hw<THRESHOLDS[lv]) dsc_add(&levels[lv],diff64,hw);
        }

        if(hw<bh){bh=hw;memcpy(bW,W,64);}
        int k=(int)(xr()%16),bit=(int)(xr()%32);
        uint32_t sv=W[k]; W[k]^=(1u<<bit);
        int nh=hw59_fast(W);
        double dE=nh-hw;
        if(dE<=0||(double)(xr()%1000000)/1e6<exp(-dE/T)){
            /* принять */
        } else W[k]=sv;
        T*=0.9999975; if(T<0.3)T=0.3;
        if(iter%300000==0){memcpy(W,bW,64);T=40.0+(double)(xr()%30);}
        if(iter%5000000==0&&bh>92){for(int i=0;i<16;i++)W[i]=(uint32_t)xr();bh=256;T=60.0;}
    }
    return NULL;
}

/* GF(2) Union-Find для DCC rank (аналог WCC_rank) */
typedef struct{int par,xr;}UF;
static int uf_find(UF *u,int i,int *xi){
    int x=0;
    while(u[i].par!=i){x^=u[i].xr;i=u[i].par;}
    *xi=x;return i;
}
static int uf_merge(UF *u,int i,int j,int sign){
    int xi,xj;int ri=uf_find(u,i,&xi),rj=uf_find(u,j,&xj);
    if(ri==rj)return(xi^xj)!=sign;
    u[ri].par=rj;u[ri].xr=xi^xj^sign;return 0;
}

/* Вычислить rank DSC кода (аналог carry_rank) */
static void analyze_dsc(DSCLevel *lv, double phi_thresh){
    long N=lv->n;
    if(N<10){printf("  [hw59<%d]: недостаточно N=%ld\n",lv->thr,N);return;}
    const char *wn[8]={"a","b","c","d","e","f","g","h"};

    printf("\n╔══════════════════════════════════════════════════╗\n");
    printf("║  DSC: hw59<%d  (N=%ld)                  \n",lv->thr,N);
    printf("╚══════════════════════════════════════════════════╝\n");
    printf("  <hw59>=%.2f\n",lv->hw_sum/N);

    /* Q-D2: зафиксированные биты delta[64] */
    int nfixed=0, fixed_w[8]={0};
    int always0[256]={0}; /* индекс w*32+b */
    for(int w=0;w<8;w++)for(int b=0;b<32;b++){
        double p=(double)lv->zero[w][b]/N;
        if(p>=phi_thresh){
            nfixed++;fixed_w[w]++;
            always0[w*32+b]=1;
        }
    }
    printf("\n  Зафиксированных бит delta[64] (P[=0]>=%.2f): %d из 256\n",phi_thresh,nfixed);
    printf("  Нефиксированных: %d из 256\n",256-nfixed);
    printf("  По словам: ");
    for(int w=0;w<8;w++) printf("%s=%d ",wn[w],fixed_w[w]);
    printf("\n  → Оставшаяся сложность: O(2^%d)\n",256-nfixed);

    /* Гистограмма P[=0] */
    int hist[11]={0};
    for(int w=0;w<8;w++)for(int b=0;b<32;b++){
        double p=(double)lv->zero[w][b]/N;
        int h=(int)(p*10); if(h>10)h=10;
        hist[h]++;
    }
    printf("  Гистограмма P[=0]: ");
    for(int h=10;h>=0;h--) if(hist[h]>0) printf("[%.1f-%.1f)=%d ", h*0.1,(h+1)*0.1,hist[h]);
    printf("\n");

    /* Union-Find rank для нефиксированных бит */
    /* Ищем пары (i,j) из нефиксированных с |phi|>0.9 */
    int free_bits[256],nfree=0;
    for(int i=0;i<256;i++) if(!always0[i]){free_bits[nfree++]=i;}
    printf("  Нефиксированных бит: %d\n",nfree);

    if(nfree>2 && N>20){
        /* Вычислим ранг через UF на нефиксированных битах */
        UF *uf=calloc(nfree,sizeof(UF));
        for(int i=0;i<nfree;i++){uf[i].par=i;uf[i].xr=0;}
        int nconstr=0;
        for(int i=0;i<nfree;i++){
            int wi=free_bits[i]/32, bi=free_bits[i]%32;
            double pi=(double)lv->sum1[wi][bi]/N;
            for(int j=i+1;j<nfree;j++){
                int wj=free_bits[j]/32, bj=free_bits[j]%32;
                double pj=(double)lv->sum1[wj][bj]/N;
                double pij=(double)(N-lv->zero[wi][bi]-lv->zero[wj][bj])/N;
                /* E[xi*xj] нет — нужно вычислить из данных... пропускаем */
                (void)pi;(void)pj;(void)pij;
                break; /* TODO: need sxy tracking */
            }
            break;
        }
        printf("  (rank анализ нефиксированных бит требует sxy tracking)\n");
        free(uf);
        (void)nconstr;
    }

    /* Топ-20 наиболее фиксированных */
    printf("\n  Топ-20 фиксированных бит delta[64] (P[=0]):\n");
    typedef struct{int w,b;double p;}FB;
    FB top[20]; int ntop=0;
    for(int w=0;w<8;w++)for(int b=0;b<32;b++){
        double p=(double)lv->zero[w][b]/N;
        if(ntop<20){top[ntop].w=w;top[ntop].b=b;top[ntop].p=p;ntop++;}
        else{
            int worst=0;
            for(int k=1;k<20;k++)if(top[k].p<top[worst].p)worst=k;
            if(p>top[worst].p){top[worst].w=w;top[worst].b=b;top[worst].p=p;}
        }
    }
    for(int k=0;k<ntop-1;k++)for(int j=k+1;j<ntop;j++)
        if(top[j].p>top[k].p){FB t=top[k];top[k]=top[j];top[j]=t;}
    for(int k=0;k<ntop;k++)
        printf("    delta[64].%s[b%2d]: P[=0]=%.4f\n",wn[top[k].w],top[k].b,top[k].p);

    /* Средний hw59 и дисперсия — для понимания что delta[64] выглядит */
    printf("\n  Полная карта P[=0] для delta[64] (нефиксированные, топ):\n");
    /* Распределение P[=0] по нефиксированным */
    int n09=0,n08=0,n07=0,n06=0,n05=0;
    for(int w=0;w<8;w++)for(int b=0;b<32;b++){
        double p=(double)lv->zero[w][b]/N;
        if(!always0[w*32+b]){
            if(p>=0.9)n09++;
            else if(p>=0.8)n08++;
            else if(p>=0.7)n07++;
            else if(p>=0.6)n06++;
            else n05++;
        }
    }
    printf("  Нефиксированных с P[=0]>=0.9: %d\n",n09);
    printf("  Нефиксированных с P[=0]>=0.8: %d\n",n08);
    printf("  Нефиксированных с P[=0]>=0.7: %d\n",n07);
    printf("  Нефиксированных с P[=0]>=0.6: %d\n",n06);
    printf("  Нефиксированных с P[=0]~0.5:  %d\n",n05);
}

int main(int argc,char**argv){
    int nthreads=6;
    long max_iters=400000000L;
    if(argc>1) max_iters=atol(argv[1]);
    if(argc>2) nthreads=atoi(argv[2]);

    printf("dsc_rank_sweep: итераций=%ld, потоков=%d\n\n",max_iters,nthreads);
    printf("Q-D2: Сколько бит delta[64] нефиксированы при hw59<85?\n");
    printf("Q-D5: Профиль P[=0] по уровням hw59\n\n");

    for(int lv=0;lv<NLEVELS;lv++){
        levels[lv].thr=THRESHOLDS[lv];
        levels[lv].target=1500;
        levels[lv].n=0;
        memset(levels[lv].zero,0,sizeof(levels[lv].zero));
        memset(levels[lv].sum1,0,sizeof(levels[lv].sum1));
        levels[lv].hw_sum=0;
        pthread_mutex_init(&levels[lv].mu,NULL);
    }

    time_t t0=time(NULL);
    TArg *args=calloc(nthreads,sizeof(TArg));
    pthread_t *th=malloc(nthreads*sizeof(pthread_t));
    for(int t=0;t<nthreads;t++){
        args[t].tid=t;
        args[t].max_iters=max_iters;
        args[t].seed=(uint64_t)time(NULL)^(uint64_t)t*0xabcdef12345LL;
        pthread_create(&th[t],NULL,sa_worker,&args[t]);
    }

    /* Мониторинг */
    for(;;){
        sleep(60);
        long ns[NLEVELS];
        for(int lv=0;lv<NLEVELS;lv++){
            pthread_mutex_lock(&levels[lv].mu);
            ns[lv]=levels[lv].n;
            pthread_mutex_unlock(&levels[lv].mu);
        }
        long elapsed=(long)(time(NULL)-t0);
        printf("[%lds] hw59<95:%ld  <90:%ld  <85:%ld  <80:%ld\n",
               elapsed,ns[0],ns[1],ns[2],ns[3]);
        fflush(stdout);
        /* Стоп когда все заполнены ИЛИ hw59<85 достигнуто */
        if(ns[2]>=1500 && ns[3]>=1500) break;
    }

    for(int t=0;t<nthreads;t++) pthread_join(th[t],NULL);
    printf("Время: %ld сек\n",(long)(time(NULL)-t0));

    for(int lv=NLEVELS-1;lv>=0;lv--)
        analyze_dsc(&levels[lv],0.90);

    printf("\n╔═══════════════════════════════════════════════════╗\n");
    printf("║  ИТОГ: РАССТОЯНИЕ ДО КОЛЛИЗИИ                    ║\n");
    printf("╚═══════════════════════════════════════════════════╝\n");
    for(int lv=NLEVELS-1;lv>=0;lv--){
        if(levels[lv].n<10)continue;
        long N=levels[lv].n;
        int nfixed=0;
        for(int w=0;w<8;w++)for(int b=0;b<32;b++){
            double p=(double)levels[lv].zero[w][b]/N;
            if(p>=0.90)nfixed++;
        }
        printf("  hw59<%2d: зафиксировано %3d бит delta[64]  сложность O(2^%d)\n",
               levels[lv].thr,nfixed,256-nfixed);
    }

    free(args);free(th);
    return 0;
}
