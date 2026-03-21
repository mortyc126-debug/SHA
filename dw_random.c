/*
 * dw_random.c — Случайный поиск DW с минимальным hw59_min
 *
 * Перебираем случайные DW с HW=1,2,3,4,8 и ищем рекорды.
 * Для каждого DW: быстрый SA (30M iter × 4 threads)
 * Сохраняем топ-50 DW с минимальным hw59_min.
 *
 * Вопрос: существует ли DW с hw59_min < 40?
 * (для любого HW)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>

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

static __thread uint64_t rng_s;
static inline uint64_t xr(){rng_s^=rng_s<<13;rng_s^=rng_s>>7;rng_s^=rng_s<<17;return rng_s;}

static int hw59_dw(const uint32_t W0[16], const uint32_t DW[16]){
    uint32_t W[2][64];
    for(int i=0;i<16;i++){W[0][i]=W0[i];W[1][i]=W0[i]^DW[i];}
    for(int i=16;i<64;i++){
        W[0][i]=s1(W[0][i-2])+W[0][i-7]+s0(W[0][i-15])+W[0][i-16];
        W[1][i]=s1(W[1][i-2])+W[1][i-7]+s0(W[1][i-15])+W[1][i-16];
    }
    uint32_t s[2][8];
    for(int j=0;j<2;j++)for(int k=0;k<8;k++)s[j][k]=IV[k];
    for(int r=0;r<64;r++)for(int j=0;j<2;j++){
        uint32_t a=s[j][0],b=s[j][1],c=s[j][2],d=s[j][3],e=s[j][4],f=s[j][5],g=s[j][6],h=s[j][7];
        uint32_t T1=h+S1(e)+CH(e,f,g)+K256[r]+W[j][r],T2=S0(a)+MAJ(a,b,c);
        s[j][7]=g;s[j][6]=f;s[j][5]=e;s[j][4]=d+T1;s[j][3]=c;s[j][2]=b;s[j][1]=a;s[j][0]=T1+T2;
    }
    int hw=0;for(int k=0;k<8;k++)hw+=__builtin_popcount(s[0][k]^s[1][k]);
    return hw;
}

/* SA worker для данного DW */
typedef struct {
    int tid; const uint32_t *DW; uint64_t seed; long max_iters;
    int result_hw; uint32_t result_W[16];
} SAArg;

static void *sa_worker(void *v){
    SAArg *a=(SAArg*)v;
    rng_s=a->seed^(uint64_t)a->tid*0xdeadbeef;
    uint32_t W[16]; for(int i=0;i<16;i++)W[i]=(uint32_t)xr();
    double T=60.; uint32_t bW[16]; int bh=256; memcpy(bW,W,64);
    for(long it=0;it<a->max_iters;it++){
        int hw=hw59_dw(W,a->DW);
        if(hw<bh){bh=hw;memcpy(bW,W,64);}
        int k=(int)(xr()%16),bit=(int)(xr()%32);
        uint32_t sv=W[k]; W[k]^=(1u<<bit);
        int nh=hw59_dw(W,a->DW); double dE=nh-hw;
        if(dE<=0||(double)(xr()%1000000)/1e6<exp(-dE/T)){}else W[k]=sv;
        T*=0.9999975;if(T<0.3)T=0.3;
        if(it%300000==0){memcpy(W,bW,64);T=40.+(double)(xr()%30);}
        if(it%5000000==0&&bh>92){for(int i=0;i<16;i++)W[i]=(uint32_t)xr();bh=256;T=60.;}
    }
    a->result_hw=bh; memcpy(a->result_W,bW,64);
    return NULL;
}

static int eval_dw_quick(const uint32_t DW[16], int nthreads, long iters, uint64_t seed){
    SAArg *args=calloc(nthreads,sizeof(SAArg));
    pthread_t *th=malloc(nthreads*sizeof(pthread_t));
    for(int t=0;t<nthreads;t++){
        args[t].tid=t;args[t].DW=DW;args[t].seed=seed^(uint64_t)t*0xabcdef;
        args[t].max_iters=iters;pthread_create(&th[t],NULL,sa_worker,&args[t]);
    }
    int best=256;
    for(int t=0;t<nthreads;t++){
        pthread_join(th[t],NULL);
        if(args[t].result_hw<best)best=args[t].result_hw;
    }
    free(args);free(th);return best;
}

/* Генератор случайного DW с заданным HW */
static void rand_dw(uint32_t DW[16], int target_hw, uint64_t *rng){
    memset(DW,0,64);
    int placed=0;
    while(placed<target_hw){
        /* случайный бит */
        uint64_t r=(*rng);(*rng)^=(*rng)<<13;(*rng)^=(*rng)>>7;(*rng)^=(*rng)<<17;
        int w=(int)(r%16); r=(*rng);(*rng)^=(*rng)<<13;(*rng)^=(*rng)>>7;(*rng)^=(*rng)<<17;
        int b=(int)(r%32);
        if(!((DW[w]>>(b))&1)){DW[w]|=(1u<<b);placed++;}
    }
}

/* Топ-50 рекордов */
#define TOPN 50
static uint32_t top_dw[TOPN][16];
static int top_hw[TOPN];
static int top_hwdw[TOPN];
static int top_n=0;
static pthread_mutex_t top_mu=PTHREAD_MUTEX_INITIALIZER;

static void record_result(const uint32_t DW[16], int hw59, int hwdw){
    pthread_mutex_lock(&top_mu);
    /* Добавить или вытеснить худший */
    if(top_n<TOPN){
        memcpy(top_dw[top_n],DW,64);
        top_hw[top_n]=hw59;top_hwdw[top_n]=hwdw;
        top_n++;
    } else {
        /* Найти максимум */
        int worst=0;
        for(int i=1;i<TOPN;i++) if(top_hw[i]>top_hw[worst])worst=i;
        if(hw59<top_hw[worst]){
            memcpy(top_dw[worst],DW,64);
            top_hw[worst]=hw59;top_hwdw[worst]=hwdw;
        }
    }
    pthread_mutex_unlock(&top_mu);
}

/* Параллельный поиск */
typedef struct {
    int tid;
    uint64_t seed;
    long sa_iters;
    int sa_threads;   /* потоков на SA */
    int n_trials;     /* сколько DW попробовать */
    int global_best;  /* best known hw59 */
} SearchArg;

static volatile int global_best_hw=256;
static pthread_mutex_t best_mu=PTHREAD_MUTEX_INITIALIZER;

/* Однопоточный SA (для встраивания в параллельный поиск) */
static int sa_single(const uint32_t DW[16], long iters, uint64_t seed){
    rng_s=seed;
    uint32_t W[16]; for(int i=0;i<16;i++)W[i]=(uint32_t)xr();
    double T=60.; uint32_t bW[16]; int bh=256; memcpy(bW,W,64);
    for(long it=0;it<iters;it++){
        int hw=hw59_dw(W,DW);
        if(hw<bh){bh=hw;memcpy(bW,W,64);}
        int k=(int)(xr()%16),bit=(int)(xr()%32);
        uint32_t sv=W[k]; W[k]^=(1u<<bit);
        int nh=hw59_dw(W,DW); double dE=nh-hw;
        if(dE<=0||(double)(xr()%1000000)/1e6<exp(-dE/T)){}else W[k]=sv;
        T*=0.9999975;if(T<0.3)T=0.3;
        if(it%300000==0){memcpy(W,bW,64);T=40.+(double)(xr()%30);}
        if(it%5000000==0&&bh>92){for(int i=0;i<16;i++)W[i]=(uint32_t)xr();bh=256;T=60.;}
    }
    return bh;
}

static void *search_worker(void *v){
    SearchArg *a=(SearchArg*)v;
    rng_s=a->seed^(uint64_t)a->tid*0xdeadbeef12LL;
    uint64_t rng=rng_s^0xfedcba9876543210ULL;

    int hw_targets[]={1,2,2,2,3,3,3,4,4,8};
    int n_targets=10;

    for(int trial=0;trial<a->n_trials;trial++){
        /* Выбор случайного HW */
        uint64_t r=rng; rng^=rng<<13;rng^=rng>>7;rng^=rng<<17;
        int hwt=hw_targets[r%n_targets];

        uint32_t DW[16];
        rand_dw(DW,hwt,&rng);

        /* Быстрый SA (однопоточный, 30M итераций) */
        uint64_t s2=rng; rng^=rng<<13;rng^=rng>>7;rng^=rng<<17;
        int hw=sa_single(DW,a->sa_iters,s2);

        int hwdw=0; for(int w=0;w<16;w++) hwdw+=__builtin_popcount(DW[w]);

        /* Проверить глобальный рекорд */
        pthread_mutex_lock(&best_mu);
        int is_record=(hw<global_best_hw);
        if(is_record){
            global_best_hw=hw;
            printf("[tid%d trial%d] ★ НОВЫЙ РЕКОРД hw59=%d (HW(DW)=%d): ",
                   a->tid,trial,hw,hwdw);
            for(int w=0;w<16;w++) if(DW[w])printf("DW[%d]=0x%08x ",w,DW[w]);
            printf("\n"); fflush(stdout);
        }
        pthread_mutex_unlock(&best_mu);

        /* Запись в топ если hw≤80 */
        if(hw<=80) record_result(DW,hw,hwdw);

        /* Прогресс каждые 100 триалов */
        if((trial+1)%100==0){
            pthread_mutex_lock(&best_mu);
            printf("[tid%d] %d trials done, global_best=%d\n",
                   a->tid,trial+1,global_best_hw);
            fflush(stdout);
            pthread_mutex_unlock(&best_mu);
        }
    }
    return NULL;
}

int main(int argc,char**argv){
    int nthreads=6;
    long sa_iters=30000000L;  /* 30M на один SA */
    int n_trials=500;         /* попыток на поток */
    if(argc>1) sa_iters=atol(argv[1]);
    if(argc>2) n_trials=atoi(argv[2]);
    if(argc>3) nthreads=atoi(argv[3]);

    printf("dw_random: случайный поиск DW с min hw59\n");
    printf("SA iters=%ldM, trials_per_thread=%d, threads=%d\n\n",
           sa_iters/1000000,n_trials,nthreads);
    printf("Стратегия: случайные DW с HW=1,2,3,4,8\n");
    printf("Цель: найти DW с hw59_min < 40 (O(2^40))\n\n");
    fflush(stdout);

    /* Базовый ориентир: Wang ΔW[0]=0x80000000 */
    {
        uint32_t DW_wang[16]={0}; DW_wang[0]=0x80000000u;
        printf("Базовый Wang ΔW[0]=0x80000000:\n");
        int hw_wang=eval_dw_quick(DW_wang,nthreads,100000000L,(uint64_t)time(NULL));
        printf("  hw59_min(Wang) = %d\n\n",hw_wang);
        if(hw_wang<global_best_hw) global_best_hw=hw_wang;
        fflush(stdout);
    }

    SearchArg *args=calloc(nthreads,sizeof(SearchArg));
    pthread_t *th=malloc(nthreads*sizeof(pthread_t));
    for(int t=0;t<nthreads;t++){
        args[t].tid=t;
        args[t].seed=(uint64_t)time(NULL)^((uint64_t)t*0xc0ffee123LL);
        args[t].sa_iters=sa_iters;
        args[t].n_trials=n_trials;
        pthread_create(&th[t],NULL,search_worker,&args[t]);
    }
    for(int t=0;t<nthreads;t++) pthread_join(th[t],NULL);

    printf("\n╔══════════════════════════════════════════════════╗\n");
    printf("║  ИТОГ RANDOM DW SEARCH                          ║\n");
    printf("╚══════════════════════════════════════════════════╝\n");
    printf("  Глобальный минимум: hw59=%d\n",global_best_hw);
    printf("  Всего записей в топ (hw59≤80): %d\n\n",top_n);

    if(top_n>0){
        /* Сортировка по hw59 */
        for(int i=0;i<top_n-1;i++)
            for(int j=i+1;j<top_n;j++)
                if(top_hw[j]<top_hw[i]){
                    uint32_t tmp[16]; memcpy(tmp,top_dw[i],64);
                    memcpy(top_dw[i],top_dw[j],64); memcpy(top_dw[j],tmp,64);
                    int ti=top_hw[i];top_hw[i]=top_hw[j];top_hw[j]=ti;
                    ti=top_hwdw[i];top_hwdw[i]=top_hwdw[j];top_hwdw[j]=ti;
                }
        printf("  Топ-%d DW (сортировано по hw59):\n",top_n);
        for(int i=0;i<top_n;i++){
            printf("  #%2d hw59=%d HW(DW)=%d: ",i+1,top_hw[i],top_hwdw[i]);
            for(int w=0;w<16;w++) if(top_dw[i][w])printf("DW[%d]=0x%08x ",w,top_dw[i][w]);
            printf("\n");
        }
    }

    if(global_best_hw<40)
        printf("\n  ★★ O(2^%d) ДОСТИГНУТА! ★★\n",global_best_hw);
    else
        printf("\n  До O(2^40): нужно снизить ещё %d\n",global_best_hw-40);

    free(args);free(th);
    return 0;
}
