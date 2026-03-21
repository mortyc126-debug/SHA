/*
 * dsc_probe2.c — DSC (Differential State Correlation) v2
 *
 * Улучшенная версия: более быстрый SA для hw59<85/90
 * Цель: измерить, сколько бит δ[64] фиксированы при hw59<85/90
 *
 * Q-D2: Оставшаяся сложность до коллизии = нефиксированные биты δ[64]
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

/* Компактный трекинг только нужного */
typedef struct {
    uint32_t diff_final[8];  /* δ[64] = финальный дифференциал */
    uint32_t diff_r[8][8];   /* δ[r] для r=0,8,16,24,32,40,48,56 */
    int hw59;
} Sample;

static int sha256_hw59(const uint32_t W0[16], uint32_t diff_out[8], uint32_t diff_r[8][8]){
    uint32_t W[2][64];
    for(int i=0;i<16;i++){W[0][i]=W0[i];W[1][i]=W0[i];}
    W[1][0]^=0x80000000u;
    for(int i=16;i<64;i++){
        W[0][i]=s1(W[0][i-2])+W[0][i-7]+s0(W[0][i-15])+W[0][i-16];
        W[1][i]=s1(W[1][i-2])+W[1][i-7]+s0(W[1][i-15])+W[1][i-16];
    }
    uint32_t s[2][8];
    for(int j=0;j<2;j++)for(int k=0;k<8;k++)s[j][k]=IV[k];

    for(int r=0;r<64;r++){
        /* сохраняем δ перед раундом r=0,8,...,56 */
        if(r%8==0 && diff_r){
            int ri=r/8;
            for(int k=0;k<8;k++) diff_r[ri][k]=s[0][k]^s[1][k];
        }
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
    for(int k=0;k<8;k++){
        uint32_t d=s[0][k]^s[1][k];
        if(diff_out) diff_out[k]=d;
        hw+=__builtin_popcount(d);
    }
    return hw;
}

/* Быстрая версия (без трекинга промежуточных раундов) */
static int sha256_fast(const uint32_t W0[16]){
    return sha256_hw59(W0,NULL,NULL);
}

/* RNG */
static __thread uint64_t rng_s;
static inline uint64_t xr(){rng_s^=rng_s<<13;rng_s^=rng_s>>7;rng_s^=rng_s<<17;return rng_s;}

/* Аккумулятор */
#define NCHECK 8  /* r=0,8,16,24,32,40,48,56 */
typedef struct {
    long N;
    /* P[bit=0] для δ[64] */
    long zero_final[8][32];
    /* среднее δ[r][w][b] */
    double mean_r[NCHECK][8][32];
    double hw_sum;
} DSCAcc;

static void acc_init(DSCAcc *a){ memset(a,0,sizeof(*a)); }

static void acc_add(DSCAcc *a, const uint32_t diff_final[8], const uint32_t diff_r[8][8], int hw){
    a->N++;
    a->hw_sum+=hw;
    for(int w=0;w<8;w++){
        for(int b=0;b<32;b++){
            if(!((diff_final[w]>>b)&1)) a->zero_final[w][b]++;
        }
    }
    for(int ri=0;ri<NCHECK;ri++)
        for(int w=0;w<8;w++)
            for(int b=0;b<32;b++)
                a->mean_r[ri][w][b]+=((diff_r[ri][w]>>b)&1);
}

/* SA worker с многоуровневым сбором */
#define HW_T1 90  /* уровень 1 */
#define HW_T2 85  /* уровень 2 (сложнее) */

typedef struct {
    int tid;
    long max_iters;
    long target1, target2;
    uint64_t seed;
    DSCAcc acc1, acc2;  /* hw59<90, hw59<85 */
    long cnt1, cnt2;
    int best_hw;
} WorkArg;

static void *sa_worker(void *v){
    WorkArg *a=(WorkArg*)v;
    rng_s=a->seed^(uint64_t)a->tid*0xdeadbeef12LL;
    acc_init(&a->acc1); acc_init(&a->acc2);
    a->cnt1=0; a->cnt2=0; a->best_hw=256;

    uint32_t W[16]; for(int i=0;i<16;i++) W[i]=(uint32_t)xr();
    uint32_t bW[16]; memcpy(bW,W,64);
    double T=80.0;
    int hw=sha256_fast(W);
    long iters=0;

    while(iters<a->max_iters){
        /* Сбор образцов */
        if(hw<HW_T1 && a->cnt1<a->target1){
            /* Отслеживаем diff */
            uint32_t df[8], dr[8][8];
            sha256_hw59(W,df,dr);
            acc_add(&a->acc1,df,dr,hw);
            a->cnt1++;
            if(hw<HW_T2 && a->cnt2<a->target2){
                acc_add(&a->acc2,df,dr,hw);
                a->cnt2++;
            }
        }

        /* Стоп если оба счётчика достигнуты */
        if(a->cnt1>=a->target1 && a->cnt2>=a->target2) break;

        iters++;
        /* SA шаг */
        int k=(int)(xr()%16),bit=(int)(xr()%32);
        uint32_t sv=W[k]; W[k]^=(1u<<bit);
        int nh=sha256_fast(W);
        double dE=nh-hw;
        if(dE<=0||(double)(xr()%1000000)/1e6<exp(-dE/T)){
            hw=nh;
            if(hw<a->best_hw){a->best_hw=hw;memcpy(bW,W,64);}
        } else W[k]=sv;

        T*=0.99999985; if(T<0.25)T=0.25;

        /* Периодические рестарты */
        if(iters%500000==0){
            memcpy(W,bW,64); T=50.0+(double)(xr()%30); hw=sha256_fast(W);
        }
        if(iters%8000000==0 && a->best_hw>95){
            for(int i=0;i<16;i++)W[i]=(uint32_t)xr();
            hw=sha256_fast(W); a->best_hw=hw; memcpy(bW,W,64); T=80.0;
        }
    }
    return NULL;
}

static void print_dsc(const DSCAcc *a, const char *label, int thresh){
    if(a->N<5){printf("  [%s]: мало образцов N=%ld\n",label,a->N);return;}
    long N=a->N;
    const char *wn[8]={"a","b","c","d","e","f","g","h"};

    printf("\n╔══════════════════════════════════════════════════╗\n");
    printf("║  DSC: %s (N=%ld, hw59<%d)\n",label,N,thresh);
    printf("╚══════════════════════════════════════════════════╝\n");
    printf("  <hw59>=%.2f\n",a->hw_sum/N);

    /* Q-D2: фиксированные биты δ[64] */
    int nfixed=0, fixed_w[8]={0};
    double best_p=0;
    for(int w=0;w<8;w++)for(int b=0;b<32;b++){
        double p=(double)a->zero_final[w][b]/N;
        if(p>0.90){nfixed++;fixed_w[w]++;}
        if(p>best_p)best_p=p;
    }
    printf("\n  Q-D2: Фиксированные биты δ[64] (P[=0]>90%%):\n");
    printf("    Всего: %d из 256  (нефиксированных: %d)\n",nfixed,256-nfixed);
    printf("    По словам: ");
    for(int w=0;w<8;w++) printf("%s=%d ",wn[w],fixed_w[w]);
    printf("\n    Max P[=0]=%.4f\n",best_p);
    printf("    → Оставшаяся сложность: O(2^%d)\n",256-nfixed);

    /* Топ фиксированных бит */
    printf("\n  Топ-15 фиксированных бит δ[64]:\n");
    typedef struct{int w,b;double p;}FB;
    FB top[20]; int ntop=0;
    for(int w=0;w<8;w++)for(int b=0;b<32;b++){
        double p=(double)a->zero_final[w][b]/N;
        if(ntop<15){top[ntop].w=w;top[ntop].b=b;top[ntop].p=p;ntop++;}
        else{
            int worst=0;
            for(int k=1;k<15;k++)if(top[k].p<top[worst].p)worst=k;
            if(p>top[worst].p){top[worst].w=w;top[worst].b=b;top[worst].p=p;}
        }
    }
    for(int k=0;k<ntop-1;k++)for(int j=k+1;j<ntop;j++)
        if(top[j].p>top[k].p){FB t=top[k];top[k]=top[j];top[j]=t;}
    for(int k=0;k<ntop;k++)
        printf("    δ[64].%s[b%2d]: P[=0]=%.4f\n",wn[top[k].w],top[k].b,top[k].p);

    /* Q-D3: профиль HW(δ[r]) по раундам */
    printf("\n  Q-D3: E[HW(δ[r])] по раундам:\n");
    int rvals[NCHECK]={0,8,16,24,32,40,48,56};
    for(int ri=0;ri<NCHECK;ri++){
        double hw_r=0;
        for(int w=0;w<8;w++)for(int b=0;b<32;b++)
            hw_r+=a->mean_r[ri][w][b]/N;
        printf("    r=%2d: E[HW(δ)]=%.1f  (%.1f%%)\n",
               rvals[ri],hw_r,100.0*hw_r/256);
    }

    /* Q-D1: proto-Wang SC — биты δ[r] почти всегда =0 */
    printf("\n  Q-D1: Proto-Wang SC биты (P[δ[r][w][b]=0] > 0.95):\n");
    int nsc=0;
    for(int ri=0;ri<NCHECK;ri++){
        int r=rvals[ri];
        for(int w=0;w<8;w++)for(int b=0;b<32;b++){
            double mean=a->mean_r[ri][w][b]/N;
            double p0=1.0-mean;
            if(p0>0.95){
                if(nsc<15)
                    printf("    r=%2d δ.%s[b%2d]: P[=0]=%.3f\n",r,wn[w],b,p0);
                nsc++;
            }
        }
    }
    if(nsc>15)printf("    ... и ещё %d бит\n",nsc-15);
    printf("  Всего proto-Wang SC (r<64): %d\n",nsc);
}

int main(int argc,char**argv){
    int nthreads=6;
    long target1=500, target2=100;
    if(argc>1) target1=atol(argv[1]);
    if(argc>2) target2=atol(argv[2]);
    if(argc>3) nthreads=atoi(argv[3]);

    printf("dsc_probe2: Differential State Correlation v2\n");
    printf("target1=%ld (hw59<%d), target2=%ld (hw59<%d), потоков=%d\n\n",
           target1,HW_T1,target2,HW_T2,nthreads);

    time_t t0=time(NULL);

    WorkArg *args=calloc(nthreads,sizeof(WorkArg));
    pthread_t *th=malloc(nthreads*sizeof(pthread_t));
    for(int t=0;t<nthreads;t++){
        args[t].tid=t;
        args[t].max_iters=500000000L;
        args[t].target1=target1/nthreads+1;
        args[t].target2=target2/nthreads+1;
        args[t].seed=(uint64_t)time(NULL)^(uint64_t)t*0xc0ffeedeadLL;
        pthread_create(&th[t],NULL,sa_worker,&args[t]);
    }

    /* Мониторинг прогресса */
    for(int sec=0;;sec++){
        sleep(10);
        long c1=0,c2=0; int bh=256;
        for(int t=0;t<nthreads;t++){c1+=args[t].cnt1;c2+=args[t].cnt2;if(args[t].best_hw<bh)bh=args[t].best_hw;}
        printf("[%3ds] hw59<%d: %ld/%ld  hw59<%d: %ld/%ld  best_hw=%d\n",
               (sec+1)*10,HW_T1,c1,target1,HW_T2,c2,target2,bh);
        fflush(stdout);
        if(c1>=target1 && c2>=target2) break;
        if(sec>60) break;  /* не более 610 сек */
    }

    for(int t=0;t<nthreads;t++) pthread_cancel(th[t]);
    for(int t=0;t<nthreads;t++) pthread_join(th[t],NULL);
    printf("Время: %ld сек\n",(long)(time(NULL)-t0));

    /* Объединяем аккумуляторы */
    DSCAcc tot1, tot2;
    acc_init(&tot1); acc_init(&tot2);
    for(int t=0;t<nthreads;t++){
        DSCAcc *a1=&args[t].acc1, *a2=&args[t].acc2;
        tot1.N+=a1->N; tot2.N+=a2->N;
        tot1.hw_sum+=a1->hw_sum; tot2.hw_sum+=a2->hw_sum;
        for(int w=0;w<8;w++)for(int b=0;b<32;b++){
            tot1.zero_final[w][b]+=a1->zero_final[w][b];
            tot2.zero_final[w][b]+=a2->zero_final[w][b];
        }
        for(int ri=0;ri<NCHECK;ri++)
            for(int w=0;w<8;w++)for(int b=0;b<32;b++){
                tot1.mean_r[ri][w][b]+=a1->mean_r[ri][w][b];
                tot2.mean_r[ri][w][b]+=a2->mean_r[ri][w][b];
            }
    }

    print_dsc(&tot1,"hw59<90",90);
    print_dsc(&tot2,"hw59<85",85);

    printf("\n╔═══════════════════════════════════════════════════╗\n");
    printf("║  ИТОГ: РАССТОЯНИЕ ДО КОЛЛИЗИИ                    ║\n");
    printf("╚═══════════════════════════════════════════════════╝\n");
    {
        long N=tot2.N; if(N<5)N=tot1.N;
        DSCAcc *a=(tot2.N>=5)?&tot2:&tot1;
        int thresh=(tot2.N>=5)?85:90;
        int nfixed=0;
        for(int w=0;w<8;w++)for(int b=0;b<32;b++){
            double p=(double)a->zero_final[w][b]/a->N;
            if(p>0.90)nfixed++;
        }
        printf("  При hw59<%d: фиксировано %d бит δ[64] из 256\n",thresh,nfixed);
        printf("  Нефиксировано: %d бит\n",256-nfixed);
        printf("  → Оставшаяся сложность: O(2^%d)\n",256-nfixed);
        if(256-nfixed<=64)
            printf("  ★ ЦЕЛЬ O(2^64) ДОСТИЖИМА!\n");
        if(256-nfixed<=40)
            printf("  ★★ ЦЕЛЬ O(2^40) ДОСТИЖИМА!\n");
        else
            printf("  До O(2^40): нужно зафиксировать ещё %d бит\n",(256-nfixed)-40);
    }

    free(args); free(th);
    return 0;
}
