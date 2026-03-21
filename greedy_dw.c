/*
 * greedy_dw.c — Жадное построение оптимального ΔW
 *
 * Алгоритм:
 *   1. Начать с DW = ΔW[0]=0x80000000 (hw59_min=76)
 *   2. Перебрать каждый возможный флип бита DW
 *   3. Для каждого: запустить SA для W и измерить hw59_min
 *   4. Взять флип с наименьшим hw59_min
 *   5. Повторить N раз
 *
 * Гипотеза: последовательное добавление бит в DW
 * снижает hw59_min: 76 → 70 → 60 → ... → 0?
 *
 * Q-G1: На сколько снижается hw59_min при добавлении каждого бита?
 * Q-G2: Сколько бит нужно добавить для hw59_min < 40?
 * Q-G3: Какова структура оптимального DW?
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

/* Параллельный SA для фиксированного DW */
typedef struct {
    int tid;
    const uint32_t *DW;
    uint64_t seed;
    long max_iters;
    int result_hw;
    uint32_t result_W[16];
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
        T*=0.9999975; if(T<0.3)T=0.3;
        if(it%300000==0){memcpy(W,bW,64);T=40.+(double)(xr()%30);}
        if(it%5000000==0&&bh>92){for(int i=0;i<16;i++)W[i]=(uint32_t)xr();bh=256;T=60.;}
    }
    a->result_hw=bh;memcpy(a->result_W,bW,64);
    return NULL;
}

/* Оценить hw59_min для данного DW (параллельный SA) */
static int eval_dw(const uint32_t DW[16], int nthreads, long iters, uint64_t seed){
    SAArg *args=calloc(nthreads,sizeof(SAArg));
    pthread_t *th=malloc(nthreads*sizeof(pthread_t));
    for(int t=0;t<nthreads;t++){
        args[t].tid=t; args[t].DW=DW; args[t].seed=seed^(uint64_t)t*0xabcd;
        args[t].max_iters=iters;
        pthread_create(&th[t],NULL,sa_worker,&args[t]);
    }
    int best=256;
    for(int t=0;t<nthreads;t++){pthread_join(th[t],NULL);if(args[t].result_hw<best)best=args[t].result_hw;}
    free(args); free(th);
    return best;
}

int main(int argc,char**argv){
    int nthreads=6;
    long eval_iters=50000000L;  /* 50M итераций для оценки */
    int greedy_steps=20;         /* число жадных шагов */
    if(argc>1) eval_iters=atol(argv[1]);
    if(argc>2) greedy_steps=atoi(argv[2]);
    if(argc>3) nthreads=atoi(argv[3]);

    printf("greedy_dw: жадное построение ΔW\n");
    printf("eval_iters=%ldM, steps=%d, nthreads=%d\n\n",eval_iters/1000000,greedy_steps,nthreads);
    printf("Q-G1: Снижение hw59_min при добавлении каждого бита?\n");
    printf("Q-G2: Сколько бит нужно для hw59_min < 40?\n\n");
    fflush(stdout);

    /* Начальный DW: Wang classical */
    uint32_t DW[16]={0}; DW[0]=0x80000000u;
    printf("=== Шаг 0: ΔW[0]=0x80000000 ===\n");
    int hw_min=eval_dw(DW,nthreads,eval_iters,(uint64_t)time(NULL));
    printf("  hw59_min = %d  HW(DW) = 1\n\n",hw_min);
    fflush(stdout);

    int current_min=hw_min;

    /* Жадные шаги */
    for(int step=1;step<=greedy_steps;step++){
        printf("=== Шаг %d: поиск лучшего флипа бита DW ===\n",step);
        fflush(stdout);

        int best_flip_w=-1, best_flip_b=-1, best_hw_after=current_min;

        for(int w=0;w<16;w++){
            for(int b=0;b<32;b++){
                uint32_t DW_new[16]; memcpy(DW_new,DW,64);
                DW_new[w]^=(1u<<b);
                /* Проверяем что DW_new ≠ 0 */
                int allz=1; for(int i=0;i<16;i++) if(DW_new[i]){allz=0;break;}
                if(allz)continue;

                int hw=eval_dw(DW_new,nthreads,eval_iters/4,(uint64_t)time(NULL)^((uint64_t)(w*32+b)*0x12345));
                if(hw<best_hw_after){
                    best_hw_after=hw; best_flip_w=w; best_flip_b=b;
                    printf("  W[%2d][bit%2d]: hw59=%d  ← новый лучший\n",w,b,hw); fflush(stdout);
                }
            }
        }

        if(best_flip_w<0){
            printf("  Нет улучшения. Остановка.\n");
            break;
        }

        /* Применяем лучший флип */
        DW[best_flip_w]^=(1u<<best_flip_b);
        int hw_confirm=eval_dw(DW,nthreads,eval_iters,(uint64_t)time(NULL));
        printf("  Применяем W[%2d][bit%2d]: hw59_min %d → %d\n",
               best_flip_w,best_flip_b,current_min,hw_confirm);

        /* Текущий DW */
        printf("  Текущий DW: "); for(int w=0;w<16;w++)if(DW[w])printf("DW[%d]=0x%08x ",w,DW[w]);
        printf("\n  HW(DW)=%d\n\n",__builtin_popcount(DW[0])+__builtin_popcount(DW[1])+__builtin_popcount(DW[2])+__builtin_popcount(DW[3])+__builtin_popcount(DW[4])+__builtin_popcount(DW[5])+__builtin_popcount(DW[6])+__builtin_popcount(DW[7])+__builtin_popcount(DW[8])+__builtin_popcount(DW[9])+__builtin_popcount(DW[10])+__builtin_popcount(DW[11])+__builtin_popcount(DW[12])+__builtin_popcount(DW[13])+__builtin_popcount(DW[14])+__builtin_popcount(DW[15]));
        fflush(stdout);

        current_min=hw_confirm;
        if(current_min<40){
            printf("★★ hw59_min=%d < 40: цель O(2^40) достигнута! ★★\n",current_min);
            break;
        }
        if(current_min==0){
            printf("★★★ КОЛЛИЗИЯ SHA-256! ★★★\n"); break;
        }
    }

    printf("\n╔══════════════════════════════════════════════════╗\n");
    printf("║  ИТОГ GREEDY DW                                  ║\n");
    printf("╚══════════════════════════════════════════════════╝\n");
    int total_hw=0; for(int w=0;w<16;w++) total_hw+=__builtin_popcount(DW[w]);
    printf("  hw59_min = %d    HW(DW) = %d\n",current_min,total_hw);
    printf("  Финальный DW:\n");
    for(int w=0;w<16;w++) if(DW[w]) printf("    DW[%2d] = 0x%08x\n",w,DW[w]);

    if(current_min<40)
        printf("\n  ★ O(2^%d) достижима!\n",current_min);
    else
        printf("\n  До O(2^40): нужно снизить ещё на %d\n",current_min-40);

    return 0;
}
