/*
 * fast_greedy.c — Быстрый жадный поиск оптимального ΔW
 *
 * КЛЮЧЕВАЯ ИДЕЯ: не нужно долгий SA для скрининга DW.
 * Достаточно вычислить delta[64](W*, DW') для фиксированного W*.
 *
 * Алгоритм:
 *   1. Найти W* с hw59(W*, DW₀)=79 (C₇₉)
 *   2. Для каждого из 512 возможных дополнительных бит:
 *      DW' = DW₀ XOR flip[k] → вычислить hw59(W*, DW')
 *   3. Выбрать flip с минимальным hw59
 *   4. Применить → новый DW, запустить SA для нового W*
 *   5. Повторить
 *
 * Скорость: 512 SHA-256 вычислений ≈ 0.001 сек (молниеносно!)
 * Затем SA для нового DW: ~60 секунд
 *
 * Q-G4: Быстрый путь к hw59_min < 40?
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

/* Параллельный SA для данного DW, возвращает {best_hw, best_W} */
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
static int find_Wstar(const uint32_t DW[16], int nthreads, long iters, uint32_t Wout[16]){
    SAArg *args=calloc(nthreads,sizeof(SAArg));
    pthread_t *th=malloc(nthreads*sizeof(pthread_t));
    for(int t=0;t<nthreads;t++){
        args[t].tid=t;args[t].DW=DW;args[t].seed=(uint64_t)time(NULL)^(uint64_t)t*0xabcdef;
        args[t].max_iters=iters;pthread_create(&th[t],NULL,sa_worker,&args[t]);
    }
    int best=256;
    for(int t=0;t<nthreads;t++){pthread_join(th[t],NULL);if(args[t].result_hw<best){best=args[t].result_hw;memcpy(Wout,args[t].result_W,64);}}
    free(args);free(th);return best;
}

/* Быстрый скрининг: для данного W* и DW₀, найти 20 лучших флипов */
static void screen_flips(const uint32_t Wstar[16], const uint32_t DW0[16],
                          int top_k, int result_flips[], int result_hw[]){
    int hw_flip[512];
    for(int di=0;di<512;di++){
        int w=di/32, b=di%32;
        uint32_t DW_new[16]; memcpy(DW_new,DW0,64);
        DW_new[w]^=(1u<<b);
        /* Пропустить DW=0 (тривиальный случай) */
        int allz=1; for(int i=0;i<16;i++) if(DW_new[i]){allz=0;break;}
        if(allz){hw_flip[di]=9999;continue;}
        hw_flip[di]=hw59_dw(Wstar,DW_new);
    }
    /* Найти top_k минимальных */
    for(int k=0;k<top_k;k++){
        int best_di=0;
        for(int di=1;di<512;di++) if(hw_flip[di]<hw_flip[best_di])best_di=di;
        result_flips[k]=best_di;
        result_hw[k]=hw_flip[best_di];
        hw_flip[best_di]=999;  /* убрать из рассмотрения */
    }
}

int main(int argc,char**argv){
    int nthreads=6;
    long sa_iters=100000000L;   /* 100M для поиска W* */
    int greedy_steps=50;
    int top_candidates=5;       /* топ-5 флипов на каждом шаге */
    if(argc>1) sa_iters=atol(argv[1]);
    if(argc>2) greedy_steps=atoi(argv[2]);
    if(argc>3) nthreads=atoi(argv[3]);

    printf("fast_greedy: быстрый жадный поиск оптимального ΔW\n");
    printf("SA iters=%ldM per DW eval, steps=%d, threads=%d\n\n",sa_iters/1000000,greedy_steps,nthreads);
    printf("Фаза 1 (скрининг): вычислить hw59(W*, DW') для 512 флипов → мгновенно\n");
    printf("Фаза 2 (SA): найти новый W* для лучшего DW → %ldM итераций\n\n",sa_iters/1000000);
    fflush(stdout);

    /* Начальный DW: Wang classical */
    uint32_t DW[16]={0}; DW[0]=0x80000000u;
    uint32_t Wstar[16]={0};

    printf("=== Инициализация: найти W* для ΔW[0]=0x80000000 ===\n");
    int hw_current=find_Wstar(DW,nthreads,sa_iters,Wstar);
    printf("  W* найден: hw59=%d\n", hw_current);
    printf("  W*: "); for(int i=0;i<16;i++)printf("%08x",Wstar[i]); printf("\n\n");
    fflush(stdout);

    /* Жадные шаги */
    int hw_history[100]={0}; int dw_hw_history[100]={0};
    hw_history[0]=hw_current;
    int total_hw_dw=1;
    dw_hw_history[0]=total_hw_dw;

    for(int step=1;step<=greedy_steps;step++){
        printf("=== Шаг %d (hw59=%d, HW(DW)=%d) ===\n",step,hw_current,total_hw_dw);
        fflush(stdout);

        /* Фаза 1: быстрый скрининг всех 512 флипов */
        int flips[20], flip_hw[20];
        screen_flips(Wstar, DW, 20, flips, flip_hw);

        printf("  Топ-10 флипов по hw59(W*, DW'):\n");
        for(int k=0;k<10;k++){
            int w=flips[k]/32, b=flips[k]%32;
            printf("    W[%2d][bit%2d]: hw59=%d\n",w,b,flip_hw[k]);
        }
        printf("  Применяем лучший флип: W[%d][bit%d] (hw59=%d → %d)\n",
               flips[0]/32,flips[0]%32,hw_current,flip_hw[0]);
        fflush(stdout);

        if(flip_hw[0]>=hw_current){
            printf("  Нет улучшения по скринингу. Пробуем SA для топ-5...\n");
            fflush(stdout);
            /* Попробуем SA для топ-5 кандидатов */
            int best_sa_hw=hw_current;
            int best_flip=-1;
            for(int k=0;k<5;k++){
                int w=flips[k]/32,b=flips[k]%32;
                uint32_t DW_try[16]; memcpy(DW_try,DW,64);
                DW_try[w]^=(1u<<b);
                uint32_t W_try[16];
                int hw_sa=find_Wstar(DW_try,nthreads,sa_iters/5,W_try);
                printf("    SA для W[%d][bit%d]: hw59=%d\n",w,b,hw_sa); fflush(stdout);
                if(hw_sa<best_sa_hw){best_sa_hw=hw_sa;best_flip=k;memcpy(Wstar,W_try,64);}
            }
            if(best_flip<0){printf("  SA тоже не улучшает. Останавливаемся.\n");break;}
            DW[flips[best_flip]/32]^=(1u<<(flips[best_flip]%32));
            hw_current=best_sa_hw;
        } else {
            /* Применяем лучший флип скрининга и запускаем SA */
            DW[flips[0]/32]^=(1u<<(flips[0]%32));
            /* Запускаем SA для нового DW */
            uint32_t W_new[16];
            int hw_sa=find_Wstar(DW,nthreads,sa_iters,W_new);
            printf("  SA для нового DW: hw59=%d\n",hw_sa); fflush(stdout);
            if(hw_sa<hw_current){
                hw_current=hw_sa;
                memcpy(Wstar,W_new,64);
            } else {
                /* Откат флипа */
                DW[flips[0]/32]^=(1u<<(flips[0]%32));
                printf("  SA не улучшила, откат.\n");
                break;
            }
        }

        total_hw_dw=0; for(int w=0;w<16;w++) total_hw_dw+=__builtin_popcount(DW[w]);
        hw_history[step]=hw_current;
        dw_hw_history[step]=total_hw_dw;

        printf("  Результат: hw59_min=%d  HW(DW)=%d\n",hw_current,total_hw_dw);
        printf("  DW: "); for(int w=0;w<16;w++)if(DW[w])printf("DW[%d]=0x%08x ",w,DW[w]);
        printf("\n\n"); fflush(stdout);

        if(hw_current==0){printf("★★★ КОЛЛИЗИЯ! ★★★\n");break;}
        if(hw_current<40){printf("★★ hw59=%d < 40: O(2^%d) достигнута! ★★\n",hw_current,hw_current);break;}
    }

    printf("\n╔══════════════════════════════════════════════════╗\n");
    printf("║  ИТОГ FAST GREEDY DW                            ║\n");
    printf("╚══════════════════════════════════════════════════╝\n");
    printf("  Шаги: "); for(int s=0;s<=greedy_steps&&hw_history[s]>0;s++) printf("%d ",hw_history[s]);
    printf("\n  HW(DW): "); for(int s=0;s<=greedy_steps&&hw_history[s]>0;s++) printf("%d ",dw_hw_history[s]);
    printf("\n  Финал: hw59_min=%d  HW(DW)=%d\n",hw_current,total_hw_dw);
    if(hw_current<40) printf("  ★ O(2^%d) достижима!\n",hw_current);
    else printf("  До O(2^40): нужно снизить ещё %d\n",hw_current-40);

    return 0;
}
