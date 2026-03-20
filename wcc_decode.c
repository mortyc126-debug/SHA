/*
 * wcc_decode.c — WCC-guided SA
 *
 * Q191: Можно ли использовать WCC ограничения для ускорения SA?
 *
 * Идея: 
 *   1. Хардкодируем 71 детерминированное ограничение carry[i]⊕carry[j]=b
 *      (из carry_rank результатов)
 *   2. Добавляем к objective: f = hw59 + λ·(число нарушений WCC)
 *   3. Измеряем: WCC-guided SA быстрее находит hw59<85?
 *
 * Дополнительно: WCC-decode — дано carry-вектор → найти W
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
static const uint32_t IV[8]={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
#define ROTR(x,n)(((x)>>(n))|((x)<<(32-(n))))
#define S0(x)(ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define S1(x)(ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define s0(x)(ROTR(x,7)^ROTR(x,18)^((x)>>3))
#define s1(x)(ROTR(x,17)^ROTR(x,19)^((x)>>10))
#define CH(e,f,g)(((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c)(((a)&(b))^((a)&(c))^((b)&(c)))

static int hw59_diff(const uint32_t W0[16]){
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
        uint32_t a=s[j][0],b=s[j][1],c=s[j][2],d=s[j][3],e=s[j][4],f=s[j][5],g=s[j][6],h=s[j][7];
        uint32_t T1=h+S1(e)+CH(e,f,g)+K256[r]+W[j][r],T2=S0(a)+MAJ(a,b,c);
        s[j][7]=g;s[j][6]=f;s[j][5]=e;s[j][4]=d+T1;s[j][3]=c;s[j][2]=b;s[j][1]=a;s[j][0]=T1+T2;
    }
    int hw=0; for(int k=0;k<8;k++) hw+=__builtin_popcount(s[0][k]^s[1][k]);
    return hw;
}

static void sha256_carry(const uint32_t W0[16], uint8_t carry[64]){
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=W0[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for(int r=0;r<64;r++){
        uint32_t T1=h+S1(e)+CH(e,f,g)+K256[r]+W[r],T2=S0(a)+MAJ(a,b,c);
        uint64_t sum=(uint64_t)T1+T2; carry[r]=(sum>>32)&1;
        h=g;g=f;f=e;e=d+T1;d=c;c=b;b=a;a=(uint32_t)sum;
    }
}

/* ===== WCC ограничения (из carry_rank при hw59<85, N=1500) ===== */
typedef struct { int i, j, sign; /* carry[i]^sign == carry[j] */ } WCCConstr;
static WCCConstr WCC[] = {
    {1,18,1},{1,32,1},{1,47,0},
    {2,5,1},{2,9,1},{2,29,0},{2,44,0},{2,61,0},
    {3,25,0},
    {5,9,0},{5,29,1},{5,44,1},{5,61,1},
    {7,57,0},
    {8,10,0},{8,17,1},{8,33,1},{8,41,1},{8,54,1},
    {9,29,1},{9,44,1},{9,61,1},
    {10,17,1},{10,33,1},{10,41,1},{10,54,1},
    {11,30,1},{11,38,0},{11,63,0},
    {12,36,0},
    {13,34,0},{13,45,0},
    {14,52,0},{14,56,0},
    {15,39,1},{15,40,1},
    {16,20,1},{16,37,0},{16,49,0},
    {17,33,0},{17,41,0},{17,54,0},
    {18,32,0},{18,47,1},
    {19,48,1},
    {20,37,1},{20,49,1},
    {22,42,0},{22,51,0},
    {24,46,1},
    {26,55,0},
    {28,31,1},{28,59,1},
    {29,44,0},{29,61,0},
    {30,38,1},{30,63,1},
    {31,59,0},
    {32,47,1},
    {33,41,0},{33,54,0},
    {34,45,0},
    {37,49,0},
    {38,63,0},
    {39,40,0},
    {41,54,0},
    {42,51,0},
    {44,61,0},
    {50,62,1},
    {52,56,0},
    {58,60,0},
};
#define NWCC (sizeof(WCC)/sizeof(WCC[0]))

static int wcc_violations(const uint8_t carry[64]){
    int v=0;
    for(int k=0;k<(int)NWCC;k++){
        int ci=carry[WCC[k].i], cj=carry[WCC[k].j];
        if((ci^cj)!=WCC[k].sign) v++;
    }
    return v;
}

static uint64_t rng_s;
static inline uint64_t xr(){rng_s^=rng_s<<13;rng_s^=rng_s>>7;rng_s^=rng_s<<17;return rng_s;}

/* ===== Сравнение: pure SA vs WCC-guided SA ===== */
typedef struct {
    int tid;
    int use_wcc;      /* 0=pure, 1=WCC-guided */
    double lambda;    /* вес WCC штрафа */
    long max_iters;
    uint64_t seed;
    /* результаты */
    int best_hw;
    long iters_to85;  /* итераций до hw59<85 (-1 если не найдено) */
    long iters_to90;
    long cnt85, cnt90;
} SAArg;

static void *sa_worker(void *v){
    SAArg *a=(SAArg*)v;
    rng_s=a->seed;
    uint32_t W[16],bestW[16];
    for(int i=0;i<16;i++) W[i]=(uint32_t)xr();
    memcpy(bestW,W,64);
    int besthw=256;
    double T=60.0;
    a->iters_to85=-1; a->iters_to90=-1;
    a->cnt85=0; a->cnt90=0;

    for(long iter=0;iter<a->max_iters;iter++){
        int hw=hw59_diff(W);
        double f=hw;
        if(a->use_wcc){
            uint8_t carry[64]; sha256_carry(W,carry);
            f+=a->lambda*wcc_violations(carry);
        }

        if(hw<90){a->cnt90++;if(a->iters_to90<0)a->iters_to90=iter;}
        if(hw<85){a->cnt85++;if(a->iters_to85<0)a->iters_to85=iter;}
        if(hw<besthw){besthw=hw;memcpy(bestW,W,64);}

        int k=(int)(xr()%16),bit=(int)(xr()%32);
        uint32_t sv=W[k]; W[k]^=(1u<<bit);
        int nh=hw59_diff(W);
        double nf=nh;
        if(a->use_wcc){
            uint8_t carry2[64]; sha256_carry(W,carry2);
            nf+=a->lambda*wcc_violations(carry2);
        }
        double dE=nf-f;
        if(dE>0&&(double)(xr()%1000000)/1e6>exp(-dE/T)) W[k]=sv;
        T*=0.9999975; if(T<0.3)T=0.3;
        if(iter%300000==0){memcpy(W,bestW,64);T=40.0+(double)(xr()%25);}
        if(iter%5000000==0&&besthw>92){for(int i=0;i<16;i++)W[i]=(uint32_t)xr();besthw=256;T=60.0;}
    }
    a->best_hw=besthw;
    return NULL;
}

int main(int argc, char **argv){
    int ntrials=4;    /* пар pure vs WCC */
    long max_iters=50000000L;
    if(argc>1) max_iters=atol(argv[1]);
    if(argc>2) ntrials=atoi(argv[2]);

    printf("wcc_decode: WCC-guided SA vs pure SA\n");
    printf("WCC ограничений: %zu (из carry_rank при hw59<85)\n\n",NWCC);
    printf("Q191: Ускоряет ли WCC штраф поиск near-collision?\n\n");

    printf("  %-10s | %-12s | %-12s | %-8s | %-8s | %s\n",
           "режим","iters→hw90","iters→hw85","cnt90","cnt85","best_hw");
    printf("  -----------+-------------+-------------+---------+---------+-------\n");

    time_t t0=time(NULL);
    /* ntrials пар: pure + wcc, запускаем по 2 потока на пару */
    int nthreads=ntrials*2;
    SAArg *args=calloc(nthreads,sizeof(SAArg));
    pthread_t *th=malloc(nthreads*sizeof(pthread_t));

    for(int t=0;t<nthreads;t++){
        args[t].tid=t;
        args[t].use_wcc=t%2;
        args[t].lambda=2.0;
        args[t].max_iters=max_iters;
        args[t].seed=(uint64_t)time(NULL)^(uint64_t)t*0x123456789abcLL;
        pthread_create(&th[t],NULL,sa_worker,&args[t]);
    }
    for(int t=0;t<nthreads;t++) pthread_join(th[t],NULL);
    printf("Время: %ld сек\n\n",(long)(time(NULL)-t0));

    long pure_to90=0,wcc_to90=0,pure_to85=0,wcc_to85=0;
    int pure_n=0,wcc_n=0;
    for(int t=0;t<nthreads;t++){
        const char *mode=args[t].use_wcc?"WCC-guided":"pure SA  ";
        printf("  %-10s | %-12ld | %-12ld | %-8ld | %-8ld | %d\n",
               mode,args[t].iters_to90,args[t].iters_to85,
               args[t].cnt90,args[t].cnt85,args[t].best_hw);
        if(args[t].use_wcc){
            if(args[t].iters_to90>0) wcc_to90+=args[t].iters_to90;
            if(args[t].iters_to85>0) wcc_to85+=args[t].iters_to85;
            wcc_n++;
        } else {
            if(args[t].iters_to90>0) pure_to90+=args[t].iters_to90;
            if(args[t].iters_to85>0) pure_to85+=args[t].iters_to85;
            pure_n++;
        }
    }

    printf("\n=== Итог: WCC-guided vs pure SA ===\n");
    printf("  Среднее iters→hw90:  pure=%-10ld  WCC=%-10ld  ускорение=%.2fx\n",
           pure_n?pure_to90/pure_n:-1, wcc_n?wcc_to90/wcc_n:-1,
           (pure_n&&wcc_n&&wcc_to90>0)?(double)(pure_to90/pure_n)/(wcc_to90/wcc_n):0.0);
    printf("  Среднее iters→hw85:  pure=%-10ld  WCC=%-10ld  ускорение=%.2fx\n",
           pure_n?pure_to85/pure_n:-1, wcc_n?wcc_to85/wcc_n:-1,
           (pure_n&&wcc_n&&wcc_to85>0)?(double)(pure_to85/pure_n)/(wcc_to85/wcc_n):0.0);

    printf("\n=== Теоретический анализ WCC decode ===\n");
    printf("Дано: 71 ограничение carry[i]⊕carry[j]=b\n");
    printf("Ранг: 24 свободных параметра из 64\n");
    printf("Подход Q191: carry-вектор C определяет 40 битов остальных через линейные уравнения\n");
    printf("Если дать SA 24 'свободных carry параметра' как подсказку → направленный поиск\n");
    printf("Статус: прямой декодер требует обратного SHA-256 — NP-hard\n");
    printf("Обходной путь: WCC как soft constraint (λ·violations) — измерено выше\n");

    free(args);free(th);
    return 0;
}
