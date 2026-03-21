/*
 * dw_deep.c — Глубокий SA для конкретных DW-кандидатов
 *
 * Берём лучшие DW из dw_random и ищем их истинный hw59_min
 * с длинным SA (500M iters × 8 threads)
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

static int eval_dw(const uint32_t DW[16], int nthreads, long iters, uint64_t seed, uint32_t Wout[16]){
    SAArg *args=calloc(nthreads,sizeof(SAArg));
    pthread_t *th=malloc(nthreads*sizeof(pthread_t));
    for(int t=0;t<nthreads;t++){
        args[t].tid=t;args[t].DW=DW;args[t].seed=seed^(uint64_t)t*0xabcdef;
        args[t].max_iters=iters;pthread_create(&th[t],NULL,sa_worker,&args[t]);
    }
    int best=256;
    for(int t=0;t<nthreads;t++){
        pthread_join(th[t],NULL);
        if(args[t].result_hw<best){best=args[t].result_hw;if(Wout)memcpy(Wout,args[t].result_W,64);}
    }
    free(args);free(th);return best;
}

/* Вычислить delta[64] для данного (W, DW) */
static int get_delta64(const uint32_t W0[16], const uint32_t DW[16], uint32_t delta[8]){
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
    int hw=0;
    for(int k=0;k<8;k++){delta[k]=s[0][k]^s[1][k];hw+=__builtin_popcount(delta[k]);}
    return hw;
}

int main(int argc,char**argv){
    setvbuf(stdout,NULL,_IONBF,0);  /* отключить буферизацию */

    int nthreads=8;
    long sa_iters=500000000L;
    if(argc>1) sa_iters=atol(argv[1]);
    if(argc>2) nthreads=atoi(argv[2]);

    printf("dw_deep: глубокий SA для лучших DW-кандидатов\n");
    printf("SA iters=%ldM, threads=%d\n\n",sa_iters/1000000,nthreads);

    /* Список DW кандидатов */
    typedef struct { const char *name; uint32_t dw[16]; } Candidate;

    Candidate cands[]={
        /* Wang classical */
        {"Wang_DW0=80000000",
         {0x80000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}},
        /* Лучший из random search: 2-бит DW с hw59=79 */
        {"Rand2_DW4+DW10",
         {0,0,0,0,0x08000000,0,0,0,0,0,0x00000002,0,0,0,0,0}},
        /* 3-бит DW с hw59=78 */
        {"Rand3_DW2+DW7+DW11",
         {0,0,0x00000001,0,0,0,0,0x01000000,0,0,0,0x00800000,0,0,0,0}},
        /* Проверяем: Wang + один доп. бит (лучший из fast_greedy) */
        {"Wang+DW14b25",
         {0x80000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0x02000000,0}},
        /* Попытка: ΔW на поздних индексах (меньше влияет на диффузию) */
        {"LateDW13b15",
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0x00008000,0,0}},
        {"LateDW14b0",
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x00000001,0}},
        {"LateDW15b0",
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x00000001}},
    };
    int ncands=sizeof(cands)/sizeof(cands[0]);

    printf("Тестируем %d кандидатов:\n\n",ncands);

    int best_hw=256;
    const char *best_name="none";
    int best_dw_idx=-1;

    for(int c=0;c<ncands;c++){
        int hwdw=0;
        for(int w=0;w<16;w++) hwdw+=__builtin_popcount(cands[c].dw[w]);
        printf("=== [%d/%d] %s (HW(DW)=%d) ===\n",c+1,ncands,cands[c].name,hwdw);

        uint32_t Wopt[16];
        int hw=eval_dw(cands[c].dw,nthreads,sa_iters,(uint64_t)time(NULL)^(uint64_t)c*0x12345,Wopt);
        printf("  hw59_min = %d\n",hw);

        if(hw<best_hw){best_hw=hw;best_name=cands[c].name;best_dw_idx=c;}

        /* Показать delta[64] для лучшего W */
        uint32_t delta[8];
        int hw2=get_delta64(Wopt,cands[c].dw,delta);
        printf("  delta[64]: %08x %08x %08x %08x  %08x %08x %08x %08x  HW=%d\n",
               delta[0],delta[1],delta[2],delta[3],delta[4],delta[5],delta[6],delta[7],hw2);

        printf("\n");
    }

    printf("╔═══════════════════════════════════════╗\n");
    printf("║  ИТОГ DEEP SEARCH                     ║\n");
    printf("╚═══════════════════════════════════════╝\n");
    printf("  Лучший: %s  hw59_min=%d\n",best_name,best_hw);
    if(best_hw<40)
        printf("  ★★ O(2^%d) ДОСТИГНУТА! ★★\n",best_hw);
    else if(best_hw<76)
        printf("  ★ Прорыв! hw59=%d < 76 (Wang Barrier!)\n",best_hw);
    else
        printf("  Wang Barrier: hw59_min=%d (все кандидаты ≥ 76)\n",best_hw);

    return 0;
}
