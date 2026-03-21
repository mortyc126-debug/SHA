/*
 * cstar_fast.c — Быстро найти C* и проверить уникальность
 *
 * Использует точно тот же SA что rank_sweep для быстрого сбора образцов hw59<80
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

static int sha256_diff(const uint32_t W0[16], uint32_t diff64[8]){
    uint32_t W[2][64];
    for(int i=0;i<16;i++){W[0][i]=W0[i];W[1][i]=W0[i];}
    W[1][0]^=0x80000000u;
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
    for(int k=0;k<8;k++){uint32_t d=s[0][k]^s[1][k];if(diff64)diff64[k]=d;hw+=__builtin_popcount(d);}
    return hw;
}

static int hw59_fast(const uint32_t W0[16]){return sha256_diff(W0,NULL);}

static __thread uint64_t rng_s;
static inline uint64_t xr(){rng_s^=rng_s<<13;rng_s^=rng_s>>7;rng_s^=rng_s<<17;return rng_s;}

/* Shared результаты */
static pthread_mutex_t mu=PTHREAD_MUTEX_INITIALIZER;
#define MAX_SAMPLES 50
static uint32_t samples_diff[MAX_SAMPLES][8];  /* delta[64] для каждого образца */
static int samples_hw[MAX_SAMPLES];
static int n_samples=0;
static int target_samples=20;

typedef struct {
    int tid; uint64_t seed; long max_iters;
} TArg;

static void *sa_worker(void *v){
    TArg *a=(TArg*)v;
    rng_s=a->seed^(uint64_t)a->tid*0xdeadbeef12LL;
    uint32_t W[16]; for(int i=0;i<16;i++)W[i]=(uint32_t)xr();
    double T=60.0; uint32_t bW[16]; int bh=256; memcpy(bW,W,64);

    for(long it=0;it<a->max_iters;it++){
        int hw=hw59_fast(W);
        if(hw<bh){bh=hw;memcpy(bW,W,64);}
        if(hw<80){
            uint32_t diff64[8];
            sha256_diff(W,diff64);
            pthread_mutex_lock(&mu);
            if(n_samples<target_samples){
                memcpy(samples_diff[n_samples],diff64,32);
                samples_hw[n_samples]=hw;
                n_samples++;
                printf("  Образец %d: hw59=%d  delta[64]=[%08x %08x %08x %08x...]\n",
                       n_samples,hw,diff64[0],diff64[1],diff64[2],diff64[3]);
                fflush(stdout);
            }
            int done=(n_samples>=target_samples);
            pthread_mutex_unlock(&mu);
            if(done)break;
        }
        int k=(int)(xr()%16),bit=(int)(xr()%32);
        uint32_t sv=W[k]; W[k]^=(1u<<bit);
        int nh=hw59_fast(W); double dE=nh-hw;
        if(dE<=0||(double)(xr()%1000000)/1e6<exp(-dE/T)){}else W[k]=sv;
        T*=0.9999975;if(T<0.3)T=0.3;
        if(it%300000==0){memcpy(W,bW,64);T=40.+(double)(xr()%30);}
        if(it%5000000==0&&bh>92){for(int i=0;i<16;i++)W[i]=(uint32_t)xr();bh=256;T=60.;}
    }
    return NULL;
}

int main(int argc,char**argv){
    int nthreads=6;
    target_samples=20;
    if(argc>1)target_samples=atoi(argv[1]);
    if(argc>2)nthreads=atoi(argv[2]);

    printf("cstar_fast: ищем C* (delta[64] при hw59<80)\n");
    printf("Сбор %d образцов hw59<80, потоков=%d\n\n",target_samples,nthreads);
    fflush(stdout);

    TArg *args=calloc(nthreads,sizeof(TArg));
    pthread_t *th=malloc(nthreads*sizeof(pthread_t));
    for(int t=0;t<nthreads;t++){
        args[t].tid=t;
        args[t].seed=(uint64_t)time(NULL)^(uint64_t)t*0xc0ffee123;
        args[t].max_iters=500000000L;
        pthread_create(&th[t],NULL,sa_worker,&args[t]);
    }
    for(int t=0;t<nthreads;t++) pthread_join(th[t],NULL);

    printf("\n=== АНАЛИЗ C* ===\n");
    printf("Собрано %d образцов hw59<80\n\n",n_samples);

    if(n_samples==0){printf("Нет образцов!\n");return 1;}

    /* Проверяем уникальность */
    int all_same=1;
    uint32_t *ref=samples_diff[0];
    for(int i=1;i<n_samples;i++){
        for(int w=0;w<8;w++){
            if(samples_diff[i][w]!=ref[w]){all_same=0;break;}
        }
        if(!all_same)break;
    }

    const char *wn[8]={"a","b","c","d","e","f","g","h"};
    printf("C* = delta[64] (первый образец hw59=%d):\n",samples_hw[0]);
    for(int w=0;w<8;w++)
        printf("  %s: 0x%08x  (HW=%d)\n",wn[w],ref[w],__builtin_popcount(ref[w]));
    printf("  HW(C*) = %d\n\n",samples_hw[0]);

    if(all_same && n_samples>1){
        printf("T_DSC_UNIQUE ПОДТВЕРЖДЕНА: все %d образцов имеют ОДИНАКОВЫЙ delta[64]!\n",n_samples);
        printf("C* = уникальный дифференциальный вывод при hw59<80\n\n");
    } else {
        printf("ВНИМАНИЕ: образцы с РАЗНЫМИ delta[64]!\n");
        for(int i=0;i<n_samples;i++){
            printf("  образец %d (hw59=%d): ",i,samples_hw[i]);
            for(int w=0;w<8;w++) printf("%08x ",samples_diff[i][w]);
            printf("\n");
        }
    }

    /* Точный список единиц в C* */
    int total_ones=0;
    printf("Биты = 1 в C*:\n");
    for(int w=0;w<8;w++){
        int cnt=0;
        for(int b=0;b<32;b++) if((ref[w]>>b)&1)cnt++;
        printf("  %s (0x%08x, HW=%d): биты ",wn[w],ref[w],cnt);
        for(int b=0;b<32;b++) if((ref[w]>>b)&1)printf("%d ",b);
        printf("\n");
        total_ones+=cnt;
    }
    printf("Итого единиц в C*: %d\n",total_ones);

    free(args);free(th);
    return 0;
}
