/*
 * rank_sweep.c — rank(WCC_k) для k = 95, 90, 85, 80, 75
 *
 * Q188: rank(WCC_k) → 0 при k→0?
 * Q189: как ранг изменяется от k=95 до k=75?
 *
 * Метод: GF(2) Union-Find для каждого уровня, как в carry_rank.c
 *        но собираем образцы всех уровней за один прогон SA.
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
        uint32_t T1=h+S1(e)+CH(e,f,g)+K256[r]+W[j][r], T2=S0(a)+MAJ(a,b,c);
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

static uint64_t rng_s;
static inline uint64_t xr(){rng_s^=rng_s<<13;rng_s^=rng_s>>7;rng_s^=rng_s<<17;return rng_s;}

/* ===== Накопитель для одного уровня ===== */
#define MAXSAMP 2000
typedef struct {
    int thr;           /* порог hw59 */
    long target;       /* сколько образцов собрать */
    long n;
    double sx[64], sxy[64*64];
} Level;

#define NLEVELS 5
static int THRESHOLDS[NLEVELS] = {95, 90, 85, 80, 75};
static Level levels[NLEVELS];
static pthread_mutex_t lmu = PTHREAD_MUTEX_INITIALIZER;

/* GF(2) ранг через Union-Find */
typedef struct { int par[64], xr[64]; } UF;
static void uf_init(UF *u){ for(int i=0;i<64;i++){u->par[i]=i;u->xr[i]=0;} }
static int uf_find(UF *u, int i, int *px){
    if(u->par[i]==i){*px=0;return i;}
    int px2; int root=uf_find(u,u->par[i],&px2);
    u->xr[i]^=px2; u->par[i]=root; *px=u->xr[i]; return root;
}
static int uf_merge(UF *u, int i, int j, int sign){
    int xi,xj; int ri=uf_find(u,i,&xi), rj=uf_find(u,j,&xj);
    if(ri==rj) return (xi^xj)!=sign;
    u->par[ri]=rj; u->xr[ri]=xi^xj^sign; return 0;
}

static int compute_rank(Level *lv, double phi_thr, int *nconstr_out, int *ncontra_out){
    long N=lv->n; if(N<50){if(nconstr_out)*nconstr_out=0;if(ncontra_out)*ncontra_out=0;return -1;}
    UF uf; uf_init(&uf);
    int nconstr=0, ncontra=0;
    for(int i=0;i<64;i++) for(int j=i+1;j<64;j++){
        double mi=lv->sx[i]/N, mj=lv->sx[j]/N;
        double vi=mi*(1-mi), vj=mj*(1-mj);
        if(vi<1e-9||vj<1e-9) continue;
        double mij=lv->sxy[i*64+j]/N;
        double phi=(mij-mi*mj)/sqrt(vi*vj);
        if(fabs(phi)<phi_thr) continue;
        int sign=(phi<0)?1:0;
        if(uf_merge(&uf,i,j,sign)) ncontra++;
        nconstr++;
    }
    /* считаем компоненты */
    int gs[64]={0};
    for(int i=0;i<64;i++){int x;gs[uf_find(&uf,i,&x)]++;}
    int comp=0; for(int i=0;i<64;i++) if(gs[i]>0) comp++;
    /* вычитаем константные */
    int nconst=0;
    for(int i=0;i<64;i++){double mi=lv->sx[i]/N;if(mi<0.05||mi>0.95)nconst++;}
    if(nconstr_out)*nconstr_out=nconstr;
    if(ncontra_out)*ncontra_out=ncontra;
    return comp-nconst;
}

/* ===== Поток SA ===== */
typedef struct { int tid; long max_iters; uint64_t seed; } TArg;
static void *sa_worker(void *v){
    TArg *a=(TArg*)v;
    rng_s=a->seed^(uint64_t)a->tid*0xdeadbeef12LL;
    uint32_t W[16]; for(int i=0;i<16;i++) W[i]=(uint32_t)xr();
    double T=60.0; uint32_t bW[16]; int bh=256; memcpy(bW,W,64);
    long iter=0;
    while(iter<a->max_iters){
        iter++;
        int hw=hw59_diff(W);
        uint8_t carry[64]; sha256_carry(W,carry);

        pthread_mutex_lock(&lmu);
        for(int lv=0;lv<NLEVELS;lv++){
            if(hw<THRESHOLDS[lv] && levels[lv].n<levels[lv].target){
                levels[lv].n++;
                for(int i=0;i<64;i++){
                    levels[lv].sx[i]+=carry[i];
                    for(int j=i+1;j<64;j++)
                        levels[lv].sxy[i*64+j]+=carry[i]*carry[j];
                }
            }
        }
        pthread_mutex_unlock(&lmu);

        if(hw<bh){bh=hw;memcpy(bW,W,64);}
        int k=(int)(xr()%16),bit=(int)(xr()%32);
        uint32_t sv=W[k]; W[k]^=(1u<<bit);
        int nh=hw59_diff(W);
        double dE=nh-hw;
        if(dE>0&&(double)(xr()%1000000)/1e6>exp(-dE/T)) W[k]=sv;
        T*=0.9999975; if(T<0.3)T=0.3;
        if(iter%300000==0){memcpy(W,bW,64);T=40.0+(double)(xr()%30);}
        if(iter%5000000==0&&bh>92){for(int i=0;i<16;i++)W[i]=(uint32_t)xr();bh=256;T=60.0;}
    }
    return NULL;
}

int main(int argc, char **argv){
    int nthreads=6;
    long max_iters=200000000L;
    if(argc>1) max_iters=atol(argv[1]);
    if(argc>2) nthreads=atoi(argv[2]);

    printf("rank_sweep: итераций=%ld, потоков=%d\n\n",max_iters,nthreads);
    printf("Q188: rank(WCC_k) → 0 при k→0?\n");
    printf("Q189: профиль rank(WCC_k) по уровням\n\n");

    for(int lv=0;lv<NLEVELS;lv++){
        levels[lv].thr=THRESHOLDS[lv];
        levels[lv].target=1500;
        levels[lv].n=0;
        memset(levels[lv].sx,0,sizeof(levels[lv].sx));
        memset(levels[lv].sxy,0,sizeof(levels[lv].sxy));
    }

    time_t t0=time(NULL);
    TArg *args=calloc(nthreads,sizeof(TArg));
    pthread_t *th=malloc(nthreads*sizeof(pthread_t));
    for(int t=0;t<nthreads;t++){
        args[t].tid=t; args[t].max_iters=max_iters;
        args[t].seed=(uint64_t)time(NULL)^(uint64_t)t*0xabcdef12345LL;
        pthread_create(&th[t],NULL,sa_worker,&args[t]);
    }
    for(int t=0;t<nthreads;t++) pthread_join(th[t],NULL);
    printf("Время: %ld сек\n\n",(long)(time(NULL)-t0));

    printf("╔══════════════════════════════════════════════════╗\n");
    printf("║         РАНГ WCC_k ПО УРОВНЯМ                   ║\n");
    printf("╚══════════════════════════════════════════════════╝\n\n");
    printf("  %-10s | %-8s | %-8s | %-8s | %-8s | %s\n",
           "уровень","N","ранг","ограничений","противор.","эффект");
    printf("  -----------+---------+---------+-------------+----------+--------\n");

    for(int lv=NLEVELS-1;lv>=0;lv--){
        int nc,nct;
        int rk=compute_rank(&levels[lv],0.95,&nc,&nct);
        printf("  hw59<%3d   | %-8ld | %-8d | %-11d | %-8d | %s\n",
               THRESHOLDS[lv],levels[lv].n,rk,nc,nct,
               (rk<0)?"недост.":
               (rk<10)?"★★★ КРИТИЧНО":
               (rk<20)?"★★ НИЗКИЙ":
               (rk<30)?"★ УМЕРЕННЫЙ":"ВЫСОКИЙ");
    }

    printf("\n📐 Подробный анализ по уровням:\n");
    for(int lv=NLEVELS-1;lv>=0;lv--){
        int nc,nct;
        int rk=compute_rank(&levels[lv],0.95,&nc,&nct);
        if(rk<0){printf("\n  [hw59<%d]: недостаточно образцов (N=%ld)\n",THRESHOLDS[lv],levels[lv].n);continue;}
        int nconst=0;
        long N=levels[lv].n;
        for(int i=0;i<64;i++){double mi=levels[lv].sx[i]/N;if(mi<0.05||mi>0.95)nconst++;}
        printf("\n  [hw59<%d]: N=%ld ранг=%d ограничений=%d константных=%d противор=%d\n",
               THRESHOLDS[lv],N,rk,nc,nconst,nct);

        /* топ-5 групп */
        /* реконструируем UF */
        UF uf; uf_init(&uf);
        for(int i=0;i<64;i++) for(int j=i+1;j<64;j++){
            double mi=levels[lv].sx[i]/N,mj=levels[lv].sx[j]/N;
            double vi=mi*(1-mi),vj=mj*(1-mj);
            if(vi<1e-9||vj<1e-9) continue;
            double mij=levels[lv].sxy[i*64+j]/N;
            double phi=(mij-mi*mj)/sqrt(vi*vj);
            if(fabs(phi)<0.95) continue;
            uf_merge(&uf,i,j,(phi<0)?1:0);
        }
        int gs[64]={0};
        for(int i=0;i<64;i++){int x;gs[uf_find(&uf,i,&x)]++;}
        /* показать группы размером >= 3 */
        for(int root=0;root<64;root++){
            if(gs[root]<3) continue;
            printf("    Группа(%d):",gs[root]);
            for(int i=0;i<64;i++){
                int xi; if(uf_find(&uf,i,&xi)==root)
                    printf(" c[%d]%s",i,(xi)?"'":"");
            }
            printf("\n");
        }
    }

    /* Сводная таблица для построения графика */
    printf("\n📈 Данные для графика rank(WCC_k) vs k:\n");
    printf("  k  | rank | N\n");
    for(int lv=0;lv<NLEVELS;lv++){
        int rk=compute_rank(&levels[lv],0.95,NULL,NULL);
        printf("  %3d| %4d | %ld\n",THRESHOLDS[lv],rk,levels[lv].n);
    }

    free(args);free(th);
    return 0;
}
