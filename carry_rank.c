/*
 * carry_rank.c — ранг carry-вектора при hw59<85
 * 
 * Q183: Полный набор детерминированных ограничений при hw59<85
 *       carry[i] = f(carry[j]) для всех пар с |φ|≈1
 *       → линейная система над GF(2)
 *       → ранг = число свободных carry битов
 *
 * Q184 (вторичный): Примерный ранг с теорией
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
    for(int i=0;i<16;i++){W[0][i]=W0[i];W[1][i]=W1[i];}
    for(int i=16;i<64;i++){
        W[0][i]=s1(W[0][i-2])+W[0][i-7]+s0(W[0][i-15])+W[0][i-16];
        W[1][i]=s1(W[1][i-2])+W[1][i-7]+s0(W[1][i-15])+W[1][i-16];
    }
    uint32_t s[2][8];
    for(int j=0;j<2;j++) for(int k=0;k<8;k++) s[j][k]=IV[k];
    for(int r=0;r<64;r++){
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

static void sha256_carry(const uint32_t W0[16], uint8_t carry[64]){
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=W0[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for(int r=0;r<64;r++){
        uint32_t T1=h+S1(e)+CH(e,f,g)+K256[r]+W[r];
        uint32_t T2=S0(a)+MAJ(a,b,c);
        uint64_t sum=(uint64_t)T1+T2;
        carry[r]=(sum>>32)&1;
        h=g;g=f;f=e;e=d+T1;d=c;c=b;b=a;a=(uint32_t)sum;
    }
}

static uint64_t rng_s;
static inline uint64_t xr(void){rng_s^=rng_s<<13;rng_s^=rng_s>>7;rng_s^=rng_s<<17;return rng_s;}

/* ===== Collect carry vectors at hw59<85 ===== */
#define MAXSAMP 4000
static uint8_t samples[MAXSAMP][64];
static int nsamp=0;
static pthread_mutex_t mu=PTHREAD_MUTEX_INITIALIZER;

typedef struct{int tid;long target;uint64_t seed;}CollArg;
static void *coll_worker(void *v){
    CollArg *a=(CollArg*)v;
    rng_s=a->seed^(uint64_t)a->tid*0xdeadbeef12345LL;
    uint32_t W[16],W1[16];
    for(int i=0;i<16;i++) W[i]=(uint32_t)xr();
    double T=60.0;
    uint32_t bestW[16]; int besthw=256;
    memcpy(bestW,W,64);
    long iters=0;
    for(;;){
        pthread_mutex_lock(&mu);
        if(nsamp>=a->target){pthread_mutex_unlock(&mu);break;}
        pthread_mutex_unlock(&mu);
        iters++;
        W1[0]=W[0]^0x80000000u;
        for(int i=1;i<16;i++) W1[i]=W[i];
        int hw=sha256_hw59(W,W1);
        if(hw<85){
            uint8_t carry[64];
            sha256_carry(W,carry);
            pthread_mutex_lock(&mu);
            if(nsamp<a->target){
                memcpy(samples[nsamp],carry,64);
                nsamp++;
            }
            pthread_mutex_unlock(&mu);
        }
        if(hw<besthw){besthw=hw;memcpy(bestW,W,64);}
        int k=(int)(xr()%16),bit=(int)(xr()%32);
        uint32_t saved=W[k];
        W[k]^=(1u<<bit);
        W1[0]=W[0]^0x80000000u;
        int nhw=sha256_hw59(W,W1);
        double dE=nhw-hw;
        if(dE>0 && (double)(xr()%1000000)/1e6>exp(-dE/T)) W[k]=saved;
        T*=0.9999975; if(T<0.3) T=0.3;
        if(iters%300000==0){memcpy(W,bestW,64);T=40.0+(double)(xr()%30);}
        if(iters%5000000==0&&besthw>92){for(int i=0;i<16;i++)W[i]=(uint32_t)xr();besthw=256;T=60.0;}
    }
    return NULL;
}

/* ===== GF(2) ранг — Union-Find для группировки ===== */
static int uf_parent[64], uf_xor[64]; /* xor[i] = XOR пути к корню */
static void uf_init(void){
    for(int i=0;i<64;i++){uf_parent[i]=i;uf_xor[i]=0;}
}
static int uf_find(int i, int *px){
    if(uf_parent[i]==i){*px=0;return i;}
    int px2; int root=uf_find(uf_parent[i],&px2);
    uf_xor[i]^=px2; uf_parent[i]=root;
    *px=uf_xor[i]; return root;
}
/* merge: установить carry[i] = sign⊕carry[j] (sign=0: равны, sign=1: анти) */
/* returns 0 if consistent, 1 if CONTRADICTION */
static int uf_merge(int i, int j, int sign){
    int xi,xj;
    int ri=uf_find(i,&xi), rj=uf_find(j,&xj);
    if(ri==rj){
        /* проверяем: xi⊕xj должно == sign */
        return (xi^xj)!=sign; /* противоречие если не совпадает */
    }
    /* merge меньший в больший */
    uf_parent[ri]=rj;
    uf_xor[ri]=xi^xj^sign;
    return 0; /* успех */
}

int main(int argc, char **argv){
    int nthreads=6;
    long target=2000;
    if(argc>1) target=atol(argv[1]);
    if(argc>2) nthreads=atoi(argv[2]);
    if(target>MAXSAMP) target=MAXSAMP;

    printf("carry_rank: собираем %ld образцов hw59<85, потоков=%d\n\n",target,nthreads);
    printf("Q183: Полный набор детерминированных ограничений carry[i]⊕carry[j]=const\n");
    printf("      → ранг carry-вектора над GF(2)\n\n");

    time_t t0=time(NULL);
    CollArg *args=calloc(nthreads,sizeof(CollArg));
    pthread_t *th=malloc(nthreads*sizeof(pthread_t));
    for(int t=0;t<nthreads;t++){
        args[t].tid=t; args[t].target=target;
        args[t].seed=(uint64_t)time(NULL)^(uint64_t)t*0xabcdef12345LL;
        pthread_create(&th[t],NULL,coll_worker,&args[t]);
    }
    for(int t=0;t<nthreads;t++) pthread_join(th[t],NULL);
    printf("Время сбора: %ld сек, образцов: %d\n\n",(long)(time(NULL)-t0),nsamp);

    int N=nsamp;
    if(N<10){printf("Недостаточно образцов.\n");return 1;}

    /* Вычислим phi для всех пар */
    double sx[64]={0}, sxy[64*64];
    memset(sxy,0,sizeof(sxy));
    for(int s=0;s<N;s++){
        for(int i=0;i<64;i++){
            sx[i]+=samples[s][i];
            for(int j=i+1;j<64;j++)
                sxy[i*64+j]+=samples[s][i]*samples[s][j];
        }
    }

    /* Найдём пары с |phi|>0.95 */
    printf("=== Детерминированные ограничения (|phi|>0.95) ===\n");
    uf_init();
    int nconstr=0, ncontradict=0;
    
    typedef struct{int i,j;double phi;}Constr;
    Constr constrs[64*64];
    int nc=0;

    for(int i=0;i<64;i++) for(int j=i+1;j<64;j++){
        double mi=sx[i]/N, mj=sx[j]/N;
        double vi=mi*(1-mi), vj=mj*(1-mj);
        if(vi<1e-9||vj<1e-9) continue;
        double mij=sxy[i*64+j]/N;
        double phi=(mij-mi*mj)/sqrt(vi*vj);
        if(fabs(phi)>0.95){
            constrs[nc].i=i; constrs[nc].j=j; constrs[nc].phi=phi; nc++;
            int sign=(phi<0)?1:0; /* φ<0: anti-correlated → XOR=1; φ>0: equal → XOR=0 */
            int contr=uf_merge(i,j,sign);
            nconstr++;
            if(contr) ncontradict++;
            printf("  c[%2d]%sc[%2d] lag=%2d: phi=%+.4f%s\n",
                   i,(phi<0)?" ⊕ ":" = ",j,j-i,phi,
                   contr?" [ПРОТИВОРЕЧИЕ!]":"");
        }
    }
    printf("\nВсего ограничений: %d, противоречий: %d\n\n",nconstr,ncontradict);

    /* Подсчёт компонент (групп) в Union-Find */
    int components=0;
    int group_size[64]={0};
    for(int i=0;i<64;i++){
        int xi; int root=uf_find(i,&xi);
        group_size[root]++;
    }
    for(int i=0;i<64;i++) if(group_size[i]>0) components++;
    int rank=components; /* число свободных carry bit групп */

    printf("=== Ранг carry-вектора (GF(2) анализ) ===\n");
    printf("Размерность {0,1}^64:              64\n");
    printf("Число детерминированных ограничений: %d\n",nconstr);
    printf("Число компонент Union-Find:          %d\n",components);
    printf("Ранг (свободных параметров):         %d\n",rank);
    printf("Степени свободы carry при hw59<85:   %d из 64\n\n",rank);

    /* Показать группы */
    printf("=== Группы (связанные carry биты) ===\n");
    for(int root=0;root<64;root++){
        if(group_size[root]<2) continue;
        printf("  Группа (корень=%d):", root);
        for(int i=0;i<64;i++){
            int xi; int r=uf_find(i,&xi);
            if(r==root) printf(" c[%d]%s",i,(xi==0)?"":"`");
        }
        printf("\n");
    }
    
    /* Константные биты */
    printf("\n=== Константные (degenerate) биты при hw59<85 ===\n");
    int nconst=0;
    for(int i=0;i<64;i++){
        double mi=sx[i]/N;
        if(mi<0.05){printf("  carry[%2d] ≈ 0 (E=%.3f)\n",i,mi);nconst++;}
        else if(mi>0.95){printf("  carry[%2d] ≈ 1 (E=%.3f)\n",i,mi);nconst++;}
    }
    printf("Константных битов: %d\n",nconst);
    printf("Полное число ограничений (det+const): %d\n",nconstr+nconst);
    printf("Оценка: carry-вектор hw59<85 лежит на многообразии dim≈%d\n",rank-nconst);

    /* Наиболее биасированные биты */
    printf("\n=== Маргинальные распределения carry[r] при hw59<85 ===\n");
    printf("  r  | E[c[r]] | смещение\n");
    for(int i=0;i<64;i++){
        double mi=sx[i]/N;
        double bias=fabs(mi-0.5);
        if(bias>0.15)
            printf("  %2d  | %.3f   | %.3f%s\n",i,mi,bias,(bias>0.35)?" ★★":"");
    }

    free(args); free(th);
    return 0;
}
