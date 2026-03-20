/*
 * tensor_rank.c — тензорный/якобианный анализ дифференциальной карты SHA-256
 *
 * Q192/Q3: Какова структура Якобиана J = ∂C(W)/∂W[i] ?
 *          (C = carry-вектор, W[0..15] = свободные параметры)
 *
 * Метод:
 *   1. Для каждого образца W (случайного и near-solution):
 *      вычисляем Якобиан J[64×16] методом конечных разностей по каждому биту W[i]
 *   2. Анализируем численный ранг J (сингулярные числа)
 *   3. Сравниваем ранг J для случайного W vs near-solution W
 *
 * Также: анализ ∂hw59/∂W[i] — градиент целевой функции
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

/*
 * Якобиан ∂carry[r]/∂bit_k(W[i]) — булевый (0 или 1)
 * J[r][i*32+bit] = carry_r(W с флипнутым битом) XOR carry_r(W)
 *
 * Размер J: 64 строки (carry биты) × 512 столбцов (биты W[0..15])
 * Вычисляем только для строк i=0..15 (свободные слова)
 */
static void compute_jacobian(const uint32_t W[16], uint8_t J[64][512]){
    uint8_t C0[64]; sha256_carry(W,C0);
    for(int i=0;i<16;i++){
        for(int bit=0;bit<32;bit++){
            uint32_t W2[16]; memcpy(W2,W,64);
            W2[i]^=(1u<<bit);
            uint8_t C1[64]; sha256_carry(W2,C1);
            for(int r=0;r<64;r++)
                J[r][i*32+bit]=C0[r]^C1[r];
        }
    }
}

/* Числовой ранг бинарной матрицы 64×512 методом Гаусса над GF(2) */
static int gf2_rank(uint8_t M[64][512]){
    /* работаем с битовыми строками */
    uint64_t rows[64][8]; /* 64 × 512 bit = 64 × 8 uint64_t */
    memset(rows,0,sizeof(rows));
    for(int r=0;r<64;r++)
        for(int c=0;c<512;c++)
            if(M[r][c]) rows[r][c/64]|=(uint64_t)1<<(c%64);

    int rank=0;
    int pivot_col[64]; memset(pivot_col,-1,sizeof(pivot_col));
    int pivot_row=0;
    for(int col=0;col<512&&pivot_row<64;col++){
        /* найти ненулевую строку в этом столбце */
        int found=-1;
        for(int r=pivot_row;r<64;r++){
            if((rows[r][col/64]>>(col%64))&1){found=r;break;}
        }
        if(found<0) continue;
        /* swap */
        if(found!=pivot_row){
            uint64_t tmp[8]; memcpy(tmp,rows[pivot_row],sizeof(tmp));
            memcpy(rows[pivot_row],rows[found],sizeof(tmp));
            memcpy(rows[found],tmp,sizeof(tmp));
        }
        /* обнулить этот столбец во всех остальных строках */
        for(int r=0;r<64;r++){
            if(r==pivot_row) continue;
            if((rows[r][col/64]>>(col%64))&1){
                for(int w=0;w<8;w++) rows[r][w]^=rows[pivot_row][w];
            }
        }
        rank++;
        pivot_row++;
    }
    return rank;
}

/* hw59-чувствительность: для каждого бита W[i] — как часто flip меняет hw59 */
static void sensitivity_analysis(int nsamples_rand, int nsamples_deep){
    printf("=== Анализ чувствительности ∂hw59/∂bit ===\n");

    /* Инициализируем RNG */
    srand(12345);
    uint64_t rng=0xdeadbeef12345678ULL;
    #define XR() (rng^=rng<<13,rng^=rng>>7,rng^=rng<<17,rng)

    /* Для случайного W */
    double sens_rand[16]={0};
    for(int s=0;s<nsamples_rand;s++){
        uint32_t W[16]; for(int i=0;i<16;i++) W[i]=(uint32_t)XR();
        int hw0=hw59_diff(W);
        for(int i=0;i<16;i++){
            int flips=0;
            for(int bit=0;bit<32;bit++){
                W[i]^=(1u<<bit); int hw1=hw59_diff(W); W[i]^=(1u<<bit);
                if(hw1!=hw0) flips++;
            }
            sens_rand[i]+=flips/32.0;
        }
    }
    printf("  Случайный W (N=%d): среднее число бит, меняющих hw59 при flip W[i]:\n",nsamples_rand);
    for(int i=0;i<16;i++) printf("    W[%2d]: %.2f/32\n",i,sens_rand[i]/nsamples_rand);

    printf("\n=== Якобиан carry-вектора ===\n");
    printf("  Вычисляем J(W) = ∂carry/∂W[0..15] (64×512 над GF(2))\n");
    printf("  Анализируем ранг J для случайных W и near-solution W\n\n");

    /* Случайные W */
    long rank_sum_rand=0; int nrand=0;
    for(int s=0;s<20;s++){
        uint32_t W[16]; for(int i=0;i<16;i++) W[i]=(uint32_t)XR();
        uint8_t J[64][512]; compute_jacobian(W,J);
        int rk=gf2_rank(J);
        rank_sum_rand+=rk; nrand++;
        if(s<5) printf("  Случайный W[%d]: rank(J)=%d\n",s,rk);
    }
    printf("  Среднее rank(J) для случайных W: %.1f\n\n",(double)rank_sum_rand/nrand);

    /* Near-solution W (SA) */
    printf("  Near-solution W (SA, hw59<90):\n");
    uint32_t bestW[16]; int besthw=256;
    uint64_t rs=0x123456789abcdefULL;
    #define XR2() (rs^=rs<<13,rs^=rs>>7,rs^=rs<<17,rs)
    uint32_t W[16]; for(int i=0;i<16;i++) W[i]=(uint32_t)XR2();
    double T=60.0; memcpy(bestW,W,64);
    for(long iter=0;iter<50000000L;iter++){
        int hw=hw59_diff(W);
        if(hw<besthw){besthw=hw;memcpy(bestW,W,64);}
        if(hw<90){
            uint8_t J[64][512]; compute_jacobian(W,J);
            int rk=gf2_rank(J);
            printf("  hw59=%d: rank(J)=%d\n",hw,rk);
            if(nsamples_deep--<=0) break;
        }
        int k=(int)(XR2()%16),bit=(int)(XR2()%32);
        uint32_t sv=W[k]; W[k]^=(1u<<bit);
        int nh=hw59_diff(W);
        double dE=nh-hw;
        if(dE>0&&(double)(XR2()%1000000)/1e6>exp(-dE/T)) W[k]=sv;
        T*=0.9999975; if(T<0.3)T=0.3;
        if(iter%300000==0){memcpy(W,bestW,64);T=40.0;}
        if(iter%5000000==0&&besthw>92){for(int i=0;i<16;i++)W[i]=(uint32_t)XR2();besthw=256;T=60.0;}
    }
    printf("  Лучшее hw59=%d за прогон\n",besthw);
    #undef XR
    #undef XR2
}

int main(void){
    printf("tensor_rank: Якобиан ∂C(W)/∂W и ранговый анализ\n\n");
    printf("Q192: Какова структура Якобиана?\n");
    printf("Q3:   Низок ли ранг дифференциальной карты N(W,ΔW)?\n\n");
    time_t t0=time(NULL);
    sensitivity_analysis(100, 5);
    printf("\nВремя: %ld сек\n",(long)(time(NULL)-t0));
    return 0;
}
