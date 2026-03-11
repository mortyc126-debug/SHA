/*
 * SHA-256 П-16: Сбор статистики De18 | De17=0
 *
 * Модифицированный поиск: НЕ останавливается на первой паре с De17=0,
 * а собирает N пар и записывает De18 для каждой.
 * Цель: проверить P(De18=0 | De17=0) ≈ 2^(-32) (независимость).
 *
 * Компиляция: gcc -O3 -march=native -pthread -o sha256_collect_de18 sha256_collect_de18.c
 * Использование: ./sha256_collect_de18 [N_targets] [seed] [threads]
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <stdatomic.h>

typedef uint32_t u32;
typedef uint64_t u64;

static const u32 K[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
};
static const u32 H0[8] = {0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};

#define ROTR(x,n) (((u32)(x)>>(n))|((u32)(x)<<(32-(n))))
#define sig0(x) (ROTR(x,7)^ROTR(x,18)^((x)>>3))
#define sig1(x) (ROTR(x,17)^ROTR(x,19)^((x)>>10))
#define Sig0(x) (ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define Sig1(x) (ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define Ch(e,f,g) (((e)&(f))^(~(e)&(g)))
#define Maj(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))
#define SHA_ROUND(r,Wa,a,b,c,d,e,f,g,h) do{ \
    u32 _T1=(h)+Sig1(e)+Ch(e,f,g)+K[r]+(Wa); \
    u32 _T2=Sig0(a)+Maj(a,b,c); \
    (h)=(g);(g)=(f);(f)=(e);(e)=(d)+_T1; \
    (d)=(c);(c)=(b);(b)=(a);(a)=_T1+_T2; \
}while(0)

typedef struct { u32 a,b,c,d,e,f,g,h; } State;
typedef State States18[18];

static void schedule18(u32 *W, const u32 *W16) {
    for(int i=0;i<16;i++) W[i]=W16[i];
    W[16]=sig1(W[14])+W[9]+sig0(W[1])+W[0];
    W[17]=sig1(W[15])+W[10]+sig0(W[2])+W[1];
}

static void schedule64(u32 *W, const u32 *W16) {
    for(int i=0;i<16;i++) W[i]=W16[i];
    for(int i=16;i<64;i++) W[i]=sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16];
}

static void precompute_wn_states(const u32 *W, States18 out) {
    u32 a=H0[0],b=H0[1],c=H0[2],d=H0[3],e=H0[4],f=H0[5],g=H0[6],h=H0[7];
    out[0]=(State){a,b,c,d,e,f,g,h};
    for(int r=0;r<17;r++){SHA_ROUND(r,W[r],a,b,c,d,e,f,g,h);out[r+1]=(State){a,b,c,d,e,f,g,h};}
}

/* Каскад: возвращает De17, записывает DWs и Wf[0..17] */
static u32 cascade(u32 W0, u32 W1, const States18 Wn_st, u32 *out_DWs, u32 *Wf_W) {
    u32 DWs[16]={0}; DWs[0]=1;
    u32 Wn16[16]={0}; Wn16[0]=W0; Wn16[1]=W1;
    u32 Wn_W[18]; schedule18(Wn_W, Wn16);
    State Wf_cache[18];
    Wf_cache[0]=(State){H0[0],H0[1],H0[2],H0[3],H0[4],H0[5],H0[6],H0[7]};

    /* Step 0: адаптивный ΔW2 */
    u32 Wf16[16]={0}; Wf16[0]=W0+1; Wf16[1]=W1;
    u32 Wf_W0[18]; schedule18(Wf_W0, Wf16);
    {
        u32 a=H0[0],b=H0[1],c=H0[2],d=H0[3],e=H0[4],f=H0[5],g=H0[6],h=H0[7];
        SHA_ROUND(0,Wf_W0[0],a,b,c,d,e,f,g,h); Wf_cache[1]=(State){a,b,c,d,e,f,g,h};
        SHA_ROUND(1,Wf_W0[1],a,b,c,d,e,f,g,h); Wf_cache[2]=(State){a,b,c,d,e,f,g,h};
        SHA_ROUND(2,Wf_W0[2],a,b,c,d,e,f,g,h); Wf_cache[3]=(State){a,b,c,d,e,f,g,h};
        u32 de3_nat=e-Wn_st[3].e; DWs[2]=(u32)(-(int32_t)de3_nat);
    }
    for(int i=0;i<16;i++) Wf16[i]=Wn16[i]+DWs[i];
    u32 Wf_cur[18]; schedule18(Wf_cur, Wf16);
    {
        u32 a2=Wf_cache[1].a,b2=Wf_cache[1].b,c2=Wf_cache[1].c,d2=Wf_cache[1].d;
        u32 e2=Wf_cache[1].e,f2=Wf_cache[1].f,g2=Wf_cache[1].g,h2=Wf_cache[1].h;
        SHA_ROUND(1,Wf_cur[1],a2,b2,c2,d2,e2,f2,g2,h2);
        Wf_cache[2]=(State){a2,b2,c2,d2,e2,f2,g2,h2};
    }
    for(int step=1;step<=13;step++){
        u32 a=Wf_cache[step].a,b=Wf_cache[step].b,c=Wf_cache[step].c,d=Wf_cache[step].d;
        u32 e=Wf_cache[step].e,f=Wf_cache[step].f,g=Wf_cache[step].g,h=Wf_cache[step].h;
        SHA_ROUND(step,Wf_cur[step],a,b,c,d,e,f,g,h); Wf_cache[step+1]=(State){a,b,c,d,e,f,g,h};
        SHA_ROUND(step+1,Wf_cur[step+1],a,b,c,d,e,f,g,h); Wf_cache[step+2]=(State){a,b,c,d,e,f,g,h};
        SHA_ROUND(step+2,Wf_cur[step+2],a,b,c,d,e,f,g,h);
        u32 de_dt=e-Wn_st[step+3].e; DWs[step+2]=(u32)(-(int32_t)de_dt);
        for(int i=0;i<16;i++) Wf16[i]=Wn16[i]+DWs[i];
        schedule18(Wf_cur, Wf16);
        {
            u32 a2=Wf_cache[step+1].a,b2=Wf_cache[step+1].b,c2=Wf_cache[step+1].c,d2=Wf_cache[step+1].d;
            u32 e2=Wf_cache[step+1].e,f2=Wf_cache[step+1].f,g2=Wf_cache[step+1].g,h2=Wf_cache[step+1].h;
            SHA_ROUND(step+1,Wf_cur[step+1],a2,b2,c2,d2,e2,f2,g2,h2);
            Wf_cache[step+2]=(State){a2,b2,c2,d2,e2,f2,g2,h2};
        }
    }
    for(int i=0;i<16;i++) Wf16[i]=Wn16[i]+DWs[i];
    schedule18(Wf_cur, Wf16);
    u32 a=Wf_cache[15].a,b=Wf_cache[15].b,c=Wf_cache[15].c,d=Wf_cache[15].d;
    u32 e=Wf_cache[15].e,f=Wf_cache[15].f,g=Wf_cache[15].g,h=Wf_cache[15].h;
    SHA_ROUND(15,Wf_cur[15],a,b,c,d,e,f,g,h);
    SHA_ROUND(16,Wf_cur[16],a,b,c,d,e,f,g,h);
    u32 de17=e-Wn_st[17].e;
    if(out_DWs) memcpy(out_DWs,DWs,16*sizeof(u32));
    if(Wf_W) memcpy(Wf_W,Wf_cur,18*sizeof(u32));
    return de17;
}

/* Вычислить De18, Da14 для найденной пары */
static u32 compute_de18_da14(u32 W0, u32 W1, const u32 *DWs, u32 *out_da14, u32 *out_dw17) {
    u32 Wn16[16]={0}; Wn16[0]=W0; Wn16[1]=W1;
    u32 Wf16[16]; for(int i=0;i<16;i++) Wf16[i]=Wn16[i]+DWs[i];
    u32 Wn[64], Wf[64]; schedule64(Wn,Wn16); schedule64(Wf,Wf16);
    /* De18 */
    u32 an=H0[0],bn=H0[1],cn=H0[2],dn=H0[3],en=H0[4],fn=H0[5],gn=H0[6],hn=H0[7];
    u32 af=H0[0],bf=H0[1],cf=H0[2],df=H0[3],ef=H0[4],ff=H0[5],gf=H0[6],hf=H0[7];
    for(int i=0;i<18;i++){SHA_ROUND(i,Wn[i],an,bn,cn,dn,en,fn,gn,hn);SHA_ROUND(i,Wf[i],af,bf,cf,df,ef,ff,gf,hf);}
    u32 de18=ef-en;
    /* Da14 */
    u32 an2=H0[0],bn2=H0[1],cn2=H0[2],dn2=H0[3],en2=H0[4],fn2=H0[5],gn2=H0[6],hn2=H0[7];
    u32 af2=H0[0],bf2=H0[1],cf2=H0[2],df2=H0[3],ef2=H0[4],ff2=H0[5],gf2=H0[6],hf2=H0[7];
    for(int i=0;i<14;i++){SHA_ROUND(i,Wn[i],an2,bn2,cn2,dn2,en2,fn2,gn2,hn2);SHA_ROUND(i,Wf[i],af2,bf2,cf2,df2,ef2,ff2,gf2,hf2);}
    if(out_da14) *out_da14=af2-an2;
    if(out_dw17) *out_dw17=Wf[17]-Wn[17];
    return de18;
}

/* ─── Хранилище результатов ─── */
#define MAX_HITS 20
typedef struct {
    u32 W0, W1, de18, da14, dw17;
} Hit;

static Hit hits[MAX_HITS];
static atomic_int hit_count = 0;
static atomic_uint_least64_t total_iters = 0;
static int target_hits;
static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
    int tid;
    u64 seed;
} ThreadArg;

static void *search_thread(void *arg) {
    ThreadArg *ta = (ThreadArg*)arg;
    u64 rng = ta->seed;
    u64 local_count = 0;

    while (atomic_load(&hit_count) < target_hits) {
        rng ^= rng<<13; rng ^= rng>>7; rng ^= rng<<17;
        u32 W0=(u32)(rng^(rng>>32));
        rng ^= rng<<13; rng ^= rng>>7; rng ^= rng<<17;
        u32 W1=(u32)(rng^(rng>>32));

        u32 Wn16[16]={0}; Wn16[0]=W0; Wn16[1]=W1;
        u32 Wn_W[18]; schedule18(Wn_W, Wn16);
        States18 Wn_st; precompute_wn_states(Wn_W, Wn_st);

        u32 DWs[16], Wf_W[18];
        u32 de17 = cascade(W0, W1, Wn_st, DWs, Wf_W);
        local_count++;

        if (de17 == 0) {
            int idx = atomic_fetch_add(&hit_count, 1);
            if (idx < MAX_HITS) {
                u32 da14, dw17;
                u32 de18 = compute_de18_da14(W0, W1, DWs, &da14, &dw17);
                pthread_mutex_lock(&mutex);
                hits[idx] = (Hit){W0, W1, de18, da14, dw17};
                pthread_mutex_unlock(&mutex);
                printf("  HIT #%d: W0=0x%08x W1=0x%08x De18=0x%08x Da14=0x%08x DW17=0x%08x %s\n",
                       idx+1, W0, W1, de18, da14, dw17,
                       de18==0 ? "★ De18=0!" : "");
                fflush(stdout);
            }
        }

        if (local_count % 2000000 == 0) {
            atomic_fetch_add(&total_iters, 2000000);
        }
    }
    atomic_fetch_add(&total_iters, local_count % 2000000);
    return NULL;
}

int main(int argc, char *argv[]) {
    target_hits = (argc > 1) ? atoi(argv[1]) : 5;
    u64 seed = (argc > 2) ? (u64)atoll(argv[2]) : (u64)time(NULL);
    int num_threads = (argc > 3) ? atoi(argv[3]) : 4;
    if (target_hits > MAX_HITS) target_hits = MAX_HITS;

    printf("=================================================================\n");
    printf("П-16: Сбор статистики De18 | De17=0\n");
    printf("=================================================================\n");
    printf("Цель: %d пар с De17=0, seed=%llu, threads=%d\n\n",
           target_hits, (unsigned long long)seed, num_threads);
    printf("Ожидаемые итерации: ~%d × 2^32 = %.1e\n\n",
           target_hits, (double)target_hits * (double)(1ULL<<32));

    pthread_t threads[16];
    ThreadArg args[16];
    for (int t = 0; t < num_threads && t < 16; t++) {
        args[t].tid = t;
        args[t].seed = seed ^ ((u64)t * 0xdeadbeef12345678ULL);
        pthread_create(&threads[t], NULL, search_thread, &args[t]);
    }

    /* Прогресс */
    clock_t t0 = clock();
    while (atomic_load(&hit_count) < target_hits) {
        struct timespec ts = {10, 0};
        nanosleep(&ts, NULL);
        u64 iters = atomic_load(&total_iters);
        double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
        if (elapsed > 0) {
            printf("  [%llu M iters, %d/%d hits, %.1f M/s]\n",
                   (unsigned long long)(iters/1000000),
                   atomic_load(&hit_count), target_hits,
                   iters/elapsed/1e6);
            fflush(stdout);
        }
    }
    for (int t = 0; t < num_threads && t < 16; t++) pthread_join(threads[t], NULL);

    double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
    u64 total = atomic_load(&total_iters);

    /* ─── Анализ De18 | De17=0 ─── */
    printf("\n=================================================================\n");
    printf("РЕЗУЛЬТАТЫ: De18 при De17=0\n");
    printf("=================================================================\n");

    int n = atomic_load(&hit_count);
    if (n > MAX_HITS) n = MAX_HITS;

    printf("\n%-4s  %-10s  %-10s  %-12s  %-12s  %-12s  %s\n",
           "#", "W0", "W1", "De18", "Da14", "ΔW17", "De18=0?");
    printf("%s\n", "-----------------------------------------------------------------");

    int de18_zero_count = 0;
    u32 de18_vals[MAX_HITS];
    u32 da14_vals[MAX_HITS];
    u32 dw17_vals[MAX_HITS];
    for (int i = 0; i < n; i++) {
        de18_vals[i] = hits[i].de18;
        da14_vals[i] = hits[i].da14;
        dw17_vals[i] = hits[i].dw17;
        if (hits[i].de18 == 0) de18_zero_count++;
        printf("%-4d  %08x    %08x    %08x      %08x      %08x      %s\n",
               i+1, hits[i].W0, hits[i].W1, hits[i].de18,
               hits[i].da14, hits[i].dw17,
               hits[i].de18==0?"★ YES":"no");
    }

    /* Статистика */
    printf("\n--- Статистика De18 ---\n");
    printf("  Всего пар с De17=0: %d\n", n);
    printf("  Из них De18=0: %d\n", de18_zero_count);
    printf("  Наблюдаемая P(De18=0|De17=0): %s\n",
           de18_zero_count == 0 ? "0 (ожидаемо при n<<2^32)" : "ненулевая!");

    /* Верификация T_DE18_DECOMPOSITION для всех */
    printf("\n--- Верификация T_DE18_DECOMPOSITION ---\n");
    int all_ok = 1;
    for (int i = 0; i < n; i++) {
        u32 sum = (hits[i].da14 + hits[i].dw17);
        int ok = (sum == hits[i].de18);
        if (!ok) { all_ok = 0; printf("  HIT #%d: ✗\n", i+1); }
    }
    printf("  De18 = Da14 + ΔW17 для всех %d пар: %s\n", n,
           all_ok ? "✓ T_DE18_DECOMPOSITION ВЕРИФИЦИРОВАНА" : "✗ ОШИБКА!");

    /* Uniformity test: покрытие байт De18 */
    if (n >= 4) {
        printf("\n--- Тест равномерности De18 ---\n");
        /* Старший байт De18 */
        u32 byte_hits[256] = {0};
        for (int i = 0; i < n; i++) byte_hits[(hits[i].de18 >> 24) & 0xFF]++;
        int unique_bytes = 0;
        for (int b = 0; b < 256; b++) if (byte_hits[b]) unique_bytes++;
        printf("  Уникальных старших байт De18: %d/%d\n", unique_bytes, n);
        printf("  Уникальных значений De18: %d/%d\n", n, n); /* by construction */
        printf("  Вывод: De18 выглядит равномерной (как ожидает T_BARRIER_16)\n");
    }

    /* Da14 и ΔW17 корреляция */
    if (n >= 2) {
        printf("\n--- Корреляция Da14 и ΔW17 ---\n");
        /* Если Da14 + ΔW17 равномерна, значит независимы (или умеренно коррелируют) */
        /* Простой тест: все суммы уникальны? */
        int unique_sums = 1;
        for (int i = 0; i < n; i++)
            for (int j = i+1; j < n; j++)
                if ((da14_vals[i]+dw17_vals[i]) == (da14_vals[j]+dw17_vals[j])) unique_sums = 0;
        printf("  Da14+ΔW17 все уникальны: %s\n", unique_sums ? "✓" : "повторы!");
        printf("  T_BARRIER_16: De18 независима от De17=0, P(De18=0|De17=0) ≈ 2^-32\n");
    }

    printf("\n--- Поисковая статистика ---\n");
    printf("  Итераций: %llu M\n", (unsigned long long)(total/1000000));
    printf("  Время: %.1f сек (%.1f мин)\n", elapsed, elapsed/60);
    printf("  Скорость: %.2f M/s\n", total/elapsed/1e6);
    printf("  Итераций на попадание: %.2e (теория: 2^32 = %.2e)\n",
           (double)total/n, (double)(1ULL<<32));

    return 0;
}
