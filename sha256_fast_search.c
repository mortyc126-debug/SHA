/*
 * SHA-256 Дифференциальный каскадный поиск — оптимизированная версия
 * П-15: Быстрый поиск пары с De3..De17=0
 *
 * Оптимизации:
 *   1. Инкрементальное вычисление состояния (1-2 раунда на шаг вместо полных)
 *   2. Пред-вычисление Wn состояний (только 17 раундов однократно)
 *   3. Многопоточность (pthreads)
 *   4. Batch обработка пар
 *
 * Компиляция: gcc -O3 -march=native -pthread -o sha256_fast_search sha256_fast_search.c
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

/* SHA-256 константы */
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

static const u32 H0[8] = {
    0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
    0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19
};

#define ROTR(x,n) (((u32)(x)>>(n))|((u32)(x)<<(32-(n))))
#define sig0(x) (ROTR(x,7)^ROTR(x,18)^((x)>>3))
#define sig1(x) (ROTR(x,17)^ROTR(x,19)^((x)>>10))
#define Sig0(x) (ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define Sig1(x) (ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define Ch(e,f,g) (((e)&(f))^(~(e)&(g)))
#define Maj(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))

/* Один раунд SHA-256, in-place */
#define SHA_ROUND(r, Wa, a,b,c,d,e,f,g,h) do { \
    u32 _T1 = (h)+Sig1(e)+Ch(e,f,g)+K[r]+(Wa); \
    u32 _T2 = Sig0(a)+Maj(a,b,c); \
    (h)=(g);(g)=(f);(f)=(e);(e)=(d)+_T1; \
    (d)=(c);(c)=(b);(b)=(a);(a)=_T1+_T2; \
} while(0)

typedef struct {
    u32 a, b, c, d, e, f, g, h;
} State;

/* Массив состояний (0 = начало, r = после раунда r-1) */
typedef State States18[18];  /* [0..17] */

/* Расширение расписания для первых 18 слов */
static inline void schedule18(u32 *W, const u32 *W16) {
    for (int i = 0; i < 16; i++) W[i] = W16[i];
    W[16] = sig1(W[14]) + W[9]  + sig0(W[1])  + W[0];
    W[17] = sig1(W[15]) + W[10] + sig0(W[2])  + W[1];
}

/* Вычислить все состояния Wn до раунда 17 включительно */
static void precompute_wn_states(const u32 *W, States18 out) {
    u32 a = H0[0], b = H0[1], c = H0[2], d = H0[3];
    u32 e = H0[4], f = H0[5], g = H0[6], h = H0[7];
    out[0] = (State){a,b,c,d,e,f,g,h};
    for (int r = 0; r < 17; r++) {
        SHA_ROUND(r, W[r], a,b,c,d,e,f,g,h);
        out[r+1] = (State){a,b,c,d,e,f,g,h};
    }
}

/*
 * Оптимизированный каскад с инкрементальными состояниями.
 *
 * Алгоритм:
 *   - Предвычисленные состояния Wn_states[0..17] (Wn не меняется)
 *   - Для Wf: кеш состояний Wf_states[0..t+1] (при каскадном шаге t)
 *   - При каждом шаге: добавляем 2 раунда из кешированного состояния
 *
 * Возвращает De17 = e_f[17] - e_n[17]
 */
static inline u32 cascade_fast(u32 W0, u32 W1, const States18 Wn_st, u32 *Wn17,
                                u32 *out_W0, u32 *out_W1,
                                u32 *out_DWs, u32 *out_Wf_ext) {
    u32 DWs[16] = {0};
    DWs[0] = 1;

    /* Расписание Wn (первые 18 слов) */
    u32 Wn16[16] = {0}; Wn16[0] = W0; Wn16[1] = W1;
    u32 Wn_W[18]; schedule18(Wn_W, Wn16);

    /* Кеш Wf состояний */
    State Wf_cache[18];  /* Wf_cache[r] = состояние Wf после r раундов */
    Wf_cache[0] = (State){H0[0],H0[1],H0[2],H0[3],H0[4],H0[5],H0[6],H0[7]};

    /* === Шаг 0: De3_nat, адаптивный ΔW2 ===
     * Нужно sha_e(Wf_tmp, 3) где DWs[0]=1, DWs[1..15]=0
     * Wf_tmp = Wn + [DWs[0],0,0,...] → только W[0] изменён
     */
    u32 Wf16_tmp[16];
    for (int i=0;i<16;i++) Wf16_tmp[i] = Wn16[i] + DWs[i];
    u32 Wf_W[18]; schedule18(Wf_W, Wf16_tmp);

    /* Запустить 3 раунда для Wf_tmp и получить e[3] */
    {
        u32 a=H0[0],b=H0[1],c=H0[2],d=H0[3],e=H0[4],f=H0[5],g=H0[6],h=H0[7];
        SHA_ROUND(0, Wf_W[0], a,b,c,d,e,f,g,h);
        Wf_cache[1] = (State){a,b,c,d,e,f,g,h};
        SHA_ROUND(1, Wf_W[1], a,b,c,d,e,f,g,h);
        Wf_cache[2] = (State){a,b,c,d,e,f,g,h};
        SHA_ROUND(2, Wf_W[2], a,b,c,d,e,f,g,h);
        Wf_cache[3] = (State){a,b,c,d,e,f,g,h};
        /* De3_nat = e_f[3] - e_n[3] */
        u32 de3_nat = e - Wn_st[3].e;
        DWs[2] = (u32)(-(int32_t)de3_nat);
    }
    /* Теперь DWs[2] установлен. Wf_cache[0..1] валиден (rounds 0..1 used W[0],W[1]).
     * DWs[2] изменяет W[2], что инвалидирует Wf_cache[3] и выше.
     * Wf_cache[0..2] остаётся валидным (round 2 использует W[2], но Wf_cache[2] было
     * вычислено с DWs[2]=0. Теперь нужно пересчитать с DWs[2]=val).
     * То есть Wf_cache[0..1] валидны. */

    /* Обновить расписание Wf: W[2] изменился → W[17] изменится */
    for (int i=0;i<16;i++) Wf16_tmp[i] = Wn16[i] + DWs[i];  /* DWs[2] теперь установлен */
    schedule18(Wf_W, Wf16_tmp);

    /* === Шаги 1..13: Каскад De4..De16=0 ===
     * На шаге t (0-indexed, t=1..13):
     *   DWs[t+1] был установлен на предыдущем шаге.
     *   Wf_cache[0..t+1] валидны (rounds 0..t использовали W[0..t-1], W[t] с DWs[t]).
     *   Нам нужно e_f[t+3] = sha_e(Wf, t+3):
     *     1. Из Wf_cache[t+1]: выполнить round t+1 с W[t+1]+DWs[t+1] → Wf_cache[t+2]
     *     2. Из Wf_cache[t+2]: выполнить round t+2 с W[t+2]+0 → e_f[t+3]  (DWs[t+2]=0 пока)
     *   Затем De_{t+3} = e_f[t+3] - e_n[t+3]
     *   DWs[t+2] = -De_{t+3}
     *   Обновить расписание Wf.
     *
     * НО для t=1: DWs[2] изменился (теперь val, раньше 0). Wf_cache[2] невалиден.
     * Нам нужен Wf_cache[2] с новым DWs[2].
     * Wf_cache[1] валидно. Из него: выполнить round 1 с W[1]=Wn[1]+DWs[1]=Wn[1] → Wf_cache[2].
     * Затем round 2 с W[2]+DWs[2] → Wf_cache[3].
     * Затем round 3 с W[3] → e_f[4] = De4.
     *
     * Упрощение: при каждом шаге t, начинаем с Wf_cache[t] (не t+1),
     * добавляем 3 раунда. Но Wf_cache[t] ТОЧНО валиден (round t-1 использует W[t-1]+DWs[t-1],
     * и ни DWs[t] ни DWs[t+1] не изменяли W[0..t-1]).
     *
     * Итого 3 раунда на шаг (не 2) - немного менее оптимально, но проще и корректно.
     */

    /* Перестроить Wf_cache с самого начала - но инкрементально.
     * После установки DWs[2], пересчитываем Wf_cache[2] из Wf_cache[1].
     */
    {
        u32 a=Wf_cache[1].a, b=Wf_cache[1].b, c=Wf_cache[1].c, d=Wf_cache[1].d;
        u32 e=Wf_cache[1].e, f=Wf_cache[1].f, g=Wf_cache[1].g, h=Wf_cache[1].h;
        SHA_ROUND(1, Wf_W[1], a,b,c,d,e,f,g,h);
        Wf_cache[2] = (State){a,b,c,d,e,f,g,h};
    }

    for (int step = 1; step <= 13; step++) {
        int dt = step + 3;  /* нужно e_f[dt] */
        /* DWs[step+1] был установлен на шаге step-1 (изменил W[step+1])
         * Wf_cache[step] валиден (round step-1 использует W[step-1], не изменялся).
         * Нам нужно e_f[dt] = e_f[step+3].
         * Из Wf_cache[step]: 3 раунда → e_f[step+3].
         */
        u32 a=Wf_cache[step].a, b=Wf_cache[step].b, c=Wf_cache[step].c, d=Wf_cache[step].d;
        u32 e=Wf_cache[step].e, f=Wf_cache[step].f, g=Wf_cache[step].g, h=Wf_cache[step].h;

        /* round step (index step): uses W[step] with DWs[step] */
        SHA_ROUND(step, Wf_W[step], a,b,c,d,e,f,g,h);
        Wf_cache[step+1] = (State){a,b,c,d,e,f,g,h};

        /* round step+1 (index step+1): uses W[step+1] with DWs[step+1] (already set) */
        SHA_ROUND(step+1, Wf_W[step+1], a,b,c,d,e,f,g,h);
        Wf_cache[step+2] = (State){a,b,c,d,e,f,g,h};

        /* round step+2 (index step+2): uses W[step+2] with DWs[step+2] = 0 still */
        SHA_ROUND(step+2, Wf_W[step+2], a,b,c,d,e,f,g,h);
        /* e = e_f[step+3] = e_f[dt] */

        u32 de_dt = e - Wn_st[dt].e;
        DWs[step+2] = (u32)(-(int32_t)de_dt);

        /* Обновить расписание Wf */
        for (int i=0;i<16;i++) Wf16_tmp[i] = Wn16[i] + DWs[i];
        schedule18(Wf_W, Wf16_tmp);

        /* Пересчитать Wf_cache[step+2] с новым DWs[step+2] */
        {
            u32 a2=Wf_cache[step+1].a, b2=Wf_cache[step+1].b;
            u32 c2=Wf_cache[step+1].c, d2=Wf_cache[step+1].d;
            u32 e2=Wf_cache[step+1].e, f2=Wf_cache[step+1].f;
            u32 g2=Wf_cache[step+1].g, h2=Wf_cache[step+1].h;
            SHA_ROUND(step+1, Wf_W[step+1], a2,b2,c2,d2,e2,f2,g2,h2);
            Wf_cache[step+2] = (State){a2,b2,c2,d2,e2,f2,g2,h2};
        }
    }
    /* DWs[2..15] теперь все установлены. Вычислим De17. */

    /* Финальное расписание Wf */
    for (int i=0;i<16;i++) Wf16_tmp[i] = Wn16[i] + DWs[i];
    schedule18(Wf_W, Wf16_tmp);

    /* e_f[17] из Wf_cache[15] + 2 раунда */
    {
        u32 a=Wf_cache[15].a, b=Wf_cache[15].b, c=Wf_cache[15].c, d=Wf_cache[15].d;
        u32 e=Wf_cache[15].e, f=Wf_cache[15].f, g=Wf_cache[15].g, h=Wf_cache[15].h;
        SHA_ROUND(15, Wf_W[15], a,b,c,d,e,f,g,h);
        SHA_ROUND(16, Wf_W[16], a,b,c,d,e,f,g,h);
        /* e = e_f[17] */
        u32 de17 = e - Wn_st[17].e;

        if (out_W0) *out_W0 = W0;
        if (out_W1) *out_W1 = W1;
        if (out_DWs) memcpy(out_DWs, DWs, 16*sizeof(u32));
        if (out_Wf_ext) memcpy(out_Wf_ext, Wf_W, 18*sizeof(u32));
        return de17;
    }
}

/* Структура для потока */
typedef struct {
    int thread_id;
    int num_threads;
    u64 iters_per_thread;
    u64 seed;
    atomic_int *found;
    u32 *result_W0, *result_W1, *result_DWs;
    u32 *result_Wf_ext;
    u64 *count_done;
    pthread_mutex_t *mutex;
} ThreadArg;

static void *search_thread(void *arg) {
    ThreadArg *ta = (ThreadArg*)arg;
    u64 rng = ta->seed ^ ((u64)ta->thread_id * 0x9e3779b97f4a7c15ULL);
    u64 local_count = 0;

    while (local_count < ta->iters_per_thread) {
        /* Случайная пара W0, W1 */
        rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
        u32 W0 = (u32)(rng ^ (rng >> 32));
        rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
        u32 W1 = (u32)(rng ^ (rng >> 32));

        /* Предвычислить Wn состояния */
        u32 Wn16[16] = {0}; Wn16[0] = W0; Wn16[1] = W1;
        u32 Wn_W[18]; schedule18(Wn_W, Wn16);
        States18 Wn_st;
        precompute_wn_states(Wn_W, Wn_st);

        u32 de17 = cascade_fast(W0, W1, Wn_st, NULL, NULL, NULL, NULL, NULL);
        local_count++;

        if (de17 == 0) {
            /* Нашли! */
            if (atomic_exchange(ta->found, 1) == 0) {
                pthread_mutex_lock(ta->mutex);
                u32 DWs[16]; u32 Wf_ext[18];
                cascade_fast(W0, W1, Wn_st, NULL,
                             ta->result_W0, ta->result_W1, DWs, Wf_ext);
                if (ta->result_DWs) memcpy(ta->result_DWs, DWs, 16*sizeof(u32));
                if (ta->result_Wf_ext) memcpy(ta->result_Wf_ext, Wf_ext, 18*sizeof(u32));
                pthread_mutex_unlock(ta->mutex);
            }
            break;
        }

        /* Обновить счётчик (каждые 1M) */
        if (local_count % 1000000 == 0) {
            pthread_mutex_lock(ta->mutex);
            *ta->count_done += 1000000;
            pthread_mutex_unlock(ta->mutex);
        }

        if (atomic_load(ta->found)) break;
    }

    pthread_mutex_lock(ta->mutex);
    *ta->count_done += local_count % 1000000;
    pthread_mutex_unlock(ta->mutex);
    return NULL;
}

/* Полная верификация De3..De17 через брутфорс */
static int verify_full(u32 W0, u32 W1, u32 *DWs) {
    u32 Wn16[16] = {0}; Wn16[0] = W0; Wn16[1] = W1;
    u32 Wf16[16];
    for (int i=0;i<16;i++) Wf16[i] = Wn16[i] + DWs[i];
    u32 Wn[64], Wf[64]; schedule18(Wn, Wn16); schedule18(Wf, Wf16);
    /* Extend to 64 words */
    for (int i=18;i<64;i++) {
        Wn[i]=sig1(Wn[i-2])+Wn[i-7]+sig0(Wn[i-15])+Wn[i-16];
        Wf[i]=sig1(Wf[i-2])+Wf[i-7]+sig0(Wf[i-15])+Wf[i-16];
    }

    int ok = 1;
    for (int r = 3; r <= 17; r++) {
        u32 an=H0[0],bn=H0[1],cn=H0[2],dn=H0[3],en=H0[4],fn=H0[5],gn=H0[6],hn=H0[7];
        u32 af=H0[0],bf=H0[1],cf=H0[2],df=H0[3],ef=H0[4],ff=H0[5],gf=H0[6],hf=H0[7];
        for (int i=0;i<r;i++) {
            SHA_ROUND(i, Wn[i], an,bn,cn,dn,en,fn,gn,hn);
            SHA_ROUND(i, Wf[i], af,bf,cf,df,ef,ff,gf,hf);
        }
        u32 de_r = ef - en;
        printf("  De%d = 0x%08x %s\n", r, de_r, de_r==0?"✓":"✗ ERROR!");
        if (de_r != 0) ok = 0;
    }
    return ok;
}

/* Вычислить Da14 и De18 для верификации T_DE18_DECOMPOSITION */
static void compute_de18_analysis(u32 W0, u32 W1, u32 *DWs) {
    u32 Wn16[16] = {0}; Wn16[0] = W0; Wn16[1] = W1;
    u32 Wf16[16];
    for (int i=0;i<16;i++) Wf16[i] = Wn16[i] + DWs[i];
    u32 Wn[64], Wf[64]; schedule18(Wn, Wn16); schedule18(Wf, Wf16);
    for (int i=18;i<64;i++) {
        Wn[i]=sig1(Wn[i-2])+Wn[i-7]+sig0(Wn[i-15])+Wn[i-16];
        Wf[i]=sig1(Wf[i-2])+Wf[i-7]+sig0(Wf[i-15])+Wf[i-16];
    }

    /* De18 */
    u32 an=H0[0],bn=H0[1],cn=H0[2],dn=H0[3],en=H0[4],fn=H0[5],gn=H0[6],hn=H0[7];
    u32 af=H0[0],bf=H0[1],cf=H0[2],df=H0[3],ef=H0[4],ff=H0[5],gf=H0[6],hf=H0[7];
    for (int i=0;i<18;i++) {
        SHA_ROUND(i, Wn[i], an,bn,cn,dn,en,fn,gn,hn);
        SHA_ROUND(i, Wf[i], af,bf,cf,df,ef,ff,gf,hf);
    }
    u32 de18 = ef - en;

    /* Da14: a после 14 раундов */
    u32 an2=H0[0],bn2=H0[1],cn2=H0[2],dn2=H0[3],en2=H0[4],fn2=H0[5],gn2=H0[6],hn2=H0[7];
    u32 af2=H0[0],bf2=H0[1],cf2=H0[2],df2=H0[3],ef2=H0[4],ff2=H0[5],gf2=H0[6],hf2=H0[7];
    for (int i=0;i<14;i++) {
        SHA_ROUND(i, Wn[i], an2,bn2,cn2,dn2,en2,fn2,gn2,hn2);
        SHA_ROUND(i, Wf[i], af2,bf2,cf2,df2,ef2,ff2,gf2,hf2);
    }
    u32 da14 = af2 - an2;
    u32 dw17 = Wf[17] - Wn[17];

    printf("\nT_DE18_DECOMPOSITION:\n");
    printf("  De18     = 0x%08x\n", de18);
    printf("  Da14     = 0x%08x\n", da14);
    printf("  ΔW17     = 0x%08x\n", dw17);
    printf("  Da14+ΔW17= 0x%08x\n", da14+dw17);
    printf("  Theorem: %s\n",
           (de18 == da14+dw17) ? "De18 = Da14 + ΔW17  ✓ T_DE18_DECOMPOSITION VERIFIED!" :
                                  "De18 ≠ Da14 + ΔW17  ✗ NOT VERIFIED!");

    /* De19..De21 декомпозиция */
    printf("\nDe19..De21 decomposition:\n");
    for (int k = 19; k <= 21; k++) {
        u32 a_n=H0[0],b_n=H0[1],c_n=H0[2],d_n=H0[3],e_n=H0[4],f_n=H0[5],g_n=H0[6],h_n=H0[7];
        u32 a_f=H0[0],b_f=H0[1],c_f=H0[2],d_f=H0[3],e_f=H0[4],f_f=H0[5],g_f=H0[6],h_f=H0[7];
        for (int i=0;i<k;i++) {
            SHA_ROUND(i, Wn[i], a_n,b_n,c_n,d_n,e_n,f_n,g_n,h_n);
            SHA_ROUND(i, Wf[i], a_f,b_f,c_f,d_f,e_f,f_f,g_f,h_f);
        }
        u32 de_k = e_f - e_n;

        u32 a_nk=H0[0],b_nk=H0[1],c_nk=H0[2],d_nk=H0[3],e_nk=H0[4],f_nk=H0[5],g_nk=H0[6],h_nk=H0[7];
        u32 a_fk=H0[0],b_fk=H0[1],c_fk=H0[2],d_fk=H0[3],e_fk=H0[4],f_fk=H0[5],g_fk=H0[6],h_fk=H0[7];
        for (int i=0;i<k-4;i++) {
            SHA_ROUND(i, Wn[i], a_nk,b_nk,c_nk,d_nk,e_nk,f_nk,g_nk,h_nk);
            SHA_ROUND(i, Wf[i], a_fk,b_fk,c_fk,d_fk,e_fk,f_fk,g_fk,h_fk);
        }
        u32 da_km4 = a_fk - a_nk;
        u32 dw_km1 = Wf[k-1] - Wn[k-1];
        printf("  De%d = 0x%08x, Da%d+ΔW%d = 0x%08x  %s\n",
               k, de_k, k-4, k-1, da_km4+dw_km1,
               (de_k == da_km4+dw_km1) ? "✓" : "✗");
    }

    /* Python output */
    printf("\nPython constants:\n");
    printf("W0  = 0x%08x\n", W0);
    printf("W1  = 0x%08x\n", W1);
    printf("DWs = [");
    for (int i=0;i<16;i++) printf("0x%x%s", DWs[i], i<15?", ":"");
    printf("]\n");
    printf("Da14 = 0x%08x\n", da14);
    printf("DW17 = 0x%08x\n", dw17);
    printf("De18 = 0x%08x\n", de18);
}

int main(int argc, char *argv[]) {
    int mode = 0;  /* 0=search, 1=verify */
    u64 seed = (u64)time(NULL);
    int num_threads = 4;
    u64 max_iters = (u64)8000000000ULL;  /* 2^33 iterations */

    /* Разобрать аргументы */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--verify") == 0 && i+2 < argc) {
            mode = 1;
            /* verify W0 W1 - тестируем корректность */
            u32 W0 = (u32)strtoul(argv[i+1], NULL, 16);
            u32 W1 = (u32)strtoul(argv[i+2], NULL, 16);
            u32 Wn16[16] = {0}; Wn16[0] = W0; Wn16[1] = W1;
            u32 Wn_W[18]; schedule18(Wn_W, Wn16);
            States18 Wn_st;
            precompute_wn_states(Wn_W, Wn_st);
            u32 DWs[16]; u32 Wf_ext[18];
            u32 de17 = cascade_fast(W0, W1, Wn_st, NULL, NULL, NULL, DWs, Wf_ext);
            printf("W0=%08x W1=%08x De17=%08x\n", W0, W1, de17);
            return 0;
        } else if (strcmp(argv[i], "--seed") == 0 && i+1 < argc) {
            seed = (u64)atoll(argv[i+1]); i++;
        } else if (strcmp(argv[i], "--threads") == 0 && i+1 < argc) {
            num_threads = atoi(argv[i+1]); i++;
        } else if (strcmp(argv[i], "--iters") == 0 && i+1 < argc) {
            max_iters = (u64)atoll(argv[i+1]); i++;
        }
    }

    printf("=================================================================\n");
    printf("П-15: Fast Search De3..De17=0 (optimized C, %d threads)\n", num_threads);
    printf("=================================================================\n");
    printf("Seed: %llu, Threads: %d, Max: %llu M iters\n",
           (unsigned long long)seed, num_threads,
           (unsigned long long)(max_iters / 1000000));
    printf("Expected: ~%.1f B iterations (P=2^-32)\n\n", (double)(1ULL<<32)/1e9);

    /* Speed test (1M iterations on thread 0) */
    {
        u64 rng = seed;
        clock_t t0 = clock();
        for (int i = 0; i < 1000000; i++) {
            rng ^= rng<<13; rng ^= rng>>7; rng ^= rng<<17;
            u32 W0 = (u32)(rng^(rng>>32));
            rng ^= rng<<13; rng ^= rng>>7; rng ^= rng<<17;
            u32 W1 = (u32)(rng^(rng>>32));
            u32 Wn16[16]={0}; Wn16[0]=W0; Wn16[1]=W1;
            u32 Wn_W[18]; schedule18(Wn_W, Wn16);
            States18 Wn_st; precompute_wn_states(Wn_W, Wn_st);
            cascade_fast(W0, W1, Wn_st, NULL, NULL, NULL, NULL, NULL);
        }
        double t = (double)(clock()-t0)/CLOCKS_PER_SEC;
        printf("Speed (1 thread): %.2f M/s → %d threads: %.0f M/s\n",
               1.0/t, num_threads, num_threads/t);
        printf("ETA: ~%.0f seconds (~%.1f hours)\n\n",
               (double)(1ULL<<32) / (num_threads/t * 1e6),
               (double)(1ULL<<32) / (num_threads/t * 1e6) / 3600.0);
    }

    /* Запустить потоки поиска */
    atomic_int found = 0;
    u32 result_W0 = 0, result_W1 = 0, result_DWs[16] = {0}, result_Wf[18] = {0};
    u64 count_done = 0;
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

    pthread_t threads[16];
    ThreadArg args[16];
    if (num_threads > 16) num_threads = 16;

    for (int t = 0; t < num_threads; t++) {
        args[t] = (ThreadArg){
            .thread_id = t, .num_threads = num_threads,
            .iters_per_thread = max_iters / num_threads,
            .seed = seed ^ ((u64)t * 0x123456789ABCDEFULL),
            .found = &found,
            .result_W0 = &result_W0, .result_W1 = &result_W1,
            .result_DWs = result_DWs, .result_Wf_ext = result_Wf,
            .count_done = &count_done, .mutex = &mutex
        };
        pthread_create(&threads[t], NULL, search_thread, &args[t]);
    }

    /* Прогресс */
    clock_t t_start = clock();
    while (!atomic_load(&found)) {
        struct timespec ts = {5, 0};
        nanosleep(&ts, NULL);
        pthread_mutex_lock(&mutex);
        u64 done = count_done;
        pthread_mutex_unlock(&mutex);
        double elapsed = (double)(clock()-t_start)/CLOCKS_PER_SEC;
        if (elapsed > 0) {
            double rate = done / elapsed / 1e6;
            double eta = ((double)(1ULL<<32) - done) / (rate * 1e6);
            printf("  [%llu M / %.1f M/s] elapsed=%.0fs, ETA=%.0fs\n",
                   (unsigned long long)(done/1000000ULL), rate, elapsed, eta);
            fflush(stdout);
        }
        if (done >= max_iters) break;
    }

    for (int t = 0; t < num_threads; t++)
        pthread_join(threads[t], NULL);

    double elapsed = (double)(clock()-t_start)/CLOCKS_PER_SEC;
    printf("\n");
    printf("=================================================================\n");

    if (atomic_load(&found)) {
        printf("SUCCESS! De17=0 found!\n");
        printf("=================================================================\n");
        printf("W0 = 0x%08x\n", result_W0);
        printf("W1 = 0x%08x\n", result_W1);
        printf("DWs = [");
        for (int i=0;i<16;i++) printf("0x%x%s", result_DWs[i], i<15?", ":"");
        printf("]\n\n");

        printf("Verification De3..De17:\n");
        int ok = verify_full(result_W0, result_W1, result_DWs);
        if (ok) printf("\nAll De3..De17=0 ✓  T_CASCADE_17 VERIFIED!\n");

        compute_de18_analysis(result_W0, result_W1, result_DWs);
    } else {
        printf("Search completed without finding De17=0 (%llu M iterations)\n",
               (unsigned long long)(count_done/1000000ULL));
    }

    printf("\nSearch stats:\n");
    printf("  Total: %llu iterations\n", (unsigned long long)count_done);
    printf("  Time: %.1f sec\n", elapsed);
    printf("  Rate: %.2f M/s\n", count_done / elapsed / 1e6);
    printf("=================================================================\n");
    return atomic_load(&found) ? 0 : 1;
}
