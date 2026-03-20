/*
 * carry_matrix.c — Полная матрица корреляций carry[64×64] по уровням hw59
 *
 * Цель: измерить phi(carry[i], carry[j]) для ВСЕХ 64×64 пар
 * по образцам hw59 от random до <85.
 *
 * Новые вопросы:
 *   Q-E: Повторяет ли лаг-структура корреляций расписание SHA (шаг 7)?
 *   Q-F: Есть ли "собственные векторы" carry-матрицы? (PCA)
 *   Q-G: Как меняется СУММА матрицы с hw59? (мера "корреляционной энергии")
 *   Q-H: Есть ли корреляции ВНУТРИ зоны r<17 при глубоком hw59?
 *
 * Компиляция: gcc -O3 -march=native -pthread -o carry_matrix carry_matrix.c -lm
 * Запуск: ./carry_matrix [target_deep] [threads]
 *         ./carry_matrix 500 6
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

/* ===== SHA-256 ===== */
static const uint32_t K[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,
    0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,
    0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,
    0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,
    0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,
    0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,
    0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,
    0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,
    0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};
static const uint32_t IV[8] = {
    0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
    0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19
};
#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))
#define S0(x) (ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define S1(x) (ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define s0(x) (ROTR(x,7)^ROTR(x,18)^((x)>>3))
#define s1(x) (ROTR(x,17)^ROTR(x,19)^((x)>>10))
#define Ch(x,y,z) (((x)&(y))^((~(x))&(z)))
#define Maj(x,y,z) (((x)&(y))^((x)&(z))^((y)&(z)))

static void sha256_carry(const uint32_t W_in[16], uint8_t carry[64],
                          uint32_t state_out59[8]) {
    uint32_t W[64];
    for (int i = 0; i < 16; i++) W[i] = W_in[i];
    for (int i = 16; i < 64; i++)
        W[i] = s1(W[i-2]) + W[i-7] + s0(W[i-15]) + W[i-16];

    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3];
    uint32_t e=IV[4],f=IV[5],g=IV[6],h=IV[7];

    for (int r = 0; r < 64; r++) {
        uint32_t T1 = h + S1(e) + Ch(e,f,g) + K[r] + W[r];
        uint32_t T2 = S0(a) + Maj(a,b,c);
        uint64_t sum = (uint64_t)T1 + T2;
        carry[r] = (sum >> 32) & 1;
        h=g; g=f; f=e; e=d+T1; d=c; c=b; b=a;
        a = (uint32_t)(sum & 0xFFFFFFFF);
        if (r == 58) {
            if (state_out59) {
                /* состояние ПОСЛЕ раунда 58 = состояние ВХОДА раунда 59 */
                state_out59[0]=a; state_out59[1]=b; state_out59[2]=c;
                state_out59[3]=d; state_out59[4]=e; state_out59[5]=f;
                state_out59[6]=g; state_out59[7]=h;
            }
        }
    }
}

static int hw59_diff(const uint32_t W[16], uint32_t dw0) {
    uint32_t W1[16];
    for (int i = 0; i < 16; i++) W1[i] = W[i];
    W1[0] ^= dw0;

    uint8_t c0[64], c1[64];
    uint32_t s59_0[8], s59_1[8];
    sha256_carry(W, c0, s59_0);
    sha256_carry(W1, c1, s59_1);

    int hw = 0;
    for (int i = 0; i < 8; i++)
        hw += __builtin_popcount(s59_0[i] ^ s59_1[i]);
    return hw;
}

static uint64_t rng_state;
static uint64_t xrand() {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

/* ===== Аккумулятор для матрицы 64×64 ===== */
typedef struct {
    double si[64], sij[64][64];
    long   n;
} CMatrix;

static void cmat_add(CMatrix *m, const uint8_t carry[64]) {
    for (int i = 0; i < 64; i++) {
        m->si[i] += carry[i];
        for (int j = i; j < 64; j++)
            m->sij[i][j] += carry[i] * carry[j];
    }
    m->n++;
}

static double cmat_phi(const CMatrix *m, int i, int j) {
    if (i > j) { int t = i; i = j; j = t; }
    if (m->n < 2) return 0;
    double n = m->n;
    double mi = m->si[i]/n, mj = m->si[j]/n;
    double mij = m->sij[i][j]/n;
    double vi = mi*(1-mi), vj = mj*(1-mj);
    if (vi<1e-10 || vj<1e-10) return 0;
    return (mij - mi*mj) / sqrt(vi*vj);
}

/* Вычислить "энергию" матрицы — сумму |phi| всех пар */
static double cmat_energy(const CMatrix *m) {
    double E = 0;
    for (int i = 0; i < 64; i++)
        for (int j = i+1; j < 64; j++)
            E += fabs(cmat_phi(m, i, j));
    return E;
}

/* Анализ лаговой структуры */
static void lag_analysis(const CMatrix *m, const char *label) {
    /* Среднее |phi(carry[r], carry[r+lag])| по всем r */
    printf("\n  [%s] Лаговая структура (avg |phi| по сдвигу):\n", label);
    printf("  lag | avg|phi| | max|phi| | топ пара\n");
    for (int lag = 1; lag <= 20; lag++) {
        double sum = 0, mx = 0;
        int mi_r = 0, mj_r = 0;
        int cnt = 0;
        for (int r = 0; r + lag < 64; r++) {
            double p = fabs(cmat_phi(m, r, r+lag));
            sum += p;
            cnt++;
            if (p > mx) { mx = p; mi_r = r; mj_r = r+lag; }
        }
        double avg = cnt > 0 ? sum/cnt : 0;
        if (lag <= 10 || mx > 0.05)
            printf("  %3d | %8.4f | %8.4f | c[%2d]↔c[%2d]\n",
                lag, avg, mx, mi_r, mj_r);
    }
}

/* ===== Поток ===== */
#define N_BUCKETS 8  /* hw59 buckets: <85,<90,<95,<100,<105,<110,<120,all */
static const int BUCKET_THR[N_BUCKETS] = {85,90,95,100,105,110,120,9999};
static const char *BUCKET_NAME[N_BUCKETS] = {
    "hw59<85","hw59<90","hw59<95","hw59<100",
    "hw59<105","hw59<110","hw59<120","ALL"
};

typedef struct {
    int       tid;
    long      target_deep;  /* цель: образцов hw59<90 */
    uint64_t  seed;
    uint32_t  dw0;
    CMatrix   bucket[N_BUCKETS];
    long      found_deep;
    pthread_mutex_t *mutex;
    CMatrix   *shared[N_BUCKETS];
} TArg;

static void *tworker(void *varg) {
    TArg *ta = (TArg*)varg;
    rng_state = ta->seed;

    CMatrix local[N_BUCKETS];
    memset(local, 0, sizeof(local));

    uint32_t W[16];
    for (int i = 0; i < 16; i++) W[i] = (uint32_t)xrand();

    double T = 40.0;
    uint32_t best_W[16];
    int best_hw = 256;
    memcpy(best_W, W, sizeof(W));

    long iters = 0;
    long max_iters = ta->target_deep * 5000;
    long deep_found = 0;

    while (iters < max_iters && deep_found < ta->target_deep) {
        iters++;

        int hw = hw59_diff(W, ta->dw0);

        /* Накапливаем для всех подходящих buckets */
        uint8_t carry[64];
        sha256_carry(W, carry, NULL);

        for (int b = 0; b < N_BUCKETS; b++)
            if (hw < BUCKET_THR[b])
                cmat_add(&local[b], carry);

        if (hw < 90) deep_found++;
        if (hw < best_hw) { best_hw = hw; memcpy(best_W, W, sizeof(W)); }

        /* SA шаг */
        int k = (int)(xrand() % 16);
        int bit = (int)(xrand() % 32);
        uint32_t saved = W[k];
        W[k] ^= (1u << bit);
        int new_hw = hw59_diff(W, ta->dw0);

        int delta = new_hw - hw;
        if (delta > 0) {
            double prob = exp(-delta / T);
            if ((xrand() & 0xFFFFFF) > (uint64_t)(prob * 0xFFFFFF))
                W[k] = saved;
        }

        T = T * 0.99998;
        if (T < 0.5) T = 0.5;

        /* Рестарт: из лучшей точки или заново */
        if (iters % 100000 == 0) {
            if (best_hw < 120) {
                memcpy(W, best_W, sizeof(W));
                /* Небольшая мутация */
                int nb = (int)(xrand()%3)+1;
                for (int m = 0; m < nb; m++) {
                    int kk = (int)(xrand()%16);
                    int bb = (int)(xrand()%32);
                    W[kk] ^= (1u<<bb);
                }
                T = 25.0;
            } else {
                for (int i = 0; i < 16; i++) W[i] = (uint32_t)xrand();
                T = 40.0;
            }
        }
    }

    /* Мерж в шаред */
    pthread_mutex_lock(ta->mutex);
    for (int b = 0; b < N_BUCKETS; b++) {
        CMatrix *dst = ta->shared[b];
        for (int i = 0; i < 64; i++) {
            dst->si[i] += local[b].si[i];
            for (int j = i; j < 64; j++)
                dst->sij[i][j] += local[b].sij[i][j];
        }
        dst->n += local[b].n;
    }
    pthread_mutex_unlock(ta->mutex);

    return NULL;
}

int main(int argc, char *argv[]) {
    long target = (argc > 1) ? atol(argv[1]) : 300;
    int  nthreads = (argc > 2) ? atoi(argv[2]) : 4;
    uint32_t dw0  = 0x80000000;

    printf("Carry Matrix: target_deep=%ld (hw59<90), threads=%d\n\n", target, nthreads);
    printf("Вопросы:\n");
    printf("  Q-E: Совпадает ли лаговая структура корреляций с шагом 7 расписания?\n");
    printf("  Q-F: Какова 'энергия' матрицы как функция hw59?\n");
    printf("  Q-G: Какие раунды наиболее коррелированы при hw59<90?\n\n");

    CMatrix shared[N_BUCKETS];
    memset(shared, 0, sizeof(shared));
    pthread_mutex_t mx = PTHREAD_MUTEX_INITIALIZER;

    pthread_t threads[32];
    TArg args[32];

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    for (int t = 0; t < nthreads; t++) {
        args[t].tid         = t;
        args[t].target_deep = target / nthreads + 1;
        args[t].seed        = 0xABCDEF01ULL*(t+1) + 0x12345678;
        args[t].dw0         = dw0;
        args[t].mutex       = &mx;
        for (int b = 0; b < N_BUCKETS; b++)
            args[t].shared[b] = &shared[b];
        pthread_create(&threads[t], NULL, tworker, &args[t]);
    }
    for (int t = 0; t < nthreads; t++) pthread_join(threads[t], NULL);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
    printf("Время: %.1f сек\n\n", elapsed);

    /* ===== РЕЗУЛЬТАТЫ ===== */

    printf("╔════════════════════════════════════════════════════════╗\n");
    printf("║         CARRY MATRIX — Полный анализ корреляций       ║\n");
    printf("╚════════════════════════════════════════════════════════╝\n\n");

    /* Образцы по bucket'ам */
    printf("📊 Образцы по уровням:\n");
    for (int b = 0; b < N_BUCKETS; b++)
        printf("  %-12s: N=%ld\n", BUCKET_NAME[b], shared[b].n);
    printf("\n");

    /* Энергия как функция hw59 — Q-F */
    printf("⚡ Корреляционная энергия (сумма |phi| всех пар):\n");
    printf("  Уровень      | N       | Энергия | Нормир. | Фактор↑\n");
    double prev_E = -1;
    for (int b = N_BUCKETS-1; b >= 0; b--) {
        if (shared[b].n < 10) continue;
        double E = cmat_energy(&shared[b]);
        double E_norm = E / (64*63/2);  /* на пару */
        double factor = (prev_E > 0) ? E/prev_E : 1.0;
        printf("  %-12s| %-7ld | %7.1f | %7.5f | %s\n",
            BUCKET_NAME[b], shared[b].n, E, E_norm,
            (factor > 1.5 && prev_E > 0) ? "*** РОСТ ***" :
            (factor > 1.2 && prev_E > 0) ? "↑" : "");
        prev_E = E;
    }
    printf("\n");

    /* Топ корреляций по каждому уровню */
    for (int b = N_BUCKETS-1; b >= 0; b -= 2) {
        if (shared[b].n < 20) continue;
        printf("🔝 Топ-15 пар при [%s] (N=%ld):\n", BUCKET_NAME[b], shared[b].n);

        /* Найти топ-15 */
        typedef struct { int i,j; double phi; } Pair;
        Pair top[2048];
        int ntop = 0;
        for (int i = 0; i < 64; i++)
            for (int j = i+1; j < 64; j++) {
                top[ntop].i = i; top[ntop].j = j;
                top[ntop].phi = cmat_phi(&shared[b], i, j);
                ntop++;
            }
        /* Сортировка по |phi| */
        for (int a = 0; a < ntop-1 && a < 15; a++)
            for (int bb = a+1; bb < ntop; bb++)
                if (fabs(top[bb].phi) > fabs(top[a].phi)) {
                    Pair tmp = top[a]; top[a] = top[bb]; top[bb] = tmp;
                }
        for (int k = 0; k < 15 && k < ntop; k++) {
            int lag = top[k].j - top[k].i;
            printf("  c[%2d]↔c[%2d] lag=%2d: phi=%+.4f%s\n",
                top[k].i, top[k].j, lag, top[k].phi,
                (fabs(top[k].phi)>0.3 ? " ★★★" :
                 fabs(top[k].phi)>0.1 ? " ★★" :
                 fabs(top[k].phi)>0.05 ? " ★" : ""));
        }
        printf("\n");
    }

    /* Лаговая структура — Q-E */
    printf("📐 Лаговая структура (ответ на Q-E: lag=7?):\n");
    lag_analysis(&shared[7], "ALL");         /* random */
    lag_analysis(&shared[3], "hw59<100");    /* глубокий */
    if (shared[1].n > 10)
        lag_analysis(&shared[1], "hw59<90");

    /* Зонный анализ — Q-H: где живут корреляции? */
    printf("\n🗺  Зонный анализ (avg |phi| в зонах):\n");
    printf("  Уровень       | Zone1(0..17) | Zone2(18..58) | Zone3(59..63)\n");
    for (int b = N_BUCKETS-1; b >= 2; b--) {
        if (shared[b].n < 20) continue;
        double z1=0,z2=0,z3=0;
        int n1=0,n2=0,n3=0;
        for (int i = 0; i < 64; i++)
            for (int j = i+1; j < 64; j++) {
                double p = fabs(cmat_phi(&shared[b],i,j));
                /* Обе точки в одной зоне */
                if (i<18 && j<18) { z1+=p; n1++; }
                else if (i>=18 && j<=58) { z2+=p; n2++; }
                else if (i>=59 && j>=59) { z3+=p; n3++; }
            }
        printf("  %-14s| %12.5f | %13.5f | %13.5f\n",
            BUCKET_NAME[b],
            n1>0?z1/n1:0, n2>0?z2/n2:0, n3>0?z3/n3:0);
    }

    /* Вывод корреляционной матрицы для глубокого уровня */
    if (shared[3].n > 30) {
        printf("\n📋 Корреляционная матрица carry[20..55] при hw59<100 (только |phi|>0.1):\n");
        printf("     ");
        for (int j = 20; j <= 55; j += 5) printf("%3d", j);
        printf("\n");
        for (int i = 20; i <= 55; i += 5) {
            printf("%3d: ", i);
            for (int j = 20; j <= 55; j += 5) {
                double p = cmat_phi(&shared[3], i, j);
                if (fabs(p) > 0.1) printf("%+.2f ", p);
                else printf(" .    ");
            }
            printf("\n");
        }
    }

    return 0;
}
