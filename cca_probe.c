/*
 * cca_probe.c — Conditional Correlation Algebra probe
 *
 * Исследование: при каком уровне hw59 включается корреляция carry?
 * Мы измеряем phi(carry[i], carry[j]) для разных популяций W:
 *   Level 0: случайные W (ожидаем phi ≈ 0, T_CARRY_NULL)
 *   Level 1: W с hw59 < 110
 *   Level 2: W с hw59 < 100
 *   Level 3: W с hw59 < 90
 *   Level 4: W с hw59 < 85
 *
 * Вопрос: есть ли ПОРОГ — уровень где phi резко вырастает?
 *
 * SHA-256 "стеклянный ящик": отслеживаем внутренние carry на каждом раунде.
 *
 * Компиляция: gcc -O3 -march=native -pthread -o cca_probe cca_probe.c -lm
 * Запуск: ./cca_probe [N_samples] [threads]
 *         ./cca_probe 50000 4
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

/* ===== SHA-256 константы ===== */
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

/* ===== Структуры для стеклянного SHA-256 ===== */

typedef struct {
    uint32_t W[64];
    uint32_t a[65], b[65], c[65], d[65];
    uint32_t e[65], f[65], g[65], h[65];
    uint32_t T1[64], T2[64];
    uint8_t  carry[64];   /* carry[r] = overflow бит T1[r]+T2[r] */
    uint8_t  carry_T1[64]; /* overflow бит при вычислении T1 */
    uint32_t H[8];
} SHA256_Glass;

static void sha256_glass(const uint32_t W_in[16], SHA256_Glass *g) {
    /* Расширение расписания */
    for (int i = 0; i < 16; i++) g->W[i] = W_in[i];
    for (int i = 16; i < 64; i++)
        g->W[i] = s1(g->W[i-2]) + g->W[i-7] + s0(g->W[i-15]) + g->W[i-16];

    /* Начальные значения */
    g->a[0]=IV[0]; g->b[0]=IV[1]; g->c[0]=IV[2]; g->d[0]=IV[3];
    g->e[0]=IV[4]; g->f[0]=IV[5]; g->g[0]=IV[6]; g->h[0]=IV[7];

    for (int r = 0; r < 64; r++) {
        /* T1 = h + S1(e) + Ch(e,f,g) + K[r] + W[r] */
        uint64_t t1_sum = (uint64_t)g->h[r] + S1(g->e[r])
                        + Ch(g->e[r],g->f[r],g->g[r])
                        + K[r] + g->W[r];
        g->T1[r]       = (uint32_t)(t1_sum & 0xFFFFFFFF);
        g->carry_T1[r] = (t1_sum >> 32) & 1;

        /* T2 = S0(a) + Maj(a,b,c) */
        g->T2[r] = S0(g->a[r]) + Maj(g->a[r],g->b[r],g->c[r]);

        /* Новые регистры */
        uint64_t sum_ae = (uint64_t)g->d[r] + g->T1[r];
        uint64_t sum_aa = (uint64_t)g->T1[r] + g->T2[r];

        g->carry[r] = (sum_aa >> 32) & 1;  /* overflow T1+T2 */

        g->h[r+1] = g->g[r];
        g->g[r+1] = g->f[r];
        g->f[r+1] = g->e[r];
        g->e[r+1] = (uint32_t)(sum_ae & 0xFFFFFFFF);
        g->d[r+1] = g->c[r];
        g->c[r+1] = g->b[r];
        g->b[r+1] = g->a[r];
        g->a[r+1] = (uint32_t)(sum_aa & 0xFFFFFFFF);
    }

    /* Финальный хэш */
    g->H[0] = IV[0]+g->a[64]; g->H[1] = IV[1]+g->b[64];
    g->H[2] = IV[2]+g->c[64]; g->H[3] = IV[3]+g->d[64];
    g->H[4] = IV[4]+g->e[64]; g->H[5] = IV[5]+g->f[64];
    g->H[6] = IV[6]+g->g[64]; g->H[7] = IV[7]+g->h[64];
}

/* ===== Вычисление hw59 для пары (W, W+ΔW) ===== */
static int compute_hw59(const uint32_t W[16], uint32_t dw0) {
    SHA256_Glass g0, g1;
    uint32_t W1[16];
    for (int i = 0; i < 16; i++) W1[i] = W[i];
    W1[0] ^= dw0;  /* XOR-дифференциал */

    sha256_glass(W, &g0);
    sha256_glass(W1, &g1);

    /* hw59 = HW(δstate на раунде 59) */
    int hw = 0;
    uint32_t regs0[8] = {g0.a[59],g0.b[59],g0.c[59],g0.d[59],
                          g0.e[59],g0.f[59],g0.g[59],g0.h[59]};
    uint32_t regs1[8] = {g1.a[59],g1.b[59],g1.c[59],g1.d[59],
                          g1.e[59],g1.f[59],g1.g[59],g1.h[59]};
    for (int i = 0; i < 8; i++)
        hw += __builtin_popcount(regs0[i] ^ regs1[i]);
    return hw;
}

/* ===== XORSHIFT64 PRNG ===== */
static uint64_t xorshift64(uint64_t *state) {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    return *state;
}

/* ===== Статистика корреляции Phi ===== */
#define N_CARRY_PAIRS 20  /* Пары раундов из Zone-2 сессии + новые */

static const int CARRY_I[N_CARRY_PAIRS] = {
    33,45,25,50, 20,28,35,40, 48,55, 10,15,20,30, 33,45,25,55, 5,60
};
static const int CARRY_J[N_CARRY_PAIRS] = {
    40,50,45,59, 28,35,40,48, 55,59, 20,25,30,40, 55,59,55,59, 15,63
};

#define N_LEVELS 5
static const char *LEVEL_NAMES[N_LEVELS] = {
    "RANDOM(all)", "hw59<110", "hw59<100", "hw59<90", "hw59<85"
};
static const int LEVEL_THRESHOLDS[N_LEVELS] = {
    9999, 110, 100, 90, 85
};

typedef struct {
    double sum_xi[N_CARRY_PAIRS];   /* Σ carry[i] */
    double sum_xj[N_CARRY_PAIRS];   /* Σ carry[j] */
    double sum_xij[N_CARRY_PAIRS];  /* Σ carry[i]*carry[j] */
    double sum_xi2[N_CARRY_PAIRS];
    double sum_xj2[N_CARRY_PAIRS];
    long   n[N_CARRY_PAIRS];
} PhiStats;

static void phi_update(PhiStats *s, const uint8_t carry[64]) {
    for (int p = 0; p < N_CARRY_PAIRS; p++) {
        int ci = carry[CARRY_I[p]];
        int cj = carry[CARRY_J[p]];
        s->sum_xi[p]  += ci;
        s->sum_xj[p]  += cj;
        s->sum_xij[p] += ci * cj;
        s->sum_xi2[p] += ci;
        s->sum_xj2[p] += cj;
        s->n[p]++;
    }
}

static double phi_compute(const PhiStats *s, int p) {
    long n = s->n[p];
    if (n < 2) return 0;
    double mean_i = s->sum_xi[p] / n;
    double mean_j = s->sum_xj[p] / n;
    double mean_ij = s->sum_xij[p] / n;
    double cov = mean_ij - mean_i * mean_j;
    double var_i = mean_i * (1 - mean_i);
    double var_j = mean_j * (1 - mean_j);
    if (var_i < 1e-12 || var_j < 1e-12) return 0;
    return cov / sqrt(var_i * var_j);
}

/* ===== Данные уровня ===== */
typedef struct {
    PhiStats phi[N_LEVELS];
    long     count[N_LEVELS];
    /* Для анализа "когда переключается" — hw59 vs max|phi| */
    double   hw59_phi_sum[256];   /* сумма max|phi| по всем парам, индекс = hw59 */
    long     hw59_count[256];
    /* Глобальная матрица корреляций carry-carry для level 0 vs level 4 */
    double   cmat_rand[64][64];   /* carry[i] x carry[j] для random */
    double   cmat_deep[64][64];   /* carry[i] x carry[j] для hw59<85 */
    long     cmat_rand_n, cmat_deep_n;
    double   cmat_rand_s[64][64], cmat_rand_si[64], cmat_rand_sj[64];
    double   cmat_deep_s[64][64], cmat_deep_si[64], cmat_deep_sj[64];
    long     cmat_rand_cnt, cmat_deep_cnt;
} LevelData;

/* ===== Поток: генерирует W, оценивает hw59, накапливает статистику ===== */
typedef struct {
    int thread_id;
    long n_samples;
    uint64_t seed;
    uint32_t dw0;
    LevelData *data;
    pthread_mutex_t *mutex;
} ThreadArg;

static void *worker(void *arg) {
    ThreadArg *ta = (ThreadArg*)arg;
    uint64_t rng = ta->seed;
    SHA256_Glass g;
    uint32_t W[16];

    LevelData local;
    memset(&local, 0, sizeof(local));

    long target = ta->n_samples;
    long collected_deep = 0;
    long iters = 0;
    long max_iters = target * 500;  /* не зависать если hw<85 редки */

    /* SA для генерации глубоких W */
    /* Сначала собираем random, потом SA-guided */
    double T_sa = 30.0;
    uint32_t best_W[16];
    int best_hw = 256;

    /* Инициализация случайного W */
    for (int i = 0; i < 16; i++)
        W[i] = (uint32_t)xorshift64(&rng);

    while (iters < max_iters) {
        iters++;

        /* Вычислить hw59 */
        int hw = compute_hw59(W, ta->dw0);

        /* Обновить статистику уровней */
        sha256_glass(W, &g);

        for (int lv = 0; lv < N_LEVELS; lv++) {
            if (hw < LEVEL_THRESHOLDS[lv]) {
                phi_update(&local.phi[lv], g.carry);
                local.count[lv]++;
            }
        }

        /* hw59 -> phi профиль (будем заполнять позже по max|phi|) */
        if (hw < 256) {
            local.hw59_count[hw]++;
        }

        /* Матрица корреляций carry[i] x carry[j] */
        if (hw >= 110) {  /* random уровень */
            for (int i = 0; i < 64; i++) {
                local.cmat_rand_si[i] += g.carry[i];
                for (int j = i; j < 64; j++) {
                    local.cmat_rand_s[i][j] += g.carry[i] * g.carry[j];
                }
            }
            local.cmat_rand_cnt++;
        }
        if (hw < 90) {  /* глубокий уровень */
            for (int i = 0; i < 64; i++) {
                local.cmat_deep_si[i] += g.carry[i];
                for (int j = i; j < 64; j++) {
                    local.cmat_deep_s[i][j] += g.carry[i] * g.carry[j];
                }
            }
            local.cmat_deep_cnt++;
            collected_deep++;
        }

        /* SA: мутация */
        if (best_hw > hw) {
            best_hw = hw;
            memcpy(best_W, W, sizeof(W));
        }

        /* SA шаг: мутировать слово W[k] бит */
        int k = (int)(xorshift64(&rng) % 12);  /* первые 12 слов */
        int b = (int)(xorshift64(&rng) % 32);
        uint32_t old_k = W[k];
        W[k] ^= (1u << b);
        int new_hw = compute_hw59(W, ta->dw0);

        double delta = new_hw - hw;
        if (delta > 0) {
            /* Принять с вероятностью exp(-delta/T) */
            double prob = exp(-delta / T_sa);
            double r01 = (xorshift64(&rng) & 0xFFFFFFF) / (double)0xFFFFFFF;
            if (r01 > prob) W[k] = old_k;  /* откат */
        }
        T_sa *= 0.99999;
        if (T_sa < 0.1) T_sa = 0.1;

        /* Перезапуск если застряли */
        if (iters % 50000 == 0) {
            if (best_hw < 130) {
                memcpy(W, best_W, sizeof(W));
            } else {
                for (int i = 0; i < 16; i++)
                    W[i] = (uint32_t)xorshift64(&rng);
            }
            T_sa = 20.0;
        }

        if (local.count[0] >= target) break;
    }

    /* Сохранить локальные данные */
    pthread_mutex_lock(ta->mutex);
    for (int lv = 0; lv < N_LEVELS; lv++) {
        ta->data->count[lv] += local.count[lv];
        for (int p = 0; p < N_CARRY_PAIRS; p++) {
            ta->data->phi[lv].sum_xi[p]  += local.phi[lv].sum_xi[p];
            ta->data->phi[lv].sum_xj[p]  += local.phi[lv].sum_xj[p];
            ta->data->phi[lv].sum_xij[p] += local.phi[lv].sum_xij[p];
            ta->data->phi[lv].n[p]       += local.phi[lv].n[p];
        }
    }
    for (int h = 0; h < 256; h++) {
        ta->data->hw59_count[h] += local.hw59_count[h];
    }
    /* Матрицы */
    for (int i = 0; i < 64; i++) {
        ta->data->cmat_rand_si[i] += local.cmat_rand_si[i];
        ta->data->cmat_deep_si[i] += local.cmat_deep_si[i];
        for (int j = i; j < 64; j++) {
            ta->data->cmat_rand_s[i][j] += local.cmat_rand_s[i][j];
            ta->data->cmat_deep_s[i][j] += local.cmat_deep_s[i][j];
        }
    }
    ta->data->cmat_rand_cnt += local.cmat_rand_cnt;
    ta->data->cmat_deep_cnt += local.cmat_deep_cnt;
    pthread_mutex_unlock(ta->mutex);

    return NULL;
}

/* ===== Вывод результатов ===== */
static void print_results(LevelData *data) {
    printf("\n");
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║          CCA PROBE — Conditional Correlation Algebra            ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");

    /* Количество образцов на каждом уровне */
    printf("📊 Образцы по уровням:\n");
    for (int lv = 0; lv < N_LEVELS; lv++) {
        printf("  Level %d [%s]: N=%ld\n", lv, LEVEL_NAMES[lv], data->count[lv]);
    }
    printf("\n");

    /* Таблица phi по уровням */
    printf("┌─────────────────────────────────────────────────────────────────────────────┐\n");
    printf("│  Пара         │ L0:RAND  │ L1:<110  │ L2:<100  │ L3:<90   │ L4:<85   │\n");
    printf("├─────────────────────────────────────────────────────────────────────────────┤\n");

    double max_phi[N_LEVELS] = {0};
    double sum_abs_phi[N_LEVELS] = {0};
    int    nonzero_phi[N_LEVELS] = {0};

    for (int p = 0; p < N_CARRY_PAIRS; p++) {
        printf("│ c[%2d]↔c[%2d] │", CARRY_I[p], CARRY_J[p]);
        for (int lv = 0; lv < N_LEVELS; lv++) {
            double ph = phi_compute(&data->phi[lv], p);
            printf(" %+7.4f  │", ph);
            double aph = fabs(ph);
            if (aph > max_phi[lv]) max_phi[lv] = aph;
            sum_abs_phi[lv] += aph;
            if (aph > 0.020) nonzero_phi[lv]++;
        }
        printf("\n");
    }
    printf("└─────────────────────────────────────────────────────────────────────────────┘\n\n");

    printf("📈 Статистика корреляций:\n");
    printf("  Уровень       │ max|phi| │ sum|phi| │ >0.020 пар\n");
    printf("  ──────────────┼──────────┼──────────┼───────────\n");
    for (int lv = 0; lv < N_LEVELS; lv++) {
        printf("  %-14s│  %6.4f  │  %6.4f  │  %d/%d\n",
            LEVEL_NAMES[lv], max_phi[lv], sum_abs_phi[lv],
            nonzero_phi[lv], N_CARRY_PAIRS);
    }
    printf("\n");

    /* Вопрос: где порог? */
    printf("🔬 Анализ порога фазового перехода:\n");
    double prev = max_phi[0];
    for (int lv = 1; lv < N_LEVELS; lv++) {
        double ratio = (prev > 1e-6) ? max_phi[lv]/prev : 0;
        printf("  L%d→L%d: max|phi| %6.4f → %6.4f  (×%.2f)%s\n",
            lv-1, lv, prev, max_phi[lv], ratio,
            (ratio > 2.0) ? "  ← ПОРОГ!" : "");
        prev = max_phi[lv];
    }
    printf("\n");

    /* Матрица корреляций: самые сильные пары при hw59<90 */
    printf("🗺  Топ-10 сильных корреляций carry[i]↔carry[j] при hw59<90:\n");

    double top_phi[64*64];
    int    top_i[64*64], top_j[64*64];
    int    ntop = 0;
    long   n = data->cmat_deep_cnt;

    if (n > 100) {
        for (int i = 0; i < 64; i++) {
            double mi = data->cmat_deep_si[i] / n;
            double vi = mi * (1-mi);
            if (vi < 1e-10) continue;
            for (int j = i+1; j < 64; j++) {
                double mj = data->cmat_deep_si[j] / n;
                double vj = mj * (1-mj);
                if (vj < 1e-10) continue;
                double mij = data->cmat_deep_s[i][j] / n;
                double phi = (mij - mi*mj) / sqrt(vi*vj);
                top_phi[ntop] = phi;
                top_i[ntop] = i;
                top_j[ntop] = j;
                ntop++;
            }
        }
        /* Простая сортировка топ-10 по |phi| */
        for (int a = 0; a < ntop-1 && a < 10; a++) {
            for (int b = a+1; b < ntop; b++) {
                if (fabs(top_phi[b]) > fabs(top_phi[a])) {
                    double tp = top_phi[a]; top_phi[a] = top_phi[b]; top_phi[b] = tp;
                    int ti = top_i[a]; top_i[a] = top_i[b]; top_i[b] = ti;
                    int tj = top_j[a]; top_j[a] = top_j[b]; top_j[b] = tj;
                }
            }
        }
        for (int k = 0; k < 10 && k < ntop; k++) {
            printf("  carry[%2d] ↔ carry[%2d]: phi=%+.4f\n",
                top_i[k], top_j[k], top_phi[k]);
        }
    } else {
        printf("  (недостаточно образцов hw59<90)\n");
    }

    /* Сравнение с random */
    printf("\n🗺  Топ-10 корреляций carry[i]↔carry[j] при RANDOM W:\n");
    ntop = 0;
    n = data->cmat_rand_cnt;
    if (n > 100) {
        for (int i = 0; i < 64; i++) {
            double mi = data->cmat_rand_si[i] / n;
            double vi = mi * (1-mi);
            if (vi < 1e-10) continue;
            for (int j = i+1; j < 64; j++) {
                double mj = data->cmat_rand_si[j] / n;
                double vj = mj * (1-mj);
                if (vj < 1e-10) continue;
                double mij = data->cmat_rand_s[i][j] / n;
                double phi = (mij - mi*mj) / sqrt(vi*vj);
                top_phi[ntop] = phi;
                top_i[ntop] = i;
                top_j[ntop] = j;
                ntop++;
            }
        }
        for (int a = 0; a < ntop-1 && a < 10; a++) {
            for (int b = a+1; b < ntop; b++) {
                if (fabs(top_phi[b]) > fabs(top_phi[a])) {
                    double tp = top_phi[a]; top_phi[a] = top_phi[b]; top_phi[b] = tp;
                    int ti = top_i[a]; top_i[a] = top_i[b]; top_i[b] = ti;
                    int tj = top_j[a]; top_j[a] = top_j[b]; top_j[b] = tj;
                }
            }
        }
        for (int k = 0; k < 10 && k < ntop; k++) {
            printf("  carry[%2d] ↔ carry[%2d]: phi=%+.4f\n",
                top_i[k], top_j[k], top_phi[k]);
        }
    }

    /* Распределение hw59 */
    printf("\n📊 Распределение hw59 (выборочно):\n");
    printf("  hw59 | count\n");
    for (int h = 75; h <= 130; h += 5) {
        long cnt = 0;
        for (int hh = h; hh < h+5 && hh < 256; hh++) cnt += data->hw59_count[hh];
        printf("  %3d-%3d | %ld\n", h, h+4, cnt);
    }
}

int main(int argc, char *argv[]) {
    long n_samples = (argc > 1) ? atol(argv[1]) : 20000;
    int n_threads  = (argc > 2) ? atoi(argv[2]) : 4;
    uint32_t dw0   = 0x80000000;  /* Wang дифференциал */

    printf("CCA Probe: N=%ld, threads=%d, ΔW[0]=0x%08x\n\n", n_samples, n_threads, dw0);
    printf("Вопросы исследования:\n");
    printf("  Q-A: При каком hw59-пороге корреляция carry ВКЛЮЧАЕТСЯ?\n");
    printf("  Q-B: Какие пары (carry[i], carry[j]) наиболее коррелированы?\n");
    printf("  Q-C: Насколько сильна корреляция? Достаточно ли для использования?\n\n");

    LevelData data;
    memset(&data, 0, sizeof(data));

    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_t threads[16];
    ThreadArg args[16];

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    for (int t = 0; t < n_threads; t++) {
        args[t].thread_id = t;
        args[t].n_samples = n_samples / n_threads;
        args[t].seed      = 12345678ULL * (t+1) + 0xDEADBEEF;
        args[t].dw0       = dw0;
        args[t].data      = &data;
        args[t].mutex     = &mutex;
        pthread_create(&threads[t], NULL, worker, &args[t]);
    }
    for (int t = 0; t < n_threads; t++) pthread_join(threads[t], NULL);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
    printf("Время: %.1f сек\n", elapsed);

    print_results(&data);

    return 0;
}
