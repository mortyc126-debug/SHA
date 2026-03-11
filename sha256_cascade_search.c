/*
 * SHA-256 Дифференциальный каскадный поиск (C версия для скорости)
 * П-15: Поиск пары с De3..De17=0 (De17=0)
 *
 * Компиляция: gcc -O3 -o sha256_cascade_search sha256_cascade_search.c
 * Использование: ./sha256_cascade_search [seed] [iters]
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef uint32_t u32;
typedef uint64_t u64;

/* SHA-256 константы */
static const u32 K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
};

static const u32 H0[8] = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
};

#define ROTR(x, n) (((x) >> (n)) | ((x) << (32 - (n))))
#define sig0(x) (ROTR(x,7) ^ ROTR(x,18) ^ ((x)>>3))
#define sig1(x) (ROTR(x,17) ^ ROTR(x,19) ^ ((x)>>10))
#define Sig0(x) (ROTR(x,2) ^ ROTR(x,13) ^ ROTR(x,22))
#define Sig1(x) (ROTR(x,6) ^ ROTR(x,11) ^ ROTR(x,25))
#define Ch(e,f,g) (((e)&(f)) ^ (~(e)&(g)))
#define Maj(a,b,c) (((a)&(b)) ^ ((a)&(c)) ^ ((b)&(c)))

/* Расширение расписания */
void schedule(u32 *W, const u32 *W16) {
    for (int i = 0; i < 16; i++) W[i] = W16[i];
    for (int i = 16; i < 64; i++) {
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]);
    }
}

/* SHA-256: R раундов, возвращает e[R] */
u32 sha_e(const u32 *W, int R) {
    u32 a = H0[0], b = H0[1], c = H0[2], d = H0[3];
    u32 e = H0[4], f = H0[5], g = H0[6], h = H0[7];
    for (int r = 0; r < R; r++) {
        u32 T1 = h + Sig1(e) + Ch(e,f,g) + K[r] + W[r];
        u32 T2 = Sig0(a) + Maj(a,b,c);
        h = g; g = f; f = e; e = d + T1;
        d = c; c = b; b = a; a = T1 + T2;
    }
    return e;
}

/* SHA-256: R раундов, возвращает всё состояние */
void sha_state(const u32 *W, int R, u32 *state) {
    u32 a = H0[0], b = H0[1], c = H0[2], d = H0[3];
    u32 e = H0[4], f = H0[5], g = H0[6], h = H0[7];
    for (int r = 0; r < R; r++) {
        u32 T1 = h + Sig1(e) + Ch(e,f,g) + K[r] + W[r];
        u32 T2 = Sig0(a) + Maj(a,b,c);
        h = g; g = f; f = e; e = d + T1;
        d = c; c = b; b = a; a = T1 + T2;
    }
    state[0] = a; state[1] = b; state[2] = c; state[3] = d;
    state[4] = e; state[5] = f; state[6] = g; state[7] = h;
}

/*
 * Каскадная функция (П-13 алгоритм):
 * 1. ΔW0=1, ΔW1=0
 * 2. ΔW2 = -De3_nat (адаптивный, De3=0 гарантировано)
 * 3. ΔW3..ΔW15 = каскад (De4..De16=0)
 * Возвращает De17 (и опционально заполняет DWs)
 */
u32 cascade_de17(u32 W0, u32 W1, u32 *DWs_out, u32 *Wn_out, u32 *Wf_out) {
    u32 Wn16[16] = {0}, Wf16[16] = {0};
    u32 DWs[16] = {0};
    u32 Wn[64], Wf[64];

    Wn16[0] = W0; Wn16[1] = W1;
    DWs[0] = 1; /* ΔW0=1 */

    /* Шаг 0: адаптивный ΔW2 → De3=0 */
    for (int i = 0; i < 16; i++) Wf16[i] = Wn16[i] + DWs[i];
    schedule(Wn, Wn16); schedule(Wf, Wf16);
    u32 de3_nat = sha_e(Wf, 3) - sha_e(Wn, 3);
    DWs[2] = (u32)(-(int32_t)de3_nat);

    /* Шаги 1..13: каскад De4..De16=0 */
    for (int step = 0; step < 13; step++) {
        int wi = step + 3;
        int dt = step + 4;
        for (int i = 0; i < 16; i++) Wf16[i] = Wn16[i] + DWs[i];
        schedule(Wn, Wn16); schedule(Wf, Wf16);
        u32 de_dt = sha_e(Wf, dt) - sha_e(Wn, dt);
        DWs[wi] = (u32)(-(int32_t)de_dt);
    }

    /* Вычислить De17 */
    for (int i = 0; i < 16; i++) Wf16[i] = Wn16[i] + DWs[i];
    schedule(Wn, Wn16); schedule(Wf, Wf16);
    u32 de17 = sha_e(Wf, 17) - sha_e(Wn, 17);

    if (DWs_out) memcpy(DWs_out, DWs, sizeof(DWs));
    if (Wn_out) memcpy(Wn_out, Wn, 64 * sizeof(u32));
    if (Wf_out) memcpy(Wf_out, Wf, 64 * sizeof(u32));

    return de17;
}

/* Быстрый xorshift64 для ГПСЧ */
static u64 rng_state;
static u32 fast_rand32() {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (u32)(rng_state ^ (rng_state >> 32));
}

int main(int argc, char *argv[]) {
    u64 seed = (argc > 1) ? (u64)atoll(argv[1]) : (u64)time(NULL);
    u64 max_iters = (argc > 2) ? (u64)atoll(argv[2]) : (u64)4000000000ULL;

    rng_state = seed ? seed : 0xdeadbeefcafe1234ULL;

    printf("=================================================================\n");
    printf("П-15: Поиск пары с De3..De17=0 (C реализация)\n");
    printf("=================================================================\n");
    printf("Seed: %llu, Max iterations: %llu\n", (unsigned long long)seed, (unsigned long long)max_iters);
    printf("Ожидаемое кол-во попыток: ~2^32 = %llu\n\n", (unsigned long long)(1ULL << 32));

    u64 count = 0;
    clock_t t_start = clock();
    clock_t t_last_report = t_start;
    u32 found_W0 = 0, found_W1 = 0, found_de17 = 1;
    u32 found_DWs[16];
    u32 found_Wn[64], found_Wf[64];

    while (count < max_iters) {
        u32 W0 = fast_rand32();
        u32 W1 = fast_rand32();

        u32 de17 = cascade_de17(W0, W1, NULL, NULL, NULL);
        count++;

        if (de17 == 0) {
            /* Нашли! */
            found_W0 = W0; found_W1 = W1; found_de17 = 0;
            /* Пересчитать с сохранением данных */
            cascade_de17(W0, W1, found_DWs, found_Wn, found_Wf);
            break;
        }

        /* Прогресс каждые 50M итераций */
        if (count % 50000000ULL == 0) {
            clock_t t_now = clock();
            double elapsed = (double)(t_now - t_start) / CLOCKS_PER_SEC;
            double rate = (double)count / elapsed / 1e6;
            double eta = ((double)(1ULL << 32) - (double)count) / (rate * 1e6);
            printf("  [%llu M / %.1f M/s] elapsed=%.1fs, ETA=%.0fs (~%.1fh)\n",
                   (unsigned long long)(count / 1000000ULL),
                   rate, elapsed, eta, eta / 3600.0);
            fflush(stdout);
        }
    }

    clock_t t_end = clock();
    double elapsed = (double)(t_end - t_start) / CLOCKS_PER_SEC;
    double rate = (double)count / elapsed / 1e6;

    printf("\n");
    printf("=================================================================\n");
    if (found_de17 == 0) {
        printf("УСПЕХ! De17=0 найдено!\n");
        printf("=================================================================\n");
        printf("W0  = 0x%08x\n", found_W0);
        printf("W1  = 0x%08x\n", found_W1);
        printf("DWs = [");
        for (int i = 0; i < 16; i++) printf("0x%x%s", found_DWs[i], i<15?", ":"");
        printf("]\n");

        /* Верификация всех De3..De17 */
        printf("\nВерификация De3..De17:\n");
        u32 Wn16[16] = {0}, Wf16[16] = {0};
        Wn16[0] = found_W0; Wn16[1] = found_W1;
        for (int i = 0; i < 16; i++) Wf16[i] = Wn16[i] + found_DWs[i];
        u32 Wn[64], Wf[64];
        schedule(Wn, Wn16); schedule(Wf, Wf16);
        int all_ok = 1;
        for (int r = 3; r <= 17; r++) {
            u32 de_r = sha_e(Wf, r) - sha_e(Wn, r);
            printf("  De%d = 0x%08x %s\n", r, de_r, de_r == 0 ? "✓" : "✗ ОШИБКА!");
            if (de_r != 0) all_ok = 0;
        }
        printf("\n%s\n", all_ok ? "Все De3..De17=0 ✓ T_CASCADE_17 ВЕРИФИЦИРОВАНА!" :
                                   "Ошибка верификации!");

        /* De18 и Da14 для П-15 */
        printf("\nТ_DE18_DECOMPOSITION анализ:\n");
        u32 de18 = sha_e(Wf, 18) - sha_e(Wn, 18);
        u32 state_n[8], state_f[8];
        sha_state(Wn, 14, state_n); sha_state(Wf, 14, state_f);
        u32 da14 = state_f[0] - state_n[0];
        u32 dw17 = Wf[17] - Wn[17];
        printf("  De18     = 0x%08x\n", de18);
        printf("  Da14     = 0x%08x\n", da14);
        printf("  ΔW17     = 0x%08x\n", dw17);
        printf("  Da14+ΔW17= 0x%08x\n", da14 + dw17);
        printf("  T_DE18_DECOMPOSITION: %s\n",
               (de18 == da14 + dw17) ? "De18 = Da14 + ΔW17 ✓ ВЕРИФИЦИРОВАНА!" :
                                       "De18 ≠ Da14 + ΔW17 ✗ НЕ выполняется!");

        /* De19..De21 и декомпозиция */
        printf("\nДекомпозиция De19..De21:\n");
        for (int k = 19; k <= 21; k++) {
            u32 de_k = sha_e(Wf, k) - sha_e(Wn, k);
            u32 state_nk[8], state_fk[8];
            sha_state(Wn, k-4, state_nk); sha_state(Wf, k-4, state_fk);
            u32 da_km4 = state_fk[0] - state_nk[0];
            u32 dw_km1 = Wf[k-1] - Wn[k-1];
            printf("  De%d = 0x%08x, Da%d+ΔW%d = 0x%08x %s\n",
                   k, de_k, k-4, k-1, da_km4+dw_km1,
                   (de_k == da_km4+dw_km1) ? "✓" : "✗");
        }

        /* Быстрый тест: P(De18=0 | De17=0) */
        printf("\nОценка P(De18=0):\n");
        printf("  De18 = 0x%08x (одна пара — не 0, ожидаемо)\n", de18);
        printf("  Для статистики нужно ≥2^32 пар с De17=0\n");

        /* Вывести полный результат в формате для Python */
        printf("\nДля Python (скопировать в p15):\n");
        printf("W0 = 0x%08x\n", found_W0);
        printf("W1 = 0x%08x\n", found_W1);
        printf("DWs = [");
        for (int i = 0; i < 16; i++) printf("0x%x%s", found_DWs[i], i<15?", ":"");
        printf("]\n");
        printf("Da14 = 0x%08x\n", da14);
        printf("ΔW17 = 0x%08x\n", dw17);
        printf("De18 = 0x%08x\n", de18);

    } else {
        printf("Поиск завершён без нахождения De17=0 (%llu итераций)\n",
               (unsigned long long)count);
        printf("P(не найти за N) = (1 - 2^-32)^N\n");
        if (count < (u64)(1ULL << 32)) {
            double p = 1.0 - (double)count / (double)(1ULL << 32);
            printf("P(не найти за %llu) ≈ %.4f%%\n",
                   (unsigned long long)count, p * 100.0);
        }
    }

    printf("\nСтатистика поиска:\n");
    printf("  Итераций: %llu\n", (unsigned long long)count);
    printf("  Время: %.1f сек\n", elapsed);
    printf("  Скорость: %.2f M/s\n", rate);
    printf("=================================================================\n");

    return (found_de17 == 0) ? 0 : 1;
}
