/*
 * П-27: Birthday search for 17-round differential pair
 * T_BIRTHDAY_COST17: δe2..δe17=0 за O(2^32) операций
 *
 * Алгоритм (T_CASCADE_17):
 *   Фиксируем W0, W1[1..15]=0.
 *   Перебираем W1[0] in [0, 2^32):
 *     1. ΔW2 = -De3_nat(W0,W1)  → De3=0
 *     2. ΔW3..ΔW15 = каскад     → De4..De16=0
 *     3. Проверяем Da13 + ΔW16 == 0  → De17=0!
 *
 * Компиляция: gcc -O3 -march=native -o birthday_search_17 birthday_search_17.c
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

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
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

static const u32 IV[8] = {
    0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
    0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19
};

#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))
#define SIG0(x)   (ROTR(x,2)  ^ ROTR(x,13) ^ ROTR(x,22))
#define SIG1(x)   (ROTR(x,6)  ^ ROTR(x,11) ^ ROTR(x,25))
#define sig0(x)   (ROTR(x,7)  ^ ROTR(x,18) ^ ((x)>>3))
#define sig1(x)   (ROTR(x,17) ^ ROTR(x,19) ^ ((x)>>10))
#define CH(e,f,g) (((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c)(((a)&(b))^((a)&(c))^((b)&(c)))

static inline void expand(u32 *W) {
    for (int i = 16; i < 18; i++)  /* только первые 2 расширенных нужны */
        W[i] = sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16];
    for (int i = 18; i < 64; i++)
        W[i] = sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16];
}

/* R раундов SHA-256, возвращает только e-регистр (индекс 4) после R раундов */
static inline u32 sha_e(const u32 *W, int R) {
    u32 a=IV[0],b=IV[1],c=IV[2],d=IV[3];
    u32 e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for (int r = 0; r < R; r++) {
        u32 T1 = h + SIG1(e) + CH(e,f,g) + K[r] + W[r];
        u32 T2 = SIG0(a) + MAJ(a,b,c);
        h=g; g=f; f=e; e=d+T1;
        d=c; c=b; b=a; a=T1+T2;
    }
    return e;
}

/* R раундов, возвращает a-регистр (индекс 0) */
static inline u32 sha_a(const u32 *W, int R) {
    u32 a=IV[0],b=IV[1],c=IV[2],d=IV[3];
    u32 e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for (int r = 0; r < R; r++) {
        u32 T1 = h + SIG1(e) + CH(e,f,g) + K[r] + W[r];
        u32 T2 = SIG0(a) + MAJ(a,b,c);
        h=g; g=f; f=e; e=d+T1;
        d=c; c=b; b=a; a=T1+T2;
    }
    return a;
}

/*
 * Вычислить De17 = Da13 + ΔW16 для пары (W0, W1) с ΔW0=1.
 * W2..W15 = 0. Адаптивный каскад ΔW2..ΔW15.
 */
static u32 compute_de17(u32 W0, u32 W1) {
    u32 Wn[64], Wf[64];
    u32 DW[16];
    memset(DW, 0, sizeof(DW));

    /* Нормальное сообщение */
    memset(Wn, 0, sizeof(Wn));
    Wn[0] = W0; Wn[1] = W1;
    expand(Wn);

    /* ΔW0=1 */
    DW[0] = 1;

    /* Шаг 0: ΔW2 → De3=0 */
    memset(Wf, 0, sizeof(Wf));
    for (int i = 0; i < 16; i++) Wf[i] = Wn[i] + DW[i];
    expand(Wf);
    u32 De3_nat = sha_e(Wf, 3) - sha_e(Wn, 3);
    DW[2] = -De3_nat;

    /* Каскад ΔW3..ΔW15: De4..De16=0 */
    for (int step = 0; step < 13; step++) {
        int wi = step+3, dt = step+4;
        memset(Wf, 0, sizeof(Wf));
        for (int i = 0; i < 16; i++) Wf[i] = Wn[i] + DW[i];
        expand(Wf);
        u32 Den_nat = sha_e(Wf, dt) - sha_e(Wn, dt);
        DW[wi] = -Den_nat;
    }

    /* Da13 и ΔW16 */
    memset(Wf, 0, sizeof(Wf));
    for (int i = 0; i < 16; i++) Wf[i] = Wn[i] + DW[i];
    expand(Wf);
    u32 Da13 = sha_a(Wf, 13) - sha_a(Wn, 13);
    u32 DW16  = Wf[16] - Wn[16];
    return Da13 + DW16;  /* = De17 */
}

/* Верификация: De3..De17 все равны 0 */
static void verify_pair(u32 W0, u32 W1) {
    u32 Wn[64], Wf[64];
    u32 DW[16];
    memset(DW, 0, sizeof(DW));
    memset(Wn, 0, sizeof(Wn));
    Wn[0]=W0; Wn[1]=W1; expand(Wn);
    DW[0]=1;

    memset(Wf, 0, sizeof(Wf));
    for (int i=0;i<16;i++) Wf[i]=Wn[i]+DW[i];
    expand(Wf);
    u32 De3_nat = sha_e(Wf,3)-sha_e(Wn,3);
    DW[2]=-De3_nat;

    for (int step=0;step<13;step++) {
        int wi=step+3, dt=step+4;
        memset(Wf,0,sizeof(Wf));
        for (int i=0;i<16;i++) Wf[i]=Wn[i]+DW[i];
        expand(Wf);
        DW[wi]=-(sha_e(Wf,dt)-sha_e(Wn,dt));
    }

    memset(Wf,0,sizeof(Wf));
    for (int i=0;i<16;i++) Wf[i]=Wn[i]+DW[i];
    expand(Wf);

    printf("Верификация De3..De17:\n");
    int all_ok = 1;
    for (int r=3;r<=17;r++) {
        u32 den=sha_e(Wn,r), def=sha_e(Wf,r);
        u32 dr=(def-den);
        printf("  De%2d = 0x%08x %s\n", r, dr, dr==0?"✓":"✗ FAIL");
        if (dr!=0) all_ok=0;
    }
    printf("Результат: %s\n", all_ok?"ВСЕ 15 НУЛЕЙ ПОДТВЕРЖДЕНЫ ✓":"ОШИБКА ✗");

    printf("\nΔW[0..15]:\n");
    for (int i=0;i<16;i++) printf("  ΔW[%2d] = 0x%08x\n", i, DW[i]);

    /* De18 */
    u32 Da14=sha_a(Wf,14)-sha_a(Wn,14);
    u32 DW17=Wf[17]-Wn[17];
    u32 De18=sha_e(Wf,18)-sha_e(Wn,18);
    printf("\nDe18 = Da14 + ΔW17 = 0x%08x + 0x%08x = 0x%08x\n",
           Da14, DW17, (Da14+DW17));
    printf("De18 (прямой): 0x%08x %s\n", De18, (Da14+DW17)==De18?"✓":"✗");
}

int main(int argc, char *argv[]) {
    setvbuf(stdout, NULL, _IONBF, 0);  /* отключить буферизацию */

    printf("=== P-27: Birthday Search de2..de17=0 ===\n");
    printf("T_BIRTHDAY_COST17: O(2^32) W1 iterations\n\n");

    u32 W0 = 0xc5bde324;
    if (argc >= 2) W0 = (u32)strtoul(argv[1], NULL, 16);

    /* Если указана пара — верифицировать */
    if (argc >= 3) {
        u32 W1 = (u32)strtoul(argv[2], NULL, 16);
        printf("Верификация пары: W0=0x%08x W1=0x%08x\n", W0, W1);
        verify_pair(W0, W1);
        return 0;
    }

    printf("W0 = 0x%08x (фиксирован)\n", W0);
    printf("W1[1..15] = 0 (фиксированы)\n");
    printf("Перебор W1[0] in [0, 2^32)...\n\n");

    /* Benchmark: скорость */
    clock_t tb = clock();
    volatile u32 dummy = 0;
    for (u32 i = 0; i < 1000; i++) dummy += compute_de17(W0, i);
    double bench_sec = (double)(clock()-tb)/CLOCKS_PER_SEC;
    double speed = 1000.0 / bench_sec;
    printf("Benchmark (1000 итераций): %.2f сек\n", bench_sec);
    printf("Скорость: %.0f iter/sec = %.2f M/sec\n", speed, speed/1e6);
    double eta_sec = 4294967296.0 / speed;
    printf("ETA для 2^32: %.0f сек = %.1f мин = %.1f ч\n\n",
           eta_sec, eta_sec/60, eta_sec/3600);
    (void)dummy;

    /* Основной поиск */
    clock_t t0 = clock();
    u64 count = 0;
    u64 found = 0;

    for (u64 w1 = 0; w1 <= 0xFFFFFFFFULL; w1++) {
        u32 de17 = compute_de17(W0, (u32)w1);
        count++;

        if (de17 == 0) {
            found++;
            double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
            printf("\n*** НАЙДЕНА ПАРА! ***\n");
            printf("W0   = 0x%08x\n", W0);
            printf("W1   = 0x%08lx\n", (unsigned long)w1);
            printf("De17 = 0x%08x\n", de17);
            printf("Итераций: %llu (2^%.2f)\n", (unsigned long long)count,
                   count>1 ? 31.9 /* log2 approx */ : 0.0);
            printf("Время: %.1f сек\n\n", elapsed);
            verify_pair(W0, (u32)w1);
            if (found >= 3) { printf("Найдено 3 пары, стоп.\n"); break; }
        }

        if (count % 10000000ULL == 0) {
            double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
            double rate = (double)count / (elapsed > 0 ? elapsed : 1);
            double eta  = (4294967296.0 - (double)count) / rate;
            printf("[%6.2f%%]  %.2f M/s  elapsed=%.0fs  ETA=%.0fs  found=%llu\n",
                   100.0*count/4294967296.0, rate/1e6, elapsed, eta,
                   (unsigned long long)found);
        }
    }

    double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
    printf("\nИтого: %llu итераций за %.1f сек, найдено: %llu\n",
           (unsigned long long)count, elapsed, (unsigned long long)found);
    return 0;
}
