/*
 * П-27B: Birthday search — многопоточная версия (pthread)
 * Параллельные диапазоны W1[0] по потокам
 *
 * Компиляция: gcc -O3 -march=native -pthread -o birthday_search_mt birthday_search_mt.c
 * Запуск:     ./birthday_search_mt [W0_hex] [N_threads]
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
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
    for (int i = 16; i < 64; i++)
        W[i] = sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16];
}

static inline u32 sha_e(const u32 *W, int R) {
    u32 a=IV[0],b=IV[1],c=IV[2],d=IV[3];
    u32 e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for (int r = 0; r < R; r++) {
        u32 T1 = h + SIG1(e) + CH(e,f,g) + K[r] + W[r];
        u32 T2 = SIG0(a) + MAJ(a,b,c);
        h=g; g=f; f=e; e=d+T1; d=c; c=b; b=a; a=T1+T2;
    }
    return e;
}

static inline u32 sha_a(const u32 *W, int R) {
    u32 a=IV[0],b=IV[1],c=IV[2],d=IV[3];
    u32 e=IV[4],f=IV[5],g=IV[6],h=IV[7];
    for (int r = 0; r < R; r++) {
        u32 T1 = h + SIG1(e) + CH(e,f,g) + K[r] + W[r];
        u32 T2 = SIG0(a) + MAJ(a,b,c);
        h=g; g=f; f=e; e=d+T1; d=c; c=b; b=a; a=T1+T2;
    }
    return a;
}

static u32 compute_de17(u32 W0, u32 W1) {
    u32 Wn[64], Wf[64];
    u32 DW[16];
    memset(DW, 0, sizeof(DW));
    memset(Wn, 0, sizeof(Wn));
    Wn[0] = W0; Wn[1] = W1;
    expand(Wn);
    DW[0] = 1;

    /* ΔW2 → De3=0 */
    memset(Wf, 0, sizeof(Wf));
    for (int i = 0; i < 16; i++) Wf[i] = Wn[i] + DW[i];
    expand(Wf);
    DW[2] = -(sha_e(Wf, 3) - sha_e(Wn, 3));

    /* Каскад ΔW3..ΔW15 */
    for (int step = 0; step < 13; step++) {
        int wi = step+3, dt = step+4;
        memset(Wf, 0, sizeof(Wf));
        for (int i = 0; i < 16; i++) Wf[i] = Wn[i] + DW[i];
        expand(Wf);
        DW[wi] = -(sha_e(Wf, dt) - sha_e(Wn, dt));
    }

    memset(Wf, 0, sizeof(Wf));
    for (int i = 0; i < 16; i++) Wf[i] = Wn[i] + DW[i];
    expand(Wf);
    return (sha_a(Wf, 13) - sha_a(Wn, 13)) + (Wf[16] - Wn[16]);
}

/* ---- Многопоточная часть ---- */

#define MAX_FOUND 16

typedef struct {
    u32 W0;
    u64 range_start;
    u64 range_end;
    int thread_id;
    int n_threads;

    /* Результаты */
    u32 found_w1[MAX_FOUND];
    int found_count;
} ThreadArgs;

/* Глобальный счётчик итераций и флаг остановки */
static atomic_ullong g_total_count = 0;
static atomic_int    g_found_total = 0;
static volatile int  g_stop = 0;
static u32 g_W0;

static void *worker(void *arg) {
    ThreadArgs *a = (ThreadArgs *)arg;
    u64 local_count = 0;

    for (u64 w1 = a->range_start; w1 < a->range_end && !g_stop; w1++) {
        u32 de17 = compute_de17(a->W0, (u32)w1);
        local_count++;

        if (de17 == 0) {
            int idx = atomic_fetch_add(&g_found_total, 1);
            if (idx < MAX_FOUND) {
                a->found_w1[a->found_count++] = (u32)w1;
                printf("\n[Thread %d] *** НАЙДЕНА ПАРА #%d ***\n",
                       a->thread_id, idx+1);
                printf("  W0 = 0x%08x\n", a->W0);
                printf("  W1 = 0x%08lx\n", (unsigned long)w1);
                printf("  iter в этом потоке: %llu\n",
                       (unsigned long long)local_count);
                fflush(stdout);
            }
            if (atomic_load(&g_found_total) >= 3)
                g_stop = 1;
        }

        /* Обновляем глобальный счётчик каждые 1M итераций */
        if ((local_count & 0xFFFFF) == 0)
            atomic_fetch_add(&g_total_count, 0x100000);
    }

    /* Дозаписываем остаток */
    atomic_fetch_add(&g_total_count, local_count & 0xFFFFF);
    return NULL;
}

/* Прогресс-поток */
static time_t g_start_time;

static void *progress_thread(void *arg) {
    (void)arg;
    while (!g_stop) {
        sleep(10);
        u64 cnt = atomic_load(&g_total_count);
        int found = atomic_load(&g_found_total);
        double elapsed = difftime(time(NULL), g_start_time);
        double rate = (elapsed > 0) ? cnt / elapsed : 0;
        double pct = 100.0 * cnt / 4294967296.0;
        double eta = (rate > 0) ? (4294967296.0 - cnt) / rate : 0;
        printf("[%5.2f%%]  %.2f M/s  elapsed=%.0fs  ETA=%.0fs  found=%d\n",
               pct, rate/1e6, elapsed, eta, found);
        fflush(stdout);
        if (cnt >= 4294967296ULL) break;
    }
    return NULL;
}

int main(int argc, char *argv[]) {
    setvbuf(stdout, NULL, _IONBF, 0);

    g_W0 = 0xc5bde324;
    int n_threads = 4;

    if (argc >= 2) g_W0 = (u32)strtoul(argv[1], NULL, 16);
    if (argc >= 3) n_threads = atoi(argv[2]);
    if (n_threads < 1 || n_threads > 64) n_threads = 4;

    printf("=== P-27 MT: Birthday Search de2..de17=0 ===\n");
    printf("W0 = 0x%08x, потоков = %d\n\n", g_W0, n_threads);

    /* Benchmark */
    {
        clock_t tb = clock();
        volatile u32 d = 0;
        for (u32 i = 0; i < 1000; i++) d += compute_de17(g_W0, i);
        double s = (double)(clock()-tb)/CLOCKS_PER_SEC;
        double spd = 1000.0/s;
        printf("Benchmark 1 поток: %.0f iter/s = %.2f M/s\n", spd, spd/1e6);
        printf("ETA 1 поток:  %.0f сек = %.1f мин\n", 4294967296.0/spd, 4294967296.0/spd/60);
        printf("ETA %d потоков: %.0f сек = %.1f мин\n\n",
               n_threads, 4294967296.0/spd/n_threads, 4294967296.0/spd/n_threads/60);
        (void)d;
    }

    pthread_t threads[64];
    ThreadArgs args[64];
    u64 chunk = (u64)0x100000000ULL / n_threads;

    g_start_time = time(NULL);

    for (int i = 0; i < n_threads; i++) {
        memset(&args[i], 0, sizeof(args[i]));
        args[i].W0          = g_W0;
        args[i].range_start = (u64)i * chunk;
        args[i].range_end   = (i == n_threads-1) ? 0x100000000ULL
                                                  : (u64)(i+1) * chunk;
        args[i].thread_id   = i;
        args[i].n_threads   = n_threads;
        pthread_create(&threads[i], NULL, worker, &args[i]);
    }

    /* Прогресс-поток */
    pthread_t pt;
    pthread_create(&pt, NULL, progress_thread, NULL);

    for (int i = 0; i < n_threads; i++)
        pthread_join(threads[i], NULL);
    g_stop = 1;
    pthread_join(pt, NULL);

    double elapsed = difftime(time(NULL), g_start_time);
    u64 total = atomic_load(&g_total_count);
    int found = atomic_load(&g_found_total);

    printf("\nИтого: %llu итераций за %.0f сек, найдено: %d\n",
           (unsigned long long)total, elapsed, found);
    if (total > 0 && elapsed > 0)
        printf("Скорость: %.2f M/s (суммарно %d потоков)\n",
               (double)total/elapsed/1e6, n_threads);

    /* Вывод найденных пар */
    for (int i = 0; i < n_threads; i++) {
        for (int j = 0; j < args[i].found_count; j++) {
            printf("\nПара: W0=0x%08x W1=0x%08x\n", g_W0, args[i].found_w1[j]);
        }
    }
    return 0;
}
