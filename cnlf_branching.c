/*
 * CARRY×NLF BRANCHING: сколько РЕАЛЬНО solutions у x + Ch(x,c1,c2) = y?
 *
 * 14 bits free → worst case 2^14 solutions.
 * But carry constraints PRUNE most branches.
 * Real branching factor = ?
 *
 * If branching ≈ 2-4 per addition → solver with small tree.
 * If branching ≈ 2^14 → hopeless.
 *
 * gcc -O3 -march=native -o cnlf_branch cnlf_branching.c
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

static inline uint32_t Ch(uint32_t e, uint32_t f, uint32_t g) {
    return (e&f)^(~e&g);
}
static inline uint32_t Maj(uint32_t a, uint32_t b, uint32_t c) {
    return (a&b)^(a&c)^(b&c);
}

/* Count ALL solutions of x + Ch(x, c1, c2) = y (mod 2^32) */
/* Bit-by-bit solver with branching */
int count_solutions_ch(uint32_t y, uint32_t c1, uint32_t c2) {
    /* Stack-based bit-by-bit solver */
    typedef struct { int bit; int carry; uint32_t x; } State;
    State stack[4096];
    int sp = 0;
    int n_solutions = 0;

    stack[sp++] = (State){0, 0, 0};

    while(sp > 0) {
        State s = stack[--sp];
        if(s.bit == 32) {
            if(s.carry == 0) n_solutions++;
            continue;
        }

        int k = s.bit;
        int yk = (y >> k) & 1;
        int c1k = (c1 >> k) & 1;
        int c2k = (c2 >> k) & 1;
        int mk = c1k ^ c2k; /* mask: does x[k] appear in Ch? */
        int vk = c2k;       /* Ch value when x[k]=0 */

        for(int xk = 0; xk <= 1; xk++) {
            int chk = (mk & xk) ^ vk; /* Ch(x,c1,c2)[k] */
            int sum = xk + chk + s.carry;
            int result = sum & 1;
            int new_carry = sum >> 1;

            if(result == yk) {
                if(sp < 4095) {
                    stack[sp++] = (State){k+1, new_carry, s.x | ((uint32_t)xk << k)};
                }
            }
        }
    }
    return n_solutions;
}

/* Same for Maj */
int count_solutions_maj(uint32_t y, uint32_t c1, uint32_t c2) {
    typedef struct { int bit; int carry; uint32_t x; } State;
    State stack[4096];
    int sp = 0;
    int n_solutions = 0;

    stack[sp++] = (State){0, 0, 0};

    while(sp > 0) {
        State s = stack[--sp];
        if(s.bit == 32) {
            if(s.carry == 0) n_solutions++;
            continue;
        }

        int k = s.bit;
        int yk = (y >> k) & 1;
        int c1k = (c1 >> k) & 1;
        int c2k = (c2 >> k) & 1;

        for(int xk = 0; xk <= 1; xk++) {
            /* Maj(x,c1,c2)[k] = xk&c1k ^ xk&c2k ^ c1k&c2k */
            int majk = (xk & c1k) ^ (xk & c2k) ^ (c1k & c2k);
            int sum = xk + majk + s.carry;
            int result = sum & 1;
            int new_carry = sum >> 1;

            if(result == yk) {
                if(sp < 4095) {
                    stack[sp++] = (State){k+1, new_carry, s.x | ((uint32_t)xk << k)};
                }
            }
        }
    }
    return n_solutions;
}

int main() {
    printf("CARRY×NLF BRANCHING FACTOR\n");
    printf("==========================\n\n");

    int N = 100000;
    srand(42);

    /* Distribution of solution count for x + Ch(x,c1,c2) = y */
    int sol_hist_ch[64] = {0}; /* sol_hist[k] = how many cases have k solutions */
    int max_sol_ch = 0;
    double avg_sol_ch = 0;

    int sol_hist_maj[64] = {0};
    int max_sol_maj = 0;
    double avg_sol_maj = 0;

    printf("Computing solution counts (N=%d)...\n", N);
    for(int trial = 0; trial < N; trial++) {
        uint32_t c1 = (uint32_t)rand() | ((uint32_t)rand()<<16);
        uint32_t c2 = (uint32_t)rand() | ((uint32_t)rand()<<16);
        uint32_t x = (uint32_t)rand() | ((uint32_t)rand()<<16);

        /* Ch */
        uint32_t y_ch = x + Ch(x, c1, c2);
        int sol_ch = count_solutions_ch(y_ch, c1, c2);
        if(sol_ch < 64) sol_hist_ch[sol_ch]++;
        if(sol_ch > max_sol_ch) max_sol_ch = sol_ch;
        avg_sol_ch += sol_ch;

        /* Maj */
        uint32_t y_maj = x + Maj(x, c1, c2);
        int sol_maj = count_solutions_maj(y_maj, c1, c2);
        if(sol_maj < 64) sol_hist_maj[sol_maj]++;
        if(sol_maj > max_sol_maj) max_sol_maj = sol_maj;
        avg_sol_maj += sol_maj;

        if(trial % 10000 == 0 && trial > 0)
            printf("  %d/%d...\r", trial, N);
    }

    printf("\nx + Ch(x, c1, c2) = y:\n");
    printf("  Average solutions: %.2f\n", avg_sol_ch / N);
    printf("  Max solutions: %d\n", max_sol_ch);
    printf("  Distribution:\n");
    for(int k = 0; k < 64; k++) {
        if(sol_hist_ch[k] > 0) {
            printf("    %2d solutions: %6d cases (%.1f%%)\n",
                   k, sol_hist_ch[k], 100.0*sol_hist_ch[k]/N);
        }
    }

    printf("\nx + Maj(x, c1, c2) = y:\n");
    printf("  Average solutions: %.2f\n", avg_sol_maj / N);
    printf("  Max solutions: %d\n", max_sol_maj);
    printf("  Distribution:\n");
    for(int k = 0; k < 64; k++) {
        if(sol_hist_maj[k] > 0) {
            printf("    %2d solutions: %6d cases (%.1f%%)\n",
                   k, sol_hist_maj[k], 100.0*sol_hist_maj[k]/N);
        }
    }

    printf("\n══════════════════════════════════════════\n");
    printf("BRANCHING FACTOR SUMMARY\n");
    printf("══════════════════════════════════════════\n\n");

    double bf_ch = avg_sol_ch / N;
    double bf_maj = avg_sol_maj / N;
    printf("  Ch:  avg %.2f solutions per equation (branching factor)\n", bf_ch);
    printf("  Maj: avg %.2f solutions per equation\n", bf_maj);
    printf("\n");

    if(bf_ch < 4 && bf_maj < 4) {
        printf("  ★ SMALL branching factor!\n");
        printf("  Per round: ~%.1f × ~%.1f = ~%.1f branches\n",
               bf_ch, bf_maj, bf_ch * bf_maj);
        printf("  Per 64 rounds: %.1f^64 ≈ 2^%.1f\n",
               bf_ch, 64 * (bf_ch > 1 ? __builtin_log2(bf_ch) : 0));
        printf("  Compare: birthday = 2^128\n");
        double total_log2 = 64 * (bf_ch > 1 ? __builtin_log2(bf_ch) : 0);
        if(total_log2 < 128)
            printf("  ★★★ BELOW BIRTHDAY: branching^64 = 2^%.0f < 2^128!\n", total_log2);
        else
            printf("  Above birthday. No improvement.\n");
    }

    return 0;
}
