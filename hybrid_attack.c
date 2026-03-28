#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

/*
 * HYBRID ATTACK: DFS partial preimage + Birthday collision
 *
 * Strategy:
 *   Phase 1: Use DFS to find N inputs with specific k-bit partial hash.
 *            Cost: N × DFS_cost(k) evaluations.
 *   Phase 2: Among N inputs with same k-bit prefix, find collision
 *            in remaining (32-k) bits via birthday.
 *            Need N ≈ 2^((32-k)/2) for birthday collision.
 *
 * Optimal k minimizes total cost = N × DFS(k) + N.
 *
 * Also: MULTI-WORD test with W[0] and W[1] both active.
 */

#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))
#define SHR(x,n)  ((x)>>(n))
#define CH(e,f,g) (((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))
#define SIG0(x) (ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define SIG1(x) (ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define sig0(x) (ROTR(x,7)^ROTR(x,18)^SHR(x,3))
#define sig1(x) (ROTR(x,17)^ROTR(x,19)^SHR(x,10))

static const uint32_t H0[8] = {
    0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
    0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19
};
static const uint32_t K_const[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

static inline uint32_t sha256_word(uint32_t w0, uint32_t w1, int tw) {
    uint32_t W[64];
    W[0]=w0; W[1]=w1;
    for(int i=2;i<16;i++) W[i]=0;
    for(int i=16;i<64;i++)
        W[i]=sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16];
    uint32_t a=H0[0],b=H0[1],c=H0[2],d=H0[3],
             e=H0[4],f=H0[5],g=H0[6],h=H0[7];
    for(int r=0;r<64;r++){
        uint32_t T1=h+SIG1(e)+CH(e,f,g)+K_const[r]+W[r];
        uint32_t T2=SIG0(a)+MAJ(a,b,c);
        h=g;g=f;f=e;e=d+T1;d=c;c=b;b=a;a=T1+T2;
    }
    uint32_t out[8]={a+H0[0],b+H0[1],c+H0[2],d+H0[3],
                     e+H0[4],f+H0[5],g+H0[6],h+H0[7]};
    return out[tw];
}

int main() {
    printf("HYBRID ATTACK: DFS + Birthday on SHA-256\n\n");

    /* ========== TEST 1: Optimal hybrid cost ========== */
    printf("=== TEST 1: OPTIMAL HYBRID COST (theoretical) ===\n\n");

    /* From honest_dfs data:
     * k=16: DFS ≈ 730 nodes
     * k=20: DFS ≈ 7384
     * k=24: DFS ≈ 58924
     *
     * For collision in H[4] (32 bits):
     * Fix k bits via DFS, birthday on (32-k) bits.
     * Need N = 2^((32-k)/2) DFS solutions for birthday.
     * Total cost = N × DFS(k)
     */

    printf("  k  | DFS(k) nodes | N=2^((32-k)/2) | Total cost  | vs birthday 2^16\n");
    printf("  ---+--------------+----------------+-------------+-----------------\n");

    struct { int k; double dfs; } data[] = {
        {0, 1}, {4, 32}, {8, 32}, {12, 56}, {16, 730},
        {20, 7384}, {24, 58924}, {28, 1000000}, {32, 853414}
    };
    int nd = 9;

    double best_cost = 1e18;
    int best_k = 0;

    for(int i=0;i<nd;i++){
        int k = data[i].k;
        double dfs_cost = data[i].dfs;
        double N_birthday = pow(2.0, (32-k)/2.0);
        double total = N_birthday * dfs_cost;
        double vs_bday = total / 65536.0;
        char marker = (total < best_cost) ? '*' : ' ';
        if (total < best_cost) { best_cost = total; best_k = k; }
        printf("  %2d | %12.0f | %14.0f | %11.0f | %15.1fx %c\n",
               k, dfs_cost, N_birthday, total, vs_bday, marker);
    }

    printf("\n  Best k=%d, total cost=%.0f (%.1f× birthday)\n\n",
           best_k, best_cost, best_cost/65536.0);

    /* ========== TEST 2: ACTUAL 2-word collision search ========== */
    printf("=== TEST 2: ACTUAL COLLISION SEARCH on H[4] ===\n");
    printf("  Brute birthday: try random W[0] until H[4] collision\n\n");

    /* Standard birthday on H[4] */
    {
        uint32_t *table_keys = calloc(100000, sizeof(uint32_t));
        uint32_t *table_vals = calloc(100000, sizeof(uint32_t));
        int table_size = 100000;
        int n_entries = 0;

        clock_t t0 = clock();
        uint64_t evals = 0;
        int found = 0;

        for(uint32_t w0 = 0; w0 < 200000 && !found; w0++) {
            evals++;
            uint32_t h4 = sha256_word(w0, 0, 4);
            uint32_t slot = h4 % table_size;

            /* Linear probing */
            for(int probe = 0; probe < 100; probe++) {
                int idx = (slot + probe) % table_size;
                if(table_vals[idx] == 0 && n_entries > 0) {
                    /* Check if this is empty or hash=0 */
                    if(table_keys[idx] == 0 && idx != 0) {
                        table_keys[idx] = w0;
                        table_vals[idx] = h4 | 1; /* mark occupied */
                        n_entries++;
                        break;
                    }
                }
                if((table_vals[idx] & ~1) == (h4 & ~1) && table_keys[idx] != w0) {
                    /* Potential collision — verify */
                    uint32_t h4_check = sha256_word(table_keys[idx], 0, 4);
                    if(h4_check == h4) {
                        double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
                        printf("  BIRTHDAY COLLISION H[4]:\n");
                        printf("    W0_a = 0x%08x, W0_b = 0x%08x\n", table_keys[idx], w0);
                        printf("    H[4] = 0x%08x\n", h4);
                        printf("    Evals: %lu (%.1fs)\n", evals, el);
                        printf("    Expected: 2^16 = 65536\n");
                        found = 1;
                        break;
                    }
                }
                if(table_keys[idx] == 0) {
                    table_keys[idx] = w0;
                    table_vals[idx] = h4;
                    n_entries++;
                    break;
                }
            }
        }
        if(!found) printf("  No collision in 200K tries\n");
        free(table_keys);
        free(table_vals);
    }

    /* ========== TEST 3: Multi-word DFS ========== */
    printf("\n=== TEST 3: MULTI-WORD DFS (W[0]+W[1]) ===\n");
    printf("  Does adding W[1] kill DFS efficiency?\n\n");

    /* Test: vary low bits of W[0] AND W[1] */
    for(int n_per_word = 8; n_per_word <= 16; n_per_word += 4) {
        int n_total = 2 * n_per_word;
        int k_target = 1; /* single bit target */
        uint32_t mask = 1;
        int trials = 10;
        uint64_t total_nodes = 0;

        clock_t t0 = clock();
        for(int trial = 0; trial < trials; trial++) {
            uint32_t w0_base = rand();
            uint32_t w1_base = rand();
            uint32_t target = sha256_word(w0_base & ((1U<<n_per_word)-1),
                                          w1_base & ((1U<<n_per_word)-1), 4) & mask;

            /* Brute DFS with propagation on last 8 bits */
            uint64_t nodes = 0;
            int found = 0;
            uint32_t total_space = 1U << n_total;
            uint32_t w0_mask = (1U << n_per_word) - 1;

            for(uint32_t x = 0; x < total_space && !found; x++) {
                nodes++;
                uint32_t w0 = x & w0_mask;
                uint32_t w1 = (x >> n_per_word) & w0_mask;
                if((sha256_word(w0, w1, 4) & mask) == target) {
                    found = 1;
                }
            }
            total_nodes += found ? nodes : total_space;
        }

        double mean = (double)total_nodes / trials;
        double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
        double eps = 1.0 - log2(mean) / n_total;
        printf("  W[0]+W[1] %d bits each (n=%d): mean=%.0f nodes, eps=%.3f (%.1fs)\n",
               n_per_word, n_total, mean, eps, el);
    }

    /* Compare: 1-word same total bits */
    printf("\n  Comparison with 1-word:\n");
    for(int n = 16; n <= 32; n += 8) {
        int trials = 10;
        uint64_t total = 0;
        for(int t = 0; t < trials; t++) {
            uint32_t target_input = rand() & ((1U<<n)-1);
            uint32_t target = sha256_word(target_input, 0, 4) & 1;
            for(uint32_t x = 0; x < (1U<<n); x++) {
                total++;
                if((sha256_word(x, 0, 4) & 1) == target) break;
            }
        }
        double mean = (double)total / trials;
        double eps = 1.0 - log2(mean) / n;
        printf("  W[0] only %d bits: mean=%.0f nodes, eps=%.3f\n", n, mean, eps);
    }

    return 0;
}
