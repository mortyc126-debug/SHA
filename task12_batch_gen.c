/*
 * ЗАДАНИЕ 12: Batch pair generator - sequential sets with 4 threads per set
 * Each set: sweep W0 range across 4 threads sharing ONE background.
 * Finds 1 pair per set (guaranteed within 2^32 sweep).
 * Outputs pairs immediately as found.
 *
 * gcc -O3 -march=native -o task12_batch_gen task12_batch_gen.c -lpthread
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#define MASK 0xFFFFFFFFU

static const uint32_t K[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
};
static const uint32_t H0[8] = {
    0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
    0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19
};

static inline uint32_t rotr(uint32_t x, int n) { return (x>>n)|(x<<(32-n)); }
static inline uint32_t Sig0(uint32_t x) { return rotr(x,2)^rotr(x,13)^rotr(x,22); }
static inline uint32_t Sig1(uint32_t x) { return rotr(x,6)^rotr(x,11)^rotr(x,25); }
static inline uint32_t sig0(uint32_t x) { return rotr(x,7)^rotr(x,18)^(x>>3); }
static inline uint32_t sig1(uint32_t x) { return rotr(x,17)^rotr(x,19)^(x>>10); }
static inline uint32_t Ch(uint32_t e,uint32_t f,uint32_t g) { return (e&f)^((~e)&g); }
static inline uint32_t Maj(uint32_t a,uint32_t b,uint32_t c) { return (a&b)^(a&c)^(b&c); }

static uint32_t wang_de17(const uint32_t Wn[16], uint32_t Wf[16]) {
    for(int i=0;i<16;i++) Wf[i]=Wn[i];
    Wf[0] = Wn[0] ^ 0x8000;
    uint32_t sn[8], sf[8];
    for(int i=0;i<8;i++){sn[i]=H0[i];sf[i]=H0[i];}

    uint32_t T1n=sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[0]+Wn[0];
    uint32_t T2n=Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
    uint32_t T1f=sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[0]+Wf[0];
    uint32_t T2f=Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);
    {uint32_t t[8]={T1n+T2n,sn[0],sn[1],sn[2],sn[3]+T1n,sn[4],sn[5],sn[6]};memcpy(sn,t,32);}
    {uint32_t t[8]={T1f+T2f,sf[0],sf[1],sf[2],sf[3]+T1f,sf[4],sf[5],sf[6]};memcpy(sf,t,32);}

    for(int r=1;r<16;r++){
        uint32_t dd=sf[3]-sn[3],dh=sf[7]-sn[7];
        uint32_t dS=Sig1(sf[4])-Sig1(sn[4]),dC=Ch(sf[4],sf[5],sf[6])-Ch(sn[4],sn[5],sn[6]);
        Wf[r]=Wn[r]+(uint32_t)(-(int32_t)(dd+dh+dS+dC));
        T1n=sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[r]+Wn[r];
        T2n=Sig0(sn[0])+Maj(sn[0],sn[1],sn[2]);
        T1f=sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[r]+Wf[r];
        T2f=Sig0(sf[0])+Maj(sf[0],sf[1],sf[2]);
        {uint32_t t[8]={T1n+T2n,sn[0],sn[1],sn[2],sn[3]+T1n,sn[4],sn[5],sn[6]};memcpy(sn,t,32);}
        {uint32_t t[8]={T1f+T2f,sf[0],sf[1],sf[2],sf[3]+T1f,sf[4],sf[5],sf[6]};memcpy(sf,t,32);}
    }
    uint32_t Wn16=sig1(Wn[14])+Wn[9]+sig0(Wn[1])+Wn[0];
    uint32_t Wf16=sig1(Wf[14])+Wf[9]+sig0(Wf[1])+Wf[0];
    T1n=sn[7]+Sig1(sn[4])+Ch(sn[4],sn[5],sn[6])+K[16]+Wn16;
    T1f=sf[7]+Sig1(sf[4])+Ch(sf[4],sf[5],sf[6])+K[16]+Wf16;
    return (sf[3]+T1f)-(sn[3]+T1n);
}

/* Per-set parallel search */
typedef struct {
    const uint32_t *bg;  /* background Wn[1..15] */
    uint64_t start, end;
    volatile int *hit;
    uint32_t hit_Wn[16], hit_Wf[16];
    uint64_t evals;
} sarg_t;

static void *set_worker(void *arg) {
    sarg_t *sa = (sarg_t *)arg;
    sa->evals = 0;
    uint32_t Wn[16], Wf[16];
    for(int i=1;i<16;i++) Wn[i]=sa->bg[i];

    for(uint64_t w0=sa->start; w0<sa->end; w0++) {
        if(*sa->hit) return NULL;
        Wn[0]=(uint32_t)w0;
        uint32_t de=wang_de17(Wn,Wf);
        sa->evals++;
        if(de==0){
            *sa->hit=1;
            memcpy(sa->hit_Wn,Wn,64);
            memcpy(sa->hit_Wf,Wf,64);
            return NULL;
        }
    }
    return NULL;
}

int main(int argc, char **argv) {
    int target = 100;
    int num_threads = 4;
    int time_limit = 600;
    if(argc>1) target=atoi(argv[1]);
    if(argc>2) num_threads=atoi(argv[2]);
    if(argc>3) time_limit=atoi(argv[3]);

    setvbuf(stdout,NULL,_IONBF,0);
    setvbuf(stderr,NULL,_IONBF,0);

    fprintf(stderr,"Task12 batch gen: target=%d, threads=%d, limit=%ds\n",
            target,num_threads,time_limit);

    FILE *rng=fopen("/dev/urandom","rb");
    time_t t_global=time(NULL);
    int total_found=0;
    uint64_t total_evals=0;

    while(total_found < target) {
        if(time(NULL)-t_global >= time_limit) break;

        uint32_t bg[16];
        fread(bg,4,16,rng);

        volatile int hit=0;
        pthread_t threads[64];
        sarg_t sargs[64];
        uint64_t chunk=0x100000000ULL/num_threads;

        for(int t=0;t<num_threads;t++){
            sargs[t].bg=bg;
            sargs[t].start=(uint64_t)t*chunk;
            sargs[t].end=(t==num_threads-1)?0x100000000ULL:(uint64_t)(t+1)*chunk;
            sargs[t].hit=&hit;
            sargs[t].evals=0;
            pthread_create(&threads[t],NULL,set_worker,&sargs[t]);
        }
        for(int t=0;t<num_threads;t++)
            pthread_join(threads[t],NULL);

        uint64_t set_evals=0;
        for(int t=0;t<num_threads;t++) set_evals+=sargs[t].evals;
        total_evals+=set_evals;

        if(hit){
            /* Find which thread got it */
            for(int t=0;t<num_threads;t++){
                /* Check if this thread found it by testing de17 */
                uint32_t Wf_test[16];
                uint32_t de=wang_de17(sargs[t].hit_Wn,Wf_test);
                if(de==0){
                    total_found++;
                    printf("PAIR\n");
                    printf("Wn =");
                    for(int w=0;w<16;w++) printf(" %08x",sargs[t].hit_Wn[w]);
                    printf("\n");
                    printf("Wf =");
                    for(int w=0;w<16;w++) printf(" %08x",sargs[t].hit_Wf[w]);
                    printf("\n");
                    double elapsed=time(NULL)-t_global;
                    fprintf(stderr,"\r  [%d/%d] %.1fs elapsed, %.1fM evals, %.1f evals/pair   \n",
                            total_found,target,elapsed,total_evals/1e6,
                            total_found>0?(double)total_evals/total_found:0);
                    break;
                }
            }
        }
    }
    fclose(rng);

    double elapsed=time(NULL)-t_global;
    fprintf(stderr,"\nDone: %d pairs in %.0fs (%.2e evals, %.1fM/s)\n",
            total_found,elapsed,(double)total_evals,total_evals/elapsed/1e6);
    return 0;
}
