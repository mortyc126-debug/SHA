/*
 * SIGNAL AMPLIFICATION: maximize Z-score of antenna.
 *
 * Method 1: OPTIMIZE antenna weights to maximize Z
 * Method 2: DIFFERENTIAL mode — compare pairs
 * Method 3: MULTI-SAMPLE averaging
 * Method 4: PER-ROUND decomposition
 *
 * gcc -O3 -march=native -o sig_amp signal_amplify.c -lm
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define ROTR(x,n) (((x)>>(n))|((x)<<(32-(n))))
#define CH(e,f,g) (((e)&(f))^((~(e))&(g)))
#define MAJ(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))
#define S0(x) (ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define S1(x) (ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define s0(x) (ROTR(x,7)^ROTR(x,18)^((x)>>3))
#define s1(x) (ROTR(x,17)^ROTR(x,19)^((x)>>10))
static const uint32_t K[64]={
0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
static const uint32_t IV[8]={
0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};

/* Per-round feature extraction: 4 features per round = 256 features total */
void extract_features(const uint32_t msg[16], int features[256]) {
    uint32_t W[64];
    for(int i=0;i<16;i++) W[i]=msg[i];
    for(int i=16;i<64;i++) W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
    uint32_t a=IV[0],b=IV[1],c=IV[2],d=IV[3],e=IV[4],f=IV[5],g=IV[6],h=IV[7];

    for(int r=0;r<64;r++){
        uint32_t erase = f ^ g;
        uint32_t s1e=S1(e), che=CH(e,f,g), s0a=S0(a), mja=MAJ(a,b,c);

        /* Feature 0: carry at erased positions (first addition h+S1(e)) */
        uint32_t v1 = h + s1e;
        features[r*4+0] = __builtin_popcount((v1 ^ (h ^ s1e)) & erase);

        /* Feature 1: carry chain length at erased positions */
        uint32_t masked = (v1 ^ (h ^ s1e)) & erase;
        int chain=0, maxc=0;
        for(int k=0;k<32;k++){if((masked>>k)&1){chain++;if(chain>maxc)maxc=chain;}else chain=0;}
        features[r*4+1] = maxc;

        /* Feature 2: generate bits (h&s1e) at erased positions */
        features[r*4+2] = __builtin_popcount((h & s1e) & erase);

        /* Feature 3: propagate chain (h^s1e) at erased */
        features[r*4+3] = __builtin_popcount((h ^ s1e) & erase);

        uint32_t t1=h+s1e+che+K[r]+W[r], t2=s0a+mja;
        uint32_t anew=t1+t2, enew=d+t1;
        h=g;g=f;f=e;e=enew;d=c;c=b;b=a;a=anew;
    }
}

int main() {
    printf("SIGNAL AMPLIFICATION\n");
    printf("====================\n\n");

    srand(42);
    uint32_t msg_correct[16];
    for(int i=0;i<16;i++) msg_correct[i]=(uint32_t)rand()|((uint32_t)rand()<<16);

    /* Extract features for correct msg */
    int feat_correct[256];
    extract_features(msg_correct, feat_correct);

    /* Collect features for N random msgs */
    int N = 10000;
    double feat_sum[256]={0}, feat_sum2[256]={0};

    for(int t=0;t<N;t++){
        uint32_t msg_w[16];
        for(int i=0;i<16;i++) msg_w[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
        int feat[256];
        extract_features(msg_w, feat);
        for(int i=0;i<256;i++){feat_sum[i]+=feat[i];feat_sum2[i]+=(double)feat[i]*feat[i];}
    }

    /* Compute Z-score per feature */
    double z_scores[256];
    for(int i=0;i<256;i++){
        double mean=feat_sum[i]/N;
        double std=sqrt(feat_sum2[i]/N-mean*mean);
        z_scores[i] = std>0 ? (feat_correct[i]-mean)/std : 0;
    }

    /* METHOD 1: Find best individual features */
    printf("METHOD 1: Best individual features (per round)\n");
    printf("───────────────────────────────────────────────\n");

    int best_idx = 0;
    double best_z = 0;
    for(int i=0;i<256;i++){
        if(fabs(z_scores[i]) > fabs(best_z)){best_z=z_scores[i];best_idx=i;}
    }
    printf("  Best feature: round %d, type %d, Z=%+.3f\n",
           best_idx/4, best_idx%4, best_z);

    /* Top 10 features */
    printf("  Top 10 features:\n");
    int sorted_idx[256];
    for(int i=0;i<256;i++) sorted_idx[i]=i;
    for(int i=0;i<256;i++) for(int j=i+1;j<256;j++)
        if(fabs(z_scores[sorted_idx[i]]) < fabs(z_scores[sorted_idx[j]]))
            {int t=sorted_idx[i];sorted_idx[i]=sorted_idx[j];sorted_idx[j]=t;}

    const char *feat_names[4]={"ch×carry","chain_len","generate","propagate"};
    for(int i=0;i<10;i++){
        int idx=sorted_idx[i];
        printf("    round %2d %-10s Z=%+.3f\n", idx/4, feat_names[idx%4], z_scores[idx]);
    }

    /* METHOD 2: OPTIMIZED WEIGHTED SUM */
    printf("\nMETHOD 2: Weighted sum (weights = Z-scores)\n");
    printf("────────────────────────────────────────────\n");
    {
        /* Use Z-scores as weights: features with high Z get high weight */
        double correct_score = 0;
        for(int i=0;i<256;i++) correct_score += z_scores[i] * feat_correct[i];

        double s=0, s2=0;
        for(int t=0;t<N;t++){
            uint32_t msg_w[16];
            for(int i=0;i<16;i++) msg_w[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
            int feat[256];
            extract_features(msg_w, feat);
            double score = 0;
            for(int i=0;i<256;i++) score += z_scores[i] * feat[i];
            s+=score; s2+=score*score;
        }
        double mean=s/N, std=sqrt(s2/N-mean*mean);
        double z_combined = std>0 ? (correct_score-mean)/std : 0;

        printf("  Correct msg score: %.1f\n", correct_score);
        printf("  Random mean: %.1f, std: %.1f\n", mean, std);
        printf("  COMBINED Z = %+.2f %s\n", z_combined,
               fabs(z_combined)>3?"★★★ STRONG":fabs(z_combined)>2?"★★ SIGNAL":"weak");
    }

    /* METHOD 3: DIFFERENTIAL — compare two msgs, predict which is "more correct" */
    printf("\nMETHOD 3: Differential (K correct vs K-1 correct)\n");
    printf("──────────────────────────────────────────────────\n");
    {
        /* For K=1..16: compute antenna for K correct words and K-1 correct.
         * Can we tell which has MORE correct words? */
        int correct_preds = 0, total = 0;
        for(int K=1;K<=16;K++){
            for(int trial=0;trial<1000;trial++){
                uint32_t msg_K[16], msg_K1[16];
                /* K correct words */
                for(int i=0;i<K;i++) msg_K[i]=msg_correct[i];
                for(int i=K;i<16;i++) msg_K[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
                /* K-1 correct words */
                for(int i=0;i<K-1;i++) msg_K1[i]=msg_correct[i];
                for(int i=K-1;i<16;i++) msg_K1[i]=(uint32_t)rand()|((uint32_t)rand()<<16);

                int feat_K[256], feat_K1[256];
                extract_features(msg_K, feat_K);
                extract_features(msg_K1, feat_K1);

                double score_K=0, score_K1=0;
                for(int i=0;i<256;i++){
                    score_K += z_scores[i]*feat_K[i];
                    score_K1 += z_scores[i]*feat_K1[i];
                }

                if(score_K > score_K1) correct_preds++;
                total++;
            }
        }
        double accuracy = (double)correct_preds/total;
        printf("  Can we predict which has MORE correct words?\n");
        printf("  Accuracy: %d/%d = %.3f (0.500 = random, >0.510 = signal)\n",
               correct_preds, total, accuracy);
        printf("  → %s\n", accuracy>0.520?"★★★ DIFFERENTIAL SIGNAL!":
               accuracy>0.510?"★★ weak signal":
               accuracy>0.505?"★ marginal":"no signal");
    }

    /* METHOD 4: MULTI-MSG verification (same msg, different measures) */
    printf("\nMETHOD 4: Z across multiple correct messages\n");
    printf("─────────────────────────────────────────────\n");
    {
        /* Is Z=1.74 specific to our one correct msg, or universal? */
        double z_values[20];
        for(int m=0;m<20;m++){
            uint32_t msg_m[16];
            srand(m*12345+67);
            for(int i=0;i<16;i++) msg_m[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
            int feat_m[256];
            extract_features(msg_m, feat_m);
            double score_m=0;
            for(int i=0;i<256;i++) score_m += z_scores[i]*feat_m[i];

            /* Compare with random baseline (already computed) */
            double s=0,s2=0;
            srand(m*99999+42);
            for(int t=0;t<2000;t++){
                uint32_t msg_w[16];
                for(int i=0;i<16;i++) msg_w[i]=(uint32_t)rand()|((uint32_t)rand()<<16);
                int feat[256];
                extract_features(msg_w, feat);
                double sc=0;
                for(int i=0;i<256;i++) sc+=z_scores[i]*feat[i];
                s+=sc;s2+=sc*sc;
            }
            double mean=s/2000, std=sqrt(s2/2000-mean*mean);
            z_values[m] = std>0?(score_m-mean)/std:0;
        }

        double z_mean=0,z_var=0;
        for(int i=0;i<20;i++) z_mean+=z_values[i]; z_mean/=20;
        for(int i=0;i<20;i++) z_var+=(z_values[i]-z_mean)*(z_values[i]-z_mean); z_var/=20;

        printf("  Z-scores for 20 different correct messages:\n  ");
        for(int i=0;i<20;i++) printf("%.1f ", z_values[i]);
        printf("\n  Mean Z = %.2f, Std = %.2f\n", z_mean, sqrt(z_var));
        printf("  → %s\n", z_mean>1.0?"★ UNIVERSAL positive signal":
               z_mean>0.5?"◆ Weak but consistent":"no universal signal");
    }

    return 0;
}
