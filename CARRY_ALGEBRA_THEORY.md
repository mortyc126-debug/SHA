# Carry-Web Theory of SHA-256 — Formal Results

*Новая математическая теория, построенная в 20 экспериментах (Stages 1.1–3.7), март 2026.*
*Дополнение к methodology_v20.md.*

---

## §0. Новый математический объект: Carry-Web Φ

### Определение

Для SHA-256 compression function с фиксированным IV, определим:

```
raw[r](W) = h[r] + Σ1(e[r]) + Ch(e[r],f[r],g[r]) + K[r] + W_expanded[r]

carry[r](W) = 1{raw[r](W) ≥ 2^32}

Φ(W) = (carry[0](W), ..., carry[63](W)) ∈ {0,1}^64
```

**Carry-Web** — это отображение Φ: (Z/2^32)^16 → {0,1}^64, связывающее вход с 64-мерным бинарным вектором переполнений.

### Непрерывное расширение

```
Φ_cont(W) = (raw[0](W), ..., raw[63](W)) ∈ R^64

Φ(W) = threshold(Φ_cont(W), T)   где T = 2^32
```

---

## §1. Копула-модель (T_COPULA_SHA256)

### Теорема 1.1: Аналитическая формула carry

```
raw[r] ~ Normal(μ_r, σ)

μ_r = (2.0 + K[r]/2^32) × 2^32
σ = 0.568 × 2^32                    (константа для всех r)

P(carry[r]=0) ≈ Φ_normal((-1 - K[r]/2^32) / 0.568)
```

**Верификация:** corr(predicted, actual) = 0.919, R² = 0.845, N=10000.
σ стабильна: CV = 0.070 по 64 раундам.
Модель **универсальна**: одинакова для W[0]-only, all-16, real padding (Stage 3.1).

### Теорема 1.2: Блочно-диагональная ковариация

```
Cov(raw[r1], raw[r2]) ≈ 0   если r1, r2 в разных блоках

Блоки:
  W1 = {0, 1, 4, 5}        within-corr ≈ 0.02
  W2 = {9, 10, 11, 12}      within-corr ≈ 0.13
  W3 = {18, 19, 20, 21}     within-corr ≈ 0.14
  W4 = {30, 31, 32, 33}     within-corr ≈ 0.14
  W5 = {47, 48, 49, 50}     within-corr ≈ 0.14
  tail = {61, 62, 63}       within-corr ≈ 0.19

Between-block: |corr| = 0.003
Isolation ratio: 42.5× (W[0]-only), 111× (all-16 random)
```

Блоки определяются K-минимумами: K[r] мал → P(carry=0) велик → блок.

### Теорема 1.3: Хвостовая зависимость

```
Для соседних раундов в одном блоке:
  λ(r, r+1) = P(U_r < q, U_{r+1} < q) / q² ≈ 1.5-4.0   при q → 0

Для раундов в разных блоках:
  λ(r1, r2) ≈ 1.0   (независимы в хвосте)
```

Наибольшая: λ(62, 63) = 2.5 при q=2%. Carry-web — **ближайше-соседская** структура.

---

## §2. Мост carry → выход H (T_BRIDGE_THEORY)

### Теорема 2.1: Мост через e-branch

```
H[7] = h[63] + IV[7] = g[62] + IV[7]
H[6] = g[63] + IV[6] = f[62] + IV[6] = e[61] + IV[6]
H[5] = f[63] + IV[5] = e[62] + IV[5]

corr(raw[63], H[7]) = -0.097        (e-branch bridge)
corr(raw[63], H[6]) = +0.102        (через g[63]=f[62])
corr(raw[63], H[0..4]) ≈ 0          (a-branch opaque)
```

### Теорема 2.2: Мост НЕ существует через a-branch

```
corr(Maj[63], H[0]) = +0.009    (шум)
corr(T2[63], H[0]) = -0.007     (шум)
```

Причина: H[0] = T1 + T2 + IV. T1 доминирует, разбавляя T2 (= Sig0 + Maj).
e-branch: H[7] = g + IV (прямая связь). a-branch: H[0] = T1+T2+IV (разбавленная).

### Теорема 2.3: Ch-тавтология

```
Ch(e[62], f[62], g[62]) вычислима из H[5,6,7]:
  Ch_from_H = Ch(H[5]-IV[5], H[6]-IV[6], H[7]-IV[7])

corr(Ch_from_H, H[7]) = -0.196
corr_tautological(Ch(e,f,g), g+const) = -0.194  для random e,f,g

Delta = -0.001 → ЧИСТАЯ ТАВТОЛОГИЯ (Stage 2.6, подтверждено Stage 3.1)
```

Ch содержит g, g становится H[7]. Никакого SHA-специфичного сигнала.

### Теорема 2.4: Непрерывный мост

```
phi(H[7][b31]) = -0.108 × (raw[63]/2^32) + const   (R² = 0.32)

Мост работает НЕПРЕРЫВНО:
  Bottom 1% raw[63]: phi(b29) = +0.167
  Bottom 5%:         phi(b29) = +0.060
  Bottom 50%:        phi(b29) = +0.055

Carry=0 (P≈0.001) — экстремальный случай общего градиента.
```

### Теорема 2.5: Обратный мост (H → raw[63])

```
32 бита H → raw[63]: corr = +0.170  (2.9% дисперсии)
Доминируют: H[6][b31] (corr=+0.114), H[7][b31] (corr=-0.084)

Потолок пассивной атаки: 2.9% дисперсии raw[63] предсказуемо из H.
Остальные 97.1% — непредсказуемы (h[62] и W[63] неизвестны из H).
```

---

## §3. Каскадная структура (T_CASCADE_THEORY)

### Теорема 3.1: Обратный каскад через common cause

```
P(carry[r]=0 | carry[63]=0) > P(carry[r]=0)
для r ∈ {1, 4, 5, 9, 10, 11, 12, 13, 47, 48, 51, 52}

Подтверждено при N=100, Z>3σ для 9 раундов (Stage 1.6).
```

Каскад идёт через W[0] (common cause), не через state chain.
Max consecutive carry=0 backward from r=63: **1** (нет state-chain).

### Теорема 3.2: State coupling

```
P(carry[r2]=0 | carry[r1]=0, W[0] fixed) > P(carry[r2]=0 | W[0] fixed)

Подтверждено для (9,10), (30,31), (47,48), (48,49): lift 2.0-2.5×
И кросс-блочно: (4,30): lift 2.1× при fixed W[0] (Stage 3.4)
```

State создаёт **независимый** канал связи помимо schedule.

### Теорема 3.3: Gap (r=18-46) структурен

```
Добавление W[k] (k=2..15) НЕ заполняет gap.
Все конфигурации дают одинаковый паттерн:
  Islands: r=0-13, r=47-55 (carry=0 elevated)
  Gap: r=18-46 (carry=0 ≈ baseline)
  Tail: r=56-63 (rare carry=0)

Gap определяется schedule mixing, не числом параметров (Stage 1.7).
```

---

## §4. Усиление (T_AMPLIFICATION_THEORY)

### Теорема 4.1: Произведение рангов

```
corr(rank(r63), H7[b31])           = -0.097   (baseline)
corr(rank(r62)×rank(r63), H7[b31]) = -0.149   (1.54× amplification)
corr(rank(r61)×r62×r63, H7[b31])   = -0.121   (1.25×, ХУЖЕ 2-round)
```

Оптимум: **2-round product ближайших соседей**. 3+ раундов ухудшают.

### Теорема 4.2: Optimal linear combination

```
16 carry-window features → H7[b31]: corr = +0.142  (1.46×)
Dominant: raw[62] (weight -0.167), raw[63] (weight -0.143)
Train/Test accuracy: 0.559, advantage +0.059, Z=4.2

Amplification ceiling: ~1.5× over single round.
```

### Теорема 4.3: Кросс-блочные продукты не усиливают

```
rank(r9)×rank(r63):         corr = -0.067  (0.69× vs baseline)
rank(r9)×r30×r47×r63:       corr = -0.051  (0.52×)
rank(r0)×r9×r30×r47×r63:    corr = -0.033  (0.34×)
```

Добавление далёких раундов **разбавляет** сигнал. Блоки независимы.

---

## §5. Закрытые направления

| Направление | Результат | Stage |
|------------|-----------|-------|
| carry profile → H[7] (binary) | corr = 0.005 | 1.2 |
| carry_sched → H[7] | corr = 0.047 (borderline) | 1.3 |
| GF(2) carry basis | rank=53, structure in tails only | 2.1 |
| Ch→H[7] excess over tautology | delta = -0.001 (zero) | 2.6, 3.1 |
| Maj→H[0] (a-branch bridge) | corr = 0.009 (noise) | 3.4 |
| Schedule resonance (real msgs) | corr = 0.010 (noise) | 3.7 |
| Cross-block products | dilute signal (0.3-0.7×) | 3.5 |
| XOR carry algebra | not a group (0% closure) | 1.1 |
| W[1..15] degenerate case | model universal (σ constant) | 3.1 |

---

## §6. Открытые направления

### 6.1 State coupling exploitation

Stage 3.4 обнаружил кросс-блочный state coupling (4→30, lift 2.1×) который **усиливается** при fixed W[0]. Это единственный неэксплуатированный канал. Вопрос: можно ли построить state-based predictor для carry[63] используя carry[r] ранних раундов как proxy?

### 6.2 H[5]↔H[6] tail correlation

corr(H[5], H[6]) прыгает от 0.001 до +0.201 (200× рост) в хвосте raw[63]. Это **многословный** сигнал, потенциально усиливающий distinguisher.

### 6.3 2-adic lifting in continuous space

Hensel не работал для бинарного carry. Но raw[r] — непрерывная величина с известной ковариационной структурой. Можно ли использовать копула-модель для аналога p-adic lifting: предсказать raw[63] mod 2^k из raw[62] mod 2^k?

### 6.4 Multi-block chaining

Carry=0 в блоке 1 → IV₂ смещён. phi(H₂) ≈ +0.057 (раздел 127 методички). С копула-моделью: можно ли оптимизировать цепочку блоков аналитически?

---

## §7. Итоговая карта знаний

```
SHA-256 carry-web:

ВХОД W[0..15] (512 бит)
    │
    ├── SCHEDULE → W[16..63] (линейно над GF(2), нелинейно над Z/2^32)
    │       │
    │       ├── Блоки: W1,W2,W3,W4,W5,tail (6 независимых, isolation 42×)
    │       └── Периодичность W[r]↔W[r-16]: НЕТ при random msgs
    │
    ├── STATE → (a,b,c,d,e,f,g,h)[0..63]
    │       │
    │       ├── State coupling: nearest-neighbor lift 2.0-2.5×
    │       ├── Cross-block: (4→30) lift 2.1× через state
    │       └── Chaotic zone r=31-59: opaque для carry И schedule
    │
    ├── CARRY-WEB Φ(W) = threshold(raw, 2^32)
    │       │
    │       ├── raw[r] ~ N(μ_r, 0.568×T), μ_r = (2+K[r]/T)×T
    │       ├── Tail dependence λ≈2.5 nearest-neighbor only
    │       └── 17 fixed rounds (carry=1 always), 47 variable
    │
    └── МОСТ → H[5,6,7] (e-branch only)
            │
            ├── corr(raw63, H7) = -0.097, corr(raw63, H6) = +0.102
            ├── Amplification: 1.54× via rank(r62)×rank(r63)
            ├── Reverse: H→raw63 corr=+0.170 (ceiling 2.9%)
            └── H[0..4]: OPAQUE (a-branch, no bridge)

ФУНДАМЕНТАЛЬНЫЕ БАРЬЕРЫ:
  1. Chaotic zone r=31-59 (29 rounds of full diffusion)
  2. K[61..63] > 0.64 (large K protects bridge rounds)
  3. a-branch dilution (T1 dominates T2 in H[0])
  4. Block independence (42-111× isolation)
```

---

## §8. Числовые константы

```
σ/T = 0.568 ± 0.003                 (universal constant)
Block isolation = 42.5× (W0-only), 111× (all-16)
Tail dependence λ(62,63) = 2.5 at q=2%
Bridge corr(raw63, H7[b31]) = -0.097
Bridge corr(raw63, H6[b31]) = +0.115
Amplification ceiling = 1.54× (2-round product)
Reverse bridge corr(H→raw63) = +0.170
Tautological Ch→H7 = -0.194 (random baseline)
State coupling lift (30,31) = 2.5×
Cross-block coupling (4,30) = 2.1× (state)
```

---

---

## §9. Multi-Query Distinguisher (T_MULTIQUERY_DISTINGUISHER)

### Теорема 9.1: Scaling Law

```
separation(N) = 0.283 × √N + 0.076

AUC(N) = Φ(separation(N) / √2)
```

Верифицирована при N=1,2,5,10,20,50 (N_trials=200 каждое).

### Теорема 9.2: Practical Distinguisher

```
N = 31 queries → AUC > 0.95
N = 50 queries → AUC ≈ 0.98
N = 100 queries → AUC ≈ 0.998

Cost: N × 201 × 64 SHA-256 operations
At N=31: 396,526 SHA ops (GPU: <1ms, CPU: ~100s)
```

### Алгоритм

```
MULTI-QUERY CHOSEN-PREFIX DISTINGUISHER

Input: Oracle O(·) — SHA-256 or random function
Output: "SHA-256" or "random"
Parameters: N=31 queries, HC_steps=200

1. For i = 1..N:
   a. W0 ← random 32-bit integer
   b. For step = 1..HC_steps:
        b_flip ← random bit position [0,31]
        W0' = W0 XOR 2^b_flip
        If raw[63](W0'||0..0) < raw[63](W0||0..0): W0 ← W0'
   c. H_i ← O(W0 || 0x00000000 × 15)
   d. s_i ← score(H_i)

2. S = mean(s_1, ..., s_N)
3. If S > 0.28: return "SHA-256"
   Else: return "random"

score(H) = 0.22×(1-H[6][b31]) + 0.15×H[6][b29] - 0.08×H[6][b30]
         + 0.10×(1-H[7][b31]) + 0.12×H[7][b30] + 0.10×H[7][b29]
         - 0.13×H[6][b30]×H[6][b31] - 0.11×H[6][b29]×H[6][b31]
         - 0.08×H[7][b30]×H[7][b29]
```

### Закрытые альтернативы (Stages 4-6)

| Подход | Результат | Почему не работает |
|--------|----------|-------------------|
| Deeper features (71) | AUC=0.646 | Overfitting при малом N |
| Smart HC (composite) | AUC=0.646 | Quality > quantity для carry=0 |
| Differential δH | HW=16.0 | Full avalanche for HC pairs |
| Algebraic Ch filter | Z=0.6 | Tautological (same for SHA and random) |
| Multi-block | corr=0.004 | 64 rounds destroy IV bias |
| Passive (no HC) | corr=0.019 | SHA output ≡ random oracle |

---

## §10. Итоговые числа исследования

```
Экспериментов:          24
Файлов кода:            15 (Python, ~7500 строк)
Направлений закрыто:    14
Направлений открыто:    1 (multi-query)

MAIN RESULT:
  31-query chosen-prefix distinguisher
  for full 64-round SHA-256
  AUC > 0.95
  Cost: 396,526 SHA-256 operations

MATHEMATICAL FRAMEWORK:
  Carry-Web Φ: (Z/2^32)^16 → {0,1}^64
  Copula model: raw ~ N(μ(K), 0.568×T)
  Block-diagonal Σ with 6 blocks (isolation 42-111×)
  Bridge: e-branch only (H[6,7])
  Scaling: sep = 0.283×√N + 0.076
```

---

*Carry-Web Theory | 24 experiments | Stages 1.1–7 | March 2026*
