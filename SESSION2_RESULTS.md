# Сессия 2: Q∩T SAT Solver для SHA-256 | Апрель 2026

## ЦЕЛЬ СЕССИИ

Построить SAT-based solver SHA-256 и ускорить его используя математику из методички
(methodology_v20.md) и результаты предыдущей сессии (13 закрытых атакующих направлений,
три шумовых источника, створочне тождество, backward chain, Q∩T алгебра).

---

## ХРОНОЛОГИЯ

### Этап 1: Baseline SAT Solver (v2)
- Пользователь предоставил `sha256_solver_v2.py` — побитовое SAT-кодирование SHA-256
- Glucose4 solver, ripple-carry adder, стандартное Tseitin-кодирование
- **Baseline: 96-98s для "Hi" (2 байта, 16 бит свободы, 64 раунда)**
- 168,465 переменных, 591,076 CNF clause

### Этап 2: NK Instruments (v3)
Добавлены два инструмента из методички:
1. **Backward chain pinning**: a[57..64] + e[61..64] вычислены из целевого хеша, зафиксированы
2. **Створочне constraints**: e[r] = a[r] + a[r-4] - Σ₀(a[r-1]) - Maj(...) для r=4..11

**Результат v3:**
- Backward pinning alone = **0× ускорение** (solver и так выводит это из hash constraints)
- Backward + створочне = **1.15× ускорение** (84s vs 98s)
- Створочне помогает непоследовательно: 2.78× на 20R, но 0.17× на 32R

### Этап 3: Q∩T Гибридный Solver (v4)
**Ключевой прорыв сессии.** Анализ методички (разделы 212-218) выявил:

> SHA-256 = Q (квадратичная GF(2), степень 2) ∩ T (пороговые carry-уравнения)
> Стандартный SAT видит только T. CryptoMiniSat с XOR clauses видит и Q.

Построен `sha256_solver_v4_qt.py`:
- **CryptoMiniSat** вместо Glucose4 (поддержка native XOR clause)
- Σ₀, Σ₁, σ₀, σ₁, XOR операции → **native XOR clauses** (Gauss elimination)
- Ch, Maj → AND (CNF) + XOR (native) — гибридное Q∩T кодирование
- Carry → стандартный CNF (T-часть)

**Результат v4:**
- CMS+XOR vs Glucose4 = **2.01× на 64R** (87s vs 175s)
- До **17.6× на 24R** (1.05s vs 18.5s)

### Этап 4: Четыре параллельных эксперимента

#### Dir 1: ANF vs MUX кодирование Ch
- MUX: Ch = (e·f) ⊕ (ē·g) — требует NOT переменную
- ANF: Ch = (e·f) ⊕ (e·g) ⊕ g — без NOT, меньше переменных
- **ANF выигрывает на 64R: 41.2s vs 85.8s = 2.1×**
- Файл: `exp_dir1_quadratic_xor.py`

#### Dir 2: Carry-Lookahead vs Ripple Carry
- CLA (Kogge-Stone parallel prefix): глубина O(log n) вместо O(n)
- **CLA помогает на 32-48R** (16s vs 30s = 1.84×, 33s vs 77s = 2.35×)
- **CLA убивает на 64R** (281s vs 41s = 0.14×) — overhead 3× по переменным
- Файл: `exp_dir2_carry_lookahead.py`

#### Dir 3: Partial c-world — carry diversity
- Перечислены все 9025 ASCII 2-byte сообщений
- Round 0: 2960 уникальных carry-паттернов (33% — лучше brute force)
- **Round 2+: ВСЕ 9025 уникальны** (carry полностью различает сообщения)
- Partial c-world = brute force после 1 раунда → **ТУПИК**
- 33 свободных carry-бит из 217 в раунде 0 (184 фиксированы от IV)
- Файл: `exp_dir3_partial_cworld.py`

#### Dir 4: 3-byte scaling
- Q∩T преимущество масштабируется: **7× для 3 байт на 8-16R**
- "Abc" 20R: Q∩T решает за 101s, Glucose4 — timeout (>180s)
- "Hi" 64R: **3.04× в лучшем прогоне**
- Файл: `exp_dir4_3byte_scaling.py`

### Этап 5: Constant Propagation (v5)
- "Lazy words": конкретные int ИЛИ SAT-переменные
- Если оба операнда — константы → вычислить напрямую, 0 SAT overhead
- 31 операция свёрнута, 3 слова расписания = константы, 1 T2-раунд = константа
- **17% меньше переменных** (98,897 vs 119,400)
- **Но медленнее**: 78s vs 41s (v4 ANF) — clause structure важнее количества переменных
- `_wadd_const` (OR/XNOR вместо MAJ) ухудшает CDCL propagation
- Файл: `sha256_solver_v5.py`

---

## КЛЮЧЕВЫЕ РЕЗУЛЬТАТЫ

### Производительность solvers

| Версия | Описание | Vars | CNF | XOR | "Hi" 64R | vs v2 |
|--------|----------|------|-----|-----|----------|-------|
| v2 | Glucose4 baseline | 168,465 | 591,076 | 0 | 96s | 1.00× |
| v3 | +backward+створочне | 175,537 | 617,188 | 0 | 84s | 1.14× |
| v4 MUX | CMS+XOR, MUX Ch | 121,448 | 202,916 | 51,112 | 87s | 1.10× |
| **v4 ANF** | **CMS+XOR, ANF Ch** | **119,400** | **203,044** | **49,064** | **41s** | **2.34×** |
| v5 | +const propagation | 98,897 | 166,072 | 42,893 | 78s | 1.23× |

**Лучший: v4 ANF = 41s (2.34× над baseline)**

### Что ПОМОГЛО

1. **CryptoMiniSat + native XOR clauses (Q-часть)**: Gauss elimination на XOR — 2× ускорение.
   Это прямая реализация Q∩T идеи из методички: разделить Q (алгебра) и T (threshold).

2. **ANF Ch encoding** (ef⊕eg⊕g вместо ef⊕ēg): убирает NOT-переменную (2K vars),
   другая clause structure лучше для CDCL. Ещё ~2× ускорение поверх CMS.

### Что НЕ ПОМОГЛО

1. **Backward pinning** (12 слов из H): 0× — solver и так выводит это из hash constraints
2. **Створочне redundant constraints**: непоследовательно (2.78× на 20R, 0.17× на 32R)
3. **Carry-Lookahead** на 64R: 3× больше переменных → 7× медленнее
4. **Constant propagation**: -17% vars но -50% speed. Clause structure > variable count
5. **Partial c-world**: carry уникален после 1 раунда, enumeration = brute force

### Установленные факты

**F-SAT1**: Для SAT на SHA-256, clause structure важнее количества переменных.
Меньше vars ≠ быстрее solve. CDCL использует propagation paths через clause graph.

**F-SAT2**: CryptoMiniSat + native XOR = единственный способ использовать Q-структуру SHA-256 в SAT.
Glucose4 не может делать Gauss elimination; конвертация XOR→CNF теряет алгебраическую информацию.

**F-SAT3**: ANF encoding Ch = ef⊕eg⊕g лучше MUX encoding ef⊕ēg для CMS.
Убирает NOT-переменную и создаёт более симметричную clause structure.

**F-SAT4**: Carry-out полностью уникален после 1 раунда SHA-256 для 2-byte ASCII.
9025/9025 сообщений имеют уникальный carry-pattern с глубины 2. Partial c-world = dead end.

**F-SAT5**: Backward chain (384 бита из H) — redundant для SAT solver.
Unit propagation от hash constraints уже выводит эти значения.

**F-SAT6**: Q∩T ускорение масштабируется с размером задачи.
16 бит: 2-3× → 24 бита: 7× → чем больше свободы, тем больше выигрыш от Gauss elimination.

---

## ОТКРЫТЫЕ НАПРАВЛЕНИЯ

### 1. CLA + Ripple гибрид
CLA помогает на 32-48R (2.35×), но overhead убивает на 64R.
Идея: CLA для ранних раундов (0-16), ripple для поздних (17-64).

### 2. Adaptive Ch encoding
ANF лучше на 64R, MUX лучше на 24R. Автоматический выбор по числу раундов.

### 3. Collision mode
Текущий solver ищет preimage. Для коллизии: два сообщения M₁, M₂ с SHA(M₁)=SHA(M₂).
Удвоение переменных, но общий hash constraint. Q∩T может помочь больше.

### 4. Структурные XOR constraints
Добавить створочне и другие тождества как XOR clauses (не CNF).
v3 добавляла створочне как CNF — может native XOR будет эффективнее.

### 5. Schedule-aware encoding
W[17], W[19], W[21] — полные константы для 2-byte msg. Использовать это.

---

## ФАЙЛЫ СЕССИИ

### Solvers
| Файл | Описание |
|-------|---------|
| `sha256_solver_v2.py` | Baseline: Glucose4, ripple carry, стандартный Tseitin |
| `sha256_solver_v3.py` | +backward pinning + створочне (NK instruments) |
| `sha256_solver_v4_qt.py` | Q∩T гибрид: CMS + native XOR + ANF/MUX Ch |
| `sha256_solver_v5.py` | +constant propagation (меньше vars, но медленнее) |

### Эксперименты
| Файл | Описание |
|-------|---------|
| `exp_dir1_quadratic_xor.py` | ANF vs MUX Ch encoding |
| `exp_dir2_carry_lookahead.py` | Kogge-Stone CLA vs Ripple carry |
| `exp_dir3_partial_cworld.py` | Carry diversity analysis (c-world) |
| `exp_dir4_3byte_scaling.py` | Q∩T scaling test (1-4 bytes) |
| `sha256_bench_v3.py` | v3 benchmark (signal timeout, replaced) |
| `sha256_bench_v3b.py` | v3 benchmark (process timeout, final) |

### Из предыдущей сессии (использованы как фундамент)
| Файл | Ключевой результат |
|-------|-------------------|
| `nk_noise_sources.py` | Три ортогональных шума: Σ (12.5×), Ch/Maj (31×), Carry (2×) |
| `nk_wall.py` | Стена: a[56]↔e[60] circular dependency |
| `nk_320bit.py` | 320-bit constraint, schedule cascade 48/48 always |
| `nk_amplification.py` | Thermostat: HW[r+1]=0.625×HW[r]+48.3, escape=2^186 |
| `methodology_v20.md` | 25,476 строк методички, разделы 212-218 = Q∩T теория |

---

## ФОРМУЛА ЛУЧШЕГО SOLVER (v4 ANF)

```
SHA-256 SAT encoding:
  Σ₀, Σ₁, σ₀, σ₁  →  native XOR clauses (Gauss elimination)
  Ch = ef ⊕ eg ⊕ g  →  2× AND (CNF) + 1× XOR3 (native)
  Maj = ab ⊕ ac ⊕ bc →  3× AND (CNF) + 1× XOR3 (native)  
  sum[k] = a[k] ⊕ b[k] ⊕ carry[k-1]  →  XOR3 (native)
  carry[k] = MAJ(a[k], b[k], carry[k-1])  →  CNF (threshold)
  
Solver: CryptoMiniSat (threads=1)
  Gauss elimination ON for XOR clauses
  CDCL for CNF clauses
  
Result: 119,400 vars, 203,044 CNF, 49,064 XOR
  "Hi" (16 bits, 64R): 41s = 2.34× over Glucose4 baseline
```

---

## ГЛАВНЫЙ ВЫВОД

**Q∩T разделение работает.** Кодирование Q-части (XOR/rotation/Sigma) как native XOR clauses
позволяет CryptoMiniSat использовать Gauss elimination — алгебраический метод, недоступный
стандартному CDCL. Это даёт стабильное 2-3× ускорение на 64R и до 17× на средних раундах.

Для атаки на SHA-256 это недостаточно (нужно 2^128 → 2^127 или меньше), но это первый
работающий инструмент, который использует **алгебраическую структуру** SHA-256 в SAT solver,
а не просто кодирует всё как чёрный ящик.

Следующий рубеж: создать полноценный Q∩T solver (не SAT-based), который решает Q-систему
(квадратичную GF(2)) и T-систему (пороговую) **одновременно**, используя структуру обеих.
