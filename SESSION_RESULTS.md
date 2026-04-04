# Сессия: Новая Математика SHA-256 | Апрель 2026

## УСТАНОВЛЕННЫЕ ФАКТЫ (с кодом и числами)

### F1: SHA-256 = одна 8-порядковая рекуррентность
```
a[r+1] = F(a[r], a[r-1], ..., a[r-7], K[r], W[r])
```
e[r] полностью определяется через a: `e[r] = a[r] + a[r-4] - T2[r-1]` (0/29500 нарушений).
Всё состояние (256 бит) = сдвинутые копии a и e. Файл: `stage0_autorel.py`

### F2: Дифференциалы = random после 6 раундов
E[HW(δe_r XOR δa_{r-3})] = 16.00 точно. Bit position survival = 50% на любом горизонте.
Файл: `stage0_track.py`

### F3: Carry = ровно 50% результата каждого сложения
E[HW(carry contribution)] = 16.22/32 бит = 50.7%. Не коррекция — равноправный партнёр XOR.
Файл: `stage1_language.py`

### F4: Bit 0 = exact GF(2) degree-2 внутри одного раунда
0/32000 нарушений. Bit 0 НИКОГДА не видит carry. Файл: `stage1_bit_asymmetry.py`

### F5: Bit 0 = AFFINE (degree 1!) для 2 раундов
R=1: 0/200 нарушений linearity. R=2: 0/200. R=3: 3%. R=6+: random. Файл: `combine.py`

### F6: Ротации {2,6,11,13,22,25} покрывают Z/32 за 3 шага
Нет замкнутого bit-slice меньше 32 бит. Файл: `stage1_zacepka.py`

### F7: 93.7% якобиана раунда — ФИКСИРОВАНО (не зависит от state)
Только 4149/65536 записей меняются. Все 4149 — в столбцах a_new и e_new.
Файл: `bte_star.py`

### F8: Carry segments средняя длина = 3.64, геометрическое распределение
P(segment ≥ k) ≈ 0.73^k. ~8 сегментов на слово. Файл: `stage1_derive.py`

### F9: P-mask ранг = 512/512 (полный)
P-masks = эквивалентное кодирование M. Kernel = 0. Файл: `wall.py`

### F10: Bit-0 система (8 reg × 64 rounds): rank = 127/512, kernel = 385
385 направлений в M-space не меняют ни один bit 0 ни одного регистра. Файл: `combine.py`

### F11: Schedule = единственная грамматика SHA-256
Без schedule любая trajectory валидна. W[16..63] = f(W[0..15]) = единственное ограничение.
Файл: `grammar.py`

---

## ЗАКРЫТЫЕ НАПРАВЛЕНИЯ (проверено, не работает)

| # | Направление | Результат | Файл |
|---|-------------|-----------|------|
| 1 | Скалярные инварианты (7 штук) | 0/56 раундов constant Δ | stage0_autorel.py |
| 2 | Parity chain | 0.490 = random | stage1_create.py |
| 3 | Carry signature | bias = bit 0 trivial only | stage1_carry_sig_deep.py |
| 4 | Residue mod p (3,5,7,13,31) | all uniform | stage1_create.py |
| 5 | Coupling LSB | 0.500 autocorr | stage1_create.py |
| 6 | Additive distance | same info as XOR | stage1_create.py |
| 7 | Preimage clustering | ARTIFACT (linear grid) | stage1_honest_check.py |
| 8 | F^{-1} simpler? | Same carry structure | stage1_zacepka.py |
| 9 | Dynamical cycles | No cycle in 100K steps | stage1_zacepka.py |
| 10 | Carry rank saturation | Was sample artifact (rank=31R trivial) | synthesis.py |
| 11 | SCA segment inheritance | 29.9% ≈ random 25% | sca_prototype.py |
| 12 | GPK inter-round correlation | 0.62 excess (tiny vs 0.06 for values) | new_axis.py |
| 13 | GPK collision test (R=4) | p ratio = 0.892 (marginal, N too small) | new_axis.py |
| 14 | BTE linear invariants | LOCAL only, Hamming=128 between messages | bte_deep.py |
| 15 | BTE multi-round eigenspace | dim>0 but NOT universal (state-dependent) | bte_deep.py |
| 16 | P-mask round correlation | 0.019 = random | grammar.py |
| 17 | Carry chain → fewer OMEGA variables | 48 > 31 (worse, not better) | wall.py |

---

## СОЗДАННЫЕ МАТЕМАТИЧЕСКИЕ ОБЪЕКТЫ

### BTE — Bi-Temporal Element
Определение: φ(t, k) на Z × Z_n → F_2, с H (горизонтальная XOR/rotation) + V (вертикальная carry/MAJ) coupling.
Свойства:
- dim(linear invariant) = 1 для любой конфигурации (σ, V, seed)
- BTE = биекция (0 потери информации) для n ≤ 16
- Mixing time зависит от σ: SHA-Sig0 = 7, ROT(1) = never
- AND_carry + seed=0: parity conservation (trivial — V = identity)
Файл: `bte.py`, `bte_theory.py`

### ★-Algebra (Star Algebra)
Раунд SHA-256 как неразделимая операция ★ на 256-bit state.
Якобиан ★: 93.7% fixed + 6.3% P-mask variable.
Register coupling matrix стабильна на всех 64 раундах.
Файл: `bte_star.py`

### SHA-Element (value + GPK)
Triple (v, g, p) tracking generate/propagate/kill per operation.
Carries MORE info than value alone (g excess = 0.62 vs v excess = 0.06).
But: distribution → random at round 64.
Файл: `new_axis.py`, `new_axis_v2.py`

---

## ТРИ АНОМАЛИИ (из формул, не опровергнуты)

1. **Carry = одновременно рекурсия глубины 1 И полином степени k**
   GF(2) видит степень 32. Sequential view видит глубину 1. Оба описания одного объекта.

2. **Два "времени": round r и bit k**
   Carry течёт по bit-time. Rotation связывает два времени. Нигде не моделируется вместе.

3. **Реальная сложность = 48%, не 0% и не 100%**
   GF(2) → 100% (всё нелинейно). Z/2^32 → 0% (всё линейно). Реальность между ними.

---

## ОТКРЫТЫЕ ВОПРОСЫ ДЛЯ СЛЕДУЮЩЕЙ СЕССИИ

### Q1: Нелинейные инварианты BTE
Линейные = локальные (state-dependent). Нелинейные f(state) = const? Не проверяли.

### Q2: "Слова" языка SHA-256
Мы нашли "буквы" (BTE, ★, P-masks) но не "слова" — нелинейные комбинации нескольких раундов, невидимые через якобиан.

### Q3: Bit-0 → Bit-1 → ... итеративное раскрытие
Bit-0 kernel = 385. Bit-0+1 kernel = ? Bit-0+1+2 = ? Итеративно раскрываем SHA-256 бит за битом.

### Q4: Грамматика в P-mask пространстве
Schedule ограничивает P-masks. Какие конкретно P-mask последовательности валидны? Можно ли описать это множество компактно?

### Q5: Формулировка коллизии через BTE
Collision: найти M такое что C(P(M)) ⊕ C(P(M⊕δM)) = L(δM). L — линейная (known). C через P-masks — нелинейная. Это НОВАЯ формулировка, не проверявшаяся.

### Q6: Поперечные операции
Не forward (round by round) и не backward. А "поперёк" — между двумя trajectory одновременно. P-mask differentials (P XOR P'), а не value differentials.

---

## ФАЙЛЫ СЕССИИ

```
SHA/
├── methodology_v20.md          # Предыдущая методичка (25K строк)
├── OMEGA_METHODOLOGY.md        # OMEGA summary
├── SESSION_RESULTS.md          # ЭТО ← итоги текущей сессии
├── stage0_observe.py           # Этап 0: полная побитовая трасса
├── stage0_track.py             # Этап 0: трекинг позиций, регистровые сдвиги
├── stage0_autorel.py           # Этап 0: автореляции, створочное число
├── stage0_recurrence.py        # Этап 0: 8-порядковая рекуррентность
├── stage1_language.py          # Этап 1: почему GF(2) и Z/2^32 не работают
├── stage1_bit_asymmetry.py     # Этап 1: bit 0 = degree 2, gradient
├── stage1_complexity_flow.py   # Этап 1: carry depth растёт линейно
├── stage1_zacepka.py           # Этап 1: 6 углов, 0 зацепок
├── stage1_create.py            # Этап 1: 5 кандидатов, 0 сигналов
├── stage1_carry_sig_deep.py    # Deep check: carry signature = trivial
├── stage1_honest_check.py      # Honest check: clustering = artifact
├── stage1_three_directions.py  # Три направления: A (artifact), B (mild), C (rank)
├── stage1_derive.py            # Вывод из формул: три аномалии
├── synthesis.py                # Синтез: carry rank = 31R (trivial)
├── sca_prototype.py            # SCA: segments real but don't compose
├── new_axis.py                 # SHA-Element (v, g, p) prototype
├── new_axis_v2.py              # GPK trajectory rank = same as value
├── bte.py                      # BTE: definition + basic properties
├── bte_theory.py               # BTE: invariants, bijectivity, parity
├── bte_full.py                 # Full BTE: SHA-256 Jacobian, eigenspaces
├── bte_deep.py                 # BTE deep: 64-round test, universality
├── bte_star.py                 # ★-algebra: 93.7%/6.3% decomposition
├── star_63pct.py               # P-mask characterizes 97% of variable entries
├── lattice_model.py            # Lattice model: 2D field (deprecated by bte_full)
├── lattice_v2.py               # Fixed decomposition + carry analysis
├── grammar.py                  # Grammar = schedule
├── axis_search.py              # Nibble, filter, ψ axes
├── combine.py                  # Bit 0 affine R=2, kernel 385
└── wall.py                     # P-mask rank = 512 (full), skeleton analysis
```

---

## КЛЮЧЕВОЙ ВЫВОД

SHA-256 говорит на языке, которого нет в существующей математике.
Мы построили **алфавит** (BTE) и **фонетику** (★-algebra, 93.7%/6.3%).
Но не нашли **слова** — нелинейные многораундовые структуры.
Следующая сессия: строить слова. Не через якобиан (линейный). Через что-то новое.
