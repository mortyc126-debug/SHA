# Комбинированная стратегия атаки SHA-256
## Structured Chaos Framework — Synthesis

---

## Все найденные инструменты

| # | Инструмент | Что даёт | Раунды | Стоимость |
|---|-----------|----------|--------|-----------|
| 1 | **Wang chain** (forward) | δe[2..16]=0 | r=2..16 (15) | O(1) adaptive |
| 2 | **Integral** (k=16) | 2 balanced бита | r=0..6 | 2^16 evaluations |
| 3 | **Schedule kernel** (R=52) | δW_xor[52..63]=0 | r=52..63 (12 слов) | dim=128 free |
| 4 | **Schedule kernel** (R=56) | δW_xor[56..63]=0 | r=56..63 (8 слов) | dim=256, **defect=3** |
| 5 | **Carry correlation** | 50× joint P for 3 words | schedule | структурная |
| 6 | **Null vectors** (defect) | HW(δH)=96 start | output | SA optimization |
| 7 | **Backward chain** | δe[r-3]=−δW[r] | от конца | нужен δW≈0 |
| 8 | **TLC decomposition** | Q:C=34:66 | все раунды | аналитическая |
| 9 | **Carry suppression** | +43 бита при r=4 | r=0..5 | state-dependent |
| 10 | **Carleman** | rank=state_dim | все раунды | теоретическая |

---

## Карта контроля (64 раунда)

```
Раунд:  0  1  2  3  4  5  6  7  ...  16 17 18 ... 51 52 53 ... 63
        ├──┤
        Integral (6R, 2 balanced при R=6)
           ├──────────────────────────┤
           Wang chain (15R, δe=0, FREE)
                                       ├──┤
                                       Барьер r=17 (O(2^32))
                                          ├──────────────┤
                                          THE GAP (~34R)
                                                         ├──────────┤
                                                         Kernel δW_xor=0
                                                         (12 слов GF2-zero)
                                                         ├──────────┤
                                                         Backward chain
                                                         (если δW_real≈0)
```

**Forward:** 17 раундов (0..16 = integral/Wang)
**Backward:** 0-12 раундов (зависит от carry в schedule)
**Gap:** 64 − 17 − backward = 47..35 раундов

---

## Комбинация 1: Wang + Kernel (без overlap)

```
Phase 1 (r=0..1):  δ входит через W[0]
Phase 2 (r=2..16): Wang chain → δe=0 (FREE, adaptive δW)
Phase 3 (r=17):    Birthday barrier → O(2^32)
Phase 4 (r=18..51): THE GAP — неконтролируемый хаос
Phase 5 (r=52..63): Schedule kernel → δW_xor=0, δW_real≈13 HW/word

Проблема: Phase 4 и Phase 5 НЕ СВЯЗАНЫ.
Kernel даёт δW_xor=0, но carry заполняет реальные δW.
Backward chain не работает при δW_real≈13.
```

**Стоимость: O(2^32) barrier + O(2^128) gap = O(2^128) (не лучше birthday)**

---

## Комбинация 2: Wang + Kernel + Null Vector SA

```
Phase 1-3: как Комбинация 1 (Wang + barrier)
Phase 4: Вместо random gap, используем NULL VECTOR стартовую точку
  - Null vector: df=0 в линейном приближении → HW(δH)≈112
  - SA от null vector: HW(δH)=96 (vs 104 от random)
  - Advantage: 8 бит

Стоимость: O(2^32) + O(2^{120}) = O(2^{120})
Экономия: 8 бит (vs birthday 2^128)
```

**Маргинальное улучшение — 8 бит.**

---

## Комбинация 3: Multiblock + Kernel + Integral

```
Block 1:
  - Используем integral: 2^16 сообщений с k=16 active bits
  - Balanced bits при R=6 → structured intermediate state H1
  - Стоимость: O(2^16)

Block 2:
  - IV = H1 (from block 1)
  - Используем schedule kernel (dim=256 at R=56)
  - Null vectors дают HW(δH)=96 стартовую точку
  - SA optimization
  - Стоимость: O(2^{120})

Total: O(2^{120})
```

**Чуть лучше, но integral не помогает block 2 напрямую.**

---

## Комбинация 4: МАКСИМАЛЬНАЯ (все инструменты)

```
Block 1: Integral set generation
  - 2^16 messages, k=16 active bits
  - Collect intermediate states at R=6 (balanced)
  - Cost: O(2^16)

Block 2: Wang + Kernel + Defect
  - IV = H1 from block 1
  - Forward: Wang chain on W[0..15] → 15 rounds δe=0
  - δW constrained to schedule kernel (dim=256, R=56)
    → δW_xor[56..63]=0
    → carry correlation 50× helps → P(low carry) enhanced
  - Barrier at r=17: O(2^32) birthday
  - Backward from end: null vector + SA → HW=96
  - Defect=3: hash output in 253-dim subspace
  - Cost: O(2^32) + O(2^{125.5})

Total: O(2^{125.5})
Savings: 2.5 bits (from defect=3)
```

---

## Главная проблема: GAP

Все комбинации упираются в **THE GAP** (r=17..51, ~34 раунда).

В этой зоне:
- Carry saturation: τ=1.8 → полная сатурация
- No δW=0 (schedule diffusion)
- No algebraic structure (degree=32, A maximally complex)
- No invariants (all biases < 0.05)

**Gap — это ЕДИНСТВЕННАЯ проблема.** Все наши инструменты покрывают r=0..16 (forward) и r=52..63 (backward/kernel). Gap r=17..51 = 34 раунда чистого хаоса.

---

## Что нужно для прорыва

| Сценарий | Что нужно | Стоимость |
|----------|----------|-----------|
| **Close gap** | Расширить forward ИЛИ backward на 34+ раундов | Невозможно текущими методами |
| **Bridge gap** | Найти коррелятор через 34 раунда хаоса | Неизвестно |
| **Reduce gap** | Расширить Wang до r=25 И backward до r=45 | Нужно: Wang-like chain для r>17 + δW_real≈0 для r<52 |
| **Bypass gap** | Алгебраическая атака, игнорирующая раундовую структуру | Carleman не работает (A max complex) |

### Конкретная цель

Для **O(2^{100})** (28 бит ниже birthday):
- Forward: 17 раундов (есть)
- Backward: нужно 19 раундов (сейчас: 0 реальных)
- Gap: 64−17−19 = 28 раундов
- Cost per gap round: ~4 бита → 28×4 = 112 > 100

Для **O(2^{64})** (64 бита ниже birthday, ПРОРЫВ):
- Forward: 17 раундов
- Backward: нужно 31 раунд
- Gap: 16 раундов
- Cost: 16×4 = 64 бит → O(2^{64}) ✓
- Нужно: backward chain на 31 раунд → δW_real≈0 для 31 слов schedule
- Стоимость δW_real=0: ~2^{16}/word × 31 words ≈ 2^{21} (если независимы)
- С carry correlation 50×: ~2^{15} (незначительная экономия)

---

## Вердикт

**Текущий best: O(2^{125.5}) = birthday − 2.5 бита.**

Для настоящего прорыва нужно решить ОДНУ из задач:
1. **Backward 31+ раундов** → нужен δW_real≈0 в schedule tail 31 слово
2. **Bridge 34 раундов** → нужен коррелятор через хаотическую зону
3. **Новая алгебра** → нужна структура, работающая ЧЕРЕЗ сатурацию

Задача 1 сводится к: найти W_base, где carry в 31 последовательном слове schedule одновременно мал. С текущей carry-корреляцией (50× на 3 слова): P ≈ 2^{-160} для 31 слова → хуже birthday.

**Вывод: SHA-256 устойчива к комбинации ВСЕХ наших техник.**
Максимальная экономия: 2.5 бита (из defect=3).

---

*Combined Strategy v1.0 | 25.03.2026*
*Файлы: 17 экспериментов в alg_attack/scf_*.py*
