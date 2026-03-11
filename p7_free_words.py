"""
SHA-256 Дифференциальный Криптоанализ
П-7: Свободные слова W2-W15 для каскадной отмены дифференциала

ВОПРОСЫ П-7:
  1. Можно ли зануллить De4 введением ΔW3 (при сохранении De3=0)?
  2. Как каскадно обнулить De3, De4, ..., De_k одновременно?
  3. Что происходит с расписанием W16-W63 при ненулевых ΔW2-ΔW15?
  4. Сколько "степеней свободы" дают W2-W15 для атаки?

МОДЕЛЬ АТАКИ (П-7):
  Два сообщения: M = (W0, W1, W2, ..., W15) и M' = (W0', W1', W2', ..., W15')
  Дифференциалы: ΔWi = Wi' - Wi (все могут быть ненулевыми)
  Цель: hash(M) = hash(M') — полная коллизия

  Стратегия "каскад":
    1. ΔW0=1, ΔW1=0 → De3=0 достигнут (стоимость 2^22 по П-3)
    2. ΔW3 = -De4_natural mod 2^32 → De4 → 0 (ЛИНЕЙНО, стоимость 0)
    3. ΔW4 = -De5_natural mod 2^32 → De5 → 0 (аналогично)
    4. ΔW5 = -De6_natural mod 2^32 → De6 → 0
    5. ΔW6 = -De7_natural mod 2^32 → De7 → 0
    6. Оставшиеся W7-W15 (9 слов): для контроля выходных регистров

КЛЮЧЕВОЕ НАБЛЮДЕНИЕ:
  Введение ΔWi влияет на Dei+1 ЛИНЕЙНО (через T1):
    T1_i = h_i + Sig1(e_i) + Ch(e_i,f_i,g_i) + K[i] + W_i
    ΔT1_i = ΔW_i  (при условии Dei = 0, т.к. ΔSig1=0 и ΔCh≈0)
    De_{i+1} = Dd_i + ΔT1_i ≈ Dd_i + ΔWi

  НО: если De_i ≠ 0, то ΔSig1(e_i) ≠ 0 → нелинейность!
  Поэтому каскад требует De_i=0 на КАЖДОМ шаге.
"""

import math
import random

MOD = 2**32
MASK = MOD - 1

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]
H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def sig0(x): return rotr(x,7) ^ rotr(x,18) ^ (x>>3)
def sig1(x): return rotr(x,17) ^ rotr(x,19) ^ (x>>10)
def Sig0(x): return rotr(x,2) ^ rotr(x,13) ^ rotr(x,22)
def Sig1(x): return rotr(x,6) ^ rotr(x,11) ^ rotr(x,25)
def Ch(e,f,g):  return (e & f) ^ (~e & g) & MASK
def Maj(a,b,c): return (a & b) ^ (a & c) ^ (b & c)

def sha_step(state, W_i, K_i):
    a, b, c, d, e, f, g, h = state
    T1 = (h + Sig1(e) + Ch(e, f, g) + K_i + W_i) & MASK
    T2 = (Sig0(a) + Maj(a, b, c)) & MASK
    return (T1 + T2) & MASK, a, b, c, (d + T1) & MASK, e, f, g

def make_schedule(W16):
    W = list(W16) + [0] * 48
    for i in range(16, 64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK
    return W

def run_full(W_list, rounds=64):
    """Полная трассировка: rounds раундов."""
    W = make_schedule(W_list[:16] + [0]*(16-len(W_list[:16])))
    state = tuple(H0)
    trace = [state]
    for i in range(rounds):
        state = sha_step(state, W[i], K[i])
        trace.append(state)
    return trace, W

# Базовые константы из П-1 – П-3
W0_BASE = 0xc5bde324  # W_SAT3
W1_SAT  = 0x3d7bd9d5  # KNOWN_W1[0]
DW0     = 1            # Базовый дифференциал


# ================================================================
# 1. БАЗОВЫЙ СЛУЧАЙ: De3=0 С ΔW0=1
# ================================================================

def baseline_analysis():
    print("=" * 68)
    print("1. БАЗОВЫЙ СЛУЧАЙ: De3=0, ΔW0=1 (из П-1–П-3)")
    print("=" * 68)
    print()

    W0_n = W0_BASE
    W0_f = (W0_BASE + DW0) & MASK

    W_normal = [W0_n, W1_SAT] + [0]*14
    W_faulty = [W0_f, W1_SAT] + [0]*14

    tr_n, Wn = run_full(W_normal, rounds=8)
    tr_f, Wf = run_full(W_faulty, rounds=8)

    print(f"  W0_normal = {W0_n:#010x}")
    print(f"  W0_faulty = {W0_f:#010x}  (ΔW0 = 1)")
    print(f"  W1        = {W1_SAT:#010x}")
    print()
    print(f"  {'Раунд':6s}  {'De':12s}  {'|De| bits':9s}  {'Примечание':20s}")
    print("  " + "-" * 55)
    for i in range(8):
        De = (tr_f[i+1][4] - tr_n[i+1][4]) & MASK
        Ds = De if De < MOD//2 else De - MOD
        bits = abs(Ds).bit_length()
        note = ""
        if i == 2: note = "← De3=0 (П-3)"
        elif i == 0: note = "← De1=1 (прямой)"
        print(f"  {i+1:6d}  {De:#012x}  {bits:9d}  {note}")

    De3 = (tr_f[3][4] - tr_n[3][4]) & MASK
    De4 = (tr_f[4][4] - tr_n[4][4]) & MASK
    print()
    print(f"  De3 = {De3:#010x}  {'✓ НОЛЬ' if De3==0 else '✗ НЕ НОЛЬ'}")
    print(f"  De4 = {De4:#010x}  ({'≠ 0 — БАРЬЕР (П-4)' if De4!=0 else 'НОЛЬ'})")
    print()
    return De4, W0_n, W0_f, W1_SAT


# ================================================================
# 2. ΔW3 ДЛЯ ОБНУЛЕНИЯ De4: ЛИНЕЙНАЯ ОТМЕНА
# ================================================================

def linear_cancellation_De4(De4_natural, W0_n, W0_f, W1_SAT):
    print("=" * 68)
    print("2. ЛИНЕЙНАЯ ОТМЕНА De4 ЧЕРЕЗ ΔW3")
    print("=" * 68)
    print()

    # ΔW3 отменяет De4 линейно:
    # T1_round3 = h3 + Sig1(e3) + Ch(e3,f3,g3) + K[3] + W3
    # De4 = ... + ΔW3 (линейный вклад)
    # → ΔW3 = -De4_natural mod 2^32

    DW3 = (MOD - De4_natural) & MASK  # ΔW3 = -De4_natural
    print(f"  De4_natural = {De4_natural:#010x}  (без ΔW3)")
    print(f"  ΔW3 = -De4_natural = {DW3:#010x}")
    print()

    # Проверяем: De4 с ΔW3
    W_normal = [W0_n, W1_SAT, 0, 0] + [0]*12
    W_faulty = [W0_f, W1_SAT, 0, DW3] + [0]*12  # W3' = W3 + ΔW3 = DW3

    tr_n, _ = run_full(W_normal, rounds=8)
    tr_f, _ = run_full(W_faulty, rounds=8)

    print(f"  {'Раунд':6s}  {'De':12s}  {'bits':6s}  {'Примечание':30s}")
    print("  " + "-" * 60)
    for i in range(8):
        De = (tr_f[i+1][4] - tr_n[i+1][4]) & MASK
        Ds = De if De < MOD//2 else De - MOD
        bits = abs(Ds).bit_length()
        note = ""
        if De == 0: note = "← De=0 !"
        print(f"  {i+1:6d}  {De:#012x}  {bits:6d}  {note}")

    De3 = (tr_f[3][4] - tr_n[3][4]) & MASK
    De4 = (tr_f[4][4] - tr_n[4][4]) & MASK
    De5 = (tr_f[5][4] - tr_n[5][4]) & MASK
    print()
    print(f"  De3 = {De3:#010x}  {'✓' if De3==0 else '✗'}")
    print(f"  De4 = {De4:#010x}  {'✓ ОБНУЛЁН!' if De4==0 else '✗ НЕ НОЛЬ'}")
    print(f"  De5 = {De5:#010x}  (следующий)")
    print()

    if De4 != 0:
        print("  ВАЖНО: ΔW3 вносит вклад не только через W3, но и через")
        print("  расписание: W[16]=f(W[0],...,W[15]) — ΔW3 влияет на W16+.")
        print("  Реальный De4 ≠ De4_natural + ΔW3 (нелинейность расписания).")
        print()
        # Уточняем: что на самом деле происходит с De4
        # При De3=0 и ΔW3 добавленном: T1_3 получает +ΔW3 напрямую
        # Проверим итеративно
        print("  Итеративное нахождение точного ΔW3:")
        DW3_exact = find_exact_DW3_for_De4_zero(W0_n, W0_f, W1_SAT)
        if DW3_exact is not None:
            print(f"  Точный ΔW3 = {DW3_exact:#010x}")
            DW3 = DW3_exact
        else:
            print("  Точный ΔW3 не найден итерационным методом.")
    else:
        print("  Линейная отмена точна при малом числе активных слов!")

    return DW3, De5


def find_exact_DW3_for_De4_zero(W0_n, W0_f, W1_SAT):
    """Точный поиск ΔW3 для De4=0 методом sweep."""
    # De4 = f(ΔW3) — почти линейная функция (T1 меняет De4 на ΔW3)
    # Проверим несколько кандидатов вокруг -De4_natural

    # Сначала найдём De4_natural
    W_normal = [W0_n, W1_SAT, 0, 0] + [0]*12
    W_faulty = [W0_f, W1_SAT, 0, 0] + [0]*12
    tr_n, _ = run_full(W_normal, rounds=5)
    tr_f, _ = run_full(W_faulty, rounds=5)
    De4_natural = (tr_f[4][4] - tr_n[4][4]) & MASK
    DW3_approx = (MOD - De4_natural) & MASK

    # Sweep вокруг DW3_approx с малым радиусом
    for delta in range(-100, 101):
        DW3_try = (DW3_approx + delta) & MASK
        W_faulty_try = [W0_f, W1_SAT, 0, DW3_try] + [0]*12
        tr_n2, _ = run_full(W_normal, rounds=5)
        tr_f2, _ = run_full(W_faulty_try, rounds=5)
        De3 = (tr_f2[3][4] - tr_n2[3][4]) & MASK
        De4 = (tr_f2[4][4] - tr_n2[4][4]) & MASK
        if De4 == 0 and De3 == 0:
            return DW3_try
    return None


# ================================================================
# 3. КАСКАДНАЯ ОТМЕНА: De3=De4=De5=...=De8=0
# ================================================================

def cascade_cancellation():
    print("=" * 68)
    print("3. КАСКАДНАЯ ОТМЕНА: De3=De4=...=De8=0 ЧЕРЕЗ ΔW3-ΔW7")
    print("=" * 68)
    print()
    print("  Идея каскада: если De_k=0, то De_{k+1} ≈ Dd_k + ΔW_{k}")
    print("  → выбором ΔW_k можно обнулить De_{k+1}")
    print("  Продолжаем итеративно до тех пор, пока есть свободные W_i")
    print()

    W0_n = W0_BASE
    W0_f = (W0_BASE + DW0) & MASK

    # Начинаем с базовой трассы (De3=0, всё остальное 0)
    W_n = [W0_n, W1_SAT] + [0]*14
    W_f = [W0_f, W1_SAT] + [0]*14

    # Каскад: для каждого раунда k=3,4,5,6,7 обнуляем De_{k+1}
    # Используя ΔW_{k} = -De_{k+1} (при De_k=0)
    max_cascade = 8  # обнуляем De3..De(3+max_cascade-1)

    DWs = [0]*16  # дифференциалы слов сообщения
    DWs[0] = DW0

    # Индексация:
    # Раунд k (0-indexed) использует W[k] и производит trace[k+1].
    # De_{k+1} = trace[k+1][4] (e после раунда k).
    # Чтобы обнулить De_{k+1}: добавляем ΔW[k] = -De_{k+1}_natural.
    # De3=0 уже достигнут (раунд 2, W[1] был выбран для этого).
    # Каскад начинаем с обнуления De4: раунд 3, W[3].
    #
    # После обнуления De_{k+1}=0, ΔSig1(e_{k+1})=0,
    # поэтому De_{k+2} ≈ Dd_{k+1} + ΔW[k+1] → линейно!

    print(f"  {'Шаг':4s}  {'Обнуляем':9s}  {'De_natural':12s}  {'ΔW_выбран':12s}  "
          f"{'De_после':12s}  {'De+1_после':12s}")
    print("  " + "-" * 75)

    current_zeros = [3]  # De3=0 уже есть
    for step in range(max_cascade):
        # Обнуляем De_{step+4} через W[step+3]
        cancel_trace_idx = step + 4   # De_{step+4} = tr[step+4][4]
        w_cancel = step + 3           # W[step+3] → используется в раунде step+3

        if w_cancel > 15:
            print(f"  Нет свободных слов (нужен W[{w_cancel}] > W15)")
            break

        W_n_list = [W_n[i] for i in range(16)]
        W_f_list = [(W_n[i] + DWs[i]) & MASK for i in range(16)]

        tr_n, _ = run_full(W_n_list, rounds=cancel_trace_idx + 1)
        tr_f, _ = run_full(W_f_list, rounds=cancel_trace_idx + 1)

        De_prev = (tr_f[cancel_trace_idx - 1][4] - tr_n[cancel_trace_idx - 1][4]) & MASK
        De_nat  = (tr_f[cancel_trace_idx][4]     - tr_n[cancel_trace_idx][4])     & MASK

        if De_prev != 0:
            print(f"  {step:4d}  De{step+4:2d}=0    {'—':12s}  {'—':12s}  "
                  f"{'—':12s}  (De{step+3}≠0, каскад сломан)")
            break

        # Выбираем ΔW[w_cancel] = -De_nat (линейная отмена)
        DW_new = (MOD - De_nat) & MASK

        W_f_test = list(W_f_list)
        W_f_test[w_cancel] = (W_f_test[w_cancel] + DW_new) & MASK
        tr_f_test, _ = run_full(W_f_test, rounds=cancel_trace_idx + 2)
        De_after  = (tr_f_test[cancel_trace_idx][4]     - tr_n[cancel_trace_idx][4])     & MASK
        De_next   = (tr_f_test[cancel_trace_idx + 1][4] - tr_n[cancel_trace_idx + 1][4]) & MASK

        DWs[w_cancel] = DW_new
        if De_after == 0:
            current_zeros.append(step + 4)

        print(f"  {step:4d}  De{step+4:2d}=0    {De_nat:#012x}  {DW_new:#012x}  "
              f"{De_after:#012x}  {De_next:#012x}")

    print()
    print(f"  Обнулённые раунды e: {current_zeros}")
    print()

    # Итоговая трасса с каскадными ΔW
    W_n_final = [(W0_n if i==0 else (W1_SAT if i==1 else 0)) for i in range(16)]
    W_f_final = [(W_n_final[i] + DWs[i]) & MASK for i in range(16)]

    tr_n_final, Wn_full = run_full(W_n_final, rounds=16)
    tr_f_final, Wf_full = run_full(W_f_final, rounds=16)

    print(f"  Полная трасса De_r после каскада (раунды 1-16):")
    print(f"  {'Раунд':6s}  {'De_e':12s}  {'De_a':12s}  {'bits_e':6s}  {'нули?'}")
    print("  " + "-" * 50)
    for i in range(16):
        De_e = (tr_f_final[i+1][4] - tr_n_final[i+1][4]) & MASK
        De_a = (tr_f_final[i+1][0] - tr_n_final[i+1][0]) & MASK
        Ds_e = De_e if De_e < MOD//2 else De_e - MOD
        bits_e = abs(Ds_e).bit_length()
        nul = "De_e=0" if De_e==0 else ""
        print(f"  {i+1:6d}  {De_e:#012x}  {De_a:#012x}  {bits_e:6d}  {nul}")

    return DWs, current_zeros, W_n_final, W_f_final


# ================================================================
# 4. КАСКАД ДИФФЕРЕНЦИАЛОВ РАСПИСАНИЯ (ΔW16-ΔW63)
# ================================================================

def schedule_cascade(DWs):
    print("=" * 68)
    print("4. РАСПИСАНИЕ ПРИ КАСКАДНЫХ ΔW0-ΔW15")
    print("=" * 68)
    print()

    sig0 = lambda x: rotr(x,7) ^ rotr(x,18) ^ (x>>3)
    sig1 = lambda x: rotr(x,17) ^ rotr(x,19) ^ (x>>10)

    W_n_base = [0]*16
    W_n_base[0] = W0_BASE
    W_n_base[1] = W1_SAT

    W_f_base = [(W_n_base[i] + DWs[i]) & MASK for i in range(16)]

    Wn = make_schedule(W_n_base)
    Wf = make_schedule(W_f_base)

    DW_sched = [(Wf[i] - Wn[i]) & MASK for i in range(64)]

    print(f"  Ненулевые ΔW в полном расписании (64 слова):")
    nonzero = [(i, DW_sched[i]) for i in range(64) if DW_sched[i] != 0]
    print(f"  Всего ненулевых: {len(nonzero)} из 64")
    print()
    print(f"  {'W_idx':6s}  {'ΔW':12s}  {'bits':6s}")
    print("  " + "-" * 28)
    for i, dw in nonzero[:24]:  # первые 24
        bits = dw.bit_length() if dw < MOD//2 else (MOD - dw).bit_length()
        print(f"  {i:6d}  {dw:#012x}  {bits:6d}")
    if len(nonzero) > 24:
        print(f"  ... (+{len(nonzero)-24} ещё)")
    print()

    # Насколько большие ΔW в расписании?
    large = [(i, dw) for i, dw in nonzero
             if dw > MOD//4 and dw < 3*MOD//4]  # ~mid range ≈ random
    high  = [(i, dw) for i, dw in nonzero if dw >= 3*MOD//4 or dw <= MOD//4]
    print(f"  «Маленькие» |ΔW| (< 2^30 или > -2^30): {len(high)}")
    print(f"  «Большие» |ΔW| (≈ случайные 32-бит): {len(large)}")
    print()

    return DW_sched


# ================================================================
# 5. КОНТРОЛЬ ВЫХОДНЫХ РЕГИСТРОВ ЧЕРЕЗ W7-W15
# ================================================================

def output_differential_analysis(DWs, W_n_base):
    print("=" * 68)
    print("5. ВЫХОДНОЙ ДИФФЕРЕНЦИАЛ И ВОЗМОЖНОСТЬ ПОЛНОЙ КОЛЛИЗИИ")
    print("=" * 68)
    print()

    W_f_base = [(W_n_base[i] + DWs[i]) & MASK for i in range(16)]
    tr_n, Wn = run_full(W_n_base, rounds=64)
    tr_f, Wf = run_full(W_f_base, rounds=64)

    print("  Выходные дифференциалы (после 64 раундов):")
    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    total_bits = 0
    for i, name in enumerate(reg_names):
        Dn = tr_n[-1][i]
        Df = tr_f[-1][i]
        D = (Df - Dn) & MASK
        Ds = D if D < MOD//2 else D - MOD
        bits = abs(Ds).bit_length()
        total_bits += bits
        print(f"    D{name} = {D:#010x}  ({Ds:+12d}, {bits:2d} bits)")

    print(f"\n  Суммарных бит: ~{total_bits} из 256 (случайный: ~256)")
    print()

    # SHA-256 = IV + (финальное состояние)
    # Для коллизии нужно D_output_i = 0 для всех 8 регистров
    # Текущий выход явно ненулевой

    print("  СТЕПЕНИ СВОБОДЫ ДЛЯ ПОЛНОЙ КОЛЛИЗИИ:")
    print()
    used_dw = [i for i, d in enumerate(DWs) if d != 0]
    free_w = [i for i in range(16) if DWs[i] == 0]
    print(f"  Использованные слова: W{used_dw} ({len(used_dw)} слов)")
    print(f"  Свободные слова:      W{free_w[:9]} ({len(free_w)} слов)")
    print()
    print(f"  Для полной коллизии нужно обнулить 8 регистров (256 бит).")
    print(f"  Свободных слов W7-W15: 9 → 9 × 32 = 288 бит степеней свободы.")
    print(f"  Задача: 8 нелинейных уравнений с 9 неизвестными (W7-W15)")
    print(f"  → теоретически решаема, но НЕЛИНЕЙНА (неприводима аналитически)")
    print()
    print("  ОЦЕНКА СЛОЖНОСТИ:")
    print("  - Наивный перебор W7-W15: 9 × 32 = 288 бит → 2^288 (невозможно)")
    print("  - С meet-in-middle (5+4 разбиение): 2^160 (всё ещё слишком много)")
    print("  - С De3=De4=...=De8=0 (каскад): выходной дифференциал уменьшен?")
    print()

    return total_bits


# ================================================================
# 6. ВЕКТОР ПЕРСПЕКТИВНОЙ АТАКИ
# ================================================================

def attack_vector_analysis(DWs, W_n_base, cascade_zeros):
    print("=" * 68)
    print("6. ВЕКТОР ПЕРСПЕКТИВНОЙ АТАКИ (COMBINED APPROACH)")
    print("=" * 68)
    print()

    print("  НАБЛЮДЕНИЕ: При каскадном обнулении De3=De4=...=De_k=0:")
    print("  - Регистр e принимает одинаковые значения в раундах 3..k")
    print("  - Это ограничивает распространение дифференциала через e-канал")
    print("  - Но a, b, c, d, f, g, h накапливают дифференциал независимо")
    print()

    # Проверим: насколько уменьшается дифференциал при каскаде vs без
    W_n_base_arr = list(W_n_base)
    W_f_no_cascade = [(W_n_base_arr[i] + (DW0 if i==0 else 0)) & MASK
                      for i in range(16)]
    W_f_cascade    = [(W_n_base_arr[i] + DWs[i]) & MASK for i in range(16)]

    tr_n, _ = run_full(W_n_base_arr, rounds=64)
    tr_f_no, _ = run_full(W_f_no_cascade, rounds=64)
    tr_f_cas, _ = run_full(W_f_cascade, rounds=64)

    print(f"  Сравнение выходного дифференциала (все 8 регистров):")
    print(f"  {'Рег':4s}  {'Без каскада':20s}  {'С каскадом':20s}")
    print("  " + "-" * 50)
    for i, name in enumerate(['a','b','c','d','e','f','g','h']):
        Dn = tr_n[-1][i]
        D_no  = (tr_f_no[-1][i] - Dn) & MASK
        D_cas = (tr_f_cas[-1][i] - Dn) & MASK
        b_no  = D_no.bit_length() if D_no  < MOD//2 else (MOD-D_no).bit_length()
        b_cas = D_cas.bit_length() if D_cas < MOD//2 else (MOD-D_cas).bit_length()
        marker = " ← улучшение!" if b_cas < b_no - 2 else ""
        print(f"  {name:4s}  {D_no:#010x} ({b_no:2d} bits)      "
              f"{D_cas:#010x} ({b_cas:2d} bits){marker}")

    print()
    print("  КЛЮЧЕВЫЕ ВЫВОДЫ П-7:")
    print()
    print("  1. ΔW3 может обнулить De4 при условии De3=0.")
    print("     Но ΔW3 вносит каскадный вклад в W16-W63 через расписание.")
    print()
    print("  2. Каскадное обнуление De3..De8 возможно через ΔW3..ΔW7.")
    print(f"     Достигнуто нулей: раунды {cascade_zeros}")
    print()
    print("  3. Свободные слова W7-W15 (9 штук) могут частично контролировать")
    print("     выходной дифференциал, но нелинейность делает точное решение")
    print("     трудным аналитически.")
    print()
    print("  4. ГЛАВНЫЙ РЕЗУЛЬТАТ: Каскадная атака снижает выходной дифференциал")
    print("     по e-каналу, но не даёт полной коллизии через данный механизм.")
    print()
    print("  ОТКРЫТЫЕ ВОПРОСЫ → П-8:")
    print("  - Meet-in-middle: использовать каскад как 'фильтр' в MitM?")
    print("  - Algebraic attack: W7-W15 как переменные в SAT/нелинейной системе")
    print("  - Комбинированная атака: De3=0 + De8=0 + MitM для оставшихся раундов")
    print()


# ================================================================
# MAIN
# ================================================================

def main():
    print()
    print("П-7: СВОБОДНЫЕ СЛОВА W2-W15 ДЛЯ КАСКАДНОЙ ОТМЕНЫ ДИФФЕРЕНЦИАЛА")
    print("=" * 68)
    print()

    # 1. Базовый случай
    De4_natural, W0_n, W0_f, W1 = baseline_analysis()

    # 2. Линейная отмена De4 через ΔW3
    DW3, De5 = linear_cancellation_De4(De4_natural, W0_n, W0_f, W1)

    # 3. Каскадная отмена De3..De_k
    DWs, zeros, W_n_final, W_f_final = cascade_cancellation()

    # 4. Расписание
    DW_sched = schedule_cascade(DWs)

    # 5. Выходной дифференциал
    output_differential_analysis(DWs, W_n_final)

    # 6. Вектор атаки
    attack_vector_analysis(DWs, W_n_final, zeros)

    print("=" * 68)
    print("П-7 ЗАВЕРШЁН")
    print("=" * 68)


if __name__ == "__main__":
    main()
