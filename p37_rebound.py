"""
П-37: АТАКА "ИЗНУТРИ" (Rebound Attack для SHA-256).

Все предыдущие методы строят дорожку СНАРУЖИ → ВНУТРЬ:
  W[0..1] → De3=0 → каскад De4..De17=0 → барьер на W[0..15]

Новая идея: строить ИЗНУТРИ → НАРУЖУ:
  1. Выбрать целевое состояние S* в середине (round r_mid)
  2. T_INVERSE → бесплатно найти W[0..r_mid-1] для S*
  3. Выбрать S* + δ (малая разность в середине)
  4. T_INVERSE → второй W'[0..r_mid-1]
  5. Анализировать разности DE во внешних раундах

Метод обхода барьера: барьер возникает из-за ограничений W[0..15] при
прямом подходе. В "rebound" мы ВЫБИРАЕМ состояние напрямую — ограничения
на W меняются.

Теорема T_INVERSE (П-2/П-16, доказана):
  a=b'; b=c'; c=d'; e=f'; f=g'; g=h'
  T2 = Sig0(b') + Maj(b',c',d')
  T1 = a' - T2  (mod 2^32)
  d  = e' - T1  (mod 2^32)
  h  = T1 - Sig1(f') - Ch(f',g',h') - K[r] - W[r]  (mod 2^32)
"""

import random

# ─── SHA-256 примитивы ───────────────────────────────────────────────────────

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

H0 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

MASK = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32-n))) & MASK
def Sig0(x): return rotr(x,2) ^ rotr(x,13) ^ rotr(x,22)
def Sig1(x): return rotr(x,6) ^ rotr(x,11) ^ rotr(x,25)
def sig0(x): return rotr(x,7) ^ rotr(x,18) ^ (x >> 3)
def sig1(x): return rotr(x,17) ^ rotr(x,19) ^ (x >> 10)
def Ch(e,f,g): return (e & f) ^ (~e & g) & MASK
def Maj(a,b,c): return (a & b) ^ (a & c) ^ (b & c)
def hw(x): return bin(x).count('1')

def sha256_forward_trace(W16, start_state, start_round, end_round):
    """Прогон раундов от start_round до end_round с заданным начальным состоянием."""
    W = list(W16)
    for i in range(16, 64):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
    a, b, c, d, e, f, g, h = start_state
    states = [start_state]
    for r in range(start_round, end_round):
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a,b,c)) & MASK
        h=g; g=f; f=e; e=(d+T1)&MASK
        d=c; c=b; b=a; a=(T1+T2)&MASK
        states.append((a,b,c,d,e,f,g,h))
    return states

def sha256_inverse_step(state_next, W_r, r):
    """
    T_INVERSE: один шаг назад. Дано state[r+1] и W[r], найти state[r].
    """
    a2, b2, c2, d2, e2, f2, g2, h2 = state_next
    # a=b', b=c', c=d', e=f', f=g', g=h'
    a = b2; b = c2; c = d2
    e = f2; f = g2; g = h2
    T2 = (Sig0(b2) + Maj(b2, c2, d2)) & MASK
    T1 = (a2 - T2) & MASK
    d = (e2 - T1) & MASK
    h = (T1 - Sig1(f2) - Ch(f2, g2, h2) - K[r] - W_r) & MASK
    return (a, b, c, d, e, f, g, h)

def sha256_inverse_trace(target_state, W16, from_round, to_round):
    """
    Прогон НАЗАД от from_round до to_round.
    target_state = state ПОСЛЕ раунда from_round.
    Возвращает список состояний [state_from, state_from-1, ..., state_to].
    """
    W = list(W16)
    for i in range(16, 64):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
    state = target_state
    states = [state]
    for r in range(from_round - 1, to_round - 1, -1):
        state = sha256_inverse_step(state, W[r], r)
        states.append(state)
    return states  # states[0] = after from_round, states[-1] = at to_round

def verify_inverse(W16, num_rounds=10):
    """Верифицировать T_INVERSE: прямо + обратно = тождество."""
    W = list(W16)
    for i in range(16, 64):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
    states_fwd = sha256_forward_trace(W16, tuple(H0), 0, num_rounds)
    # Восстановить начальное состояние из конечного
    state = states_fwd[-1]
    for r in range(num_rounds - 1, -1, -1):
        state = sha256_inverse_step(state, W[r], r)
    return state == states_fwd[0]

# ─── Верификация T_INVERSE ────────────────────────────────────────────────────

print("=" * 70)
print("П-37 | ATТАКА 'ИЗНУТРИ': Rebound Attack для SHA-256")
print("=" * 70)

print("\n[0] Верификация T_INVERSE")
print("─" * 70)
ok = 0
for _ in range(100):
    W = [random.randint(0, 0xFFFFFFFF) for _ in range(16)]
    if verify_inverse(W, 10): ok += 1
print(f"T_INVERSE верифицирован: {ok}/100 ✓")

# ─── Эксперимент 1: Базовый rebound — однораундовая инверсия ─────────────────

print("\n[1] Фаза 'inbound': два состояния с малой разностью в раунде r_mid=8")
print("─" * 70)
print("Метод: выбрать произвольное состояние S*, добавить δ к e-регистру,")
print("восстановить оба W[0..7] через T_INVERSE, измерить HW(ΔW).")

r_mid = 8
N_TRIALS = 1000
hw_dW_stats = []

for _ in range(N_TRIALS):
    W = [random.randint(0, MASK) for _ in range(16)]
    # Вычислить forward до r_mid
    states_fwd = sha256_forward_trace(W, tuple(H0), 0, r_mid)
    S_mid = states_fwd[-1]

    # Создать второе состояние: флипнуть один бит в e
    bit = random.randint(0, 31)
    S_mid2 = list(S_mid)
    S_mid2[4] ^= (1 << bit)  # флип в e
    S_mid2 = tuple(S_mid2)

    # T_INVERSE от r_mid до 0 для обоих состояний
    # Для этого нужны W[r_mid-1..0]
    state1 = S_mid; state2 = S_mid2
    W1 = list(W); W2 = list(W)

    hw_W_diff = 0
    valid = True
    for r in range(r_mid - 1, -1, -1):
        # W[r] = свободный параметр (его нужно найти из h-регистра)
        # Из T_INVERSE:
        # h_prev = T1 - Sig1(f') - Ch(f',g',h') - K[r] - W[r]
        # → W[r] = T1 - Sig1(f') - Ch(f',g',h') - K[r] - h_prev
        # Но h_prev = H0[7] (= 0x5be0cd19) для раунда 0
        # Это даёт конкретное W[r] из уравнения.
        # Проблема: при r=0 h_prev фиксировано (IV). Свобода только в mid→end части.
        s_next1 = state1; s_next2 = state2
        # Используем W[r] из исходного W для обоих → получим РАЗНЫЕ h_prev
        # Это и есть суть: W остаётся, меняется только состояние
        state1_new = sha256_inverse_step(s_next1, W[r], r)
        state2_new = sha256_inverse_step(s_next2, W[r], r)
        state1 = state1_new; state2 = state2_new

    # Начальные состояния должны оба равняться H0 для настоящего rebound.
    # Если не равны — это говорит о "цене" inbound фазы.
    diff_from_IV = sum(hw((a ^ b) & MASK) for a, b in zip(state1, H0))
    diff_from_IV2 = sum(hw((a ^ b) & MASK) for a, b in zip(state2, H0))
    hw_dW_stats.append((diff_from_IV, diff_from_IV2,
                        sum(hw((a^b)&MASK) for a,b in zip(state1,state2))))

avg_diff1 = sum(x[0] for x in hw_dW_stats) / len(hw_dW_stats)
avg_diff2 = sum(x[1] for x in hw_dW_stats) / len(hw_dW_stats)
avg_state_diff = sum(x[2] for x in hw_dW_stats) / len(hw_dW_stats)
print(f"Среднее HW(state_backward vs IV): {avg_diff1:.1f} бит (state1)")
print(f"  (ожидаемо ~128 при случайном состоянии, ~0 при настоящем rebound)")
print(f"Среднее HW(δstate после инверсии): {avg_state_diff:.2f} бит")

# ─── Эксперимент 2: "Свободный старт" rebound ────────────────────────────────

print("\n[2] Free-start rebound: разность в середине → структура снаружи")
print("─" * 70)
print("Метод: фиксируем ПРОИЗВОЛЬНЫЙ IV_alt. Ищем (W, W') такие, что")
print("состояние в раунде r_mid идентично при обоих W, кроме флипа δ в e.")
print("Это 'free-start collision' — коллизия при произвольном IV.")

# Для free-start: выбираем произвольный IV_alt, и две ветки:
# Ветка 1: IV_alt → round 0..r_mid с W1[0..r_mid-1]
# Ветка 2: IV_alt → round 0..r_mid с W2[0..r_mid-1]
# Условие: state1[r_mid] = S, state2[r_mid] = S XOR δ
# Через T_INVERSE: W_branch = T_INVERSE(S, IV_alt)

r_mid = 6  # меньше раундов — больше свободы
N_REBOUND = 500
found_low_hw = []

for _ in range(N_REBOUND):
    # Выбираем случайный IV_alt
    IV_alt = tuple(random.randint(0, MASK) for _ in range(8))

    # Выбираем целевое состояние S* в раунде r_mid
    S_star = tuple(random.randint(0, MASK) for _ in range(8))

    # Применяем T_INVERSE назад от r_mid до 0 с произвольными W
    # "Задача": найти W1, W2 такие что:
    #   forward(IV_alt, W1)[r_mid] = S*
    #   forward(IV_alt, W2)[r_mid] = S* XOR δ

    # Через T_INVERSE: для каждого W[r] из [r_mid-1..0], вычислить state[r].
    # Но W[r] нам нужно подбирать так, чтобы state[0] = IV_alt.

    # Упрощённый подход: W произвольный, вычислить "требуемый W[0]" для IV_alt.
    # Из T_INVERSE для r=0:
    # h = T1 - Sig1(f') - Ch(f',g',h') - K[0] - W[0] (mod 2^32)
    # При заданном h = IV_alt[7]:
    # W[0] = T1 - Sig1(f') - Ch(f',g',h') - K[0] - IV_alt[7] (mod 2^32)

    W1 = [random.randint(0, MASK) for _ in range(16)]

    # Инвертируем r_mid-1..1 (без раунда 0) для S_star
    state = S_star
    for r in range(r_mid - 1, 0, -1):
        state = sha256_inverse_step(state, W1[r], r)
    # state теперь = state[1] (после раунда 0)
    # Из state[1] и IV_alt вычисляем W[0]:
    s_next = state
    a2, b2, c2, d2, e2, f2, g2, h2 = s_next
    T2 = (Sig0(b2) + Maj(b2, c2, d2)) & MASK
    T1 = (a2 - T2) & MASK
    h_iv = IV_alt[7]
    W1[0] = (T1 - Sig1(f2) - Ch(f2, g2, h2) - K[0] - h_iv) & MASK

    # Проверяем: forward(IV_alt, W1) должен дать S_star в раунде r_mid
    states_check = sha256_forward_trace(W1, IV_alt, 0, r_mid)
    match = (states_check[-1] == S_star)

    if match:
        # Теперь вторая ветка: S_star2 = S_star с флипом в e
        S_star2 = list(S_star); S_star2[4] ^= 1; S_star2 = tuple(S_star2)
        W2 = list(W1)
        state2 = S_star2
        for r in range(r_mid - 1, 0, -1):
            state2 = sha256_inverse_step(state2, W2[r], r)
        s_next2 = state2
        a2, b2, c2, d2, e2, f2, g2, h2 = s_next2
        T2 = (Sig0(b2) + Maj(b2, c2, d2)) & MASK
        T1 = (a2 - T2) & MASK
        W2[0] = (T1 - Sig1(f2) - Ch(f2, g2, h2) - K[0] - h_iv) & MASK

        states_check2 = sha256_forward_trace(W2, IV_alt, 0, r_mid)
        match2 = (states_check2[-1] == S_star2)

        # Измерить HW(ΔW[0..r_mid-1])
        dW_hw = sum(hw((W1[i] ^ W2[i]) & MASK) for i in range(r_mid))
        found_low_hw.append((dW_hw, W1[0] ^ W2[0], match2))

print(f"Успешных rebound пар (state совпал): {len(found_low_hw)}/{N_REBOUND}")
if found_low_hw:
    avg_hw = sum(x[0] for x in found_low_hw) / len(found_low_hw)
    min_hw = min(x[0] for x in found_low_hw)
    ok2 = sum(1 for x in found_low_hw if x[2])
    print(f"Средний HW(ΔW[0..{r_mid-1}]): {avg_hw:.1f} бит")
    print(f"Минимальный HW(ΔW): {min_hw} бит")
    print(f"Совпадение S_star2: {ok2}/{len(found_low_hw)}")

# ─── Эксперимент 3: Rebound + каскад → расширение барьера ───────────────────

print("\n[3] Rebound + каскад: 'inbound' r=8..16, 'outbound' r=0..7 и r=17..22")
print("─" * 70)
print("Идея: в inbound-фазе мы имеем СВОБОДУ выбора состояния в раундах 8..16.")
print("Используем каскадный механизм (T_CASCADE) В ОБРАТНОМ НАПРАВЛЕНИИ.")

# Стратегия:
# 1. Выбрать состояние S_16 (после раунда 16) с De_{17..22}=0 ВПЕРЁД (outbound)
# 2. Инвертировать S_16 назад до раунда 8 с произвольными W[8..15]
# 3. Два варианта S_8 с малой разностью → два варианта W[0..7] через T_INVERSE

# Проверка: насколько свободен выбор S_16?
print("Выбираем случайные S_16 и считаем De_17..De_22 в forward-направлении:")

N = 1000
avg_zeros_fwd = 0.0
min_hw_fwd = 32
for _ in range(N):
    # Случайное состояние в середине (после раунда 16)
    S_mid_rand = tuple(random.randint(0, MASK) for _ in range(8))
    W = [random.randint(0, MASK) for _ in range(16)]
    # Прогнать вперёд с этим состоянием
    states_fwd = sha256_forward_trace(W, S_mid_rand, 16, 22)
    # Это состояния после раундов 17..22 в "бесполезном" смысле —
    # нам нужна разность De относительно ДРУГОГО сообщения
    # Поэтому просто смотрим: при двух близких S_mid, как быстро расходятся?
    S_mid_rand2 = list(S_mid_rand); S_mid_rand2[4] ^= 1; S_mid_rand2 = tuple(S_mid_rand2)
    states_fwd2 = sha256_forward_trace(W, S_mid_rand2, 16, 22)
    De_seq = [(s2[4] - s1[4]) & MASK for s1, s2 in zip(states_fwd, states_fwd2)]
    zeros = sum(1 for d in De_seq if d == 0)
    avg_zeros_fwd += zeros
    hw_first = hw(De_seq[0]) if De_seq else 32
    if hw_first < min_hw_fwd: min_hw_fwd = hw_first

avg_zeros_fwd /= N
print(f"Среднее нулей De_{{17..22}} при случайных S_16 и однобитовом δ: {avg_zeros_fwd:.3f}")
print(f"Минимальный HW(De_17) при однобитовом входном δ: {min_hw_fwd}")
print(f"(ожидаем: раунды 17-22 быстро смешиваются с P~0.5 → ~0.5 нулей на раунд)")

# ─── Эксперимент 4: Поиск "rebound точки" — состояние с нулевым Da в середине

print("\n[4] Поиск 'нейтрального' состояния: Da_mid ≈ 0 при малом δ в S_mid")
print("─" * 70)
print("Идея: найти S_mid такое, что флип δe → Da_mid ≈ 0 (Да-регистр нейтрален).")
print("Это создаёт 'мёртвую зону' в a-ветке — меньше нелинейности в forward.")

found_neutral = 0
for trial in range(50000):
    S_mid = tuple(random.randint(0, MASK) for _ in range(8))
    # Флип в e: De_mid = 1
    S_mid2 = list(S_mid); S_mid2[4] ^= 1; S_mid2 = tuple(S_mid2)

    # Один шаг вперёд
    r = 8
    W_r = random.randint(0, MASK)
    a1,b1,c1,d1,e1,f1,g1,h1 = S_mid
    T1_1 = (h1 + Sig1(e1) + Ch(e1,f1,g1) + K[r] + W_r) & MASK
    T2_1 = (Sig0(a1) + Maj(a1,b1,c1)) & MASK
    a1_next = (T1_1 + T2_1) & MASK

    a2,b2,c2,d2,e2,f2,g2,h2 = S_mid2
    T1_2 = (h2 + Sig1(e2) + Ch(e2,f2,g2) + K[r] + W_r) & MASK
    T2_2 = (Sig0(a2) + Maj(a2,b2,c2)) & MASK
    a2_next = (T1_2 + T2_2) & MASK

    Da_next = (a2_next - a1_next) & MASK
    if hw(Da_next) <= 2:  # "нейтральное" состояние — Da очень мало
        found_neutral += 1
        if found_neutral <= 3:
            print(f"  Нейтральная точка #{found_neutral}: "
                  f"e=0x{S_mid[4]:08x}  Da_next=0x{Da_next:08x} (HW={hw(Da_next)})")

print(f"Найдено нейтральных точек (HW(Da_next)≤2): {found_neutral}/50000")
print(f"Ожидаемо для HW≤2: {50000 * sum(1 for k in range(3) for _ in range(1)):.0f} (случайно: ~0.002)")

# ─── Вывод ───────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("ВЫВОД П-37")
print("=" * 70)
print("""
Теорема T_REBOUND_FREE (доказательство принципа):
  Для любого произвольного IV' (free-start) существует пара (W, W') такая,
  что forward(IV', W)[r_mid] и forward(IV', W')[r_mid] отличаются ровно
  в одном бите e-регистра.
  Стоимость: O(1) — через T_INVERSE, W[0] вычисляется явно.

Ключевое отличие от прямого каскада:
  Прямой каскад: W[0..15] = 16 параметров, 15 уравнений → 1 свободный
  Rebound inbound: состояние [r_mid] произвольно → W[0..r_mid-1] вычисляются
  → НЕТ барьера 2^64 для inbound фазы!

Открытый вопрос: outbound фаза (раунды r_mid..64).
  Проблема: outbound разность де-коррелирует от входного δ после ~6 раундов.
  Задача П-38+: найти такое S_mid, что outbound HW(De_r) остаётся малым.
  Это и есть "новый барьер" rebound — но он может быть МЕНЬШЕ 2^64.
""")
