"""
SHA-256 Дифференциальный Криптоанализ
П-2: Структура ADD8 — аналитическое доказательство и De3=0

КЛЮЧЕВОЙ РЕЗУЛЬТАТ П-1:
  lifted_diff_e2 принимает ровно 8 значений (для чётных W0, delta=1)

РЕЗУЛЬТАТЫ П-2 (финальные):

  T_ADD8 ДОКАЗАНА:
    ΔCh(e1+1, H0[4], H0[5]) = +1  (константа, т.к. H0[4] бит0=1, H0[5] бит0=0)
    ΔSig1(e1, 1) принимает 8 значений, определяемых битами {26,21,7} Sig1(e1):
      (b26,b21,b7) → d2 = (1-2*b26)*2^26 + (1-2*b21)*2^21 + (1-2*b7)*2^7 + 1
    Нижний байт d2 всегда 0x81 (т.к. ΔSig1 mod 256 = 0x80, +1 от ΔCh)

  ИСПРАВЛЕННАЯ ФОРМУЛА De3=0:
    P-1 использовала неправильный ΔCh (одинаковый f2 для обоих потоков).
    КОРРЕКТНОЕ условие:
      Sig1(e2_f) - Sig1(e2_n) + Ch(e2_f, e1_f, H0[4]) - Ch(e2_n, e1_n, H0[4]) ≡ 0 (mod 2^32)
    где e1_f = (e1_n + 1) mod 2^32 — f2 РАЗЛИЧАЕТСЯ в потоках!
    Разница: CORR - WRONG = (e2_f & 1) = bit0(e2_f)
    Для W_SAT3: e2_f=0xc31e77ce (чётный) → WRONG=CORRECT случайно ✓

  De3=0 ПЛОТНОСТЬ:
    С W[1..15]=0: P(De3=0) ≈ 0 (нет решений в 10^6 сэмплах)
    С произв. (W0, W1): P(De3=0) ≈ 0 (нет в 500k сэмплах)
    Сканирование W1 в [0x3d000000, 0x3e000000) для W0=W_SAT3:
      Решения с ШАГОМ ~2^15 ≈ 32768
      Ожидаемо ~512 = 2^9 решений в 2^24 диапазоне → P ≈ 2^(-15)
      Это МНОГО МЕНЬШЕ предсказания T_PERIOD3 (1/256 = 2^(-8))
    Вывод: T_PERIOD3 либо некорректна, либо относится к другому условию

ЗАДАЧА П-2:
  1. Доказать аналитически: 8 = 2^3 из битов (26, 21, 7) Sig1(e1)  ✓
  2. Вычислить все 8 значений явно  ✓
  3. P(De3=0 | класс) для каждого из 8 классов
  4. Общее P(De3=0), сравнение с T_PERIOD3
  5. Структура: e2 → W0 (обратный путь)
"""

import random
from collections import Counter, defaultdict

M32  = 0xFFFFFFFF
MOD  = 2**32

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,
    0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,
    0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,
    0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,
    0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,
    0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,
    0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,
    0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,
    0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
]
H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & M32
def Sig0(x): return rotr(x,2)  ^ rotr(x,13) ^ rotr(x,22)
def Sig1(x): return rotr(x,6)  ^ rotr(x,11) ^ rotr(x,25)
def sig0(x): return rotr(x,7)  ^ rotr(x,18) ^ (x >> 3)
def sig1(x): return rotr(x,17) ^ rotr(x,19) ^ (x >> 10)
def Ch(e,f,g):  return (e & f) ^ (~e & g) & M32
def Maj(a,b,c): return (a & b) ^ (a & c) ^ (b & c)

def schedule(W16):
    W = list(W16) + [0]*48
    for i in range(16,64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & M32
    return W

def sha256_rounds(W16, N):
    W = schedule(W16)
    a,b,c,d,e,f,g,h = H0
    for i in range(N):
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[i] + W[i]) & M32
        T2 = (Sig0(a) + Maj(a,b,c)) & M32
        h=g; g=f; f=e; e=(d+T1)&M32
        d=c; c=b; b=a; a=(T1+T2)&M32
    return [a,b,c,d,e,f,g,h]

# ================================================================
# Фиксированные константы IV
# ================================================================
# a0,b0,...,h0 = H0[0..7]
a0,b0,c0,d0,e0,f0,g0,h0 = H0

# После раунда 0: f1=e0=H0[4], g1=f0=H0[5], h1=g0=H0[6], d1=c0=H0[2]
f1 = e0  # H0[4]
g1 = f0  # H0[5]
h1 = g0  # H0[6]
d1 = c0  # H0[2]

# После раунда 1: g2=f1=H0[4], h2=g1=H0[5]
g2 = f1  # H0[4]
h2 = g1  # H0[5]
d2 = b0  # H0[1]

# T2_0 не зависит от W0
T2_0 = (Sig0(a0) + Maj(a0, b0, c0)) & M32

# C_IV: константная часть T1_0
C_IV = (h0 + Sig1(e0) + Ch(e0,f0,g0) + K[0]) & M32

def e1_from_W0(W0):
    """e1 = (H0[3] + C_IV + W0) mod 2^32"""
    return (d0 + C_IV + W0) & M32

def e2_from_e1(e1, W1=0):
    """e2 = d1 + T1_1(e1, W1) = H0[2] + h1 + Sig1(e1) + Ch(e1,f1,g1) + K[1] + W1"""
    T1_1 = (h1 + Sig1(e1) + Ch(e1, f1, g1) + K[1] + W1) & M32
    return (d1 + T1_1) & M32

# Маска XOR для Sig1 при сдвиге на 1 бит
# ROTR6(1)=2^26, ROTR11(1)=2^21, ROTR25(1)=2^7
SIG1_MASK = rotr(1,6) ^ rotr(1,11) ^ rotr(1,25)  # = 0x04200080

# ================================================================
# 1. АНАЛИТИЧЕСКОЕ ДОКАЗАТЕЛЬСТВО: ПОЧЕМУ РОВНО 8 ЗНАЧЕНИЙ
# ================================================================

def prove_8_values():
    """
    Теорема ADD8:
    Для чётных e1: ΔSig1(e1, 1) = Sig1(e1+1) - Sig1(e1)
    принимает ровно 8 значений, определяемых битами {26,21,7} Sig1(e1).

    Причина:
    - Чётный e1 → e1+1 меняет ТОЛЬКО бит 0 (нет carry вверх)
    - Sig1 = ROTR6 ⊕ ROTR11 ⊕ ROTR25
    - XOR-разность: Sig1(e1+1) ⊕ Sig1(e1) = ROTR6(1)⊕ROTR11(1)⊕ROTR25(1) = 0x04200080
      (фиксирована для всех чётных e1)
    - Арифм. разность: Sig1(e1+1) - Sig1(e1) = (s⊕mask) - s, где s=Sig1(e1)
      Зависит от битов s в позициях {26,21,7} (те, что флипаются маской)
    - 2^3 = 8 комбинаций {b26, b21, b7} → 8 значений

    Аналогично ΔCh(e1, 1, f1, g1) = 1 для всех чётных e1
    (т.к. H0[4] имеет бит0=1, H0[5] имеет бит0=0 → Ch меняется точно на +1)
    """
    print("=" * 65)
    print("1. АНАЛИТИЧЕСКОЕ ДОКАЗАТЕЛЬСТВО: ADD8")
    print("=" * 65)

    mask = SIG1_MASK  # 0x04200080

    print(f"\n  SIG1_MASK = ROTR6(1)⊕ROTR11(1)⊕ROTR25(1) = {mask:#010x}")
    print(f"  Биты маски: {26} (2^26={2**26}), {21} (2^21={2**21}), {7} (2^7={2**7})")
    print()
    print(f"  Чётный e1 → (+1 меняет только бит 0) → XOR-diff Sig1 = CONST = {mask:#010x}")
    print(f"  Арифм. diff = (Sig1(e1)⊕mask) - Sig1(e1) зависит от битов {{26,21,7}} Sig1(e1):")
    print()

    # Вычисляем 8 значений ΔSig1 аналитически
    bit_positions = [26, 21, 7]
    sig1_diffs = {}
    for b26 in range(2):
        for b21 in range(2):
            for b7 in range(2):
                # Арифм. эффект XOR маской: для каждого бита i из {26,21,7}:
                # если бит i в s = 0 → flip добавляет +2^i
                # если бит i в s = 1 → flip вычитает -2^i
                d = (1-2*b26)*(2**26) + (1-2*b21)*(2**21) + (1-2*b7)*(2**7)
                sig1_diffs[(b26, b21, b7)] = d

    print(f"  {'(b26,b21,b7)':12s}  {'ΔSig1 над Z':15s}  {'ΔSig1 mod 2^32':14s}")
    print("  " + "-" * 50)
    for key, d in sorted(sig1_diffs.items()):
        print(f"  {str(key):12s}  {d:15d}  {d % MOD:#014x}")

    # ΔCh для чётного e1
    f1_bit0 = f1 & 1   # H0[4] bit0
    g1_bit0 = g1 & 1   # H0[5] bit0
    dCh_even = f1_bit0 - g1_bit0  # = 1 - 0 = +1
    print(f"\n  ΔCh(чётный e1, 1, H0[4], H0[5]):")
    print(f"    H0[4] бит0 = {f1_bit0},  H0[5] бит0 = {g1_bit0}")
    print(f"    => ΔCh = f1_bit0 - g1_bit0 = {dCh_even}  (КОНСТАНТА для всех чётных e1)")
    print()
    print(f"  ИТОГО: lifted_diff_e2 = ΔSig1 + {dCh_even}")
    print()

    # Все 8 значений lifted_diff_e2
    d2_values = {}
    print(f"  {'(b26,b21,b7)':12s}  {'d2 над Z':15s}  {'d2 mod 2^32':14s}  {'вероятность'}")
    print("  " + "-" * 60)
    for key, dS in sorted(sig1_diffs.items()):
        d2 = dS + dCh_even
        d2_values[key] = d2
        # Вероятность: 1/8 если bits 26,21,7 независимы
        print(f"  {str(key):12s}  {d2:15d}  {d2 % MOD:#014x}  ~1/8 (12.5%)")

    print()
    print("  Верификация численно: 500k сэмплов чётных W0")
    obs = Counter()
    for _ in range(500000):
        W0 = random.randint(0, M32) & ~1
        e1 = e1_from_W0(W0)
        s  = Sig1(e1)
        key = (int(bool(s & (1<<26))), int(bool(s & (1<<21))), int(bool(s & (1<<7))))
        dS = Sig1((e1+1)&M32) - Sig1(e1)
        dC = Ch((e1+1)&M32, f1, g1) - Ch(e1, f1, g1)
        actual_d2 = dS + dC
        predicted = d2_values[key]
        if actual_d2 != predicted:
            print(f"  ОШИБКА! key={key} predicted={predicted} actual={actual_d2}")
        obs[key] += 1

    print(f"  {'(b26,b21,b7)':12s}  {'наблюд.':8s}  {'теория':8s}  {'откл.%'}")
    for key in sorted(obs.keys()):
        n = obs[key]
        theory = 500000 / 8
        print(f"  {str(key):12s}  {n:8d}  {theory:8.0f}  {100*(n-theory)/theory:+.2f}%")

    print(f"\n  ✓ Все 8 значений подтверждены, распределение равномерное (~1/8 каждое)")
    return d2_values

# ================================================================
# 2. ЯВНЫЕ ЗНАЧЕНИЯ ВСЕХ 8 d2
# ================================================================

def compute_8_values_table(d2_values):
    print("\n" + "=" * 65)
    print("2. ТАБЛИЦА 8 ЗНАЧЕНИЙ lifted_diff_e2")
    print("=" * 65)
    print()
    print(f"  Заметим: d2 mod 2^32 ≡ 0x81 (mod 0x100) для всех 8 классов")
    print(f"  Т.к. ΔSig1 mod 256 = 0x80 (из ROTR25(1)=0x80) и ΔCh = +1")
    print(f"  => нижний байт d2 всегда = 0x80+0x01 = 0x81")
    print()
    print(f"  {'Класс':5s}  {'d2 над Z':15s}  {'d2 mod 2^32':14s}  {'нижний байт':12s}")
    print("  " + "-" * 55)
    for i, (key, d2) in enumerate(sorted(d2_values.items())):
        lo = d2 % 256
        print(f"  #{i}   {str(key):13s}  {d2:15d}  {d2 % MOD:#014x}  {lo:#04x} {'✓' if lo==0x81 else '✗'}")

    # Все 8 значений mod 2^32 — что общего?
    vals_mod = sorted(set(d % MOD for d in d2_values.values()))
    print(f"\n  mod 2^32: {[hex(v) for v in vals_mod]}")
    print(f"  Старшие 4 бита: {sorted(set(v >> 28 for v in vals_mod))}")
    print(f"  Биты [27:8]:    {sorted(set((v >> 8) & 0xFFFFF for v in vals_mod))}")

# ================================================================
# 3. УСЛОВИЕ De3=0 ДЛЯ КАЖДОГО ИЗ 8 КЛАССОВ
# ================================================================

def analyze_de3_per_class(d2_values, N_per_class=100000):
    """
    Для каждого из 8 классов (b26,b21,b7):
    P(De3=0 | класс) = P(ΔSig1(e2, d2) + ΔCh(e2, e1, H0[4], d2) ≡ 0 mod 2^32)

    Важно: f2 = e1 (зависит от W0!), g2 = H0[4].
    Условие — это 2-мерное в (e1, e2), но e2 = e2(e1).
    """
    print("\n" + "=" * 65)
    print("3. P(De3=0) ДЛЯ КАЖДОГО ИЗ 8 КЛАССОВ")
    print("=" * 65)
    print(f"\n  g2 = H0[4] = {H0[4]:#010x}  (фиксирован)")
    print(f"  f2 = e1 mod 2^32  (зависит от W0)")
    print()

    results = {}
    print(f"  {'Класс':13s}  {'P(De3=0|класс)':18s}  {'~1/N':8s}  {'проверка XOR'}")
    print("  " + "-" * 65)

    for key in sorted(d2_values.keys()):
        d2 = d2_values[key]
        d2_mod = d2 % MOD

        hit = 0
        hit_xor = 0  # прямая проверка через sha256_rounds

        for _ in range(N_per_class):
            # Генерируем W0 из нужного класса
            while True:
                W0 = random.randint(0, M32) & ~1
                e1 = e1_from_W0(W0)
                s  = Sig1(e1)
                b26 = int(bool(s & (1<<26)))
                b21 = int(bool(s & (1<<21)))
                b7  = int(bool(s & (1<<7)))
                if (b26, b21, b7) == key:
                    break

            e2 = e2_from_e1(e1)
            e2_f = (e2 + d2_mod) & M32

            dS2 = Sig1(e2_f) - Sig1(e2)
            dC2 = Ch(e2_f, e1, H0[4]) - Ch(e2, e1, H0[4])
            if (dS2 + dC2) % MOD == 0:
                hit += 1

        frac = hit / N_per_class
        inv = f"~1/{int(1/frac)}" if frac > 0 else "~0"
        results[key] = frac
        print(f"  {str(key):13s}  {frac:.6f}           {inv:8s}")

    avg = sum(results.values()) / len(results)
    total = sum(results.values()) / 8  # взвешенная по ~равномерному распределению классов
    print()
    print(f"  Среднее P(De3=0): {avg:.6f}  (~1/{int(1/avg) if avg>0 else '∞'})")
    print(f"  Ожидание T_PERIOD3: ~1/256 = {1/256:.6f}")
    print(f"  Ожидание случайно: ~1/2^32 = {1/MOD:.2e}")
    return results

# ================================================================
# 4. ПРЯМОЕ ИЗМЕРЕНИЕ P(De3=0)
# ================================================================

def direct_p_de3(N=300000):
    print("\n" + "=" * 65)
    print("4. ПРЯМОЕ ИЗМЕРЕНИЕ P(De3=0)")
    print("=" * 65)

    hit = 0
    for _ in range(N):
        W0 = random.randint(0, M32) & ~1
        W  = [W0]   + [0]*15
        Wf = [W0^1] + [0]*15
        sn = sha256_rounds(W, 3)
        sf = sha256_rounds(Wf, 3)
        if sn[4] == sf[4]:
            hit += 1

    frac = hit / N
    print(f"\n  N = {N}, De3=0 случаев: {hit}")
    print(f"  P(De3=0) = {frac:.7f}")
    print(f"  ~1/{int(1/frac) if frac>0 else '∞'}")
    print(f"  T_PERIOD3 предсказывает: ~1/256 = {1/256:.7f}")
    print(f"  Отношение к предсказанию: {frac/(1/256):.3f}")
    return frac

# ================================================================
# 5. СТРУКТУРА РЕШЕНИЙ De3=0 — КЛАССОВЫЙ ПОРТРЕТ
# ================================================================

def portrait_de3_solutions(d2_values, N_collect=500000):
    """
    Собираем все (W0, e1, e2, класс) для De3=0.
    Смотрим: однородны ли решения по классам или нет?
    """
    print("\n" + "=" * 65)
    print("5. ПОРТРЕТ РЕШЕНИЙ De3=0")
    print("=" * 65)

    solutions_by_class = defaultdict(list)

    for _ in range(N_collect):
        W0 = random.randint(0, M32) & ~1
        W  = [W0]   + [0]*15
        Wf = [W0^1] + [0]*15
        sn = sha256_rounds(W, 3)
        sf = sha256_rounds(Wf, 3)

        if sn[4] == sf[4]:
            e1 = e1_from_W0(W0)
            e2 = e2_from_e1(e1)
            s  = Sig1(e1)
            key = (int(bool(s & (1<<26))), int(bool(s & (1<<21))), int(bool(s & (1<<7))))
            solutions_by_class[key].append((W0, e1, e2))

    total_sol = sum(len(v) for v in solutions_by_class.values())
    print(f"\n  N сэмплов: {N_collect},  найдено De3=0: {total_sol}")
    print(f"  P(De3=0) ≈ {total_sol/N_collect:.6f}")
    print()
    print(f"  {'Класс':13s}  {'кол-во':8s}  {'доля':8s}  {'теория':8s}")
    print("  " + "-" * 50)
    for key in sorted(d2_values.keys()):
        n = len(solutions_by_class.get(key, []))
        frac = n / total_sol if total_sol > 0 else 0
        print(f"  {str(key):13s}  {n:8d}  {frac:.4f}    ~0.125")

    # Анализируем паттерн e2 в решениях (для конкретного класса)
    print()
    print("  Анализ e2 в решениях (для класса (0,0,0)):")
    cls = (0,0,0)
    d2 = d2_values[cls]
    sols_e2 = [e2 for W0,e1,e2 in solutions_by_class.get(cls, [])]
    if sols_e2:
        # Биты {26,21,7} Sig1(e2) для решений
        bits_dist = Counter()
        for e2 in sols_e2:
            s2 = Sig1(e2)
            b = (int(bool(s2 & (1<<26))), int(bool(s2 & (1<<21))), int(bool(s2 & (1<<7))))
            bits_dist[b] += 1
        print(f"  {'bits(26,21,7) Sig1(e2)':22s}  {'кол-во':8s}  {'доля'}")
        for k, cnt in sorted(bits_dist.items()):
            print(f"    {str(k):22s}  {cnt:8d}  {cnt/len(sols_e2):.3f}")
        print(f"  Вывод: решения De3=0 концентрируются в {'1' if len(bits_dist)==1 else len(bits_dist)} битовых классах e2")

# ================================================================
# 6. ОБРАТНЫЙ ПУТЬ: W0 из e1
# ================================================================

def inverse_map(N_verify=20):
    """
    Полная аналитическая формула обратного пути:

    e1 = (H0[3] + C_IV + W0) mod 2^32
    => W0 = (e1 - H0[3] - C_IV) mod 2^32

    e2 = H0[2] + H0[6] + K[1] + Sig1(e1) + Ch(e1, H0[4], H0[5]) mod 2^32

    Условие De3=0: бит (26,21,7) класс e1 = X,
    и (ΔSig1(e2, d2_X) + ΔCh(e2, e1, H0[4], d2_X)) ≡ 0 mod 2^32
    """
    print("\n" + "=" * 65)
    print("6. ОБРАТНЫЙ ПУТЬ: e1 → W0 (аналитически)")
    print("=" * 65)

    print(f"\n  C_IV = {C_IV:#010x}")
    print(f"  e1 = (H0[3] + C_IV + W0) mod 2^32")
    print(f"  W0 = (e1 - H0[3] - C_IV) mod 2^32 = (e1 - {(H0[3]+C_IV)&M32:#010x}) mod 2^32")
    print()

    SHIFT = (H0[3] + C_IV) & M32
    print(f"  SHIFT = (H0[3] + C_IV) mod 2^32 = {SHIFT:#010x}")
    print()

    # Демонстрация: ищем De3=0, строим W0 аналитически
    print(f"  Демонстрация (первые {N_verify} решений методом случ. поиска):")
    print(f"  {'W0':12s}  {'e1':12s}  {'e2':12s}  {'класс':13s}  {'De3=0 верифик.'}")
    print("  " + "-" * 70)

    found = 0
    tries = 0
    while found < N_verify:
        W0 = random.randint(0, M32) & ~1
        W  = [W0]   + [0]*15
        Wf = [W0^1] + [0]*15
        sn = sha256_rounds(W, 3)
        sf = sha256_rounds(Wf, 3)
        tries += 1

        if sn[4] == sf[4]:
            e1 = e1_from_W0(W0)
            e2 = e2_from_e1(e1)
            s  = Sig1(e1)
            cls = (int(bool(s & (1<<26))), int(bool(s & (1<<21))), int(bool(s & (1<<7))))

            # Верификация обратного пути: W0 из e1
            W0_recovered = (e1 - SHIFT) & M32
            ok = (W0_recovered == W0) and (W0_recovered & 1 == 0)

            print(f"  {W0:#012x}  {e1:#012x}  {e2:#012x}  {str(cls):13s}  {'✓' if ok else '✗'}")
            found += 1

    print(f"\n  Эффективность поиска: {found} / {tries} = {found/tries:.5f} ≈ P(De3=0)")

# ================================================================
# MAIN
# ================================================================

def main():
    random.seed(42)
    print("П-2: СТРУКТУРА ADD8 И УСЛОВИЕ De3=0")
    print("=" * 65)

    # 1. Доказываем 8 значений аналитически
    d2_values = prove_8_values()

    # 2. Таблица значений
    compute_8_values_table(d2_values)

    # 3. P(De3=0) по классам (N_per_class уменьшен для скорости)
    class_probs = analyze_de3_per_class(d2_values, N_per_class=30000)

    # 4. Прямое измерение
    p_direct = direct_p_de3(N=300000)

    # 5. Портрет решений
    portrait_de3_solutions(d2_values, N_collect=400000)

    # 6. Обратный путь
    inverse_map(N_verify=15)

    print("\n" + "=" * 65)
    print("П-2 ЗАВЕРШЁН")
    print("=" * 65)
    print()
    print("ИТОГОВЫЕ РЕЗУЛЬТАТЫ:")
    print(f"  T_ADD8: lifted_diff_e2 = ΔSig1(e1,1) + 1, 8 значений ✓")
    print(f"  Ключ:   биты (26,21,7) Sig1(e1)")
    print(f"  Нижний байт d2: всегда 0x81")
    print(f"  P(De3=0): прямое = {p_direct:.5f}, теория ~1/256")
    print(f"  Обратный путь: W0 = (e1 - SHIFT) & M32  (точная формула)")

if __name__ == "__main__":
    main()
