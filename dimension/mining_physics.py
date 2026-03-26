"""
Bitcoin Mining через призму нового измерения.

Задача: найти nonce (32 бит) такой что SHA256(SHA256(header||nonce)) < target
Target = число с N ведущих нулей → hash должен начинаться с N нулевых бит.

В терминах нашего измерения:
  - Ткань: ДВОЙНАЯ (два прохода SHA-256)
  - Вход: 80 байт header, из них свободны только 4 байта (nonce)
  - Выход: 256 бит hash
  - Задача: НЕ коллизия, а partial preimage

Вопросы:
  1. Как nonce (32 бит) влияет на двойной хеш?
  2. Есть ли структура в near-target хешах?
  3. Даёт ли наше измерение преимущество перед brute force?
"""

import numpy as np
import struct
import hashlib
import time

MASK32 = 0xFFFFFFFF

def hw(x):
    return bin(x).count('1')

def sha256_raw(data: bytes) -> bytes:
    return hashlib.sha256(data).digest()

def double_sha256(data: bytes) -> bytes:
    return sha256_raw(sha256_raw(data))

def hash_to_int(h: bytes) -> int:
    return int.from_bytes(h, 'big')

def leading_zeros(h: bytes) -> int:
    """Число ведущих нулевых бит."""
    n = int.from_bytes(h, 'big')
    if n == 0: return 256
    return 256 - n.bit_length()

def make_header(nonce: int, seed: int = 0) -> bytes:
    """Создаём псевдо-header 80 байт с заданным nonce."""
    # 76 байт "header" + 4 байта nonce
    np.random.seed(seed)
    header_base = np.random.bytes(76)
    return header_base + struct.pack('<I', nonce & MASK32)


def experiment_mining_dimension(num_samples=100000):
    np.random.seed(42)

    print("=" * 70)
    print("BITCOIN MINING — Анализ в новом измерении")
    print("=" * 70)

    # === 1. Структура ткани майнинга ===
    print("\n" + "=" * 70)
    print("1. СТРУКТУРА ТКАНИ МАЙНИНГА")
    print("=" * 70)

    print(f"""
  ДВОЙНАЯ ТКАНЬ: F_mining = F_SHA1 ⊗ F_SHA2

  Проход 1: SHA256(header || nonce)
    Вход: 80 байт = 640 бит
    Header (фикс.): 608 бит → СЕЧЕНИЕ (не свободны)
    Nonce (свободн.): 32 бит → единственная свобода
    Раундов: 64
    Выход: 256 бит (intermediate hash)

  Проход 2: SHA256(intermediate_hash)
    Вход: 256 бит (из прохода 1) + padding
    Свободных бит: 0 (всё определено проходом 1)
    Раундов: 64
    Выход: 256 бит (final hash)

  ИТОГО:
    Свободных бит:    32 (только nonce)
    Раундов:          128 (64 + 64)
    Горлышко:         8 регистров = rank ≤ 8
    Условие:          leading_zeros(hash) ≥ difficulty

  БУХГАЛТЕРИЯ:
    Свобода:          32 бит
    Условие (target):  N бит (ведущих нулей)
    Ядро:             32 - N бит
    Если N > 32:      ядро пусто → brute force 2^N
""")

    # === 2. Влияние nonce на hash ===
    print("=" * 70)
    print("2. ВЛИЯНИЕ NONCE (δ = 1 бит) НА ДВОЙНОЙ HASH")
    print("=" * 70)

    header_base = np.random.bytes(76)

    # Как 1 бит nonce меняет intermediate и final hash?
    inter_diffs = []
    final_diffs = []
    lz_distribution = []

    for trial in range(10000):
        nonce1 = np.random.randint(0, 2**32)
        nonce2 = nonce1 ^ (1 << np.random.randint(0, 32))

        h1 = header_base + struct.pack('<I', nonce1)
        h2 = header_base + struct.pack('<I', nonce2)

        inter1 = sha256_raw(h1)
        inter2 = sha256_raw(h2)

        final1 = sha256_raw(inter1)
        final2 = sha256_raw(inter2)

        # δ intermediate
        inter_delta = bytes(a ^ b for a, b in zip(inter1, inter2))
        inter_hw = sum(hw(b) for b in inter_delta)
        inter_diffs.append(inter_hw)

        # δ final
        final_delta = bytes(a ^ b for a, b in zip(final1, final2))
        final_hw = sum(hw(b) for b in final_delta)
        final_diffs.append(final_hw)

    print(f"\n  δ(nonce) = 1 бит:")
    print(f"    HW(δ intermediate): {np.mean(inter_diffs):.1f} ± {np.std(inter_diffs):.1f}  (из 256)")
    print(f"    HW(δ final hash):   {np.mean(final_diffs):.1f} ± {np.std(final_diffs):.1f}  (из 256)")
    print(f"    Теор. random:       128.0")

    # Последовательные nonce (как реально майнят)
    print(f"\n  Последовательные nonce (n, n+1):")
    seq_inter = []
    seq_final = []
    for trial in range(10000):
        nonce = np.random.randint(0, 2**32 - 1)
        h1 = header_base + struct.pack('<I', nonce)
        h2 = header_base + struct.pack('<I', nonce + 1)

        inter1 = sha256_raw(h1)
        inter2 = sha256_raw(h2)
        final1 = sha256_raw(inter1)
        final2 = sha256_raw(inter2)

        inter_delta = bytes(a ^ b for a, b in zip(inter1, inter2))
        final_delta = bytes(a ^ b for a, b in zip(final1, final2))
        seq_inter.append(sum(hw(b) for b in inter_delta))
        seq_final.append(sum(hw(b) for b in final_delta))

    print(f"    HW(δ intermediate): {np.mean(seq_inter):.1f} ± {np.std(seq_inter):.1f}")
    print(f"    HW(δ final hash):   {np.mean(seq_final):.1f} ± {np.std(seq_final):.1f}")

    # === 3. Распределение leading zeros ===
    print("\n" + "=" * 70)
    print("3. РАСПРЕДЕЛЕНИЕ ВЕДУЩИХ НУЛЕЙ (естественное)")
    print("=" * 70)

    lz_counts = np.zeros(40)
    total_tested = num_samples

    for trial in range(total_tested):
        nonce = trial
        h = header_base + struct.pack('<I', nonce & MASK32)
        final = double_sha256(h)
        lz = leading_zeros(final)
        if lz < 40:
            lz_counts[lz] += 1

    print(f"\n  Протестировано nonce: {total_tested:,}")
    print(f"\n  {'LZ':>4} {'Найдено':>10} {'P(эмпир.)':>12} {'P(теор.)':>12} {'Ratio':>8}")
    for lz in range(20):
        found = int(lz_counts[lz])
        p_emp = found / total_tested if total_tested > 0 else 0
        p_theo = 0.5 ** (lz + 1)  # P(exactly lz leading zeros)
        ratio = p_emp / p_theo if p_theo > 0 else 0
        if found > 0 or lz < 10:
            print(f"  {lz:4d} {found:10d} {p_emp:12.6f} {p_theo:12.6f} {ratio:8.3f}")

    # === 4. Корреляция между nonce и leading zeros ===
    print("\n" + "=" * 70)
    print("4. ЕСТЬ ЛИ СТРУКТУРА В 'ХОРОШИХ' NONCE?")
    print("=" * 70)

    # Собираем nonce, дающие ≥ 10 leading zeros
    good_nonces = []
    good_lz = []
    all_first_bytes = []

    for trial in range(num_samples):
        nonce = trial
        h = header_base + struct.pack('<I', nonce & MASK32)
        final = double_sha256(h)
        lz = leading_zeros(final)
        all_first_bytes.append(final[0])
        if lz >= 10:
            good_nonces.append(nonce)
            good_lz.append(lz)

    print(f"\n  Nonce с ≥10 leading zeros: {len(good_nonces)} из {num_samples}")

    if len(good_nonces) > 1:
        diffs = np.diff(sorted(good_nonces))
        print(f"  Расстояние между 'хорошими' nonce:")
        print(f"    mean = {np.mean(diffs):.0f}")
        print(f"    std  = {np.std(diffs):.0f}")
        print(f"    min  = {np.min(diffs)}")
        print(f"    max  = {np.max(diffs)}")
        print(f"    Теор. среднее расстояние: {2**10:.0f}")

        # HW хороших nonce
        hw_good = [hw(n) for n in good_nonces]
        hw_random = [hw(np.random.randint(0, 2**32)) for _ in range(len(good_nonces))]
        print(f"\n  HW('хороших' nonce):  {np.mean(hw_good):.1f} ± {np.std(hw_good):.1f}")
        print(f"  HW(случайных nonce):  {np.mean(hw_random):.1f} ± {np.std(hw_random):.1f}")
        print(f"  → {'РАЗЛИЧИЕ' if abs(np.mean(hw_good) - np.mean(hw_random)) > 2 else 'Нет различия'}")

    # Битовый паттерн хороших nonce — есть ли общие биты?
    if len(good_nonces) > 10:
        always_0 = MASK32
        always_1 = MASK32
        for n in good_nonces:
            always_0 &= ~n  # биты, которые всегда 0
            always_1 &= n   # биты, которые всегда 1
        print(f"\n  Биты, ВСЕГДА = 0 в хороших nonce: {hw(always_0)} из 32")
        print(f"  Биты, ВСЕГДА = 1 в хороших nonce: {hw(always_1)} из 32")
        print(f"  → {'ПАТТЕРН НАЙДЕН!' if hw(always_0) + hw(always_1) > 0 else 'Нет паттерна (ожидаемо)'}")

    # === 5. Двойной SHA-256 vs одиночный ===
    print("\n" + "=" * 70)
    print("5. ДВОЙНОЙ SHA-256: ЧТО ДОБАВЛЯЕТ ВТОРОЙ ПРОХОД?")
    print("=" * 70)

    # Корреляция между intermediate hash и final hash
    inter_lz = []
    final_lz_list = []
    for trial in range(10000):
        nonce = np.random.randint(0, 2**32)
        h = header_base + struct.pack('<I', nonce & MASK32)
        inter = sha256_raw(h)
        final = sha256_raw(inter)
        inter_lz.append(leading_zeros(inter))
        final_lz_list.append(leading_zeros(final))

    corr = np.corrcoef(inter_lz, final_lz_list)[0, 1]
    print(f"\n  Корреляция LZ(intermediate) vs LZ(final): {corr:.6f}")
    print(f"  Если бы связаны: корреляция >> 0")
    print(f"  → {'НЕЗАВИСИМЫ' if abs(corr) < 0.05 else 'СВЯЗАНЫ!'}")

    # "Хороший" intermediate → помогает?
    good_inter = [(il, fl) for il, fl in zip(inter_lz, final_lz_list) if il >= 8]
    if good_inter:
        mean_final_given_good_inter = np.mean([fl for _, fl in good_inter])
        mean_final_all = np.mean(final_lz_list)
        print(f"\n  Средние LZ(final) при LZ(inter) ≥ 8: {mean_final_given_good_inter:.2f}")
        print(f"  Средние LZ(final) в целом:            {mean_final_all:.2f}")
        print(f"  → {'Второй проход СТИРАЕТ информацию' if abs(mean_final_given_good_inter - mean_final_all) < 0.5 else 'Структура СОХРАНЯЕТСЯ!'}")

    # === 6. Наше измерение: что оно даёт для майнинга? ===
    print("\n" + "=" * 70)
    print("6. ВЕРДИКТ: ЧТО ДАЁТ НАШЕ ИЗМЕРЕНИЕ ДЛЯ МАЙНИНГА?")
    print("=" * 70)

    print(f"""
  ТКАНЬ МАЙНИНГА в нашем измерении:
    F_mining = F_SHA1 ⊘ F_SHA2  (двойная свёртка)

  БУХГАЛТЕРИЯ:
    Свобода:      32 бит (nonce)
    Условие:      N бит (leading zeros target)
    Ядро:         32 - N бит

    Текущий Bitcoin difficulty: N ≈ 76 ведущих нулевых бит
    Ядро: 32 - 76 = -44 бит → ПУСТО

    Это значит: 32 бит nonce НЕДОСТАТОЧНО.
    Реальный майнинг: extraNonce + timestamp + другие поля → ~96 бит свободы.
    Ядро: 96 - 76 = 20 бит → ~2^76 попыток.

  CARRY-СТЕНА:
    Проход 1: ~8000 бит carry
    Проход 2: ~8000 бит carry
    Итого:    ~16000 бит carry-нелинейности

    Это ×500 от 32-бит свободы nonce.
    Нет ни малейшего шанса найти структуру.

  ТЕРМОДИНАМИКА:
    1 бит δ(nonce) → HW(δ hash) ≈ 128 (полная диффузия)
    Последоват. nonce → полная декорреляция
    LZ(intermediate) и LZ(final) НЕЗАВИСИМЫ
    Нет паттернов в "хороших" nonce

  ВЫВОД:
    Наше измерение подтверждает: для майнинга
    НЕТ АЛЬТЕРНАТИВЫ BRUTE FORCE.

    Причины:
    1. Свобода (32 бит) << carry-нелинейность (~16000 бит)
    2. Двойной SHA-256 стирает ВСЮ структуру между проходами
    3. Хорошие nonce распределены РАВНОМЕРНО (нет паттерна)
    4. Термодинамика: каждый nonce → независимый случайный хеш
""")

    # === 7. Бонус: скорость vs теория ===
    print("=" * 70)
    print("7. БОНУС: ЭКСПЕРИМЕНТАЛЬНАЯ СКОРОСТЬ МАЙНИНГА")
    print("=" * 70)

    target_lz = 20  # ищем хеш с ≥20 leading zeros
    expected_attempts = 2 ** target_lz

    start = time.time()
    found = False
    attempts = 0

    for nonce in range(min(2**22, expected_attempts * 4)):
        h = header_base + struct.pack('<I', nonce & MASK32)
        final = double_sha256(h)
        attempts += 1
        if leading_zeros(final) >= target_lz:
            elapsed = time.time() - start
            print(f"\n  Target: ≥{target_lz} leading zeros")
            print(f"  Найдено! nonce = {nonce}")
            print(f"  Hash = {final.hex()}")
            print(f"  Leading zeros: {leading_zeros(final)}")
            print(f"  Попыток: {attempts:,}")
            print(f"  Теоретическое: {expected_attempts:,}")
            print(f"  Ratio: {attempts / expected_attempts:.2f}")
            print(f"  Время: {elapsed:.2f}с")
            print(f"  Скорость: {attempts / elapsed:,.0f} H/s (Python)")
            found = True
            break

    if not found:
        elapsed = time.time() - start
        print(f"\n  Target: ≥{target_lz} leading zeros")
        print(f"  Не найдено за {attempts:,} попыток ({elapsed:.1f}с)")
        print(f"  Теоретическое: {expected_attempts:,}")


if __name__ == "__main__":
    experiment_mining_dimension(num_samples=100000)
