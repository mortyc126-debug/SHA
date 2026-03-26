"""
Проверка: есть ли паттерн в хороших nonce?
Используем СЛУЧАЙНЫЕ nonce вместо последовательных, чтобы убрать артефакт.
"""
import numpy as np
import struct
import hashlib

def sha256_raw(data): return hashlib.sha256(data).digest()
def double_sha256(data): return sha256_raw(sha256_raw(data))
def leading_zeros(h):
    n = int.from_bytes(h, 'big')
    if n == 0: return 256
    return 256 - n.bit_length()
def hw(x): return bin(x).count('1')

MASK32 = 0xFFFFFFFF
np.random.seed(42)
header_base = np.random.bytes(76)

NUM = 500000
target_lz = 10

good_nonces = []
all_hw = []

for _ in range(NUM):
    nonce = np.random.randint(0, 2**31)  # случайный nonce
    h = header_base + struct.pack('<I', nonce & MASK32)
    final = double_sha256(h)
    lz = leading_zeros(final)
    all_hw.append(hw(nonce))
    if lz >= target_lz:
        good_nonces.append(nonce)

good_hw = [hw(n) for n in good_nonces]

print(f"Протестировано: {NUM:,}")
print(f"Хороших (≥{target_lz} LZ): {len(good_nonces)}")
print()
print(f"HW(хороших nonce):  {np.mean(good_hw):.1f} ± {np.std(good_hw):.1f}")
print(f"HW(всех nonce):     {np.mean(all_hw):.1f} ± {np.std(all_hw):.1f}")
print(f"→ {'РАЗЛИЧИЕ!' if abs(np.mean(good_hw) - np.mean(all_hw)) > 1.0 else 'Нет различия (ожидаемо)'}")

# Битовые позиции
if len(good_nonces) > 20:
    bit_freq_good = np.zeros(32)
    bit_freq_all = np.zeros(32)
    for n in good_nonces:
        for b in range(32):
            if n & (1 << b): bit_freq_good[b] += 1
    bit_freq_good /= len(good_nonces)

    for n in range(min(NUM, 100000)):
        nonce = np.random.randint(0, 2**31)
        for b in range(32):
            if nonce & (1 << b): bit_freq_all[b] += 1
    bit_freq_all /= min(NUM, 100000)

    print(f"\nЧастота бит = 1 по позициям:")
    print(f"{'Бит':>4} {'Хорошие':>10} {'Все':>10} {'Δ':>8}")
    max_delta = 0
    for b in range(32):
        delta = abs(bit_freq_good[b] - bit_freq_all[b])
        max_delta = max(max_delta, delta)
        if delta > 0.03 or b < 5 or b > 27:
            print(f"{b:4d} {bit_freq_good[b]:10.4f} {bit_freq_all[b]:10.4f} {delta:8.4f}")

    print(f"\nМакс. отклонение: {max_delta:.4f}")
    print(f"→ {'ПАТТЕРН!' if max_delta > 0.05 else 'Нет паттерна — хорошие nonce РАВНОМЕРНЫ'}")
