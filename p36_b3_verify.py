"""
B3 уточнение: Верификация чередующегося каскада.
Предыдущий скрипт имел ошибку в проверке (смещение индекса).
Проверяем правильно: δe_r и δa_r при чередующемся каскаде.
"""
import random
MASK = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32-n))) & MASK
def sig0(x):  return rotr(x,7) ^ rotr(x,18) ^ (x>>3)
def sig1(x):  return rotr(x,17) ^ rotr(x,19) ^ (x>>10)
def Sig0(x):  return rotr(x,2) ^ rotr(x,13) ^ rotr(x,22)
def Sig1(x):  return rotr(x,6) ^ rotr(x,11) ^ rotr(x,25)
def Ch(e,f,g):  return ((e&f) ^ (~e&g)) & MASK
def Maj(a,b,c): return (a&b) ^ (a&c) ^ (b&c)
def hw(x): return bin(x).count('1')

K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
     0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
     0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
     0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
     0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
     0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def make_schedule(W16):
    W = list(W16) + [0]*48
    for i in range(16, 64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK
    return W

def sha_rounds(W, R):
    a,b,c,d,e,f,g,h = IV
    states = [[a,b,c,d,e,f,g,h]]
    for r in range(R):
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a,b,c)) & MASK
        h=g; g=f; f=e; e=(d+T1)&MASK
        d=c; c=b; b=a; a=(T1+T2)&MASK
        states.append([a,b,c,d,e,f,g,h])
    return states

def alternating_cascade_correct(W0, W1, DW0=1):
    """
    Чередующийся каскад (ИСПРАВЛЕННЫЙ):
    ΔW_i выбирается для обнуления цели на раунде i+1:
      i=2: De3=0
      i=3: Da4=0
      i=4: De5=0
      i=5: Da6=0
      ...
      i=14: De15=0
      i=15: Da16=0
    Итого: δe∈{3,5,7,9,11,13,15} и δa∈{4,6,8,10,12,14,16}
    """
    Wn = [W0, W1] + [0]*14
    DWs = [0]*16; DWs[0] = DW0

    for i in range(2, 16):  # i = W-индекс и текущий раунд
        target_r = i + 1    # обнуляем на раунде i+1
        Wfc = [(Wn[k]+DWs[k])&MASK for k in range(16)]
        sn = sha_rounds(make_schedule(Wn), target_r)
        sf = sha_rounds(make_schedule(Wfc), target_r)

        if i % 2 == 0:  # чётный i → обнуляем δe на раунде i+1
            nat = (sf[target_r][4] - sn[target_r][4]) & MASK
        else:            # нечётный i → обнуляем δa на раунде i+1
            nat = (sf[target_r][0] - sn[target_r][0]) & MASK

        DWs[i] = (-nat) & MASK

    Wf = [(Wn[k]+DWs[k])&MASK for k in range(16)]
    sn_full = sha_rounds(make_schedule(Wn), 18)
    sf_full = sha_rounds(make_schedule(Wf), 18)
    return DWs, sn_full, sf_full

print("=" * 70)
print("B3 ВЕРИФИКАЦИЯ: Чередующийся каскад δe∈{3,5,...,15}, δa∈{4,6,...,16}")
print("=" * 70)

# Тест на 5 случайных парах
for trial in range(5):
    W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
    DWs, sn, sf = alternating_cascade_correct(W0, W1)
    de_zeros = []; da_zeros = []
    row = []
    for r in range(1, 19):
        if r >= len(sn): break
        de_r = (sf[r][4] - sn[r][4]) & MASK
        da_r = (sf[r][0] - sn[r][0]) & MASK
        row.append((r, de_r, da_r))
        if r in range(3, 16, 2) and de_r == 0: de_zeros.append(r)
        if r in range(4, 17, 2) and da_r == 0: da_zeros.append(r)
    print(f"\nПара {trial+1}: W0=0x{W0:08x} W1=0x{W1:08x}")
    print(f"{'r':>3} | {'δe_r HW':>9} | {'δa_r HW':>9} | Цель      | Достигнуто")
    print("-" * 55)
    for r, de_r, da_r in row:
        if r < 3 or r > 17: continue
        if r % 2 == 1 and r >= 3:   # нечётные r → δe должно быть 0 (поставлено ΔW_{r-1})
            ok = "δe=0 ✓" if de_r==0 else f"δe≠0 ✗ HW={hw(de_r)}"
        elif r % 2 == 0 and r >= 4: # чётные r → δa должно быть 0 (поставлено ΔW_{r-1})
            ok = "δa=0 ✓" if da_r==0 else f"δa≠0 ✗ HW={hw(da_r)}"
        else:
            ok = "—"
        print(f"  {r:>2} | {hw(de_r):>9} | {hw(da_r):>9} | {'δe→0' if r%2==1 else 'δa→0':>9} | {ok}")
    print(f"  δe=0 на нечётных {list(range(3,16,2))}: {len(de_zeros)}/7")
    print(f"  δa=0 на чётных   {list(range(4,17,2))}: {len(da_zeros)}/7")

# Сводная статистика по N=200 парам
print("\n\nСводная статистика (N=200 пар):")
import statistics as st
de_odd_rates = []; da_even_rates = []
for _ in range(200):
    W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
    _, sn, sf = alternating_cascade_correct(W0, W1)
    de_count = sum(1 for r in range(3,16,2) if (sf[r][4]-sn[r][4])&MASK==0)
    da_count = sum(1 for r in range(4,17,2) if (sf[r][0]-sn[r][0])&MASK==0)
    de_odd_rates.append(de_count)
    da_even_rates.append(da_count)

print(f"  δe=0 на нечётных r (из 7): avg={st.mean(de_odd_rates):.3f}  "
      f"min={min(de_odd_rates)}  max={max(de_odd_rates)}")
print(f"  δa=0 на чётных   r (из 7): avg={st.mean(da_even_rates):.3f}  "
      f"min={min(da_even_rates)}  max={max(da_even_rates)}")
print(f"\n  Ожидание при P=1 (детерм.): avg=7.000")
print(f"  Если avg=7: чередующийся каскад детерминирован (P=1) ✓")
print(f"  Если avg<7: частичный успех, есть взаимные зависимости ✗")
