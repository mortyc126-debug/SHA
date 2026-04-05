#!/usr/bin/env python3
"""
SHA-256: полная система из 64 раундов.
Для каждого раунда: все уравнения, зависимости, неизвестные.
Генерируется автоматически из структуры SHA-256.
"""

MASK32 = 0xFFFFFFFF
K=[0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV=[0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

lines = []

def add(s):
    lines.append(s)

add("# SHA-256: ПОЛНАЯ СИСТЕМА 64 РАУНДОВ")
add("# Все уравнения, все зависимости, все неизвестные")
add("#")
add("# Обозначения:")
add("#   a[r], e[r] — активные регистры после раунда r")
add("#   W[r] — слово расписания (W[0..15] = вход M, W[16..63] = schedule)")
add("#   K[r] — раундовая константа (известна)")
add("#   Σ₀(x) = ROTR(x,2) ⊕ ROTR(x,13) ⊕ ROTR(x,22)")
add("#   Σ₁(x) = ROTR(x,6) ⊕ ROTR(x,11) ⊕ ROTR(x,25)")
add("#   Ch(e,f,g) = (e∧f) ⊕ (¬e∧g)")
add("#   Maj(a,b,c) = (a∧b) ⊕ (a∧c) ⊕ (b∧c)")
add("#   σ₀(x) = ROTR(x,7) ⊕ ROTR(x,18) ⊕ SHR(x,3)")
add("#   σ₁(x) = ROTR(x,17) ⊕ ROTR(x,19) ⊕ SHR(x,10)")
add("#   Все операции + mod 2³²")
add("#")
add("# Shift register:")
add("#   b[r] = a[r-1], c[r] = a[r-2], d[r] = a[r-3]")
add("#   f[r] = e[r-1], g[r] = e[r-2], h[r] = e[r-3]")
add("#")
add("# Створочне тождество (доказано 600K/600K):")
add("#   e[r] = a[r] + a[r-4] - Σ₀(a[r-1]) - Maj(a[r-1],a[r-2],a[r-3])")
add("")

# ================================================================
# НАЧАЛЬНОЕ СОСТОЯНИЕ
# ================================================================

add("=" * 80)
add("НАЧАЛЬНОЕ СОСТОЯНИЕ (r=0)")
add("=" * 80)
add(f"  a[0] = IV[0] = 0x{IV[0]:08x}  (ИЗВЕСТНО, константа)")
add(f"  b[0] = IV[1] = 0x{IV[1]:08x}")
add(f"  c[0] = IV[2] = 0x{IV[2]:08x}")
add(f"  d[0] = IV[3] = 0x{IV[3]:08x}")
add(f"  e[0] = IV[4] = 0x{IV[4]:08x}")
add(f"  f[0] = IV[5] = 0x{IV[5]:08x}")
add(f"  g[0] = IV[6] = 0x{IV[6]:08x}")
add(f"  h[0] = IV[7] = 0x{IV[7]:08x}")
add("")

# ================================================================
# РАСПИСАНИЕ
# ================================================================

add("=" * 80)
add("РАСПИСАНИЕ (W[0..63])")
add("=" * 80)
add("")
for r in range(16):
    add(f"  W[{r:>2}] = M[{r}]   (ВХОД, неизвестное, 32 бита)")
add("")
for r in range(16, 64):
    add(f"  W[{r:>2}] = σ₁(W[{r-2:>2}]) + W[{r-7:>2}] + σ₀(W[{r-15:>2}]) + W[{r-16:>2}]")
add("")
add("  Зависимости schedule:")
add("    W[16] зависит от: W[14], W[9], W[1], W[0]")
add("    W[17] зависит от: W[15], W[10], W[2], W[1]")
add("    ...")
add("    W[63] зависит от: W[61], W[56], W[48], W[47]")
add("    → ВСЕ W[16..63] = функции W[0..15] = функции M")
add("")

# ================================================================
# 64 РАУНДА
# ================================================================

for r in range(64):
    add("=" * 80)
    add(f"РАУНД {r}")
    add("=" * 80)

    # Shift register identities
    if r == 0:
        b_src = f"IV[1] = 0x{IV[1]:08x}"
        c_src = f"IV[2] = 0x{IV[2]:08x}"
        d_src = f"IV[3] = 0x{IV[3]:08x}"
        f_src = f"IV[5] = 0x{IV[5]:08x}"
        g_src = f"IV[6] = 0x{IV[6]:08x}"
        h_src = f"IV[7] = 0x{IV[7]:08x}"
        a_src = f"IV[0] = 0x{IV[0]:08x}"
        e_src = f"IV[4] = 0x{IV[4]:08x}"
    else:
        b_src = f"a[{r-1}]"
        c_src = f"a[{r-2}]" if r >= 2 else f"IV[{1}]"
        d_src = f"a[{r-3}]" if r >= 3 else f"IV[{3-r+r}]"
        if r == 1: d_src = f"IV[2] = 0x{IV[2]:08x}"
        elif r == 2: d_src = f"IV[1] = 0x{IV[1]:08x}"
        elif r == 3: d_src = f"IV[0] = 0x{IV[0]:08x}"
        else: d_src = f"a[{r-3}]"

        f_src = f"e[{r-1}]"
        g_src = f"e[{r-2}]" if r >= 2 else f"IV[{5}]"
        if r == 1: h_src = f"IV[6] = 0x{IV[6]:08x}"
        elif r == 2: h_src = f"IV[5] = 0x{IV[5]:08x}"
        elif r == 3: h_src = f"IV[4] = 0x{IV[4]:08x}"
        else: h_src = f"e[{r-3}]"

        a_src = f"a[{r}]"
        e_src = f"e[{r}]"

    add(f"")
    add(f"  Входное состояние (через shift register):")
    add(f"    a = {a_src}")
    add(f"    b = {b_src}")
    add(f"    c = {c_src}")
    add(f"    d = {d_src}")
    add(f"    e = {e_src}")
    add(f"    f = {f_src}")
    add(f"    g = {g_src}")
    add(f"    h = {h_src}")

    # Determine what's known/unknown at this round
    w_status = "ВХОД (неизвестное)" if r < 16 else "schedule (функция M)"

    add(f"")
    add(f"  Слово расписания:")
    add(f"    W[{r}] = {w_status}")
    if r >= 16:
        add(f"    W[{r}] = σ₁(W[{r-2}]) + W[{r-7}] + σ₀(W[{r-15}]) + W[{r-16}]")
    add(f"    K[{r}] = 0x{K[r]:08x} (константа)")

    add(f"")
    add(f"  Уравнения раунда:")

    # Equation 1: T1
    if r == 0:
        add(f"    T1[{r}] = IV[7] + Σ₁(IV[4]) + Ch(IV[4],IV[5],IV[6]) + K[{r}] + W[{r}]")
        add(f"           = 0x5be0cd19 + Σ₁(0x510e527f) + Ch(0x510e527f,0x9b05688c,0x1f83d9ab) + 0x{K[r]:08x} + W[{r}]")
        add(f"           = 0xf377ed68 + W[{r}]   ← ЛИНЕЙНО по W[{r}]")
    else:
        add(f"    T1[{r}] = {h_src} + Σ₁({e_src}) + Ch({e_src},{f_src},{g_src}) + K[{r}] + W[{r}]")

    # Equation 2: T2
    add(f"    T2[{r}] = Σ₀({a_src}) + Maj({a_src},{b_src},{c_src})")
    if r == 0:
        add(f"           = Σ₀(0x6a09e667) + Maj(0x6a09e667,0xbb67ae85,0x3c6ef372)")
        add(f"           = 0x08909ae5   ← КОНСТАНТА (не зависит от M)")
    add(f"    → T2[{r}] НЕ зависит от W[{r}]")

    # Equation 3: a_new
    add(f"    a[{r+1}] = T1[{r}] + T2[{r}]")
    if r == 0:
        add(f"            = 0xfc08884d + W[0]   ← ЛИНЕЙНО по W[0]")

    # Equation 4: e_new
    add(f"    e[{r+1}] = {d_src} + T1[{r}]")

    # Створочне
    if r >= 3:
        add(f"")
        add(f"  Створочне тождество (e через a-chain, без W):")
        add(f"    e[{r+1}] = a[{r+1}] + a[{r-3}] - Σ₀(a[{r}]) - Maj(a[{r}],a[{r-1}],a[{r-2}])")

    # W recovery
    add(f"")
    add(f"  Восстановление W[{r}] (если знать state):")
    add(f"    T2[{r}] = Σ₀(a[{r}]) + Maj(a[{r}],a[{r-1 if r>0 else 'IV1'}],a[{r-2 if r>1 else 'IV2'}])")
    add(f"    T1[{r}] = a[{r+1}] - T2[{r}]")
    add(f"    W[{r}]  = T1[{r}] - h[{r}] - Σ₁(e[{r}]) - Ch(e[{r}],f[{r}],g[{r}]) - K[{r}]")

    # Dependency summary
    deps_a = set()
    deps_e = set()
    if r == 0:
        deps_a.add("W[0]")
        deps_e.add("W[0]")
    else:
        # a[r+1] depends on T1+T2, which depends on state[r] + W[r]
        # state[r] = (a[r],...,a[r-3], e[r],...,e[r-3])
        for i in range(4):
            if r-i >= 1: deps_a.add(f"a[{r-i}]")
            if r-i >= 1: deps_e.add(f"a[{r-i}]")
        for i in range(4):
            if r-i >= 1: deps_a.add(f"e[{r-i}]")
            if r-i >= 1: deps_e.add(f"e[{r-i}]")
        deps_a.add(f"W[{r}]")
        deps_e.add(f"W[{r}]")

    add(f"")
    add(f"  Зависимости:")
    add(f"    a[{r+1}] зависит от: {', '.join(sorted(deps_a))}")
    add(f"    e[{r+1}] зависит от: {', '.join(sorted(deps_e))}")

    # Count unknowns
    if r < 16:
        add(f"    Новое неизвестное: W[{r}] (32 бита)")
    else:
        add(f"    Новых неизвестных: 0 (W[{r}] определён schedule)")

    add("")

# ================================================================
# ФИНАЛИЗАЦИЯ
# ================================================================

add("=" * 80)
add("ФИНАЛИЗАЦИЯ (feedforward)")
add("=" * 80)
add("")
add("  H[0] = a[64] + IV[0]   mod 2³²")
add("  H[1] = a[63] + IV[1]   (= b[64] + IV[1])")
add("  H[2] = a[62] + IV[2]   (= c[64] + IV[2])")
add("  H[3] = a[61] + IV[3]   (= d[64] + IV[3])")
add("  H[4] = e[64] + IV[4]")
add("  H[5] = e[63] + IV[5]   (= f[64] + IV[5])")
add("  H[6] = e[62] + IV[6]   (= g[64] + IV[6])")
add("  H[7] = e[61] + IV[7]   (= h[64] + IV[7])")
add("")
add("  Из H извлекаются:")
add("    a[61..64] = H[3]-IV[3], H[2]-IV[2], H[1]-IV[1], H[0]-IV[0]")
add("    e[61..64] = H[7]-IV[7], H[6]-IV[6], H[5]-IV[5], H[4]-IV[4]")
add("    a[57..60] через backward (створочне), БЕЗ W")
add("    e[57..60] НЕИЗВЕСТНЫ без W[60..63]")
add("")

# ================================================================
# ПОДСЧЁТ
# ================================================================

add("=" * 80)
add("ИТОГОВЫЙ ПОДСЧЁТ")
add("=" * 80)
add("")
add("  Всего уравнений: 64 раунда × 2 (a,e) = 128 уравнений по 32 бита")
add("  Всего неизвестных: W[0..15] = 16 слов = 512 бит")
add("  Schedule: W[16..63] = 48 слов, определены через W[0..15]")
add("  Число независимых переменных: 512 бит (= M)")
add("  Число ограничений (H = target): 256 бит")
add("  Свобода: 512 - 256 = 256 бит")
add("")
add("  Для КОЛЛИЗИИ: два M с одинаковым H")
add("    Переменные: M₁ (512 бит) + M₂ (512 бит) = 1024 бит")
add("    Ограничения: SHA(M₁) = SHA(M₂) = 256 бит")
add("    Свобода: 1024 - 256 = 768 бит")
add("    Но: birthday на 256-бит выходе = 2^128")
add("")
add("  СТЕНА: schedule связывает W[0..15] с W[16..63].")
add("  Изменение ЛЮБОГО W[k] для k<16 каскадирует через ВСЕ 48 слов W[16..63].")
add("  Backward chain (a[57..64] из H) не обходит schedule.")
add("  Birthday = 2^128 — оптимален.")

# Write output
output = '\n'.join(lines)
print(output[:3000])
print(f"\n... [{len(lines)} total lines]")

with open('/home/user/SHA/SHA256_64_EQUATIONS.md', 'w') as f:
    f.write(output)

print(f"\nWritten to SHA256_64_EQUATIONS.md ({len(output)} bytes)")
