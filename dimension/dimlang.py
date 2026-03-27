"""
DIMENSION LANGUAGE (DimLang) — язык программирования нашего измерения.

Это НЕ обычный язык. Его примитивы = объекты нашего измерения.
Его операции = физика нашего измерения.
Программа на DimLang = путешествие через fabric SHA-256.

ПРИМИТИВЫ (типы):
  State   — 8 регистров × 32 бита (точка в пространстве)
  Word    — 32-битное слово (единица информации)
  Message — 16 слов (входной вектор)
  Delta   — разница между двумя State/Message (вектор отклонения)
  Path    — последовательность State через раунды (траектория)
  Fabric  — множество Path (ткань пространства)
  Mark    — метка на Fabric (точка интереса)

ОПЕРАЦИИ (встроенные):
  repair(path, node, round)   — починить node (a или e) в раунде
  propagate(path, rounds)     — протолкнуть path вперёд на N раундов
  section(fabric, round)      — срез ткани на раунде r
  overlay(path1, path2)       — наложить два пути, измерить δ
  fold(fabric)                — свернуть ткань (найти пересечения)
  absorb(word, rounds)        — измерить поглощение слова
  split(path, round)          — разделить path на два отрезка в точке r
  merge(path1, path2, round)  — соединить два пути в точке r
  cost(operation)             — стоимость операции в единицах измерения

НОВЫЕ (неисследованные):
  warp(path, transform)       — деформировать path нестандартно
  resonate(paths)             — найти резонанс между путями
  tunnel(path, from_r, to_r)  — "пройти" через раунды без вычисления
  invert(path, round)         — обратить path от раунда назад

Мы НЕ ЗНАЕМ пределов этих операций. Тестируем.
"""

import numpy as np
import struct, hashlib
from collections import namedtuple
from typing import List, Tuple, Optional, Dict

MASK32 = 0xFFFFFFFF

# ═══════════════════════════════════════════════════════════════
# SHA-256 PRIMITIVES (hardware layer)
# ═══════════════════════════════════════════════════════════════

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(*args):
    r = 0
    for a in args: r = (r + a) & MASK32
    return r
def sub32(x, y): return (x - y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
]
IV_CONST = (0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
            0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19)


# ═══════════════════════════════════════════════════════════════
# DIMENSION TYPES
# ═══════════════════════════════════════════════════════════════

class State:
    """Точка в пространстве SHA-256: 8 регистров."""
    __slots__ = ['regs']
    def __init__(self, regs=None):
        self.regs = tuple(regs) if regs else IV_CONST

    def __getitem__(self, i): return self.regs[i]
    def __eq__(self, other): return self.regs == other.regs
    def __repr__(self): return f"State({','.join(hex(r)[:6] for r in self.regs)})"

    @property
    def a(self): return self.regs[0]
    @property
    def e(self): return self.regs[4]

    def delta(self, other):
        """Расстояние до другого State."""
        return Delta([self.regs[i] ^ other.regs[i] for i in range(8)])


class Delta:
    """Вектор отклонения между двумя State."""
    __slots__ = ['bits']
    def __init__(self, bits):
        self.bits = tuple(bits)

    @property
    def weight(self):
        return sum(hw(b) for b in self.bits)

    @property
    def is_zero(self):
        return all(b == 0 for b in self.bits)

    def __repr__(self): return f"Delta(w={self.weight})"


class Word:
    """32-битное слово — единица информации."""
    __slots__ = ['val']
    def __init__(self, val=0): self.val = val & MASK32
    def __repr__(self): return f"W({hex(self.val)})"
    def __xor__(self, other): return Word(self.val ^ other.val)
    def __eq__(self, other): return self.val == other.val


class Message:
    """16 слов — входной вектор."""
    __slots__ = ['words']
    def __init__(self, words=None):
        if words is None:
            self.words = [Word(np.random.randint(0, 2**32)) for _ in range(16)]
        else:
            self.words = [Word(w) if isinstance(w, int) else w for w in words]

    def __getitem__(self, i): return self.words[i]
    def __setitem__(self, i, v):
        self.words[i] = Word(v) if isinstance(v, int) else v

    def expand(self):
        """Расширить до 64 слов через schedule."""
        W = [w.val for w in self.words]
        for r in range(16, 64):
            W.append(add32(sigma1(W[r-2]), W[r-7], sigma0(W[r-15]), W[r-16]))
        return [Word(w) for w in W]

    def raw(self):
        return [w.val for w in self.words]

    def delta(self, other):
        return [Word(self.words[i].val ^ other.words[i].val) for i in range(16)]

    def copy(self):
        return Message([w.val for w in self.words])


class Path:
    """Траектория через раунды: последовательность State + Words."""
    __slots__ = ['states', 'schedule', 'message', 'n_rounds']

    def __init__(self, message, n_rounds=64):
        self.message = message
        self.n_rounds = n_rounds
        self.schedule = message.expand()

        # Compute all states
        self.states = [State()]
        a,b,c,d,e,f,g,h = IV_CONST
        for r in range(n_rounds):
            w = self.schedule[r].val
            T1 = add32(h, Sigma1(e), Ch(e,f,g), K[r], w)
            T2 = add32(Sigma0(a), Maj(a,b,c))
            h,g,f,e = g,f,e,add32(d,T1)
            d,c,b,a = c,b,a,add32(T1,T2)
            self.states.append(State((a,b,c,d,e,f,g,h)))

    @property
    def final_hash(self):
        s = self.states[self.n_rounds]
        return State([add32(IV_CONST[i], s[i]) for i in range(8)])

    def at(self, r):
        """State в раунде r."""
        return self.states[r]

    def section(self, r):
        """Срез: state + word на раунде r."""
        return self.states[r], self.schedule[r] if r < len(self.schedule) else None


class Mark:
    """Метка на fabric — точка интереса."""
    __slots__ = ['path_idx', 'round', 'note', 'value']
    def __init__(self, path_idx, round, note="", value=0):
        self.path_idx = path_idx
        self.round = round
        self.note = note
        self.value = value
    def __repr__(self): return f"Mark(p={self.path_idx},r={self.round},{self.note})"


class Fabric:
    """Множество Path — ткань пространства."""
    def __init__(self):
        self.paths: List[Path] = []
        self.marks: List[Mark] = []

    def weave(self, message, n_rounds=64):
        """Добавить path в ткань."""
        p = Path(message, n_rounds)
        idx = len(self.paths)
        self.paths.append(p)
        return idx

    def mark(self, path_idx, round, note="", value=0):
        m = Mark(path_idx, round, note, value)
        self.marks.append(m)
        return m


# ═══════════════════════════════════════════════════════════════
# DIMENSION OPERATIONS
# ═══════════════════════════════════════════════════════════════

def overlay(path1: Path, path2: Path, r=None) -> Delta:
    """Наложить два пути, измерить δ."""
    if r is not None:
        return path1.at(r).delta(path2.at(r))
    return path1.final_hash.delta(path2.final_hash)


def repair_a(state1: State, state2: State, target_state1_next: State, r: int) -> Word:
    """A-repair: найти W[r] для path2, чтобы a2_new = a1_new."""
    a2,b2,c2,d2,e2,f2,g2,h2 = state2.regs
    a1_new = target_state1_next.a
    T2_2 = add32(Sigma0(a2), Maj(a2,b2,c2))
    T1_needed = sub32(a1_new, T2_2)
    W_needed = sub32(T1_needed, add32(h2, Sigma1(e2), Ch(e2,f2,g2), K[r]))
    return Word(W_needed)


def repair_e(state1: State, state2: State, target_e_next, r: int) -> Word:
    """E-repair: найти W[r] для path2, чтобы e2_new = target."""
    a2,b2,c2,d2,e2,f2,g2,h2 = state2.regs
    T1_needed = sub32(target_e_next, d2)
    W_needed = sub32(T1_needed, add32(h2, Sigma1(e2), Ch(e2,f2,g2), K[r]))
    return Word(W_needed)


def propagate_one(state: State, w: Word, r: int) -> State:
    """Один раунд: State + Word → новый State."""
    a,b,c,d,e,f,g,h = state.regs
    T1 = add32(h, Sigma1(e), Ch(e,f,g), K[r], w.val)
    T2 = add32(Sigma0(a), Maj(a,b,c))
    return State((add32(T1,T2), a, b, c, add32(d,T1), e, f, g))


def absorb(message: Message, word_idx: int, bit: int, n_rounds: int) -> int:
    """Измерить поглощение: сколько бит выхода меняет 1 бит входа."""
    p1 = Path(message, n_rounds)
    m2 = message.copy()
    m2[word_idx] = m2[word_idx].val ^ (1 << bit)
    p2 = Path(m2, n_rounds)
    return overlay(p1, p2).weight


def schedule_delta(msg1: Message, msg2: Message) -> List[int]:
    """δ в expanded schedule."""
    s1 = msg1.expand()
    s2 = msg2.expand()
    return [hw(s1[i].val ^ s2[i].val) for i in range(64)]


# ═══════════════════════════════════════════════════════════════
# NEW OPERATIONS: неисследованные
# ═══════════════════════════════════════════════════════════════

def warp(path: Path, transform_fn, rounds=range(64)):
    """Деформировать path: применить transform к state на каждом раунде.
    transform_fn(state, round) → modified_state.
    Что это делает? Мы не знаем. Тестируем."""
    warped_states = [path.at(0)]
    for r in rounds:
        s = path.at(r + 1)
        s_warped = transform_fn(s, r)
        warped_states.append(s_warped)
    return warped_states


def resonate(paths: List[Path], r: int) -> Dict:
    """Найти резонанс: какие пары путей максимально БЛИЗКИ на раунде r?
    Возвращает карту расстояний."""
    n = len(paths)
    distances = {}
    min_dist = 256
    min_pair = None
    for i in range(n):
        for j in range(i+1, n):
            d = overlay(paths[i], paths[j], r).weight
            distances[(i,j)] = d
            if d < min_dist:
                min_dist = d
                min_pair = (i,j)
    return {'distances': distances, 'min_dist': min_dist, 'min_pair': min_pair}


def tunnel(path: Path, from_r: int, to_r: int) -> Delta:
    """Туннель: как меняется δ между from_r и to_r для ПАРЫ путей?
    Не вычисляя промежуточные раунды — прямое отображение.
    Это НОВАЯ операция: skip rounds."""
    # В стандартной физике нельзя "пропустить" раунды.
    # Но что если между from_r и to_r есть ПРЯМОЕ отображение?
    # Тестируем: можно ли ПРЕДСКАЗАТЬ state[to_r] из state[from_r] без промежуточных?
    return path.at(from_r).delta(path.at(to_r))


def invert_round(state_after: State, w: Word, r: int) -> State:
    """Обратить 1 раунд: из state[r+1] и W[r] получить state[r]."""
    a_new, b_new, c_new, d_new, e_new, f_new, g_new, h_new = state_after.regs
    # b_new = a_old, c_new = b_old, d_new = c_old
    # f_new = e_old, g_new = f_old, h_new = g_old
    a_old = b_new
    b_old = c_new
    c_old = d_new
    e_old = f_new
    f_old = g_new
    g_old = h_new

    # T1 = a_new - T2, where T2 = Sigma0(a_old) + Maj(a_old,b_old,c_old)
    T2 = add32(Sigma0(a_old), Maj(a_old, b_old, c_old))
    T1 = sub32(a_new, T2)

    # d_old: e_new = d_old + T1 → d_old = e_new - T1
    d_old = sub32(e_new, T1)

    # h_old: T1 = h_old + Sigma1(e_old) + Ch(e_old,f_old,g_old) + K[r] + W[r]
    h_old = sub32(T1, add32(Sigma1(e_old), Ch(e_old,f_old,g_old), K[r], w.val))

    return State((a_old, b_old, c_old, d_old, e_old, f_old, g_old, h_old))


# ═══════════════════════════════════════════════════════════════
# COMPOSITE OPERATIONS (programs in DimLang)
# ═══════════════════════════════════════════════════════════════

def dual_repair(msg_base: Message, delta_w0: int, strategy='a') -> Tuple[Message, int]:
    """
    Dual repair: начать с δW[0], repair через r=1..15.
    strategy='a': repair a-node (pipe cascade → e converges)
    strategy='e': repair e-node (a diverges)
    strategy='alternate': чередовать a и e repair
    strategy='smart': repair тот node, у которого больше δ
    Returns: (modified message, final delta weight)
    """
    msg_prime = msg_base.copy()
    msg_prime[0] = msg_base[0].val ^ delta_w0

    p1 = Path(msg_base)
    s1 = State()
    s2 = State()

    # Round 0: keep original
    s1 = propagate_one(s1, Word(msg_base[0].val), 0)
    s2 = propagate_one(s2, Word(msg_prime[0].val), 0)

    for r in range(1, 16):
        s1_next = propagate_one(s1, Word(msg_base[r].val), r)

        if strategy == 'a':
            w_needed = repair_a(s1, s2, s1_next, r)
        elif strategy == 'e':
            w_needed = repair_e(s1, s2, s1_next.e, r)
        elif strategy == 'alternate':
            if r % 2 == 0:
                w_needed = repair_a(s1, s2, s1_next, r)
            else:
                w_needed = repair_e(s1, s2, s1_next.e, r)
        elif strategy == 'smart':
            da = hw(s1.a ^ s2.a)
            de = hw(s1.e ^ s2.e)
            if da >= de:
                w_needed = repair_a(s1, s2, s1_next, r)
            else:
                w_needed = repair_e(s1, s2, s1_next.e, r)
        else:
            w_needed = Word(msg_base[r].val)

        msg_prime[r] = w_needed.val
        s2 = propagate_one(s2, w_needed, r)
        s1 = s1_next

    p_prime = Path(msg_prime)
    return msg_prime, overlay(p1, p_prime).weight


def fabric_search(n_paths: int, n_rounds: int = 64) -> Fabric:
    """Создать ткань из N случайных путей, искать близкие пары."""
    fab = Fabric()
    for _ in range(n_paths):
        msg = Message()
        fab.weave(msg, n_rounds)
    return fab


# ═══════════════════════════════════════════════════════════════
# EXPORT for testing
# ═══════════════════════════════════════════════════════════════

__all__ = [
    'State', 'Delta', 'Word', 'Message', 'Path', 'Fabric', 'Mark',
    'overlay', 'repair_a', 'repair_e', 'propagate_one', 'absorb',
    'schedule_delta', 'warp', 'resonate', 'tunnel', 'invert_round',
    'dual_repair', 'fabric_search',
]
