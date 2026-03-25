#!/usr/bin/env python3
"""
SCF: ТЕОРИЯ ПРЕДСТАВЛЕНИЙ — SHA-256 как действие группы.

Идея: SHA round function действует на state Z_{2^32}^8.
Множество всех round functions (параметризованных W[r], K[r])
образует ПОЛУГРУППУ. Коллизия = два элемента полугруппы
с одинаковым действием на IV.

Конкретнее: SHA(W) = f_{W[63]} ∘ ... ∘ f_{W[0]} (IV)
где f_w(state) = round(state, w, K[r]).

Collision: f_{W1}(IV) = f_{W2}(IV), W1 ≠ W2.

Вопросы теории представлений:
1. Какова ГРУППА СИММЕТРИЙ round function?
   Есть ли g: state→state с f(g(s)) = g(f(s)) (коммутирует)?
2. Есть ли ИНВАРИАНТЫ: I(state) такие что I(f(state)) = I(state)?
3. Есть ли ОРБИТЫ: подмножества state, замкнутые под f?

Если найдём инвариант I: collision → I(SHA(W1)) = I(SHA(W2))
автоматически. Поиск в пространстве инвариантов ДЕШЕВЛЕ.
"""
import os, sys, math

MASK32 = 0xFFFFFFFF
K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK32
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def add32(*args):
    s=0
    for x in args: s=(s+x)&MASK32
    return s
def sub32(a,b): return (a-b)&MASK32
def hw(x): return bin(x).count('1')

def expand_real(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def sha_compress(W16, iv=None):
    if iv is None: iv=IV
    W=expand_real(W16); s=list(iv)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return [add32(iv[i],s[i]) for i in range(8)]

def sha_round_one(state, w, k):
    a,b,c,d,e,f,g,h = state
    T1=add32(h,Sig1(e),Ch(e,f,g),k,w)
    T2=add32(Sig0(a),Maj(a,b,c))
    return [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]


# ============================================================
# EXP 1: ADDITIVE INVARIANTS — I(state) = Σ state[i] mod 2^k
# ============================================================
def exp1_additive_invariants(N):
    print("="*70)
    print("EXP 1: ADDITIVE INVARIANTS")
    print("  Test: I(state) = Σ state[i] (mod 2^k)")
    print("  Does I(round(state)) = I(state) + f(W,K)?")
    print("  If f(W,K) depends ONLY on W,K (not state) → invariant!")
    print("="*70)

    for k in [1, 2, 4, 8, 32]:
        mod = 1 << k
        mask = mod - 1

        # Test: is Σstate[i] mod 2^k preserved?
        # I_after - I_before should be constant for fixed W,K
        diffs = set()
        consistent = True

        W_fixed = int.from_bytes(os.urandom(4),'big')
        K_fixed = K[0]

        for trial in range(min(N, 200)):
            state = [int.from_bytes(os.urandom(4),'big') for _ in range(8)]
            I_before = sum(state) & mask

            state_after = sha_round_one(state, W_fixed, K_fixed)
            I_after = sum(state_after) & mask

            diff = (I_after - I_before) & mask
            diffs.add(diff)

            if len(diffs) > 1:
                consistent = False
                break

        if consistent and len(diffs) == 1:
            print(f"  mod 2^{k}: Σstate INVARIANT! diff = {diffs.pop()} (constant)")
        else:
            print(f"  mod 2^{k}: NOT invariant ({len(diffs)} distinct diffs)")

    # Test weighted sums
    print(f"\n  Weighted sums: I = Σ c_i·state[i]")
    for weights in [(1,1,1,1,1,1,1,1), (1,0,0,0,1,0,0,0), (1,-1,1,-1,1,-1,1,-1),
                    (1,1,1,1,-1,-1,-1,-1), (0,0,0,1,0,0,0,1)]:
        diffs = set()
        for trial in range(200):
            state = [int.from_bytes(os.urandom(4),'big') for _ in range(8)]
            I_before = sum(w*s for w,s in zip(weights, state)) & MASK32
            state_after = sha_round_one(state, int.from_bytes(os.urandom(4),'big'), K[0])
            I_after = sum(w*s for w,s in zip(weights, state_after)) & MASK32
            diffs.add((I_after - I_before) & MASK32)
            if len(diffs) > 3: break

        w_str = ",".join(f"{w:+d}" for w in weights)
        is_inv = len(diffs) == 1
        print(f"    [{w_str}]: {'INVARIANT ★' if is_inv else f'{len(diffs)} diffs'}")


# ============================================================
# EXP 2: XOR INVARIANTS — I(state) = ⊕ state[i] or subsets
# ============================================================
def exp2_xor_invariants(N):
    print("\n" + "="*70)
    print("EXP 2: XOR INVARIANTS")
    print("  I = state[i] ⊕ state[j] ⊕ ... for subsets")
    print("="*70)

    W_fixed = int.from_bytes(os.urandom(4),'big')

    # Test all single and pair XORs
    print(f"  Single register: I = state[i]")
    for i in range(8):
        diffs = set()
        for trial in range(200):
            state = [int.from_bytes(os.urandom(4),'big') for _ in range(8)]
            I_before = state[i]
            after = sha_round_one(state, W_fixed, K[0])
            I_after = after[i]
            diffs.add(I_before ^ I_after)
            if len(diffs) > 3: break
        print(f"    state[{i}]: {'INVARIANT ★' if len(diffs)==1 else f'{len(diffs)} diffs'}")

    # Shift register: state[1] = state_prev[0]
    # So: state_after[1] ⊕ state_before[0] should be 0!
    print(f"\n  Shift register invariants:")
    shifts = [(1,0,'b=a'), (2,1,'c=b'), (3,2,'d=c'),
              (5,4,'f=e'), (6,5,'g=f'), (7,6,'h=g')]
    for after_reg, before_reg, name in shifts:
        all_zero = True
        for trial in range(200):
            state = [int.from_bytes(os.urandom(4),'big') for _ in range(8)]
            after = sha_round_one(state, W_fixed, K[0])
            if after[after_reg] != state[before_reg]:
                all_zero = False; break
        print(f"    {name}: {'EXACT INVARIANT ★★★' if all_zero else 'broken'}")


# ============================================================
# EXP 3: FUNCTIONAL INVARIANTS — nonlinear functions of state
# ============================================================
def exp3_functional(N):
    print("\n" + "="*70)
    print("EXP 3: FUNCTIONAL INVARIANTS")
    print("  Test nonlinear functions: I = f(state) preserved by round")
    print("="*70)

    # I1: state[0] + state[4] (a + e)
    # After round: a_new = T1+T2, e_new = d+T1
    # a_new + e_new = T1+T2+d+T1 = 2·T1+T2+d
    # Before: a+e
    # Diff = 2·T1+T2+d - a - e = 2·T1 + (Σ0(a)+Maj) + d - a - e
    # This depends on state → not invariant. But mod 2?
    # Mod 2: 2·T1 = 0, so diff mod 2 = T2+d-a-e = Σ0(a)+Maj(a,b,c)+d-a-e

    for name, func in [
        ("a+e mod 2^32", lambda s: add32(s[0], s[4])),
        ("a⊕e", lambda s: s[0] ^ s[4]),
        ("a+b+c+d mod 2^32", lambda s: add32(s[0],s[1],s[2],s[3])),
        ("e+f+g+h mod 2^32", lambda s: add32(s[4],s[5],s[6],s[7])),
        ("Σ_all mod 2^32", lambda s: add32(*s)),
        ("a·e mod 2^32", lambda s: (s[0]*s[4])&MASK32),
        ("parity(a)⊕parity(e)", lambda s: (hw(s[0])^hw(s[4]))&1),
    ]:
        diffs = set()
        for trial in range(200):
            state = [int.from_bytes(os.urandom(4),'big') for _ in range(8)]
            W_r = int.from_bytes(os.urandom(4),'big')
            I_before = func(state)
            after = sha_round_one(state, W_r, K[0])
            I_after = func(after)
            diff = sub32(I_after, I_before) if isinstance(I_before, int) and I_before > 1 else I_after ^ I_before
            diffs.add(diff & MASK32)
            if len(diffs) > 5: break

        is_inv = len(diffs) == 1
        print(f"  {name:30s}: {'INVARIANT ★' if is_inv else f'{min(len(diffs),6)}+ diffs'}")


# ============================================================
# EXP 4: COLLISION INVARIANT — properties shared by collisions
# ============================================================
def exp4_collision_properties(N):
    print("\n" + "="*70)
    print("EXP 4: COLLISION INVARIANT SEARCH")
    print("  If collisions have property P → search in P-space")
    print("  Test on REDUCED-ROUND collisions (findable)")
    print("="*70)

    # Find reduced-round collisions (R=16 via SAT or birthday)
    # For R=1: trivial — T1 determines everything
    # Use R=2-3 where birthday is feasible

    for R in [1, 2]:
        print(f"\n  R={R} rounds:")

        def sha_R(W16, R_rounds=R):
            W=expand_real(W16); s=list(IV)
            for r in range(R_rounds):
                a,b,c,d,e,f,g,h=s
                T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
                T2=add32(Sig0(a),Maj(a,b,c))
                s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
            return [add32(IV[i],s[i]) for i in range(8)]

        # Birthday search for R-round collision
        seen = {}
        collision_pairs = []

        for trial in range(min(N * 1000, 100000)):
            W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            H = sha_R(W)
            H_key = tuple(H)

            if H_key in seen:
                W_prev = seen[H_key]
                if W_prev != W:
                    collision_pairs.append((W_prev, W))
                    if len(collision_pairs) >= 5:
                        break
            else:
                seen[H_key] = W

        print(f"    Found {len(collision_pairs)} collisions in {trial+1} tries")

        if collision_pairs:
            # Analyze properties of collision pairs
            print(f"    Properties of collision pairs:")
            for idx, (Wa, Wb) in enumerate(collision_pairs[:3]):
                dW = [Wa[i] ^ Wb[i] for i in range(16)]
                n_diff = sum(1 for d in dW if d != 0)
                hw_total = sum(hw(d) for d in dW)
                sum_a = add32(*Wa) & MASK32
                sum_b = add32(*Wb) & MASK32
                sum_diff = sub32(sum_a, sum_b)

                print(f"      Pair {idx}: {n_diff} words differ, "
                      f"HW(δW)={hw_total}, "
                      f"Σ(W) diff=0x{sum_diff:08x}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    exp1_additive_invariants(N)
    exp2_xor_invariants(N)
    exp3_functional(N)
    exp4_collision_properties(N)
