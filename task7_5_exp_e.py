#!/usr/bin/env python3
"""
ЗАДАНИЕ 7.5: Bit-Slice Complexity
Experiment E: Carry structure of SHA-256 output — sensitivity bit 0 vs bit 31
Priority: HIGHEST
"""
import struct, os, time

MASK = 0xFFFFFFFF
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
H0 = (0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19)

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def Ch(e,f,g): return (e&f)^((~e)&g)&MASK
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def add(*args):
    s=0
    for x in args: s=(s+x)&MASK
    return s
def rw(n): return list(struct.unpack(f'>{n}I', os.urandom(4*n)))

def sha256_compress(M):
    """Full SHA-256 compression: returns (state64, hash)."""
    W = list(M[:16])
    for i in range(16,64):
        W.append(add(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    a,b,c,d,e,f,g,h = H0
    for r in range(64):
        T1=add(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add(Sig0(a),Maj(a,b,c))
        h,g,f,e,d,c,b,a = g,f,e,add(d,T1),c,b,a,add(T1,T2)
    state64 = (a,b,c,d,e,f,g,h)
    hash_out = tuple(add(state64[i],H0[i]) for i in range(8))
    return state64, hash_out

def getbit(word, bit):
    return (word >> bit) & 1

def flipbit_msg(M, word_idx, bit_idx):
    """Flip bit bit_idx of word word_idx in message M."""
    M2 = list(M)
    M2[word_idx] ^= (1 << bit_idx)
    return M2

def main():
    print("="*70)
    print("EXPERIMENT E: Carry structure — sensitivity bit 0 vs bit 31")
    print("="*70)
    t0 = time.time()

    # E1: Verify analytical carry structure
    print("\n  E1: Verify H[i][0] = state64[i][0] XOR IV[i][0] (no carry)")
    ok = 0
    for _ in range(10000):
        M = rw(16)
        state64, hash_out = sha256_compress(M)
        all_match = True
        for i in range(8):
            h_bit0 = getbit(hash_out[i], 0)
            s_bit0 = getbit(state64[i], 0)
            iv_bit0 = getbit(H0[i], 0)
            expected = s_bit0 ^ iv_bit0
            if h_bit0 != expected:
                all_match = False
                break
        if all_match:
            ok += 1
    print(f"    H[i][0] = state64[i][0] XOR IV[i][0]: {ok}/10000")

    # E2: Verify H[i][1] has carry dependency (degree 2)
    print("\n  E2: Verify H[i][1] involves carry from bit 0")
    ok = 0
    for _ in range(10000):
        M = rw(16)
        state64, hash_out = sha256_compress(M)
        all_match = True
        for i in range(8):
            s0 = getbit(state64[i], 0)
            s1 = getbit(state64[i], 1)
            iv0 = getbit(H0[i], 0)
            iv1 = getbit(H0[i], 1)
            carry0 = s0 & iv0  # carry from bit 0
            expected = s1 ^ iv1 ^ carry0
            if getbit(hash_out[i], 1) != expected:
                all_match = False
                break
        if all_match:
            ok += 1
    print(f"    H[i][1] = s[1] XOR iv[1] XOR (s[0] AND iv[0]): {ok}/10000")

    # E3: Sensitivity of state64[i][bit_j] to each input bit
    # For each output word i and bit position j, flip each input bit
    # and measure how often state64[i][j] changes.
    print("\n  E3: Sensitivity of state64[i][j] to input bits")
    print("      (how many of 512 input bits affect each output bit)")

    n_samples = 200  # per input bit
    # Test representative output bits: word 0, bits 0,1,7,15,31
    test_bits = [(0,0),(0,1),(0,7),(0,15),(0,31),
                 (3,0),(3,31),(7,0),(7,31)]

    for out_word, out_bit in test_bits:
        sensitive_count = 0
        for in_word in range(16):
            for in_bit in range(32):
                flips = 0
                for _ in range(n_samples):
                    M = rw(16)
                    _, h1 = sha256_compress(M)
                    M2 = flipbit_msg(M, in_word, in_bit)
                    _, h2 = sha256_compress(M2)
                    if getbit(h1[out_word], out_bit) != getbit(h2[out_word], out_bit):
                        flips += 1
                if flips > n_samples * 0.1:  # >10% flip rate → sensitive
                    sensitive_count += 1
        print(f"    H[{out_word}][{out_bit:2d}]: sensitivity = {sensitive_count}/512"
              f"  ({'FULL' if sensitive_count == 512 else 'REDUCED!'})")

    # E4: More detailed — flip rate distribution for bit 0 vs bit 31
    print("\n  E4: Flip rate distribution (expected ~50% for each input bit)")
    n_detail = 500
    for out_word, out_bit in [(0,0), (0,31)]:
        rates = []
        for in_word in range(16):
            for in_bit in range(32):
                flips = 0
                for _ in range(n_detail):
                    M = rw(16)
                    _, h1 = sha256_compress(M)
                    M2 = flipbit_msg(M, in_word, in_bit)
                    _, h2 = sha256_compress(M2)
                    if getbit(h1[out_word], out_bit) != getbit(h2[out_word], out_bit):
                        flips += 1
                rates.append(flips / n_detail)
        avg_rate = sum(rates) / len(rates)
        min_rate = min(rates)
        max_rate = max(rates)
        # Count how many are significantly different from 0.5
        biased = sum(1 for r in rates if abs(r - 0.5) > 0.05)
        print(f"    H[{out_word}][{out_bit:2d}]: avg flip rate = {avg_rate:.4f}, "
              f"min={min_rate:.4f}, max={max_rate:.4f}, "
              f"biased (>5% from 0.5): {biased}/512")

    elapsed = time.time() - t0
    print(f"\n  Total time: {elapsed:.1f}s")

if __name__ == "__main__":
    main()
