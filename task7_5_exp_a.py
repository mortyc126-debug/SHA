#!/usr/bin/env python3
"""
ЗАДАНИЕ 7.5 Experiment A: Algebraic degree estimation via cube sums
"""
import struct, os, time, random

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

def sha256_hash(M):
    W = list(M[:16])
    for i in range(16,64):
        W.append(add(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    a,b,c,d,e,f,g,h = H0
    for r in range(64):
        T1=add(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add(Sig0(a),Maj(a,b,c))
        h,g,f,e,d,c,b,a = g,f,e,add(d,T1),c,b,a,add(T1,T2)
    return tuple(add(s,iv) for s,iv in zip((a,b,c,d,e,f,g,h), H0))


def cube_sum(base_msg, cube_indices, out_word, out_bit):
    """Compute XOR sum of output bit over a cube (affine subspace).
    cube_indices: list of (word_idx, bit_idx) pairs — cube variables.
    base_msg: the "offset" point.
    Returns 0 or 1."""
    dim = len(cube_indices)
    xor_sum = 0
    base = list(base_msg)
    for mask in range(1 << dim):
        M = list(base)
        for j in range(dim):
            if mask & (1 << j):
                w, b = cube_indices[j]
                M[w] ^= (1 << b)
        H = sha256_hash(M)
        xor_sum ^= (H[out_word] >> out_bit) & 1
    return xor_sum


def main():
    print("="*70)
    print("EXPERIMENT A: Algebraic degree estimation (cube sums)")
    print("="*70)
    t0 = time.time()

    # For each dimension, test if cube sum = 0 for all cubes
    # If all zero at dim d → degree < d
    # If any nonzero at dim d → degree >= d

    # Test output bits: H[0][0], H[0][1], H[0][15], H[0][31], H[7][0], H[7][31]
    test_outputs = [(0,0), (0,1), (0,15), (0,31), (7,0), (7,31)]

    # All 512 input bit indices
    all_inputs = [(w, b) for w in range(16) for b in range(32)]

    # Dimensions to test (limit by computation: 2^dim evaluations per cube)
    dims = [1, 2, 4, 8, 10, 12, 14, 16]
    n_cubes = 50  # cubes per dimension

    for out_word, out_bit in test_outputs:
        print(f"\n  Output bit H[{out_word}][{out_bit}]:")
        degree_lower = 0  # degree >= degree_lower

        for dim in dims:
            t1 = time.time()
            nonzero_found = False
            for trial in range(n_cubes):
                base = list(struct.unpack('>16I', os.urandom(64)))
                cube_idx = random.sample(all_inputs, dim)
                s = cube_sum(base, cube_idx, out_word, out_bit)
                if s != 0:
                    nonzero_found = True
                    break
            elapsed_dim = time.time() - t1
            if nonzero_found:
                degree_lower = dim
                status = f"degree >= {dim} (nonzero at trial {trial+1})"
            else:
                status = f"all {n_cubes} cubes zero → degree < {dim}?"
            print(f"    dim={dim:2d}: {status}  [{elapsed_dim:.1f}s]")

            if not nonzero_found:
                # If all cubes are zero, degree < dim, stop
                print(f"    → Estimated degree: [{degree_lower}, {dim})")
                break

            # Time guard
            if elapsed_dim > 30:
                print(f"    (stopping: too slow for higher dims)")
                print(f"    → Degree >= {degree_lower} (lower bound only)")
                break
        else:
            print(f"    → Degree >= {degree_lower} (didn't find upper bound in tested dims)")

    elapsed = time.time() - t0
    print(f"\n  Total time: {elapsed:.1f}s")

if __name__ == "__main__":
    main()
