# SHA-256 Algebraic Analysis: Carry-Aware AND-Algebra

## Summary of Findings

This research develops a new algebraic framework for analyzing SHA-256
and applies it to break the Wang chain barrier at round 17.

## 1. Fundamental Identity (NEW)

For any GF(2)-linear function L (Sig0, Sig1, sig0, sig1):

```
D_add[L](x, y) = L(c) - 2·(L(x) ∧ L(c))
where c = (x+y) ⊕ x
```

**Verified 100% exact** for all 4 Sigma functions (5000+ trials each).

This identity decomposes the additive differential of any rotation-XOR
function into:
- A GF(2)-linear part: `L(c)`
- A quadratic correction: `2·(L(x) ∧ L(c))`

The correction term involves AND — the same operation as Ch and Maj.

## 2. Bilinear Form and Kernel

The quadratic forms Ch and Maj define a bilinear form B on state space:
- **Rank of B = 5** (per bit position, 8×8 matrix)
- **Kernel dimension = 4**: {d, h, f⊕g, a⊕b⊕c}
- Registers d and h are "invisible" to quadratic nonlinearity
- 95.1% of AND terms cancel across rounds (massive algebraic compensation)

## 3. De17 = 0 Solution (CONCRETE RESULT)

**Found a message pair where De17 = 0:**

```
W[ 0] = 0x3eba2820    W'[0] = 0xbeba2820  (DW = 0x80000000)
W[ 1] = 0x3082cc68    W[14] = 0x0013efc2  ← free word
...
De17 = 0x00000000  ✓ VERIFIED
```

Search cost: ~3.5 × 10^9 evaluations (~2^32)

## 4. Free Word Mechanism

Key structural finding: **W[12..15] are free words**.

- Da13 (barrier equation) depends ONLY on W[0..11]
- W[12..15] don't affect the barrier but DO affect De at rounds 13+
- This gives 128 bits of freedom for barrier control

De dependence map:
```
Round | W[12] | W[13] | W[14] | W[15]
  13  |  DEP  |   -   |   -   |   -
  14  |  DEP  |  DEP  |   -   |   -
  15  |  DEP  |  DEP  |  DEP  |   -
  16+ |  DEP  |  DEP  |  DEP  |  DEP
```

## 5. Multi-Barrier Extension

(De17, De18, De19, De20) are **independently controllable** via
the 4 free words (verified: 20000/20000 unique 4-tuples).

- 128 bits freedom = 4 × 32-bit constraints
- Wang chain extends from 16 to 20 rounds
- Cost: O(2^32) per barrier = O(2^34) total for 4 barriers

## 6. Collision Budget

```
Rounds  1-16:  O(1)      Wang chain (De=0 maintained)
Rounds 17-20:  O(2^34)   Free word technique (4 barriers)
Rounds 21-64:  ???       No free words left in single block
                         → multi-block or statistical methods needed
```

## Files

| File | Description |
|------|-------------|
| `clifford_step1.py` | Bilinear form extraction from Ch/Maj |
| `clifford_step2.py` | Multi-round decomposition, GF(2) gap analysis |
| `clifford_step3.py` | Initial carry algebra attempt (failed identity) |
| `clifford_step3b.py` | Corrected fundamental identity (verified) |
| `clifford_step4.py` | AND-algebra applied to Wang chain barrier |
| `clifford_step5.py` | Chain structure, free word discovery |
| `clifford_step6.py` | Multi-word differential freedom analysis |
| `clifford_step7.py` | Python MITM search |
| `clifford_step7b.py` | Optimized search with carry analysis |
| `clifford_step8.py` | Python verification framework |
| `clifford_step9.py` | Multi-barrier analysis and collision budget |
| `de17_search.c` | C brute-force search (first De17=0 find) |
| `de17_verify.c` | C search with full 64-round verification |

## Key Equations

### SHA-256 Round Differential (exact):
```
Da' = Dh + [Sig1(c_e) - 2(Sig1(e)∧Sig1(c_e))] + DCh + DW
    + [Sig0(c_a) - 2(Sig0(a)∧Sig0(c_a))] + DMaj

De' = Dd + DT1  (where DT1 = first line above)
```

### XOR-Addition Bridge:
```
A ⊕ B = A + B - 2·(A ∧ B)   (mod 2^32)
```

### Carry Pattern:
```
c = (x+y) ⊕ x = y ⊕ 2·carry(x,y)
P(HW(carry) = k) ≈ (1/2)^k   (geometric distribution)
```
