# CARRY THEORY — Новый раздел математики

## Аксиоматика
Carry system CS(B) = (G, B, φ, c) where:
- G = digit group, B = base, φ = digit addition, c = carry function
- C1: φ(a,b) + B·c(a,b,0) = a+b in Z
- C2: c_i = c(digit_i(a), digit_i(b), c_{i-1}), c_{-1} = 0
- C3: c uniquely determined by (G, B)

## Fundamental Constants

η(B) = (3·log_B(3) - 4)/4   for B=2: η = (3ln3-4ln2)/(4ln2) = 0.18872
Δ(B) = Σ η_i               for B=2: Δ = 1/(6ln2) = 0.24044
ξ(B) = 1/ln(B/(B-1))        for B=2: ξ = 1/ln2 = 1.44270

## Carry Algebra (CA1-CA5)
- Carry ≠ group (context-dependent)
- Carry ≠ module (nonlinear in both arguments)
- Carry-path non-unique: C_{k-1} Catalan paths for sum of k operands
  T1 of SHA-256: 14 carry-paths for 5 operands
- Carry-path entropy: log₂ C_{k-1} = 3.81 bits for T1

## Carry Geometry (CG1-CG2)
- Carry ball: B_r(a) = {b : d_carry(a,b) ≤ r}
- B_0(a) = {b : a AND b = 0}, |B_0| = 2^{n-HW(a)}
- Carry space has NEGATIVE curvature (balls smaller than flat)

## Carry Number Theory (CN1-CN3)
- Carry-multiplicity μ(a) = E_b[d_carry(a,b)]
- Exact formula: μ(a) = Σ_j (1 - 2^{-(g_j+1)}) where g_j = gap after j-th one
- μ_max = k (all ones at start), μ_min = k(1-2^{-n/k}) (uniform)
- For HW=16, n=32: μ ≈ 12

## Carry Transform (CTr1-CTr3)
- CT[f](a) = Σ_x f(x)·(-1)^{overflow(a,x)}
- Nonlinear analog of Walsh-Hadamard
- CT-spectrum of SHA-256: FLAT (carry not correlated through 64 rounds)
- CT useful for shallow functions (few addition layers)

## Carry Complexity (CC1-CC2, CC-Open)
- CARRY-Search ∈ P (bit-by-bit reconstruction, O(n))
- CARRY-Collision ∈ NP (verify in poly time)
- OPEN: CARRY-Collision NP-hard?

## SHA-512 spectral invariant: RETRACTED as attack
144 common eigenvalues exist (structural fact), but monodromy = linear
approximation, carry = independent nonlinear part. Spectral filtering
doesn't reduce carry work. Collision cost = 2^256 (birthday).
