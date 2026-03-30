# Carry Operator, Carry Calculus, Carry Homotopy

## I. Carry Operator C(a,b) = (a+b - a⊕b)/2

C(a,0) = 0, C(a,a) = a, C(a,~a) = 0, C(a,b) = C(b,a)
|Z(C)| = 3^n (carry-free pairs, P = (3/4)^n)
|Im(C_a)| = 2^{HW(a)}, |fiber| = 2^{n-HW(a)}
C is idempotent: C(C(a,b), C(a,b)) = C(a,b)

## II. Carry Calculus

∂C/∂a_k: lower-triangular matrix
  [k] = b_k ⊕ c_{k-1} (diagonal)
  [k+m] decays as (1/2)^m (sub-diagonal, carry self-healing)
det(∇_a C) = ∏(b_i ⊕ c_{i-1}), P(det ≠ 0) = 2^{-n}

Carry-Laplacian: Δ_C f(a) = Σ_k f(a⊕e_k) · ∂C/∂a_k
Carry-harmonic: Δ_C f = 0
|Harm_C| ≈ 2^{n/2} — COINCIDES with birthday bound!

## III. Carry Homotopy

Carry-graph G_C: vertices = {0,1}^n, edges = carry-free pairs
Connected, diameter 2 (via vertex 0)

Carry-complex K_C: simplicial complex of carry-free cliques
Max simplex: {0, e_0, ..., e_{n-1}} (dimension n)
K_C CONTRACTIBLE (0 is cone point: C(0,a) = 0 ∀a)
H_k(K_C) = 0 for all k ≥ 1 (trivial topology)
