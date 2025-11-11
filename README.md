\
# Stabilizer LDPC Randomizer — Paper-Style README (English)

**Author:** Kenta Kasai

---

## Abstract

This program implements a randomization algorithm for stabilizer LDPC matrices \(H_X, H_Z\). It repeatedly applies **2×2 cross-swaps** to \(H_X\) while maintaining the **stabilizer commutation condition**

\[
H_X H_Z^\top \oplus H_Z H_X^\top = 0 \quad (\text{over GF(2)}).
\]

Every proposed change of \(H_X\) is followed by a **local symmetric repair** of \(H_Z\) obtained by solving a **general stabilizer linear system** on a restricted domain. The implementation is fully sparse (set-based) and preserves both row and column weights.

---

## 1. Background and Objective

In quantum LDPC code design, one must construct sparse matrices \((H_X, H_Z)\) that satisfy the stabilizer commutation constraint. This program introduces controlled randomness into \(H_X\) via weight-preserving cross-swaps and then performs localized algebraic repairs of \(H_Z\) so that the pair continues to commute. The method directly enforces the general stabilizer constraint without special-structure assumptions.

### Contributions
- Local repair formulated and solved as a **general (symmetric) GF(2) linear system** on a small domain.
- Adoption of solutions that satisfy **row/column weight balance** constraints within the domain.
- **Exact DFS** fallback for very small systems.
- Fully sparse (set-of-indices) data structure for clarity and extensibility.

---

## 2. Preliminaries and Notation

- Field: GF(2). Addition is XOR; set operations use symmetric differences.
- A sparse binary matrix is represented by **row-support sets** \(\mathrm{rows}[i]\subseteq[n)\) and **column-support sets** \(\mathrm{cols}[j]\subseteq[m)\).
- Parity inner product:
  \[
  \langle a,b\rangle := |a\cap b| \bmod 2
  \]
  (implemented by `dot_parity`).
- Initial matrices use a tiled block-ID construction with block size \(P\), block rows \(d_c\), and block columns \(d_r\):
  \[
  m=d_c P,\quad n=d_r P,\quad \text{row weight}=d_r,\quad \text{column weight}=d_c.
  \]

**Stabilizer commutation.** A pair \((H_X,H_Z)\) is commuting iff
\[
H_X H_Z^\top \oplus H_Z H_X^\top = O_{m\times m},
\]
equivalently, for all \((k,i)\),
\[
\langle H_X[k,:], H_Z[i,:]\rangle \oplus \langle H_Z[k,:], H_X[i,:]\rangle = 0.
\]

---

## 3. Cross-Swap (2×2)

Given rows \(i_1,i_2\) and columns \(j_1,j_2\), the \(2\times2\) submatrix
\[
\begin{bmatrix} a & b \\ c & d \end{bmatrix}
\]
is *swappable* iff \((a,d)=(1,1)\) and \((b,c)=(0,0)\), or vice versa. Toggling all four entries then **preserves row/column weights** (`is_valid_cross_pattern`, `cross_swap_in_place`). We apply this only to \(H_X\) to diversify structure while keeping degrees fixed.

**Lemma (Weight preservation).** Each affected row/column loses two 1s and gains two 1s (or the reverse), hence degrees remain unchanged.

---

## 4. Local Repair Formulation (General Stabilizer System)

Let \(H'_X\) denote the swapped matrix. We seek a correction \(\Delta\) such that \(H'_Z := H_Z \oplus \Delta\) restores commutation. Over GF(2),
\[
H'_X (H_Z \oplus \Delta)^\top \oplus (H_Z \oplus \Delta) H'_X{}^\top = O,
\]
which linearizes to
\[
H'_X H_Z^\top \oplus H_Z H'_X{}^\top \oplus H'_X \Delta^\top \oplus \Delta H'_X{}^\top = O.
\]
Thus, \(\Delta\) must satisfy
\[
H'_X \Delta^\top \oplus \Delta H'_X{}^\top = H'_X H_Z^\top \oplus H_Z H'_X{}^\top.
\tag{1}
\]

We restrict \(\Delta\) to a **local domain** \(I\times J\) determined by the touched columns:
- \(I\): rows reachable via \(H_Z\),
- \(J\): columns reachable from \(I\) via \(H_Z\) or \(H'_X\),
- \(K\): rows reachable from \(J\) via \(H'_X\) or \(H_Z\) (allowing \(K\cap I\)).

The linear system (1) is assembled on \(I\times J\) (`assemble_general`).

**Balance constraints.** To preserve row/column weights within \(I\times J\),
\[
\sum_{j\in J} \sigma_{ij}=0,\qquad \sum_{i\in I} \sigma_{ij}=0,
\]
where \(\sigma_{ij}=+1\) for \(0\to1\) and \(-1\) for \(1\to0\). This is checked by `balanced_solution_from_vector`.

---

## 5. Solvers

- Solve (1) by sparse Gaussian elimination over GF(2) (`solve_gf2_sparse`).  
- Apply \(\Delta\) only if **balance** constraints are satisfied.  
- Accept the candidate \((H'_X,H'_Z)\) **only after global stabilizer validation**.  
- If the system is tiny (by default, \(\le 24\) variables), run an **exact DFS** (`exact_balanced_solver`) with balance constraints.

---

## 6. Algorithm (Pseudo-Code)

```text
Input: P, d_c, d_r, steps, seed
H_X ← block_id(P, d_c, d_r)
H_Z ← block_id(P, d_c, d_r)

for t = 1..steps:
  for trial = 1..max_trials:
    pick (i1,i2,j1,j2) with valid cross pattern on H_X
    H'_X ← swap(H_X, i1,i2,j1,j2)
    assert deg(H'_X) == deg(H_X)  // weight preservation

    (I,J,K) ← build_local_sets_general(H'_X, H_Z, touched={j1,j2})
    (A,b)   ← assemble_general(H'_X, H_Z, I,J,K)  // Eq. (1) on I×J
    x       ← solve_gf2_sparse(A,b)

    if x is None or not balanced(x):
       if |x| small: x ← exact_balanced_solver(A,b, balance)
       else: continue

    H'_Z ← apply_delta(H_Z, x)
    if stabilizer_ok(H'_X, H'_Z):
       (H_X, H_Z) ← (H'_X, H'_Z)   // commit & print
       break  // next t
  if not committed: report "no feasible swap"
```

**Correctness sketch.** (i) Cross-swaps preserve degrees of \(H_X\). (ii) Any \(\Delta\) solving (1) recovers commutation for the updated pair. (iii) Balance preserves local row/column weights of \(H_Z\). (iv) The final global check prevents inconsistent commits.

---

## 7. Complexity and Parameters

- **Global validation:** worst-case \(O(m^2 \bar{w})\) where \(\bar{w}\) is the average row weight.  
- **Local assembly:** \(O(|I|\cdot|J|)\).  
- **Sparse elimination:** practical up to a few hundred variables.  
- **DFS fallback:** exponential; default cap is 24 variables.

**Parameter guidelines.**
- `max_domain_size = |I|\cdot|J|` around 200–600 is usually practical.
- `dfs_var_cap` around 20–26 is reasonable; larger values may explode runtime.
- Start with small `steps` to limit output size.

---

## 8. Implementation Map (Key Functions)

- `SparseBinMat`: sparse matrix container (`add`, `toggle`, `copy`, `nnz`).
- `dot_parity`: parity inner product (iterate over the smaller set).
- `is_valid_cross_pattern` / `cross_swap_in_place` / `find_valid_cross_swap`: cross-swap utilities.
- `make_block_id_sparse`: block-ID initializer.
- `stabilizer_violation_positions`: global commutation checker.
- `build_local_sets_general`: extracts local domain \(I,J,K\).
- `assemble_general`: sparse assembly of (1).
- `solve_gf2_sparse`: sparse Gaussian elimination over GF(2).
- `balanced_solution_from_vector`: balance check for row/column weights.
- `exact_balanced_solver`: exact DFS with balance constraints.
- `adaptive_local_repair`: end-to-end local repair pipeline (general system only).
- `repeat_swaps_and_repairs_safe`: main loop (propose → repair → verify → commit).
- `print_sparse_01`: full 0/1 matrix printing for debugging.

---

## 9. Usage

```bash
python3 stab-randomizer.py
```
Default parameters in the script:
- `P=25, d_c=3, d_r=4, steps=2000, seed=42, max_trials=300`
- In `adaptive_local_repair`: `max_domain_size=400, dfs_var_cap=24`

Use **Python 3.10+**. If you include `from __future__ import annotations`, place it at the **top of the file**.

---

## 10. Limitations and Future Work

- Acceptance and convergence rates are empirical at present.
- The neighborhood design (\(I,J,K\)) can be further optimized.
- Code distance, rank, and degeneracy evaluation are not yet implemented.
- Parallelizing swap search and validation is left as future work.

---


## License
Intended for research and educational release. Please choose and state an open-source license (e.g., MIT, BSD-3, Apache-2.0) when publishing.

## Acknowledgment
Prepared in a paper-style format with equations and a brief correctness sketch for academic reuse.
