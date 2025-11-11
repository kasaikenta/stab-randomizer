# Stabilizer LDPC Randomizer — Technical README (English, No LaTeX)

**Author:** Kenta Kasai

---

## Abstract

This program implements a randomization algorithm for stabilizer LDPC matrices `HX` and `HZ`. It repeatedly applies 2x2 cross‑swaps to `HX` and performs a local symmetric repair of `HZ` so that the global stabilizer commutation constraint remains satisfied:

```
HX HZ^T XOR HZ HX^T = 0    (over GF(2))
```

All data structures are fully sparse (set‑based). The procedure preserves row and column weights by construction and validates the stabilizer constraint before committing any change.

---

## 1. Background and Objective

The goal is to construct sparse binary matrices (`HX`, `HZ`) that commute in the stabilizer sense. The algorithm adds controlled randomness to `HX` via weight‑preserving 2x2 cross‑swaps and then repairs `HZ` locally by solving a small general stabilizer linear system on a restricted domain. The approach is algebraic, uses only GF(2) arithmetic, and does not assume special structure beyond the stabilizer constraint itself.

### Contributions

- Local repair framed as a general (symmetric) GF(2) linear system on a small domain
- Adoption of solutions that preserve row/column weights inside the domain
- Exact DFS fallback for tiny systems
- Purely sparse (set‑of‑indices) implementation using only the Python standard library

---

## 2. Preliminaries and Notation

- Field: GF(2). Addition is XOR; we use set symmetric difference for sparse operations.
- Sparse binary matrices are represented by two lists of sets:
  - `rows[i]` = set of column indices `j` with `A[i,j] = 1`
  - `cols[j]` = set of row indices `i` with `A[i,j] = 1`
- Parity inner product between two supports `a` and `b` is `parity(|a ∩ b|)`.
- Initial construction: a tiled block‑ID pattern with parameters `P` (block size), `dc` (block rows) and `dr` (block columns). This yields
  - number of rows `m = dc * P`
  - number of columns `n = dr * P`
  - uniform row weight `dr` and column weight `dc`

**Stabilizer commutation (global check).** For all row pairs `(k, i)`,
```
parity( HX[k] ∩ HZ[i] ) XOR parity( HZ[k] ∩ HX[i] ) = 0
```
The function `stabilizer_violation_positions(HX, HZ)` enumerates all violating pairs.

---

## 3. Cross‑Swap (2x2) on HX

Pick two distinct rows `i1`, `i2` and two distinct columns `j1`, `j2`. Consider the 2x2 submatrix
```
[ a  b ]
[ c  d ]
```
It is swappable iff either `(a, d) = (1, 1)` and `(b, c) = (0, 0)` or the opposite. Toggling all four entries preserves both row and column weights. We apply this operation only to `HX`.

**Weight preservation.** Each affected row/column loses two 1s and gains two 1s (or vice versa), so degrees remain unchanged.

---

## 4. Local Repair (General Stabilizer System)

Let `HX'` be the swapped matrix. We seek a binary correction `Δ` such that `HZ' = HZ XOR Δ` restores commutation. Over GF(2) the required condition is equivalent to the linear system:
```
HX' Δ^T XOR Δ HX'^T = HX' HZ^T XOR HZ HX'^T    (Eq. 1)
```

We restrict `Δ` to a local domain `I × J` determined by the columns touched in the swap:
- `I`: rows reachable via `HZ`
- `J`: columns reachable from `I` via `HZ` or `HX'`
- `K`: rows reachable from `J` via `HX'` or `HZ` (K may intersect I)

The system in Eq. 1 is assembled only on `I × J` and only for rows in `K` (function `assemble_general`).

**Balance constraints (weight preservation inside I × J).**
- For every `i` in `I`, the number of `0→1` flips equals the number of `1→0` flips across columns in `J`.
- For every `j` in `J`, the number of `0→1` flips equals the number of `1→0` flips across rows in `I`.
The function `balanced_solution_from_vector` checks these constraints for a candidate solution.

---

## 5. Solvers

1. Solve Eq. 1 using sparse Gaussian elimination over GF(2) (`solve_gf2_sparse`).
2. Accept a candidate `Δ` only if balance constraints hold.
3. Apply `Δ` to obtain `HZ'` and run the global stabilizer validation. Only if it passes do we commit the change.
4. If the system is tiny (default: at most 24 variables), an exact DFS solver (`exact_balanced_solver`) is available as a fallback.

---

## 6. Algorithm Outline

```
Input: P, dc, dr, steps, seed
HX ← block_id(P, dc, dr)
HZ ← block_id(P, dc, dr)

for t in 1..steps:
  committed ← false
  for trial in 1..max_trials:
    (i1,i2,j1,j2) ← find_valid_cross_swap(HX)
    if none: break

    HX' ← HX with cross-swap at (i1,i2,j1,j2)
    assert row/col weights preserved

    (I,J,K) ← build_local_sets_general(HX', HZ, touched={j1,j2})
    (A,b)   ← assemble_general(HX', HZ, I,J,K)
    x       ← solve_gf2_sparse(A,b)

    if x is none or not balanced(x):
       if num_vars(x) small: x ← exact_balanced_solver(A,b)
       else: continue

    HZ' ← apply_delta(HZ, x)
    if stabilizer_ok(HX', HZ'):
       HX ← HX'; HZ ← HZ'; committed ← true; print matrices; break

  if not committed:
    print "no feasible swap"
```

**Correctness sketch.** Cross‑swaps preserve degrees of `HX`. Any `Δ` that satisfies Eq. 1 recovers commutation for the updated pair. Balance constraints preserve row/column weights in the local domain of `HZ`. A final global check prevents inconsistent commits.

---

## 7. Complexity and Parameters

- Global validation: roughly `O(m^2 * average_row_weight)` in the naive form.
- Local assembly: scales with `|I| * |J|`.
- Sparse elimination: practical up to a few hundred variables (depends on sparsity).
- DFS fallback: exponential; capped by a variable limit (default 24).

Recommended settings:
- `max_domain_size = |I| * |J|` around 200–600
- `dfs_var_cap` around 20–26
- Keep `steps` small for interactive runs to reduce output size

---

## 8. Key Functions (Map)

- `SparseBinMat`: sparse matrix container (`add`, `toggle`, `copy`, `nnz`)
- `dot_parity`: parity of set intersection
- `is_valid_cross_pattern`, `cross_swap_in_place`, `find_valid_cross_swap`: cross‑swap utilities
- `make_block_id_sparse`: block‑ID initializer
- `stabilizer_violation_positions`: global commutation checker
- `build_local_sets_general`: extracts local domain `(I, J, K)`
- `assemble_general`: builds the local linear system (Eq. 1)
- `solve_gf2_sparse`: sparse Gaussian elimination over GF(2)
- `balanced_solution_from_vector`: weight‑balance check
- `exact_balanced_solver`: exact DFS for tiny systems
- `adaptive_local_repair`: orchestrates the local repair pipeline
- `repeat_swaps_and_repairs_safe`: main loop (swap → repair → verify → commit)
- `print_sparse_01`: full 0/1 matrix printer (debug)

---

## 9. Usage

```
python3 stab-randomizer.py
```
Default parameters in the script:
- `P=25`, `dc=3`, `dr=4`, `steps=2000`, `seed=42`, `max_trials=300`
- Local repair defaults: `max_domain_size=400`, `dfs_var_cap=24`

Python 3.10 or newer is recommended. If you use `from __future__ import annotations`, place it at the top of the file.

---

## 10. Limitations and Future Work

- Acceptance and convergence rates are currently empirical.
- Neighborhood design `(I, J, K)` can be further optimized.
- Code distance, rank, and degeneracy evaluation are not included.
- Parallelization for swap search and validation is left for future work.


---

## License

Intended for research and educational release. Please include an open‑source license file (e.g., MIT, BSD‑3, Apache‑2.0) when publishing.
