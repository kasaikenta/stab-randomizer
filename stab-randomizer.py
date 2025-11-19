from __future__ import annotations
# -*- coding: utf-8 -*-
"""
Stabilizer LDPC: repeated cross-swap + local repair (general stabilizer case)

- Fully sparse (set-based) implementation suitable for larger P (e.g., P=15+)
- Preserves H_X row/col weights (via valid 2×2 cross-swap only)
- Preserves H_Z row/col weights on commit (balanced Δ)
- Prints full matrices (no truncation) on each commit; default run is large, so
  prefer smaller 'steps' when experimenting.

Debug-oriented design choices:
- If the CSS-like local system fails, fall back to the general stabilizer system.
- If kernel search is insufficient on tiny domains, try an exact balanced solver.
- Logic is written to be test-friendly; add tests rather than changing existing ones.
"""
from dataclasses import dataclass
from typing import List, Set, Dict, Tuple, Optional
import random

# ========================
# Sparse binary matrix API
# ========================

@dataclass
class SparseBinMat:
    """
    Sparse binary matrix represented as sets:
      - rows[i] = set of column indices j where A[i, j] = 1
      - cols[j] = set of row indices i where A[i, j] = 1

    Attributes
    ----------
    m : int
        Number of rows.
    n : int
        Number of columns.
    rows : List[Set[int]]
        Per-row sets of active column indices.
    cols : List[Set[int]]
        Per-column sets of active row indices.
    """
    m: int
    n: int
    rows: List[Set[int]]
    cols: List[Set[int]]

    @staticmethod
    def zeros(m: int, n: int) -> "SparseBinMat":
        """
        Create an m×n all-zero sparse binary matrix.

        Parameters
        ----------
        m : int
            Number of rows.
        n : int
            Number of columns.

        Returns
        -------
        SparseBinMat
        """
        return SparseBinMat(m, n, [set() for _ in range(m)], [set() for _ in range(n)])

    def copy(self) -> "SparseBinMat":
        """
        Deep copy of the sparse matrix (row/col sets are cloned).

        Returns
        -------
        SparseBinMat
        """
        return SparseBinMat(self.m, self.n, [set(s) for s in self.rows], [set(s) for s in self.cols])

    def add(self, i: int, j: int) -> None:
        """
        Set A[i, j] ← 1 (idempotent).

        Parameters
        ----------
        i : int
            Row index.
        j : int
            Column index.
        """
        if j not in self.rows[i]:
            self.rows[i].add(j)
            self.cols[j].add(i)

    def toggle(self, i: int, j: int) -> None:
        """
        Flip A[i, j] over GF(2): 0↔1.

        Parameters
        ----------
        i : int
            Row index.
        j : int
            Column index.
        """
        if j in self.rows[i]:
            self.rows[i].remove(j)
            self.cols[j].remove(i)
        else:
            self.rows[i].add(j)
            self.cols[j].add(i)

    def nnz(self) -> int:
        """
        Number of nonzeros (Hamming weight).

        Returns
        -------
        int
        """
        return sum(len(r) for r in self.rows)


# ================
# GF(2) primitives
# ================

def dot_parity(a: Set[int], b: Set[int]) -> int:
    """
    Parity of the set intersection: (|a ∩ b| mod 2).

    Parameters
    ----------
    a, b : Set[int]
        Supports of two binary vectors (as sets of indices).

    Returns
    -------
    int
        0 or 1 (parity).
    """
    # Iterate over the smaller set for speed
    if len(a) > len(b):
        a, b = b, a
    s = 0
    for x in a:
        if x in b:
            s ^= 1
    return s


# ===============
# Constructions
# ===============

def is_valid_cross_pattern(H: SparseBinMat, i1: int, i2: int, j1: int, j2: int) -> bool:
    """
    Check whether (i1,i2,j1,j2) forms a swappable 2×2 "cross" pattern:

        [a  b]
        [c  d]

    A valid cross-swap toggles the diagonal pair (a,d) ↔ (b,c) while preserving
    row/column weights, i.e., either (a,d)=1 and (b,c)=0 or vice versa.

    Returns
    -------
    bool
    """
    a = (j1 in H.rows[i1]); b = (j2 in H.rows[i1])
    c = (j1 in H.rows[i2]); d = (j2 in H.rows[i2])
    return (a and d and not b and not c) or (b and c and not a and not d)


def cross_swap_in_place(HX: SparseBinMat, i1: int, i2: int, j1: int, j2: int) -> None:
    """
    Perform the 2×2 cross-swap in-place by toggling the four entries.

    Notes
    -----
    Callers are responsible for ensuring the pattern is valid so that
    row/column weights are preserved.
    """
    HX.toggle(i1, j1)
    HX.toggle(i1, j2)
    HX.toggle(i2, j1)
    HX.toggle(i2, j2)


def find_valid_cross_swap(
    HX: SparseBinMat,
    rng: random.Random,
    max_trials: int = 800
) -> Optional[Tuple[int, int, int, int]]:
    """
    Randomly search for a valid 2×2 cross-swap on H_X.

    Algorithm
    ---------
    1. Uniformly sample two distinct rows (i1,i2).
    2. Form the set of candidate columns touched by either row and, if there
       are at least two, uniformly sample two columns (j1,j2).
    3. Check whether the 2×2 submatrix induced by (i1,i2,j1,j2) is a
       cross-pattern whose diagonals can be toggled without changing row/column
       weights via ``is_valid_cross_pattern``.
    4. Repeat up to ``max_trials`` times, returning the first valid tuple.

    Parameters
    ----------
    HX : SparseBinMat
        Current X-check matrix.
    rng : random.Random
        RNG instance for reproducibility.
    max_trials : int, optional
        Maximum random trials, by default 800.

    Returns
    -------
    Optional[Tuple[int,int,int,int]]
        (i1, i2, j1, j2) if found; otherwise None.
    """
    m, n = HX.m, HX.n
    for _ in range(max_trials):
        i1, i2 = rng.sample(range(m), 2)
        # Candidate columns are those touched by either row
        cand = list(HX.rows[i1] | HX.rows[i2])
        if len(cand) < 2:
            continue
        j1, j2 = rng.sample(cand, 2)
        if is_valid_cross_pattern(HX, i1, i2, j1, j2):
            return i1, i2, j1, j2
    return None


def make_block_id_sparse(P: int, dc: int, dr: int) -> SparseBinMat:
    """
    Build a (dc·P)×(dr·P) block-diagonal "ID" pattern that yields
    uniform row/column weights:
      - Each row has exactly dr ones.
      - Each column has exactly dc ones.

    For every block-row br∈[0, dc) and block-col bc∈[0, dr),
    we place a P×P identity along the shared diagonal.

    Parameters
    ----------
    P : int
        Circulant size per block.
    dc : int
        Column weight (per column).
    dr : int
        Row weight (per row).

    Returns
    -------
    SparseBinMat
    """
    H = SparseBinMat.zeros(dc * P, dr * P)
    for br in range(dc):
        for bc in range(dr):
            for t in range(P):
                H.add(br * P + t, bc * P + t)
    return H


# =========================
# Stabilizer global check
# =========================

def stabilizer_violation_positions(HX: SparseBinMat, HZ: SparseBinMat) -> List[Tuple[int, int]]:
    """
    Return all (k, i) where the stabilizer commutation constraint fails:

        H_X H_Z^T ⊕ H_Z H_X^T = 0  (over GF(2))

    Concretely, we flag (k, i) when
        dot_parity(HX[k, :], HZ[i, :]) ⊕ dot_parity(HZ[k, :], HX[i, :]) = 1.

    Assumes HX.m == HZ.m.

    Parameters
    ----------
    HX, HZ : SparseBinMat
        X- and Z-check matrices.

    Returns
    -------
    List[Tuple[int,int]]
        Pairs of row indices indicating violations.

    Complexity
    ----------
    O(m^2 · avg_row_weight) in the worst case (naïve double loop).
    """
    bad = []
    for k in range(HX.m):
        HXk = HX.rows[k]
        HZk = HZ.rows[k]
        for i in range(HX.m):
            if dot_parity(HXk, HZ.rows[i]) ^ dot_parity(HZk, HX.rows[i]):
                bad.append((k, i))
    return bad


# ==============================
# Neighborhood builders (sparse)
# ==============================

def build_local_sets_general(
    HXp: SparseBinMat,
    HZ: SparseBinMat,
    touched_cols: Set[int],
    max_I: int = 16,
    max_J: int = 32,
    max_K: int = 32
) -> Tuple[Set[int], Set[int], Set[int]]:
    """
    Build local index sets (I, J, K) for the general stabilizer case (K may
    intersect I). Use this when the CSS-like builder fails.

    Algorithm
    ---------
    1. Seed ``I`` with rows of H_Z incident to the touched columns until the
       ``max_I`` cap is met.
    2. Expand to ``J`` by collecting every column touched by rows in ``I`` from
       both H'_X and H_Z (capped by ``max_J``).
    3. Expand again to ``K`` via rows touching the provisional ``J`` columns in
       either matrix (capped by ``max_K``).
    4. Abort early if any frontier is empty so the caller can try another swap.

    Parameters
    ----------
    HXp : SparseBinMat
    HZ : SparseBinMat
    touched_cols : Set[int]
    max_I, max_J, max_K : int

    Returns
    -------
    Tuple[Set[int], Set[int], Set[int]]
    """
    I: Set[int] = set()
    for j in touched_cols:
        for i in HZ.cols[j]:
            I.add(i)
            if len(I) >= max_I:
                break
        if len(I) >= max_I:
            break
    if not I:
        return set(), set(), set()

    J: Set[int] = set()
    for i in list(I):
        for j in (HZ.rows[i] | HXp.rows[i]):
            J.add(j)
            if len(J) >= max_J:
                break
        if len(J) >= max_J:
            break
    if not J:
        return set(), set(), set()

    K: Set[int] = set()
    for j in list(J):
        for k in (HXp.cols[j] | HZ.cols[j]):
            K.add(k)
            if len(K) >= max_K:
                break
        if len(K) >= max_K:
            break
    return I, J, K


# ==========================
# Local linear systems (A x = b)
# ==========================

def assemble_general(
    HXp: SparseBinMat,
    HZ: SparseBinMat,
    I: Set[int],
    J: Set[int],
    K: Set[int]
):
    """
    Assemble the general local linear system over GF(2):

        For k∈K and i∈I:
          (H'_X)_{k,J} Δ_{i,J}^T ⊕ Δ_{k,J} (H'_X)_{i,J}^T
          = (H'_X H_Z^T ⊕ H_Z H'_X^T)_{k,i}

    When k∉I, the Δ_{k,J} term is implicitly zero.

    Algorithm
    ---------
    * Sort I,J,K to obtain deterministic index lists and assign each (i,j)∈I×J
      a column in the linear system.
    * For every k∈K build one equation per i∈I by:
        a) Adding the variables corresponding to overlaps between row k of H'_X
           and the J-domain (first bilinear term).
        b) When k∈I, adding the variables representing row i overlaps so the
           second bilinear term shares coefficients with Δ_{k,J}.
    * Emit the commutator parity on the right-hand side using the original HX
      and HZ rows.

    Returns
    -------
    row_eqs : List[Set[int]]
        Sparse equations representing the linear system.
    rhs : List[int]
        Right-hand side bits.
    var_index : Dict[Tuple[int, int], int]
        Mapping (i, j) → column indices in the linear system.
    I_list, J_list : List[int]
        Ordered lists describing the row/column domains used above.
    """
    I_list = sorted(I); J_list = sorted(J); K_list = sorted(K)
    var_index: Dict[Tuple[int,int], int] = {
        (i, j): idx
        for idx, (i, j) in enumerate((ii, jj) for ii in I_list for jj in J_list)
    }

    row_eqs: List[Set[int]] = []
    rhs: List[int] = []
    Jset = set(J_list)

    for k in K_list:
        HXk_on_J = HXp.rows[k] & Jset
        for i in I_list:
            eq = set()
            # (H'_X)_{k,J} Δ_{i,J}^T
            for j in HXk_on_J:
                eq.add(var_index[(i, j)])
            # Δ_{k,J} (H'_X)_{i,J}^T when k∈I
            if k in I:
                for j in (HXp.rows[i] & Jset):
                    eq.add(var_index[(k, j)])
            e = dot_parity(HXp.rows[k], HZ.rows[i]) ^ dot_parity(HZ.rows[k], HXp.rows[i])
            row_eqs.append(eq)
            rhs.append(e)

    return row_eqs, rhs, var_index, I_list, J_list


# ======================
# GF(2) sparse elimination
# ======================

def solve_gf2_sparse(row_eqs: List[Set[int]], rhs: List[int], nvars: int) -> Optional[List[int]]:
    """
    Solve A x = b over GF(2) where each row of A is given as a set of active
    column indices (sparse support). Uses a simple sparse Gaussian elimination
    with symmetric-difference operations on sets.

    Algorithm
    ---------
    1. Clone each row support so mutations do not leak to callers.
    2. Sweep the rows from top to bottom, picking the smallest column index in
       the current support as the candidate pivot. If that column already has a
       pivot row, symmetric-difference the two rows (elimination); otherwise
       record the pivot and eliminate the column from all lower rows.
    3. Abort immediately when a zero row has RHS=1 (inconsistent system).
    4. Perform a backward sweep that assigns pivot variables first and treats
       unpivoted variables as free zeros.

    Parameters
    ----------
    row_eqs : List[Set[int]]
        Sparse rows (supports).
    rhs : List[int]
        Right-hand side.
    nvars : int
        Number of variables (columns).

    Returns
    -------
    Optional[List[int]]
        A 0/1 solution vector of length nvars if consistent; otherwise None.
    """
    A = [set(r) for r in row_eqs]
    b = rhs[:]
    pivot_of_col: Dict[int, int] = {}
    col_of_row: Dict[int, int] = {}

    # Forward elimination
    for r in range(len(A)):
        if not A[r]:
            if b[r]:
                return None  # inconsistent: 0 = 1
            continue
        while A[r]:
            c = min(A[r])
            if c in pivot_of_col:
                pr = pivot_of_col[c]
                A[r] ^= A[pr]; b[r] ^= b[pr]
            else:
                pivot_of_col[c] = r
                col_of_row[r] = c
                # Eliminate below
                for rr in range(r + 1, len(A)):
                    if c in A[rr]:
                        A[rr] ^= A[r]; b[rr] ^= b[r]
                break
        if not A[r] and b[r]:
            return None

    # Back-substitution (free vars default to 0)
    x = [0] * nvars
    for r in range(len(A) - 1, -1, -1):
        if r not in col_of_row:
            continue
        c = col_of_row[r]
        s = 0
        for j in A[r]:
            if j != c:
                s ^= x[j]
        x[c] = s ^ b[r]
    return x


# =====================================
# Balanced solution (weights preserved)
# =====================================

def balanced_solution_from_vector(
    x: List[int],
    var_index: Dict[Tuple[int, int], int],
    HZ: SparseBinMat,
    I_list: List[int],
    J_list: List[int]
) -> bool:
    """
    Check whether the proposed Δ (encoded by x) preserves the row and column
    weights of H_Z within the local domain (I×J). Toggling from 0→1 counts +1,
    and 1→0 counts −1; all row and column totals must net to zero.

    Returns
    -------
    bool
        True if row/col balances are preserved; False otherwise.
    """
    R = {i: 0 for i in I_list}  # row balances
    C = {j: 0 for j in J_list}  # column balances
    inv = {idx: ij for ij, idx in var_index.items()}

    for idx, val in enumerate(x):
        if val:
            i, j = inv[idx]
            # +1 if HZ[i,j] was 0 (becomes 1), −1 if it was 1 (becomes 0)
            s = 1 - 2 * (1 if j in HZ.rows[i] else 0)
            R[i] += s
            C[j] += s

    return all(v == 0 for v in R.values()) and all(v == 0 for v in C.values())


def exact_balanced_solver(
    row_eqs: List[Set[int]],
    rhs: List[int],
    var_index: Dict[Tuple[int, int], int],
    HZ: SparseBinMat,
    I_list: List[int],
    J_list: List[int],
    max_nodes: int = 200_000
) -> Optional[List[int]]:
    """
    Exact brute-force (DFS) solver for tiny local systems with balance constraints.
    Useful when elimination under-constrains the solution and we need an exact,
    weight-preserving Δ.

    Algorithm
    ---------
    * Precompute equation adjacencies and the ±1 contribution that each
      variable would apply to every affected row/column balance in H_Z.
    * Order the variables by decreasing equation degree for stronger pruning.
    * Run depth-first search that first tries assigning 0, then 1. When testing
      1 the solver flips cached parities and accumulates row/column balances,
      unwinding on backtrack.
    * Terminate when every equation parity matches ``rhs`` and all balances are
      zero, or when the node budget ``max_nodes`` is exceeded.

    Parameters
    ----------
    row_eqs, rhs, var_index : see assemblers
    HZ : SparseBinMat
        Current H_Z to evaluate 0→1 / 1→0 signs.
    I_list, J_list : List[int]
        Ordered local row/column domains.
    max_nodes : int
        Early cutoff for the DFS search tree.

    Returns
    -------
    Optional[List[int]]
        A balanced solution vector if found; otherwise None.
    """
    v = len(var_index)

    # Build variable→equations adjacency
    var_eqs: List[List[int]] = [[] for _ in range(v)]
    for r, cols in enumerate(row_eqs):
        for c in cols:
            var_eqs[c].append(r)

    inv = {idx: ij for ij, idx in var_index.items()}

    # Precompute +1 / −1 signs per variable for row/col balances
    sgn_row = [0] * v
    sgn_col = [0] * v
    i_to_pos = {i: p for p, i in enumerate(I_list)}
    j_to_pos = {j: p for p, j in enumerate(J_list)}
    for idx, (i, j) in inv.items():
        sgn = 1 - 2 * (1 if j in HZ.rows[i] else 0)
        sgn_row[idx] = (i_to_pos[i], sgn)
        sgn_col[idx] = (j_to_pos[j], sgn)

    # DFS state
    eqp = [0] * len(row_eqs)  # current LHS parity per equation
    R = [0] * len(I_list)     # row balances
    C = [0] * len(J_list)     # column balances
    deg = [0] * v             # heuristic: higher degree first
    for cols in row_eqs:
        for c in cols:
            deg[c] += 1
    order = sorted(range(v), key=lambda x: (-deg[x], x))  # static ordering

    x = [0] * v
    nodes = 0

    def dfs(t: int) -> Optional[List[int]]:
        nonlocal nodes
        nodes += 1
        if nodes > max_nodes:
            return None
        if t == v:
            # All equations and balances must be satisfied
            for r in range(len(row_eqs)):
                if eqp[r] != rhs[r]:
                    return None
            if any(R) or any(C):
                return None
            return x[:]

        j = order[t]

        # Try 0
        sol = dfs(t + 1)
        if sol is not None:
            return sol

        # Try 1 (toggle this variable)
        x[j] = 1
        changed = []
        for r in var_eqs[j]:
            eqp[r] ^= 1
            changed.append(r)
        ir, s = sgn_row[j]; jr, s2 = sgn_col[j]
        R[ir] += s; C[jr] += s2

        sol = dfs(t + 1)

        # Revert
        for r in changed:
            eqp[r] ^= 1
        R[ir] -= s; C[jr] -= s2
        x[j] = 0

        return sol

    return dfs(0)


# ==============================
# Adaptive local repair pipeline
# ==============================

def adaptive_local_repair(
    HX: SparseBinMat,
    HZ: SparseBinMat,
    HXp: SparseBinMat,
    touched_cols: Set[int],
    *,
    max_domain_size: int = 400,
    dfs_var_cap: int = 24
) -> Tuple[bool, SparseBinMat]:
    """
    Locally repair H_Z after a proposed swap on H_X, using ONLY the general
    stabilizer linearization (no CSS-like fast path).

    Strategy
    --------
    1) Build local neighborhoods (I, J, K) via `build_local_sets_general`.
    2) Assemble the general GF(2) system (symmetrized):
         (H'_X)_{k,J} Δ_{i,J}^T ⊕ Δ_{k,J} (H'_X)_{i,J}^T
         = (H'_X H_Z^T ⊕ H_Z H'_X^T)_{k,i}
    3) Solve with sparse elimination. Enforce row/col weight balance of H_Z.
    4) If underdetermined and tiny, try an exact DFS with balance constraints.

    On success: return (True, repaired_HZ); else: (False, original HZ).
    """
    # 1) General local domain
    I, J, K = build_local_sets_general(HXp, HZ, touched_cols)
    if not (I and J and K):
        return False, HZ
    if len(I) * len(J) > max_domain_size:
        # Domain is too large; skip this proposal (caller will try another swap)
        return False, HZ

    # 2) Assemble symmetric (general) system
    row_eqs, rhs, var_index, I_list, J_list = assemble_general(HXp, HZ, I, J, K)

    # 3) Solve by sparse elimination and check balance
    x = solve_gf2_sparse(row_eqs, rhs, len(var_index))
    if x is not None and balanced_solution_from_vector(x, var_index, HZ, I_list, J_list):
        HZ_new = HZ.copy()
        inv = {idx: ij for ij, idx in var_index.items()}
        for idx, val in enumerate(x):
            if val:
                i, j = inv[idx]
                HZ_new.toggle(i, j)
        # Global stabilizer check
        if not stabilizer_violation_positions(HXp, HZ_new):
            return True, HZ_new

    # 4) Exact DFS fallback (only if the number of variables is small enough)
    if len(var_index) <= dfs_var_cap:
        x = exact_balanced_solver(row_eqs, rhs, var_index, HZ, I_list, J_list)
        if x is not None:
            HZ_new = HZ.copy()
            inv = {idx: ij for ij, idx in var_index.items()}
            for idx, val in enumerate(x):
                if val:
                    i, j = inv[idx]
                    HZ_new.toggle(i, j)
            if not stabilizer_violation_positions(HXp, HZ_new):
                return True, HZ_new

    return False, HZ


# ======================
# Pretty printer (full)
# ======================

def print_sparse_01(A: SparseBinMat, name: str) -> None:
    """
    Print a full 0/1 view of a sparse matrix (no truncation).

    Parameters
    ----------
    A : SparseBinMat
    name : str
        Label to print before the matrix.
    """
    print(f"{name} (shape={A.m}x{A.n}, nnz={A.nnz()}):")
    for i in range(A.m):
        row = ["1" if j in A.rows[i] else "0" for j in range(A.n)]
        print("".join(row))


# ===============
# Main runner
# ===============

def repeat_swaps_and_repairs_safe(
    P: int = 10,
    dc: int = 3,
    dr: int = 4,
    steps: int = 1000,
    seed: int = 42,
    max_trials: int = 300
) -> Tuple[SparseBinMat, SparseBinMat]:
    """
    Repeatedly apply random cross-swaps to H_X and repair H_Z locally.
    Validate the stabilizer constraint (H_X H_Z^T ⊕ H_Z H_X^T = 0) at each step.

    Algorithm
    ---------
    1. Initialize H_X and H_Z with the balanced block-diagonal pattern produced
       by ``make_block_id_sparse``.
    2. For every outer iteration, repeatedly sample a valid cross-swap via
       ``find_valid_cross_swap`` and simulate it on a copy of H_X.
    3. Invoke ``adaptive_local_repair`` to obtain a candidate H_Z update, and
       discard the attempt unless the commutator violations disappear.
    4. Commit the new pair (H_X, H_Z) only when weights and commutation checks
       succeed, otherwise try another swap until ``max_trials`` is exhausted.
    5. Emit verbose diagnostics—including the full matrices—after each commit
       and re-check the constraint at the very end.

    Parameters
    ----------
    P : int
        Block period for the initial construction.
    dc : int
        Column weight.
    dr : int
        Row weight.
    steps : int
        Number of outer iterations.
    seed : int
        Random seed.
    max_trials : int
        Maximum attempts for finding a valid swap per step.

    Returns
    -------
    Tuple[SparseBinMat, SparseBinMat]
        The final (H_X, H_Z).
    """
    rng = random.Random(seed)

    HX = make_block_id_sparse(P, dc, dr)
    HZ = make_block_id_sparse(P, dc, dr)

    success_count = 0

    for t in range(1, steps + 1):
        committed = False

        for trial in range(1, max_trials + 1):
            pick = find_valid_cross_swap(HX, rng, max_trials=800)
            if pick is None:
                break

            i1, i2, j1, j2 = pick
            HXp = HX.copy()

            # Check row/column weights before and after the swap
            row_w_before = [len(HX.rows[i]) for i in range(HX.m)]
            col_w_before = [len(HX.cols[j]) for j in range(HX.n)]

            cross_swap_in_place(HXp, i1, i2, j1, j2)

            row_w_after = [len(HXp.rows[i]) for i in range(HXp.m)]
            col_w_after = [len(HXp.cols[j]) for j in range(HXp.n)]
            assert row_w_before == row_w_after and col_w_before == col_w_after, \
                "H_X weights changed during the swap"

            # Attempt a local repair
            ok, HZ_new = adaptive_local_repair(HX, HZ, HXp, {j1, j2})

            if ok:
                # Validate the stabilizer constraint after repair
                viol = stabilizer_violation_positions(HXp, HZ_new)
                if viol:
                    print(f"[step {t}] ⚠️ {len(viol)} commutation violations — discarding repair")
                    continue  # Discard and try again

                # No issues → commit
                HX, HZ = HXp, HZ_new
                success_count += 1
                print(f"[step {t}] committed at trial {trial} (Global E=0, weights preserved)")
                print_sparse_01(HX, "H_X")
                print_sparse_01(HZ, "H_Z")

                # Double-check: final stabilizer condition
                final_viol = stabilizer_violation_positions(HX, HZ)
                if final_viol:
                    print(f"[check] ⚠️ {len(final_viol)} commutation violations remain after commit")
                else:
                    print(f"[check] ✅ Stabilizer constraint OK")

                #print(f"[summary] commits so far: {success_count}/{t}")
                committed = True
                break

        if not committed:
            print(f"[step {t}] no feasible swap found. (commits so far: {success_count}/{t})")

    print(f"[final] total commits: {success_count}/{steps}")

    # Final commutation check
    final_viol = stabilizer_violation_positions(HX, HZ)
    if final_viol:
        print(f"[final check] ⚠️ {len(final_viol)} commutation violations remain in the final state")
    else:
        print("[final check] ✅ Final H_X and H_Z satisfy the commutation constraint")

    return HX, HZ


if __name__ == "__main__":
    # Consider reducing 'steps' for interactive runs to avoid huge output.
    repeat_swaps_and_repairs_safe()
