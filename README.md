\
# Stabilizer LDPC Randomizer — 論文調 README

**著者:** 笠井健太  
**概要:** 本プログラムは，スタビライザ LDPC 行列対 \(H_X, H_Z\) に対して，2×2 交差スワップ（cross-swap）により \(H_X\) をランダム化し，それに伴う可換条件の破れを **一般スタビライザ系（対称形）** の局所解法で修復するアルゴリズムを提供する．設計上，(i) 行・列重みの保存，(ii) スタビライザ可換条件
\[
H_X H_Z^\top \oplus H_Z H_X^\top = 0 \quad (\text{GF}(2))
\]
の維持をコミット条件とし，(iii) 小規模領域に対しては厳密 DFS による完全探索を備える．実装は完全疎（集合型）であり，Python 標準ライブラリのみで動作する．

---

## 1. 背景と目的

量子 LDPC 符号の構成では，スタビライザ可換条件を満たす疎行列対 \((H_X, H_Z)\) の設計が根幹となる．本プログラムは，**行・列重みを保存する微小操作（交差スワップ）** を繰り返し適用しつつ，可換条件を満たすよう **局所的に \(H_Z\) を修復** することで，構造を保ちながら多様な実例を生成することを目的とする．CSS 構造に依らず，一般のスタビライザ条件を直接扱う点が特徴である．

### 貢献
- **一般スタビライザ系のみ**を用いた局所修復（CSS-like の近似は用いない）．
- 局所領域 \(I\times J\) に制限した **対称形線形系** の組立と疎ガウス消去．
- **行・列重み保存（balance）制約** を満たす解のみを採用．
- 微小規模では **厳密 DFS** による完全探索フォールバック．
- 完全疎データ構造（行・列のインデックス集合）のみで高可読・高拡張性．

---

## 2. 予備知識と記法

- 体は GF(2)，加法は XOR，集合の対称差 \(\triangle\) を用いる．
- 疎二値行列は **行サポート集合** \(\text{rows}[i]\subseteq[n)\) と **列サポート集合** \(\text{cols}[j]\subseteq[m)\) で表現する．
- パリティ内積は
  \[
  \langle a,b\rangle := |a\cap b|\bmod 2
  \]
  で与える（実装では `dot_parity`）．
- 初期行列はブロック ID 構成：サイズ \(P\)，ブロック行数 \(d_c\)，ブロック列数 \(d_r\) に対し，
  \[
  m=d_c P,\quad n=d_r P,\quad \text{row-weight}=d_r,\quad \text{col-weight}=d_c.
  \]

**スタビライザ可換条件.** 行列対 \((H_X,H_Z)\) は
\[
H_X H_Z^\top \oplus H_Z H_X^\top = O_{m\times m}
\]
を満たすとき可換である．本実装では，全対 \((k,i)\) に対し
\[
\langle H_X[k,:], H_Z[i,:]\rangle \oplus \langle H_Z[k,:], H_X[i,:]\rangle = 0
\]
を検証（`stabilizer_violation_positions`）する．

---

## 3. 交差スワップ（Cross-Swap）

2 行 \(i_1,i_2\) と 2 列 \(j_1,j_2\) の 2×2 部分
\[
\begin{bmatrix}
a & b\\
c & d
\end{bmatrix}
\]
が \((a,d)=(1,1),(b,c)=(0,0)\) **または** その逆であれば，4 要素を同時にトグルすることで **各行・各列の重みは不変** となる（`is_valid_cross_pattern`，`cross_swap_in_place`）．本操作は \(H_X\) にのみ適用し，行・列次数を保ったまま構造を撹乱する．

**補題（重み保存）.** 上記条件下でのトグルは，各行・各列で 1 が 2 つ消えて 2 つ生じる（またはその逆）ため，次数は不変．

---

## 4. 局所修復問題の定式化（一般スタビライザ系）

スワップ後の \(H'_X\) に対し，\(H_Z\) を \(\Delta\)（0/1 のトグル行列）で更新して \(H'_Z:=H_Z\oplus\Delta\) とする．可換条件は
\[
H'_X (H_Z\oplus\Delta)^\top \oplus (H_Z\oplus\Delta) H'_X{}^\top = O
\]
であり，GF(2) では線形化して
\[
H'_X H_Z^\top \oplus H_Z H'_X{}^\top \oplus H'_X \Delta^\top \oplus \Delta H'_X{}^\top = O.
\]
したがって \(\Delta\) は
\[
H'_X \Delta^\top \oplus \Delta H'_X{}^\top = H'_X H_Z^\top \oplus H_Z H'_X{}^\top
\tag{1}
\]
を満たさねばならない．本実装では，スワップで影響を受けやすい **局所領域** を
- \(I\): 触れた列から \(H_Z\) で辿れる行集合，
- \(J\): \(I\) から \(H_Z\) または \(H'_X\) で辿れる列集合，
- \(K\): \(J\) から \(H'_X\) または \(H_Z\) で辿れる行集合（\(K\cap I\) を許す）
として抽出し，\(\Delta\) のサポートを \(I\times J\) に **制限** して (1) を行ごとに組み立てる（`assemble_general`）．

**平衡制約（balance）.** \(\Delta\) 適用により，領域 \(I\times J\) 内の **各行・各列の 1 の個数が不変** となるよう
\[
\sum_{j\in J} \sigma_{i j} = 0,\quad \sum_{i\in I} \sigma_{i j} = 0
\]
を課す（\(\sigma_{ij}=+1\) if \(0\to1\), \(-1\) if \(1\to0\)）．`balanced_solution_from_vector` が検査する．

---

## 5. 解法

- まず疎ガウス消去（GF(2), `solve_gf2_sparse`）で (1) の連立を解く．
- 得られた解が平衡制約を満たせば \(\Delta\) を適用して候補 \(H'_Z\) を構成．
- 候補 \((H'_X,H'_Z)\) はグローバル可換検証（全 \((k,i)\)）を通過した場合にのみ採用．
- なお未知数が少数（既定 24 以下）の場合は **厳密 DFS**（`exact_balanced_solver`）で完全探索を行う（平衡制約込み）．

---

## 6. アルゴリズム（擬似コード）

```text
Input: P, d_c, d_r, steps, seed
H_X ← block_id(P, d_c, d_r)
H_Z ← block_id(P, d_c, d_r)

for t = 1..steps:
  for trial = 1..max_trials:
    pick (i1,i2,j1,j2) with valid cross pattern on H_X
    H'_X ← swap(H_X, i1,i2,j1,j2)
    assert deg(H'_X) == deg(H_X)  // 行列重み保存

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

**正当性スケッチ.** (i) 交差スワップは \(H_X\) の重みを保存．(ii) (1) を満たす \(\Delta\) を領域 \(I\times J\) に適用すると可換条件が回復．(iii) 平衡制約により \(H_Z\) の（少なくとも局所領域内の）行・列重みが保たれる．(iv) コミット前のグローバル検証により不整合は排除される．

---

## 7. 計算量と実装注意

- **大域検証**：最悪 \(O(m^2 \cdot \bar{w})\)（\(\bar{w}\) は平均行重み）．
- **局所系の組立**：\(|I|\cdot|J|\) にほぼ比例．
- **疎消去**：非零の疎度に依存（実測では数百変数規模まで実用的）．
- **DFS**：指数的（既定では変数 24 程度までに制限）．

**パラメータ指針.**
- `max_domain_size = |I|\cdot|J|` は 200–600 程度が実用的な目安．
- `dfs_var_cap` は 20–26 程度．大きくしすぎると急激に遅くなる．
- ログが大きいので `steps` は小さく始めること．

---

## 8. 実装詳細（対応関数）

- `SparseBinMat`：疎行列（行・列サポート）／`add`／`toggle`／`copy`／`nnz`．
- `dot_parity`：パリティ内積（小さい集合で走査）．
- `is_valid_cross_pattern`／`cross_swap_in_place`／`find_valid_cross_swap`：交差スワップ．
- `make_block_id_sparse`：初期ブロック ID 生成．
- `stabilizer_violation_positions`：全ペア検証（違反リスト）．
- `build_local_sets_general`：局所領域 \(I,J,K\) の抽出．
- `assemble_general`：式 (1) の疎行列化．
- `solve_gf2_sparse`：疎ガウス消去（GF(2)）．
- `balanced_solution_from_vector`：行・列平衡制約の検査．
- `exact_balanced_solver`：小規模厳密 DFS．
- `adaptive_local_repair`：上記を束ねた局所修復パイプライン（一般系のみ）．
- `repeat_swaps_and_repairs_safe`：メインループ（提案→修復→検証→コミット）．
- `print_sparse_01`：0/1 全量表示（デバッグ向け）．

---

## 9. 使用方法

```bash
python3 stab-randomizer.py
```
主な引数（関数内デフォルト）:
- `P=25, d_c=3, d_r=4, steps=2000, seed=42, max_trials=300`
- `adaptive_local_repair`: `max_domain_size=400, dfs_var_cap=24`

**Python 3.10+** を推奨．`from __future__ import annotations` を用いる場合は **ファイル先頭** に置くこと．

---

## 10. 限界と今後の課題

- 受理率・収束性の解析はヒューリスティックに留まる．
- 局所領域の設計（\(I,J,K\) の拡張規則）は改良の余地がある．
- 行列距離・ランク・位数などの符号パラメータ評価は未実装．
- 並列化（スワップ候補探索・検証の並列）は将来的改善点．

---

## 参考文献（例）
- D. Gottesman, *Stabilizer Codes and Quantum Error Correction*, PhD thesis, 1997.  
- A. R. Calderbank and P. W. Shor, *Good quantum error-correcting codes exist*, Phys. Rev. A, 1996.  
- A. M. Steane, *Multiple-Particle Interference and Quantum Error Correction*, Proc. R. Soc. A, 1996.  
- J.-P. Tillich and G. Zémor, *Quantum LDPC codes with positive rate and minimum distance proportional to n*, ISIT 2009.  

> **注**：本実装は CSS-like の近似系を用いず，(1) の一般対称形のみを解く．

---

## ライセンス
研究・教育目的での公開を想定．正式な公開に際してはライセンスを明記すること（例：MIT, BSD-3, Apache-2.0 など）．

## 謝辞
本 README は，学術資料としての再利用を想定し，数式と正当性スケッチを併記した．
