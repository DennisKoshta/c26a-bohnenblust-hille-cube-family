# An explicit family for the Boolean Bohnenblust–Hille inequality: $2(1-2^{-d})^{(d+1)/(2d)}$ at every degree

*Note, version of 2026-10-26 (first written 2026-10-23). Not peer reviewed. See the AI-use disclosure at the end.*

## 1. Statement

For $f:\{\pm1\}^n\to\mathbb R$ write $f(x)=\sum_{S\subseteq[n]}\widehat f(S)\,x^S$, with $x^S=\prod_{i\in S}x_i$ and
$\|f\|_\infty=\max_{x\in\{\pm1\}^n}|f(x)|$. For $d\ge1$ put $p_d=\frac{2d}{d+1}$, and let
$\mathrm{BH}^{\le d}_{\{\pm1\}}$ be the least $C$ such that

$$\Big(\sum_{|S|\le d}|\widehat f(S)|^{p_d}\Big)^{1/p_d}\le C\,\|f\|_\infty$$

for all $n$ and all $f$ of degree at most $d$. Every single $f$ gives the lower bound
$\mathrm{BH}^{\le d}_{\{\pm1\}}\ge R(f):=\|\widehat f\|_{p_d}/\|f\|_\infty$.

**Theorem.** For every $d\ge1$ there is $f_d:\{\pm1\}^{3d}\to\{0,\pm1\}$ with the following properties:

- its degree is at most $d$;
- it has exactly $(2^d-1)2^d$ nonzero Fourier coefficients, all equal to $\pm2^{-d}$;
- $\Pr[f_d=0]=2^{-d}$.

Consequently

$$\mathrm{BH}^{\le d}_{\{\pm1\}}\;\ge\;R_d:=2\,(1-2^{-d})^{\frac{d+1}{2d}} .$$

For $d\ge2$ this is strictly larger than $2^{(d-1)/d}$, the value of the degree-$d$ address function. Also
$R_d<2$ and $R_d\to2$, with $2-R_d\sim\frac{d+1}{d}2^{-d}$.

| $d$ | $R_d$ | address $2^{(d-1)/d}$ | $d(1-\log_2R_d)$ |
|---|---|---|---|
| 2 | 1.611854897735 | 1.414214 | 0.6226 |
| 3 | 1.829652855011 | 1.587401 | 0.3853 |
| 4 | 1.920932230356 | 1.681793 | 0.2328 |
| 5 | 1.962262141550 | 1.741101 | 0.1374 |
| 6 | 1.981711051583 | 1.781797 | 0.0795 |
| 7 | 1.991056425313 | 1.811447 | 0.0453 |

(The digits are truncated lower ends of exact interval enclosures; see §4.)

**Relation to Tao's table.** The entry is constant 26a of github.com/teorth/optimizationproblems,
`constants/26a.md`, read on 2026-10-23 at repository HEAD `46dbf03`, committed 2026-10-17. It defines
$C_{26a}=\sup_{d\ge1}\mathrm{BH}^{\le d}_{\{\pm1\}}$ and lists the lower bound $2$, citing ADGP2025. That bound is
the limit of the address-function values $2^{(d-1)/d}$. The page adds that "chasing incremental improvements of
lower bounds seems less compelling than establishing any finite uniform upper bound", and that a construction
forcing the constant above $100$ "would already be a genuinely interesting result". Since every $R_d<2$, **the
theorem does not change the table's lower bound for $C_{26a}$.** It improves the best lower bound we know of at
each fixed degree $d\ge2$.

**Literature.** The relevant papers, with dates and what each says:

- **Arunachalam, Dutt, Escudero Gutiérrez, Palazuelos.** Math. Ann. 392 (2025); arXiv:2405.10933. This is the
  table's source for the bound 2, via the address function. Fact 7 there shows that Boolean ($\pm1$-valued) $f$
  satisfy $R(f)\le2^{(d-1)/d}$; Proposition 2(b) below reproves this.
- **Defant, Mastyło, Pérez.** Math. Ann. 374 (2019); arXiv:1706.03670. Subexponential upper bound.
- **Slote.** arXiv:2608.01424, 2026-08-02. The paper says it is not aware of any lower bound for the classical cube
  constant that is asymptotically greater than 2.
- **Upper-bound papers, autumn 2026.** None contains a lower-bound construction:
  - Ivanisvili, arXiv:2609.12427, 2026-09-11 (polynomial growth);
  - Slote–Tseng–Volberg, arXiv:2609.21144, 2026-09-17;
  - Barbosa–Núñez-Alarcón–Pellegrino–Raposo, arXiv:2610.15363, 2026-10-11.
- **Pellegrino–Teixeira, "Typical Bohnenblust–Hille ratios".** arXiv:2609.31779, 2026-09-24. Random homogeneous
  polynomials; not a lower-bound construction.
- **Pellegrino–Teixeira, arXiv:1604.07595 (2016).** Conjectures that $2^{1-1/m}$ is the sharp real constant for
  **$m$-linear forms**. Our $f_d$ are not $m$-linear forms (see §5), so the theorem says nothing about that
  conjecture.

**Literature check today (2026-10-23).**

- Short arXiv keyword searches ("Bohnenblust–Hille", "Bohnenblust Hille Hamming cube", "Bohnenblust Hille lower
  bound") returned only the papers above.
- In Tao's repository at `46dbf03`, `26a.md` still lists the lower bound 2. The pull requests touching it (#49,
  typo fixes, 2026-03-11; #137) do not report a new lower bound.
- Re-check 2026-10-26: `26a.md` is unchanged at main `e95af44` (last edit to the file 2026-03-09). No open pull
  request touches it. An arXiv search for submissions from 2026-10-20 to 2026-10-27 mentioning "Bohnenblust" found
  nothing.
- I found nothing that supersedes or anticipates the theorem. The search was a short keyword search, not a
  systematic review. In particular, the literature on real homogeneous polynomial BH constants was not searched for
  square-free examples (see §5).

## 2. Notation and two lemmas

Identify $x\in\{\pm1\}^n$ with $\xi\in\mathbb F_2^n$ via $x_i=(-1)^{\xi_i}$. For $\alpha\in\mathbb F_2^n$ put
$\chi_\alpha(\xi)=(-1)^{\alpha\cdot\xi}$, which is $x^S$ with $S=\operatorname{supp}\alpha$, and let $|\alpha|$ be its
weight.

**Lemma 1 (blocks).** Let $V\le\mathbb F_2^n$ have dimension $k$, let
$A=\{\xi:\alpha\cdot\xi=\alpha\cdot\xi_0\ \forall\alpha\in V\}$ be the affine subspace through $\xi_0$, and let
$\beta\in\mathbb F_2^n$. Then

$$\chi_\beta\mathbf 1_A=2^{-k}\sum_{\alpha\in V}\chi_\alpha(\xi_0)\,\chi_{\beta+\alpha}.$$

So $\chi_\beta\mathbf 1_A$ has $2^k$ Fourier coefficients $\pm2^{-k}$, supported on the coset $\beta+V$, and its
degree is $\max_{\alpha\in V}|\beta+\alpha|$.

*Proof.* $\mathbf 1_A(\xi)=2^{-k}\sum_{\alpha\in V}(-1)^{\alpha\cdot(\xi+\xi_0)}$. Multiply by $\chi_\beta$. $\square$

**Lemma 2 (path cosets).** Let $P$ be a set of $m+1$ coordinates, and let $E_P$ be the space of even-weight vectors
supported on $P$; it has dimension $m$ and is spanned by $e_u+e_v$ over the edges $uv$ of any path through $P$. Let
$C_P$ be the set of vectors supported on $P$ whose weight is congruent to $m$ mod 2. Then $C_P$ is a coset of $E_P$,
and every element of $C_P$ has weight at most $m$.

*Proof.* $C_P$ is $E_P$ itself if $m$ is even, and $e_u+E_P$ for any $u\in P$ if $m$ is odd. The weight is at most
$m+1$, and $m+1$ has the wrong parity. $\square$

## 3. The construction and the proof

For $d=1$ take $f_1=(a+b)/2$: here $R=1=R_1$. Let $d\ge2$, index $i\in\mathbb Z_d$, and use $3d$ coordinates
$z_i,a_i,b_i$. Define the edge forms $\varepsilon_i=e_{z_i}+e_{z_{i+1}}$ (the edges of a $d$-cycle on the $z$'s; for
$d=2$ the two edges coincide as vectors) and the pair forms $\pi_i=e_{a_i}+e_{b_i}$. Fix $c\in\mathbb F_2^{d}$: $c=0$
if $d$ is odd, $c=e_0$ if $d$ is even. In both cases $d+|c|$ is odd.

For each proper subset $J\subsetneq\mathbb Z_d$ (there are $2^d-1$ of them) define

$$w_J(i)=\mathbf 1_J(i-1)+c_i\in\mathbb F_2,\qquad V_J=\operatorname{span}\big(\{\varepsilon_i\}_{i\in J}\cup\{\pi_i\}_{i\notin J}\big),$$

$$A_J=\{\xi:\ \varepsilon_i\cdot\xi=w_J(i)\ (i\in J),\ \ \pi_i\cdot\xi=w_J(i)\ (i\notin J)\}.$$

**Structure of $V_J$.** Because $J\ne\mathbb Z_d$, the edges $\{\varepsilon_i\}_{i\in J}$ form a disjoint union of
paths on the cycle. A path with $m$ edges spans $E_P$ on its $m+1$ vertices. Together with the pair forms, the $d$
defining vectors are independent and supported on disjoint coordinate sets, so $\dim V_J=d$ and $A_J$ is a nonempty
affine subspace of codimension $d$.

**The coset.** Let $\beta_J+V_J$ be the product of the following factors:

- $\{e_{a_i},e_{b_i}\}$ for each $i\notin J$, i.e. the coset $e_{a_i}+\langle\pi_i\rangle$;
- $C_P$ for each path $P$ of $J$;
- $\{0\}$ on all remaining coordinates.

By Lemma 2 every element has weight at most $(d-|J|)+\sum_P m_P=(d-|J|)+|J|=d$.

Put

$$f_d=\sum_{J\subsetneq\mathbb Z_d}\chi_{\beta_J}\mathbf 1_{A_J}.$$

By Lemma 1 each term has degree at most $d$.

**Claim 1 (disjoint cosets).** If $J\ne K$, pick $i\in J\triangle K$, say $i\notin J$ and $i\in K$. Every element of
$\beta_J+V_J$ restricted to $\{a_i,b_i\}$ is $e_{a_i}$ or $e_{b_i}$. Every element of $\beta_K+V_K$ restricted there
is $0$. So the cosets are disjoint.

**Claim 2 (disjoint live sets).** Let $J\ne K$ and $y=\mathbf 1_J+\mathbf 1_K\in\mathbb F_2^{\mathbb Z_d}$, so
$y\ne0$.

- *Case $y$ not all-ones.* The cyclic word $y$ has some $i$ with $y_{i-1}=1$ and $y_i=0$. Either $i\in J\cap K$,
  in which case $\varepsilon_i\in V_J\cap V_K$, or $i\notin J\cup K$, in which case $\pi_i\in V_J\cap V_K$. That form
  equals $w_J(i)$ on $A_J$ and $w_K(i)$ on $A_K$, and $w_J(i)+w_K(i)=y_{i-1}=1$. Hence $A_J\cap A_K=\emptyset$.
- *Case $y$ all-ones.* Then $K=\mathbb Z_d\setminus J$, and both are nonempty. Let
  $\alpha=\sum_{i\in J}\varepsilon_i$. Since $\sum_{i\in\mathbb Z_d}\varepsilon_i=0$, also
  $\alpha=\sum_{i\in K}\varepsilon_i$, so $\alpha\in V_J\cap V_K$. Its value on $A_J$ plus its value on $A_K$ is
  $$\sum_{i\in J}w_J(i)+\sum_{i\in K}w_K(i)=\#\{i:\ i-1,i\ \text{on the same side}\}+|c|=d-t+|c|,$$
  where $t$ is the number of side changes around the cycle, which is even. So the sum is $\equiv d+|c|\equiv1$, and
  again $A_J\cap A_K=\emptyset$.

**Conclusion.**

- By Claim 2, at each point at most one term is nonzero, and each term is $\pm1$ on its live set. So
  $f_d\in\{0,\pm1\}$, $\|f_d\|_\infty=1$, and $\Pr[f_d\ne0]=(2^d-1)2^{-d}$.
- By Lemma 1 and Claim 1, $\widehat{f_d}$ is $\pm2^{-d}$ on a disjoint union of $2^d-1$ cosets of size $2^d$, and
  zero elsewhere.

Hence

$$R(f_d)=\big((2^d-1)2^d\,2^{-dp_d}\big)^{1/p_d}=\big(4^d(1-2^{-d})\big)^{\frac{d+1}{2d}}2^{-d}=2(1-2^{-d})^{\frac{d+1}{2d}}.$$

**Comparison with the address function.** $R_d>2^{(d-1)/d}$ is equivalent to $(1-2^{-d})^{(d+1)/2}>\tfrac12$. By
Bernoulli's inequality, with exponent $(d+1)/2\ge1$, the left side is at least $1-(d+1)2^{-d-1}$. That exceeds
$\tfrac12$ exactly when $d+1<2^d$, which holds for every $d\ge2$. $\square$

**Example ($d=2$, after relabelling).** With $u_i=(a_i+b_i)/2$, $v_i=(a_i-b_i)/2$ and $Z_\pm=(z_0\pm z_1)/2$, each
pair of the form $(u,v)$ has exactly one nonzero entry, and that entry is $\pm1$. The function is

$$4f_2=(a_0-b_0)(a_1+b_1)+(z_0+z_1)(a_0+b_0)+(z_0-z_1)(a_1-b_1),$$

that is, $f_2=v_0u_1+Z_+u_0+Z_-v_1$.

- The three products are never nonzero together: the first two clash on pair 0, the first and third on pair 1, and
  the last two on $Z_\pm$.
- So $f_2\in\{0,\pm1\}$, it vanishes on $1/4$ of the cube, and it has 12 coefficients $\pm\frac14$.
- Therefore $R=12^{3/4}/4=1.6118\ldots>\sqrt2$.

## 4. Certificates and checking

Certificates use the format `bh_cube/v1`, specified in `problems/bh_cube/FORMAT.md`. Each one lists the integer
coefficients of $2^df$; the claimed lower bound is a 12-digit truncation.

| file | construction | $n$ | $d$ | terms | sup | SHA-256 |
|---|---|---|---|---|---|---|
| `instances/d2_n5.json` | search | 5 | 2 | 12 | 4 | `13691d36694b3bc1fcc2ed8ba62d785f97c3bb0900c0e9f5f2d91f0a7444385d` |
| `instances/d3_n8.json` | search | 8 | 3 | 56 | 8 | `af475381723e51c82340c9d618f99cf9769af2d0c0062e1ae4cb3eaa91fb2e53` |
| `instances/d4_n11.json` | clique search | 11 | 4 | 240 | 16 | `4d303d60ce61ddea5435648bec1cf7410c8bea7992127c906727c7f9122277ad` |
| `instances/cycle_d5_n15.json` | §3 | 15 | 5 | 992 | 32 | `bbbd7e89949ffa937ad05f2dd57842b6b5c608d82573121e5a6024d6c2dcb8c5` |
| `instances/cycle_d6_n18.json` | §3 | 18 | 6 | 4032 | 64 | `f9b01bcf183adb0aabcead0de136f0512d29f470da06f6dc58eda369c9637459` |
| `instances/cycle_d7_n21.json` | §3 | 21 | 7 | 16256 | 128 | `5b22ba5c5b60dd5f6074ee6547f26f3727b2cf5064dd40998d578eb201a929d3` |

The $d\le4$ files are the functions first found by search. They live on $3d-1$ variables (one pair is replaced by a
single variable) and have the same ratio $R_d$. They are certified numerically but are not literally the §3
functions. `python -m problems.bh_cube.cycle d` writes the §3 function for any $d$; for $d=5,6,7$ its output is
byte-identical to the committed files (re-checked 2026-10-23).

Code hashes (SHA-256):

- `problems/bh_cube/cycle.py`: `a197e7463d6e70eb55a5d6c0988a08a803f5eaef0a857201079ac04c074ccc65`
- `problems/bh_cube/core.py`: `ea580c1271c2a7994d8b0dc666a424f7d7ea1fa2b917067eda385338daee47c1`
- `verify/bh_cube.py`: `b496c595a3bd697c05c3f20404b214be38c7359c81bed5bc2bca2f19960c3e47`

**How the certificates were checked.**

1. **Search-side scorer** (`problems/bh_cube/core.py`). It computes the sup norm exactly:
   an integer fast Walsh–Hadamard transform over all $2^n$ points. It encloses $\|\widehat f\|_{p_d}$ with mpmath
   interval arithmetic at 120 digits, and compares exact rational endpoints with the claim.

   An internal code review found that an earlier version rounded the endpoint through a 53-bit float. That let
   forged claims up to $\sim10^{-16}$ above $R$ pass. It is fixed in the version included here; the committed 12-digit claims were
   never affected.
2. **Brute force.** For the $d\le4$ files, $\|f\|_\infty$ and $R$ were recomputed by naive evaluation at every
   point, and they agree.
3. **Independent verifier** (`verify/bh_cube.py`; red-control tests in
   `verify/test_bh_cube.py`). It was written from `FORMAT.md`
   alone by a different AI agent (model gemini-4-argon) and does not import from `problems/`. It
   uses:
   - an exact WHT with overflow-safe sums, cross-checked in pure Python for $n\le14$;
   - exact rational interval endpoints;
   - red controls: tampered coefficients, wrong degree, claims forged at $+10^{-20}$, malformed JSON, and runs
     under `python -O`.

   Rerun on 2026-10-23, and again on 2026-10-26 from a scratch copy of this repository, it reports **VALID** for all six files. Example output for
   `cycle_d7_n21.json`: $R\in[1.9910564253132065783420630,\,1.9910564253132065783420631]$, sup $=128$.
   `python -I verify/test_bh_cube.py` passes (6 tests).
4. **Proof check.** The same independent agent read the proof (in an internal description equivalent to §3, and the `cycle.py`
   docstring) and reported it complete and sound. It covered coset disjointness, both live-set cases, and the
   parity of $d-t+|c|$ for both parities of $d$. Using only that description, it rebuilt $f_d$ for
   $d=2,\dots,6$ and verified each one.

   This was a reading by an AI agent. **The proof has not been formalised (no Lean or other proof assistant), and
   no human has checked it.**

## 5. What is not claimed

- **No record for $C_{26a}$.** Every $R_d<2$, so $\sup_d R_d=2$, equal to the listed bound. Nothing here beats 2.
- **These {0,±1}-valued functions cannot beat 2.**

  **Proposition 2.**
  - (a) If $F$ is integer-valued of degree at most $d$, then $\widehat F\in2^{-d}\mathbb Z$.
  - (b) Hence every $\{0,\pm1\}$-valued $f$ of degree at most $d$ has $R(f)\le2\sigma^{(d+1)/d}$, where
    $\sigma^2=\Pr[f\ne0]$. If $\sigma<1$ this is less than 2; if $\sigma=1$ ($f$ Boolean) then
    $R(f)\le2^{(d-1)/d}$.

  *Proof.*
  - For (a), write $F$ in the variables $y_i=(1-x_i)/2\in\{0,1\}$:
    $F=\sum_{|T|\le d}c_T\prod_{i\in T}y_i$ with $c_T\in\mathbb Z$ by Möbius inversion. Each
    $\prod_{i\in T}y_i=2^{-|T|}\prod_{i\in T}(1-x_i)$ has coefficients in $2^{-d}\mathbb Z$.
  - For (b), every nonzero $|\widehat f(S)|$ is at least $2^{-d}$, and $p_d<2$. So
    $\sum|\widehat f|^{p_d}\le2^{d(2-p_d)}\sum|\widehat f|^2=2^{d(2-p_d)}\sigma^2$. Taking the $1/p_d$ power, with
    $(2-p_d)/p_d=1/d$, gives $R\le2\sigma^{2/p_d}=2\sigma^{(d+1)/d}$.
  - For Boolean $f$, apply (a) to $(1-f)/2$. This gives $\widehat f\in2^{1-d}\mathbb Z$, and the same computation
    gives $R\le2^{(d-1)/d}$ (ADGP2025, Fact 7). $\square$

  Within sums of codimension-$d$ blocks with disjoint live sets and disjoint cosets, $2^d-1$ blocks is the maximum.
  With $2^d$ blocks $f$ would be Boolean with all coefficients $\pm2^{-d}$, contradicting (b). So $R_d$ is optimal
  for that class.

  I have **not** proved that $R_d$ is optimal among all flat $\{0,\pm1\}$ functions. That would follow from
  $\Pr[f=0]\ge2^{-d}$, which I have not proved.
- **No claim about $\mathrm{BH}^{\le d}_{\{\pm1\}}$ itself.** $R_d$ is a lower bound only; no matching upper bound
  is claimed. Beating 2 would need $f$ with values outside $\{0,\pm1\}$ after scaling. For values in
  $\frac1K\mathbb Z$, the same argument gives $R\le2\sigma^{(d+1)/d}K^{1/d}$. Numerical searches over arbitrary real
  coefficients at $n\le11$ (LP ascent and smooth ascent) found nothing above $R_d$. That is evidence, not a theorem.
- **No claim about homogeneous polynomials or multilinear forms.**
  - For $d\ge3$, $f_d$ mixes the degrees $d,d-2,\dots$, so it says nothing about Pellegrino–Teixeira's conjecture
    for $m$-linear forms.
  - $f_2$ is a square-free 2-homogeneous polynomial. I did not search the literature on real homogeneous polynomial
    BH constants for square-free examples, so I make no claim relative to it. (The real-polynomial lower bounds I
    know of use squared variables, which have no analogue on the cube.)
- **No novelty claim beyond the search in §1.** A construction like this may exist in a source I did not find.

## 6. Reproduction

See `README.md` for the commands, which have been run from a scratch copy, and `SHA256SUMS` for the file hashes.
Running everything takes well under a minute on one core; $d=7$ ($2^{21}$ points) dominates. The search code that led
to the construction (clique search over tight block options, LP and smooth ascent) is not included. It is not needed
to check the theorem.

## 7. AI-use disclosure

The construction, the proof and the code were found and written with Claude (Anthropic, Claude Opus 5.5), directed
by Dennis Koshta. The route was numerical search first, then structure extraction, then a brute-force search over
sign rules, which produced the closed form in §3. The independent verifier and the proof reading in §4 were done by
a separate AI agent (gemini-4-argon). **No human has checked the proof or the certificates.**
