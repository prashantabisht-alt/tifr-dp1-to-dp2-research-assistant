# TCRW Fig 4 — Knowledge brief before Fortran port

*Scope:* compile everything needed to reproduce Fig 4 in Fortran. What the
paper shows, what the Python does, whether Python has any issues, and how to
structure the Fortran implementation panel by panel.

---

## 1. Paper physics — what Fig 4 is actually doing

Fig 4 shifts the paper from Markov-chain dynamics (Fig 1–3) to the **spectral
/ non-Hermitian band-theory description** of the one-walker transition matrix
$P$. The central object is the 4-band Bloch matrix

$$
P(\mathbf{k}) =
\begin{pmatrix}
0 & R_1 + C_1 e^{+ik_x} & 0 & R_2 + C_2 e^{-ik_x} \\
R_2 + C_2 e^{+ik_y} & 0 & R_1 + C_1 e^{-ik_y} & 0 \\
0 & R_2 + C_2 e^{+ik_x} & 0 & R_1 + C_1 e^{-ik_x} \\
R_1 + C_1 e^{+ik_y} & 0 & R_2 + C_2 e^{-ik_y} & 0
\end{pmatrix},
\tag{paper Eq. 1}
$$

with
$$
C_1=(1-\omega)(1-D_r),\quad C_2=\omega(1-D_r),\quad
R_1=\omega D_r,\quad R_2=(1-\omega) D_r.
$$

Interpretation:
- row/column index $d\in\{0,1,2,3\}\equiv\{\uparrow,\rightarrow,\downarrow,\leftarrow\}$
- only $d\to d\pm1$ mod 4 matrix elements are non-zero (every step is exactly one rotation)
- $C$-entries carry a lattice phase $e^{\pm i k_\alpha}$ because the chiral step translates; $R$-entries have no phase because the noise step is pure rotation.
- rotation sense: the chiral rule in the paper rotates with one chirality
  (say CW with probability $\omega$) and the noise rule rotates with the
  **opposite** chirality — this is what breaks degeneracy and gives a
  topological gap.

$P(\mathbf{k})$ is **non-Hermitian** (not even normal) because the two rules
rotate oppositely, so the left and right eigenvectors differ and spectra can
go complex even with real rates. That is the whole reason Fig 4 is interesting:
TCRW is a stochastic model whose generator has the structure of a non-Hermitian
Bloch Hamiltonian with a topological gap.

### Symmetries (paper supp. / needed for Fig 4's interpretation)
Let $\sigma_\alpha$ be Pauli matrices in the sublattice $\mathbb{Z}_2$ (even/odd $d$)
and in a second $\mathbb{Z}_2$ (direction within a sublattice).

- **Inversion**: $I P(\mathbf{k}) I^{-1}=P(-\mathbf{k})$ with $I=\sigma_x\otimes\mathbb{1}$.
- **Time reversal** (classical, anti-unitary): $P(\mathbf{k})^*=P(-\mathbf{k})$.
- **Sublattice / chiral**: $S P(\mathbf{k}) S^{-1}=-P(\mathbf{k})$ with $S=\mathbb{1}\otimes\sigma_z$.

  → the spectrum comes in $\pm\lambda$ pairs, and because $P$ is real at $\mathbf{k}=0$ they also come in complex-conjugate pairs. So the four eigenvalues at every $\mathbf{k}$ fall into two classes: a real pair $\pm a(\mathbf{k})$ and a pure-imaginary pair $\pm i b(\mathbf{k})$, or more generally $\{\lambda,-\lambda,\lambda^*,-\lambda^*\}$.

The gap closes exactly at $\omega=\tfrac12$ (Fig 4b/e) because at
$\omega=\tfrac12$ the chiral-CW and chiral-CCW rates are equal, so the
non-Hermitian splitting vanishes and the real $\pm a$ and imaginary $\pm ib$
pairs meet.

### Panel-by-panel content

| panel | what it shows | x-axis | y-axis |
|-------|----------------|--------|--------|
| 4a | cartoon: at $\omega=1$ the CCW edge current of internal state $d$ is illustrated | — | — |
| 4b | PBC band structure along $\Gamma\text{–}X\text{–}M\text{–}\Gamma$. Four bands; gap closes at $\omega=0.5$. | $\mathbf{k}$ along BZ path | $\mathrm{Re}\lambda,\ \mathrm{Im}\lambda$ |
| 4c | OBC spectrum on tiny $L=2$ lattice: $\mathrm{Re}\lambda$ vs $\mathrm{Im}\lambda$. Introduces edge modes. | $\mathrm{Re}\lambda$ | $\mathrm{Im}\lambda$ |
| 4d | OBC at fixed $\omega=1$, $L=10$: $\mathrm{Re}\lambda$ vs $D_r$. Shows coalescence / EP formation as $D_r\to0$. | $D_r$ | $\mathrm{Re}\lambda$ |
| 4e | OBC at fixed $D_r=0.1$, $L=10$: $\mathrm{Re}\lambda$ vs $\omega$. Gap-closing at $\omega=0.5$. | $\omega$ | $\mathrm{Re}\lambda$ |
| 4f | OBC complex plane at $\omega=1$, colored by edge-localization $\rho_\partial=\sum_{\text{boundary}}\|v\|^2$. Shows ring of edge states. | $\mathrm{Re}\lambda$ | $\mathrm{Im}\lambda$ |
| 4g | Same but at fixed $D_r$, varying $\omega$. | $\mathrm{Re}\lambda$ | $\mathrm{Im}\lambda$ |
| 4h | Hybrid BC (periodic in $y$, open in $x$): spectrum as function of $k_y$; edge bands visible. | $k_y$ | $\mathrm{Re}\lambda$ |
| 4i | PBC bands re-plotted in the $(\cos k_x,\cos k_y)$ plane ("band circle"), colored by Re$\lambda$ to visualize how bands wrap. | $\cos k_x$ | $\cos k_y$ |

---

## 2. Python implementation audit

### 2.1 File → panel map

| file | functions used | panels | status |
|------|-----------------|--------|--------|
| `tcrw_spectrum.py` | `build_Pk`, `pbc_band_structure` (Hungarian continuity), `pbc_full_bz`, `obc_spectrum`, `classify_eigenvalues`, `fig4b_band_structure`, `fig4_obc_complex_plane`, `fig4_gap_closing`, `fig4_pbc_vs_obc` | 4b, 4c, 4f, 4g | ✅ audited |
| `tcrw_obc.py` | `build_transition_matrix(ω,D_r,L)`, `build_transition_matrix_generic(ω,D_r,mask)` | 4c–4g (OBC) | ✅ audited |
| `tcrw_fig4_extra.py` | `compute_edge_weight`, `obc_spectrum_with_weights`, `fig4d_spectrum_vs_Dr`, `fig4e_spectrum_vs_omega`, `fig4i_band_circle`, `fig8_spectrum_grid` | 4d, 4e, 4i, Fig 8 grid | ✅ audited |
| `tcrw_fig4h_hpbc.py` | `build_hpbc_matrix(ω,D_r,L,k_y)`, `fig4h_hpbc_spectrum` | 4h only | ✅ audited |

### 2.2 Central ingredient: `build_Pk`

```python
def build_Pk(omega, D_r, kx, ky):
    C1 = (1 - omega) * (1 - D_r)   # chiral CCW
    C2 = omega       * (1 - D_r)   # chiral CW
    R1 = omega       * D_r         # noise  CCW
    R2 = (1 - omega) * D_r         # noise  CW
    ekx = np.exp(1j*kx); eky = np.exp(1j*ky)
    P = np.array([
        [0,             R1 + C1*ekx,   0,             R2 + C2/ekx],
        [R2 + C2*eky,   0,             R1 + C1/eky,   0         ],
        [0,             R2 + C2*ekx,   0,             R1 + C1/ekx],
        [R1 + C1*eky,   0,             R2 + C2/eky,   0         ],
    ], dtype=complex)
    return P
```

Row-by-row comparison with paper Eq (1) → exact match. Conventions:
- CCW: $d\to(d-1)\bmod 4 \equiv (d+3)\bmod 4$
- CW:  $d\to(d+1)\bmod 4$
- chiral step uses CW with prob $\omega$ → gives $C_2$ amplitude in "$+1$ direction";
  noise step uses CCW with prob $\omega$ → gives $R_1$ amplitude in "$-1$ direction".

### 2.3 OBC matrix (`tcrw_obc.py::build_transition_matrix`)

State index: `s = d*L^2 + y*L + x`, dimension $4L^2$. Matrix is **column-stochastic**
(each column sums to 1). For each $(x,y,d)$ the column contains:
- noise term: weight $R_1$ into $(x,y,(d-1)\bmod4)$ and weight $R_2$ into $(x,y,(d+1)\bmod4)$;
- chiral term if $(x+DX[d],\,y+DY[d])$ is inside the $L\times L$ box:
  weight $C_1$ into $((x+DX,y+DY),(d-1)\bmod4)$ and $C_2$ into $((x+DX,y+DY),(d+1)\bmod4)$;
- if the chiral step is **blocked** (boundary): total weight $(1-D_r)$ stays at
  $(x,y,d)$ — i.e. a diagonal element, no rotation.

This "stay" term is what breaks the sublattice symmetry at the boundary and
lets edge modes appear. It is the same convention the Fortran `tcrw_step.f90`
kernel uses (`step_type=2` branch: no move, no rotation).

### 2.4 HPBC matrix (`tcrw_fig4h_hpbc.py::build_hpbc_matrix`)

$(4L)\times(4L)$ complex matrix parametrized by $k_y\in[-\pi,\pi]$:
- state = $(x,d)$ with $x\in\{0,\dots,L-1\}$, $d\in\{0,1,2,3\}$;
- $y$-periodicity is enforced by multiplying chiral transitions whose step has
  $DY[d]\neq0$ by $e^{+i k_y DY[d]}$;
- $x$-boundaries are open → blocked chiral steps become a $(1-D_r)$ stay-term
  just like in OBC.

### 2.5 Numerical methods & resolution

Per `TCRW_ACCURACY_REPORT.md` the final Python resolutions are:
- 4(b): PBC bands — $N_k=100$ per BZ segment with Hungarian band-continuity tracking.
- 4(d/e): OBC sweeps — 50 parameter points, $L=12$.
- 4(f/g): OBC complex plane — $L=12$, edge weights from right eigenvectors.
- 4(h): HPBC — $N_{k_y}=200$, $L=10$.
- 4(i): PBC on $(\cos k_x,\cos k_y)$ — $N_k=150$ per segment.

All solves use `numpy.linalg.eig` on a **general complex** matrix (returns
both right eigenvectors). OBC matrix is stored dense (up to $\sim 576\times576$
at $L=12$ — trivial).

### 2.6 Issues found

- **Python, Fig 4 panels:** nothing incorrect. `TCRW_ACCURACY_REPORT.md`
  logged each panel as verified, and a direct row-by-row comparison of
  `build_Pk` with paper Eq (1) reproduces every entry including the
  direction of each $e^{\pm ik_\alpha}$ phase.
- **Cosmetic inconsistency in existing Fortran `tcrw_fig3hij.f90`
  header comment:** it says "Noise rotates CCW with prob $(1-\omega)$, CW
  with prob $\omega$". That is backwards — `tcrw_step.f90` does the opposite
  (`r_rot<omega` → CCW for the noise branch). The *physics worked example*
  further down in the same comment block ("$\omega=1$ noise from $\leftarrow$
  gives $\downarrow$") is actually consistent with the kernel (CCW of
  $\leftarrow$ is $\downarrow$). So only the single-sentence description is
  wrong; no numbers change. Worth fixing when we next touch that file so
  the repo's internal documentation stays clean.

---

## 3. Fortran implementation plan

### 3.1 Ordering

Go in order of increasing LAPACK pain:

1. **4(b) PBC bands** — $4\times4$ complex matrix only. Can even be done
   analytically via the sublattice block structure, but cleanest is `zgeev`
   on 4×4, which costs nothing. Produces 4 bands along $\Gamma$–$X$–$M$–$\Gamma$.
   Use this as the LAPACK smoke-test.
2. **4(c) OBC, $L=2$** — 32×32 complex. `zgeev` trivial. Good first OBC run.
3. **4(d,e) OBC sweeps, $L=10$** — 400×400 complex, parameter sweep of
   25–50 points. Each `zgeev` is ~$10^8$ flops; negligible at this size.
4. **4(f,g) OBC complex plane, $L=10$** — same matrix as (d,e) but we also
   need **right eigenvectors** to compute the edge-localization weight
   $\rho_\partial=\sum_{(x,y)\in\partial}|v_{(x,y,d)}|^2$.
5. **4(h) HPBC, $L=10$** — 40×40 complex matrix parametrized by $k_y$;
   $N_{k_y}=200$ sweeps. Right eigenvectors needed for edge coloring.
6. **4(i) PBC band circle** — only requires `build_Pk` + $N_k^2$ eigensolves of
   the 4×4; plot in $(\cos k_x,\cos k_y)$ plane. No LAPACK beyond the 4×4.

### 3.2 Files I'll add under `fortran_reproduction/`

Single shared matrix-build module to avoid duplication:

```
tcrw_bloch.f90       ! build_Pk (4x4 complex), column-stochastic check
tcrw_obc_matrix.f90  ! build_P_obc(omega,Dr,L)  returns dense 4L^2 x 4L^2
tcrw_hpbc_matrix.f90 ! build_P_hpbc(omega,Dr,L,ky) returns 4L x 4L
tcrw_eig.f90         ! thin wrapper around zgeev (eigvals + right eigvecs)

tcrw_fig4b.f90       ! PBC bands Γ-X-M-Γ
tcrw_fig4c.f90       ! OBC L=2 complex plane
tcrw_fig4de.f90      ! OBC sweep vs D_r and vs omega (combined)
tcrw_fig4fg.f90      ! OBC complex plane w/ edge weights
tcrw_fig4h.f90       ! HPBC sweep vs ky
tcrw_fig4i.f90       ! Band circle
```

Each dumps `*.txt` data; gnuplot files match the existing `tcrw_figXY.gnu`
pattern (qt/pdf/both modes), added as `tcrw-fig4X-build/run/plot/all`
targets in `run.sh`. LAPACK linkage already auto-detected by `run.sh`
(Accelerate on macOS, `-llapack -lblas` on Linux).

### 3.3 LAPACK calls we need

Only two calls are required for the whole of Fig 4:

- `ZGEEV('N','V', n, A, lda, W, VL, 1, VR, ldvr, WORK, lwork, RWORK, info)`
  — right eigenvectors only (`VL` side disabled). Needed for 4(c)–4(h).
- `ZGEEV('N','N', ...)` — eigenvalues only. Cheap variant for 4(b), 4(d), 4(e), 4(i).

Workspace query: call `ZGEEV` once with `lwork=-1` to have LAPACK return the
optimal `lwork` in `WORK(1)`, then allocate and call again. Standard.

### 3.4 Sanity checks (must pass before any plot is trusted)

For every matrix we build:

1. **Column-stochastic OBC:** $\max_j |1-\sum_i P_{ij}|$ must be $\lesssim 10^{-14}$.
2. **Perron–Frobenius:** exactly one eigenvalue equals 1 (to
   $\sim 10^{-12}$), all others $|\lambda|<1$.
3. **Sublattice-symmetry:** eigenvalues must come in $\pm\lambda$ pairs. A
   practical check: for every $\lambda$ in the spectrum, $-\lambda$ is also
   in the spectrum to tolerance. (Boundary "stay" terms break this slightly
   in OBC — expect the deviation to scale like $1/L$, which is itself a good
   diagnostic.)
4. **Real-at-$\mathbf{k}=0$:** `build_Pk(ω,D_r,0,0)` is real; its spectrum is
   $\{\pm(1-D_r),\pm(1-2\omega)D_r\}$ or similar depending on sign
   convention — easy cross-check against Python.
5. **Bipartite pair classification:** after diagonalization, classify each
   eigenvalue into the real-pair / imaginary-pair sector
   (like `classify_eigenvalues` in Python) to reproduce Fig 4b/e colorings.

### 3.5 Output-format contract with gnuplot

Match existing style so we don't redo plotting plumbing:

- **4(b)**: columns `k_index  k_path  Re(λ_1..4)  Im(λ_1..4)`.
- **4(c/f/g)**: columns `Re(λ)  Im(λ)  ρ_edge` (one point per eigenvalue per parameter value).
- **4(d)**: columns `D_r  Re(λ_1..4L²)` — transpose-friendly wide file.
- **4(e)**: columns `ω   Re(λ_1..4L²)`.
- **4(h)**: columns `k_y  Re(λ_1..4L)  ρ_edge_1..4L`.
- **4(i)**: columns `cos(k_x)  cos(k_y)  Re(λ_1..4)`.

### 3.6 Risks / what could bite us

- **LAPACK not installed** on the target machine. `run.sh` has a lapack_check
  hook; if it fails we stop cleanly and print the install one-liner.
- **Eigenvalue ordering is non-deterministic across LAPACK builds.** For 4(b)
  we need to track bands across $\mathbf{k}$. Python uses Hungarian matching on
  $|\lambda_i(k_n)-\lambda_j(k_{n+1})|$; I'll port the same algorithm (it's ~30
  lines of Fortran, just an $O(n^3)$ assignment on $n=4$).
- **Right eigenvectors of a non-Hermitian matrix are not orthogonal.** For the
  edge-weight coloring this is fine (we only need the $L^2$ norm of the
  eigenvector's real-space profile), but it matters if we later compute a
  Zak phase — that needs the **biorthogonal** basis, i.e. also the left
  eigenvectors from `ZGEEV(..., 'V','V', ...)`. Keep this in mind for Fig 8,
  but Fig 4 does not need it.
- **Mass-conservation of the "stay" term**: easy to double-count when coding
  OBC. Build rule: first build the bulk (noise + chiral) columns, then for
  each blocked chiral move add a diagonal $(1-D_r)$ at the source state.
  Column sums must be 1 exactly (integer arithmetic in the rate definitions,
  then cast to real).

### 3.7 First concrete deliverable for the next step

Panel 4(b). It needs only:
- `tcrw_bloch.f90::build_Pk(ω,D_r,kx,ky,P)`  (returns 4×4 complex matrix)
- a 4×4 `ZGEEV` call
- a $\Gamma$–$X$–$M$–$\Gamma$ path generator with `Nk=100` per segment
- a Hungarian band-continuity tracker on 4 bands
- `tcrw_fig4b.gnu` that draws 4 colored curves with a twin-axis Re/Im split,
  matching `tcrw_spectrum.py::fig4b_band_structure`.

Once 4(b) reproduces the Python curves to eyeball precision at
$(\omega,D_r)\in\{(1,0.1),(0.5,0.1),(0.2,0.1)\}$, the rest of Fig 4 is
mechanical.

---

## 4. Takeaways

- Python is clean for all of Fig 4. No physics fixes required before porting.
- One cosmetic comment error in `tcrw_fig3hij.f90` header, worth fixing when we touch it next; does not affect any output.
- The Fortran port needs `ZGEEV` only. The largest matrix is $4L^2\times4L^2=576\times576$ at $L=12$, trivial for LAPACK.
- Natural first panel to implement: 4(b) — $4\times4$ matrix, Hungarian tracking, no new LAPACK complications. Gives us the smoke test for the whole toolchain.
