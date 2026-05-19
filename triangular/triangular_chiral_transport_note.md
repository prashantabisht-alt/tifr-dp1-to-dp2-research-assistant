# Triangular chiral RTW: transport checkpoint

**Date:** 2026-05-19  
**Status:** first single-particle transport layer verified.

This note records the clean result before moving on to real-space propagators,
KMC, first-passage, and edge modes.

---

## 1. Model

The walker lives on a triangular lattice with axial site coordinates

\[
\mathbf n=(n_1,n_2)
\]

and six internal director states

\[
m=0,1,2,3,4,5.
\]

The six lattice steps are

\[
\Delta_0=(1,0),\quad
\Delta_1=(0,1),\quad
\Delta_2=(-1,1),
\]

\[
\Delta_3=(-1,0),\quad
\Delta_4=(0,-1),\quad
\Delta_5=(1,-1).
\]

One run event moves only in the current director:

\[
(\mathbf n,m)\to(\mathbf n+\Delta_m,m)
\]

at rate \(v\).

Tumble events change the director at the same site:

\[
m\to m+1 \quad \text{at rate } \gamma_+,
\]

\[
m\to m-1 \quad \text{at rate } \gamma_-,
\]

\[
m\to m+3 \quad \text{at rate } \gamma_r.
\]

The chirality bias is parameterized by

\[
\gamma_+=\frac{\gamma}{2}(1+b),\qquad
\gamma_-=\frac{\gamma}{2}(1-b),
\qquad b\in[-1,1].
\]

Achiral:

\[
b=0,\qquad \gamma_+=\gamma_-.
\]

Chiral:

\[
b\ne0,\qquad \gamma_+\ne\gamma_-.
\]

This is the triangular-lattice analogue of the Mallikarjun-Pal chiral
run-and-tumble walker: ref. 44 has four orientations; here we have six.

---

## 2. Master Equation

Let \(P_m(\mathbf n,t)\) be the probability of being at site \(\mathbf n\)
with director \(m\). Then

\[
\partial_t P_m(\mathbf n,t)
=
vP_m(\mathbf n-\Delta_m,t)
+\gamma_+P_{m-1}(\mathbf n,t)
+\gamma_-P_{m+1}(\mathbf n,t)
+\gamma_rP_{m+3}(\mathbf n,t)
-(v+\gamma_++\gamma_-+\gamma_r)P_m(\mathbf n,t).
\]

All director indices are modulo 6.

The equation is just bookkeeping:

\[
\text{rate of change}=\text{arrivals}-\text{departures}.
\]

---

## 3. Bloch Matrix

Use the Fourier convention

\[
\widetilde P_m(\mathbf k,t)=
\sum_{\mathbf n} e^{i\mathbf k\cdot\mathbf n}P_m(\mathbf n,t).
\]

Then

\[
\partial_t \widetilde{\mathbf P}(\mathbf k,t)
=
M(\mathbf k)\widetilde{\mathbf P}(\mathbf k,t).
\]

The diagonal entries are

\[
M_{m,m}(\mathbf k)
=
v e^{i\mathbf k\cdot\Delta_m}
-
(v+\gamma_++\gamma_-+\gamma_r).
\]

The off-diagonal entries are

\[
M_{m+1,m}=\gamma_+,\qquad
M_{m-1,m}=\gamma_-,\qquad
M_{m+3,m}=\gamma_r.
\]

The convention is:

\[
\frac{dP}{dt}=MP,
\]

so **columns are source states** and **rows are destination states**.

At \(\mathbf k=0\), the run contribution cancels:

\[
v e^{i\mathbf k\cdot\Delta_m}-v=0.
\]

Therefore probability conservation requires

\[
\sum_m M_{m,j}(0)=0.
\]

This is verified numerically.

---

## 4. Symmetry Checks

The model has triangular \(C_6\) rotational symmetry even when chiral.
Numerically,

\[
\mathrm{spec}\,M(R_{60}\mathbf k)
=
\mathrm{spec}\,M(\mathbf k)
\]

to machine precision.

Chirality breaks mirror symmetry, not \(C_6\) symmetry. Under mirror reflection,
clockwise and counter-clockwise turns swap:

\[
\gamma_+\leftrightarrow\gamma_-.
\]

Numerical check at \(b=0.35\):

\[
\text{mirror covariance error with same rates}=0.28,
\]

\[
\text{mirror covariance error with swapped rates}=0.
\]

At \(b=0\), mirror symmetry is restored:

\[
\text{mirror covariance error}=0.
\]

---

## 5. Diffusion Tensor

The long-time transport is described by a diffusion tensor

\[
D=
\begin{pmatrix}
D_{xx} & D_{xy}\\
D_{yx} & D_{yy}
\end{pmatrix}.
\]

The symmetric/even part controls ordinary spreading:

\[
D_{\rm even}=\frac{D_{xx}+D_{yy}}{2}.
\]

The antisymmetric/odd part is the chiral response:

\[
D_{\rm odd}=\frac{D_{xy}-D_{yx}}{2}.
\]

For \(C_6\) symmetry,

\[
D_{xx}=D_{yy},\qquad D_{xy}=-D_{yx}.
\]

The eigenvalue curvature gives only the symmetric part:

\[
\lambda_0(\mathbf k)
=
-D^{\rm sym}_{ij}k_i k_j+O(k^3).
\]

The odd part cannot be seen from this scalar dispersion because

\[
k_iD^{\rm odd}_{ij}k_j=0.
\]

The full tensor is computed by a Green-Kubo / pseudoinverse calculation. For a
Poisson jump process, the total diffusion has two parts:

\[
D_{\rm total}=D_{\rm persistent}+D_{\rm jump}.
\]

The jump-noise term is essential; without it the Green-Kubo even diffusion does
not match the eigenvalue curvature.

---

## 6. Verified Numerical Values

For the test parameters

\[
v=1,\qquad \gamma=0.8,\qquad \gamma_r=0.2,
\]

the code gives:

| \(b\) | \(D_{\rm even}\) | \(D_{\rm odd}\) |
|---:|---:|---:|
| 0.0 | 0.875000 | 0.000000 |
| 0.35 | 0.822410 | 0.173503 |
| 0.8 | 0.672297 | 0.292576 |

Interpretation:

- at \(b=0\), the walker is achiral and \(D_{\rm odd}=0\);
- as \(|b|\) grows, \(D_{\rm odd}\) grows;
- chirality reduces ordinary spreading but increases the transverse/chiral
  response.

At the fully chiral point

\[
v=1,\qquad \gamma=1,\qquad b=1,\qquad \gamma_r=0,
\]

the exact values are

\[
D_{\rm even}=\frac12,
\qquad
D_{\rm odd}=\frac{\sqrt3}{4},
\qquad
\frac{D_{\rm odd}}{D_{\rm even}}=\frac{\sqrt3}{2}.
\]

These exact values are reproduced numerically.

---

## 7. Verification Status

Verification script:

```text
triangular/verify_diffusion_tensor.py
```

Current result:

```text
47/47 checks passed
```

The checks include:

- probability conservation at \(\mathbf k=0\);
- eigenvalue curvature \(D_{\rm even}\) versus Green-Kubo \(D_{\rm even}\);
- \(C_6\) isotropy:
  \[
  D_{xx}=D_{yy};
  \]
- antisymmetric odd part:
  \[
  D_{xy}=-D_{yx};
  \]
- \(D_{\rm odd}=0\) at \(b=0\);
- antisymmetry under chirality reversal:
  \[
  D_{\rm odd}(b)=-D_{\rm odd}(-b);
  \]
- analytical DFT calculation of the circulant rotation generator;
- exact fully chiral values at \(b=1,\gamma_r=0\).

The DFT cross-check is important because it independently diagonalizes the
director-tumbling generator. A sign error in the DFT pseudoinverse originally
flipped \(D_{\rm odd}\), but the corrected exponent

\[
\omega^{k(m'-m)}
\]

now agrees with Green-Kubo to machine precision.

---

## 8. Current Figures

Bias scan:

```text
triangular/D_vs_bias.png
```

Reversal-rate scan:

```text
triangular/D_vs_bias_gamma_r_scan.png
```

Square-vs-triangular comparison:

```text
triangular/outputs/square_vs_triangular_diffusion.png
```

The square-vs-triangular comparison is the first direct indication of new
triangular-lattice physics: at matched rates, the \(C_6\) model has a stronger
odd-diffusive response than the \(C_4\) square-lattice analogue.

---

## 9. Next Step

The next physics target is the exact real-space propagator:

\[
P(\mathbf n,t)
=
\frac{1}{L^2}
\sum_{\mathbf k}
e^{-i\mathbf k\cdot\mathbf n}
\mathbf 1^T e^{M(\mathbf k)t}\widetilde{\mathbf P}(\mathbf k,0).
\]

This will produce the first real-space probability cloud for the triangular
chiral RTW, before writing the Fortran KMC.

