# Ref. 44 companion-code audit: Mallikarjun-Pal CRTP

**Source repo:** `RMStatMech/Chiral-Run-and-Tumble-Particle`  
**Local copy:** `triangular/papers/RMStatMech_companion_code/`  
**Remote HEAD checked:** `9de10b900248b4f05d53c00318053d4aa6cf38f3`  
**Audit date:** 2026-05-19

This is the companion code for the square/4-orientation chiral
run-and-tumble particle paper that Kabir pointed us toward as ref. 44.

The key point:

\[
\text{Ref. 44} = 4\text{-state square/continuous-space CRTP},
\]

while our current project is

\[
\text{ours} = 6\text{-state triangular-lattice chiral RTW}.
\]

So this repo is a **template and benchmark**, not code to copy directly.

---

## 1. What the Python notebook contains

File:

`triangular/papers/RMStatMech_companion_code/Chiral_RTP.ipynb`

The notebook has two main simulation blocks.

### 1.1 Marginal position distribution

It starts a CRTP at the origin with isotropic initial orientation. The four
internal states are

\[
\{0,1,2,3\}
\quad\leftrightarrow\quad
\left\{0,\frac{\pi}{2},\pi,\frac{3\pi}{2}\right\}.
\]

The velocity direction is

\[
\hat{\mathbf e}_m=(\cos\theta_m,\sin\theta_m).
\]

For the marginal \(x\)-distribution they only track the \(x\)-component:

\[
f_x(0)=1,\qquad f_x(1)=0,\qquad f_x(2)=-1,\qquad f_x(3)=0.
\]

In that first block they set

\[
\Gamma_2=0,
\]

and use only two competing exponential clocks:

\[
\Gamma_1 = \text{left-turn rate},\qquad
\Gamma_3 = \text{right-turn rate}.
\]

Example in the notebook:

\[
(\Gamma_1,\Gamma_3)=(0.4,1.0).
\]

The code uses the first-reaction method: draw one exponential waiting time for
each possible tumble channel, take the minimum, move ballistically for that
time, then update the orientation according to which clock won.

### 1.2 MFPT between absorbing vertical boundaries

The second Python block computes the mean first-passage time between absorbing
boundaries at

\[
x=\pm L.
\]

Here all three rates are present:

\[
(\Gamma_1,\Gamma_2,\Gamma_3)=(0.4,0.2,1.0)
\]

in the notebook example.

The event channels are:

\[
m\to m+1 \quad \text{at rate } \Gamma_1,
\]

\[
m\to m+2 \quad \text{at rate } \Gamma_2,
\]

\[
m\to m-1 \quad \text{at rate } \Gamma_3.
\]

This is exactly the square analogue of our triangular event structure.

---

## 2. What the Mathematica notebook contains

File:

`triangular/papers/RMStatMech_companion_code/CRTP.nb`

The Mathematica notebook is not a trajectory simulator. It is the analytic
MFPT supplement.

Important pieces:

- It is explicitly titled **Mean First-Passage Time (MFPT)**.
- It defines the rates \(\Gamma_1,\Gamma_2,\Gamma_3\).
- It discusses a chirality parameter \(\epsilon\) as the bias between left and
  right turning rates.
- It plots MFPT versus starting position \(x\).
- It scans \(\Gamma_1\) and plots the maximum MFPT versus \(\Gamma_1\).

The relevant text says the bias is reflected in the absolute difference
between

\[
\Gamma_1
\quad\text{and}\quad
\Gamma_3.
\]

They often choose a scale where one of the two rates is fixed to one and the
other is varied, e.g.

\[
\Gamma_1=\epsilon,\qquad \Gamma_3=1.
\]

That is a different parameterization from ours, but the physics is the same:
chirality means left and right turning are unequal.

---

## 3. Translation dictionary

| Ref. 44 square CRTP | Meaning | Our triangular RTW |
|---|---|---|
| \(m=0,1,2,3\) | four internal directions | \(m=0,\ldots,5\) |
| \(\theta_m=m\pi/2\) | square orientation angle | triangular director \(\Delta_m\) |
| speed \(v\) | ballistic/run speed | hop/run rate \(v\) |
| \(\Gamma_1\) | left / CCW turn | \(\gamma_+\) |
| \(\Gamma_2\) | reversal | \(\gamma_r\) |
| \(\Gamma_3\) | right / CW turn | \(\gamma_-\) |
| \(\Gamma_1\ne\Gamma_3\) | chirality | \(\gamma_+\ne\gamma_-\) |
| \(x=\pm L\) absorbing walls | 1D projected MFPT | triangular strip/domain MFPT |

Our preferred symmetric parameterization is

\[
\gamma_+=\frac{\gamma}{2}(1+b),
\qquad
\gamma_-=\frac{\gamma}{2}(1-b),
\qquad
b\in[-1,1].
\]

The inverse map is

\[
\gamma=\gamma_++\gamma_-,
\qquad
b=\frac{\gamma_+-\gamma_-}{\gamma_++\gamma_-}.
\]

So if we want to reproduce their example

\[
\Gamma_1=0.4,\qquad \Gamma_3=1.0,
\]

then in our notation

\[
\gamma_+=0.4,\qquad
\gamma_-=1.0,
\]

\[
\gamma=1.4,\qquad
b=\frac{0.4-1.0}{1.4}\approx -0.428571.
\]

The sign of \(b\) depends only on which direction we call left/CCW.

---

## 4. What we should reuse

### Reuse conceptually

1. **First-reaction KMC structure.**  
   Their simulation draws one exponential waiting time per event channel and
   picks the smallest. This is exactly what our Fortran/Python KMC should do.

2. **Observable list.**  
   They compute:
   - marginal position distributions,
   - MFPT to absorbing walls,
   - MFPT versus starting position,
   - optimal chirality/minimum-search-time curves.

3. **Presentation logic.**  
   Their paper asks: how does chirality change spreading and search? We can ask
   the same question, but with 6 directions and triangular geometry.

### Do not reuse blindly

1. Their Python code is **trajectory code**, not a Bloch-matrix code.
2. Their model has **4 orientations**, not 6.
3. Their spatial motion is continuous ballistic motion; our current code is a
   lattice Poisson-jump RTW.
4. Their MFPT boundaries are \(x=\pm L\); our triangular project has richer
   boundary choices: vertical strip, armchair/zigzag edges, rhombus, hexagon.

---

## 5. Immediate next tasks for our project

1. Keep `triangular_chiral_rtw.py` as the canonical triangular Bloch-matrix
   model.
2. Keep `triangular_chiral_propagator.py` as the exact Fourier propagator.
3. Add current/quiver plots from orientation-resolved probabilities:

\[
\mathbf J(\mathbf n,t)
=
v\sum_{m=0}^{5}P_m(\mathbf n,t)\Delta_m.
\]

4. Build a clean KMC file for the triangular chiral RTW using the first-reaction
   method, then compare KMC against the exact propagator.
5. After that, port the ref. 44 MFPT question:

\[
\text{square: } x=\pm L
\quad\longrightarrow\quad
\text{triangular: strip / zigzag / armchair absorbing boundaries}.
\]

That MFPT section is the clearest direct continuation of ref. 44.
