# 2026-05-19 — Kabir discussion note: Ref. 44 is continuous-space

Today I told Kabir that the Mallikarjun-Pal chiral run-and-tumble paper
(TCRW ref. 44) is **not** a lattice model in the spatial coordinate.
Kabir had been thinking of it as closer to a square-lattice model, but the
paper/code uses continuous ballistic motion with four discrete orientation
states.

Important correction for our framing:

\[
\text{Mallikarjun-Pal / ref. 44}
=
\text{continuous space}
+
\text{four internal directions}.
\]

Our model is different:

\[
\text{our Phase 1 model}
=
\text{triangular lattice}
+
\text{six internal directions}.
\]

So the relation is:

- Ref. 44 gives the **chirality convention** and observable list.
- Ref. 44 does **not** give the triangular-lattice transition matrix.
- We should not call our work a direct triangular-lattice version of their
  lattice model, because their spatial coordinate is not lattice-based.
- The careful phrase is:

> We build a triangular-lattice, six-director analogue of the
> Mallikarjun-Pal chiral run-and-tumble particle, preserving their
> left/right/reversal reorientation structure but replacing continuous
> ballistic motion by Poisson hops on a triangular lattice.

This distinction matters for:

1. **Short-time MSD.**  
   Their continuous CRTP has ballistic short-time behavior

   \[
   \langle r^2(t)\rangle \sim v^2t^2.
   \]

   Our lattice Poisson-jump RTW has jump-noise short-time behavior

   \[
   \langle r^2(t)\rangle \sim vt.
   \]

2. **MFPT comparison.**  
   Their absorbing boundaries are continuous vertical walls \(x=\pm L\).
   Our triangular problem should use lattice strip/domain boundaries
   such as zigzag, armchair, rhombus, or hexagonal domains.

3. **Novelty.**  
   Our novelty is not merely "same paper on triangular." It is:

   \[
   4\text{-orientation continuous CRTP}
   \longrightarrow
   6\text{-orientation triangular-lattice chiral RTW}.
   \]

This note should be kept in mind before writing the intro or explaining the
project to Kabir again.
