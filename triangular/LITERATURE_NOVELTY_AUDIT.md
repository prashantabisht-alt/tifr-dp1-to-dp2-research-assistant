# Literature novelty audit — triangular chiral RTW

**Date:** 2026-05-19  
**Question:** Is our Phase 1 model already done somewhere?

## Target model being checked

Our Phase 1 target is deliberately narrow:

- continuous-time Poisson run-and-tumble walker,
- fixed triangular / six-neighbor lattice,
- six internal director states,
- deterministic hop in the current director at rate `v`,
- left/right turning rates
  \[
  \gamma_+ = \frac{\gamma(1+b)}{2}, \qquad
  \gamma_- = \frac{\gamma(1-b)}{2},
  \]
- optional reversal rate \(\gamma_r\),
- observables: diffusion tensor including \(D_{\rm odd}\), propagator, first
  passage/search, and later OBC edge/localization diagnostics.

## Strict Verdict

No exact duplicate was found.

The closest papers split the ingredients:

- Mallikarjun-Pal: continuous-time chiral RTW rates, but continuum space and
  four orientations.
- Wójcik-Kalz: lattice odd diffusion and edge-current/topology framing, but
  square lattice and discrete-time coin-step dynamics.
- Osat-Speck TCRW: square-lattice edge/topology machinery, but discrete-time
  Markov transition matrix and different update rule.
- Marris/Giuggioli, Batchelor-Henry, Gilbert-Sanders: lattice geometry,
  persistence, and first passage, but not our chiral six-director RTW.
- Oropesa-de Castro-Löwen-Liarte: triangular run-and-tumble-like motion on
  adaptive networks, but not our fixed-lattice chiral single-walker generator.

## Safe Novelty Claim

Use this wording:

> We study a continuous-time six-director rotation-chiral run-and-tumble
> walker on the fixed triangular lattice, and compute its odd diffusivity,
> propagator, first-passage/search observables, and boundary signatures.

Even sharper:

> Existing topological/chiral lattice-walker models are mostly square-lattice
> discrete-time coin or Markov walks. Our contribution is the continuous-time
> \(C_6\) triangular run-and-tumble generator and its transport/search/edge
> consequences.

## Unsafe Claims

Avoid these:

- "first chiral random walk"
- "first triangular random walk"
- "first triangular active walker"
- "first first-passage theory on triangular lattices"
- "Wójcik-Kalz is the same model on square"
- "Mallikarjun-Pal is the square-lattice version"

## Highest-Risk / Must-Cite Papers

| Paper | Why it matters | Exact relation to us |
|---|---|---|
| Wójcik & Kalz 2026, *The chiral random walk* | Closest lattice odd-diffusion/topological edge-current reference | Square, discrete-time coin-step; not our continuous-time RTW |
| Mallikarjun & Pal 2023, *Chiral run-and-tumble walker* | Same left/right/reversal rate mechanism and MFPT/search style | Continuum space, four orientations; not lattice/C6 |
| Osat, Meyberg, Metson & Speck 2026, *Topological chiral random walker* | Strongest square-lattice OBC/topology method template | Discrete-time square model, different update rule |
| Hargus, Epstein & Mandadapu 2021, *Odd Diffusivity of Chiral Random Motion* | Owns the odd-diffusivity language and Green-Kubo framing | Conceptual/method reference |
| Oropesa, de Castro, Löwen & Liarte 2026 | Triangular run-and-tumble-adjacent active particles on networks | Adaptive/trail-mediated network; not fixed-lattice chiral RTW |
| Marris & Giuggioli 2024 | Persistent/anti-persistent lattice first-passage theory | Important FPT background; not chiral C6 RTW |
| Marris, Sarvaharman & Giuggioli 2023 | Exact dynamics in hexagonal/honeycomb domains | Geometry/boundary reference; terminology trap |
| Ouvry & Polychronakos 2019/2020 | Contains "triangular chiral walk" language | Algebraic area / Hofstadter / combinatorics; title trap |
| Larralde 1997 | Historical chiral persistent random walk transport | Not triangular/C6, but important for chiral persistent-walk history |
| Batchelor & Henry 2002 | Exact triangular random walk with absorbing boundaries | Non-active/non-persistent; method/background |

## Papers To Add Or Skim Next

These are not currently all in `triangular/papers/` and should be added or
skimmed before writing the intro/FPT section.

1. Larralde 1997, *Transport properties of a two-dimensional "chiral"
   persistent random walk*. This protects us from overclaiming chiral
   persistent-walk transport.
2. Ouvry & Polychronakos 2020, *Lattice walk area combinatorics...*. This is
   the follow-up title trap to the 2019 triangular chiral-walk combinatorics
   paper.
3. Sarmiento et al. 2025, *First-Passage-Time Asymmetry for Biased
   Run-and-Tumble Processes*. This is relevant to our FPT/search section.
4. Batchelor & Henry 2002, *Exact solution for random walks on the triangular
   lattice with absorbing boundaries*. This is a baseline for triangular
   absorbing-boundary geometry.
5. Gilbert & Sanders 2010, *Diffusion coefficients for multi-step persistent
   random walks on lattices*. We already have the 2009 Lorentz-gas paper; the
   2010 diffusion-coefficient paper is also relevant.

## Separation Of Branches

Phase 1 single-walker branch:

- chiral RTW,
- one particle,
- transport/FPT/search/edge signatures,
- no pressure, no equation of state, no MIPS.

Later many-particle branch:

- hard hexagon / hard-core triangular gas,
- activity-induced deviations from equilibrium equation of state,
- pressure and MIPS,
- Whitelam-Klymko-Mandal, Solon et al., Jaleel-Mandal-Thomas-Rajesh,
  Panagiotopoulos, Baxter hard hexagon.

Do not mix these branches in the same derivation.

## Intro Skeleton

Random walks on triangular/hexagonal lattices and persistent lattice walks
have a long history, including exact results for bounded domains and
first-passage statistics. Recent work has also introduced lattice chiral
random walks with odd diffusion and topological edge currents on square
lattices, and continuous-space chiral run-and-tumble particles with
first-passage/search optimization. Here we combine these threads in a
different setting: a continuous-time six-director run-and-tumble process on
the fixed triangular lattice, where the \(C_6\) geometry changes the transport
tensor and boundary response.

