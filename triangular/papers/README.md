# Papers folder — guide

Reading list for the triangular-lattice active-walker project. Organised
by *phase of the project* rather than by topic.

## Phase 1 — chiral run-and-tumble walker on triangular (now)

### Closest prior art — MUST READ AND MUST CITE

These two define what's already known and what we're competing against
for priority. Read both before any code.

* **`Wojcik_Kalz_2026_chiral_random_walk_odd_diffusion.pdf`** —
  arXiv:2602.09920, **published Feb 2026**. Square-lattice chiral random
  walker with internal director and tunable chirality parameter $p$,
  interpolating from diffusive random walk to deterministic quantum walk.
  They derive odd diffusion and edge currents on the square lattice.
  Important correction: this is a **discrete-time coin-step model**, not
  our continuous-time RTW. We cite them prominently as the closest lattice
  odd-diffusion/topological-edge reference; our Phase 1 is the
  continuous-time triangular counterpart, not a literal extension of the
  same dynamics.
* **`ChiralActive_obstacle_lattices_NatComm_2024.pdf`** —
  Nature Communications, Feb 2024. Chiral active particle in
  *continuum*, navigating around square vs triangular obstacle lattices.
  Independent evidence that triangular geometry produces qualitatively
  different behaviour for chiral particles (the cage effect). Supports
  Kabir's "triangular is different" prediction.

### Novelty-risk / title-near papers — audit before claiming priority

These are not the exact project, but they are close enough that a referee
could ask about them. The safe claim is **not** "no chiral triangular random
walk exists." The safe claim is: we study an **active six-state
rotation-chiral run-and-tumble walker on the triangular lattice**, with
transport, first passage, and edge/localisation observables in the spirit
of refs 42-44.

* **`Gilbert_Sanders_2009_persistent_triangular_lorentz.pdf`** —
  persistent random-walk description of deterministic diffusion on a
  triangular Lorentz gas. This is a near miss because it uses direction
  memory on triangular geometry, but it is not the Mallikarjun-Pal-style
  active chiral RTW with tunable clockwise/counter-clockwise tumbling.
* **`Marris_Giuggioli_2024_persistent_antipersistent_first_passage.pdf`** —
  persistent/anti-persistent lattice walks and first-passage theory in
  bounded and unbounded domains. Useful for our FPT section and for checking
  what is already known about persistence-controlled search.
* **`Marris_Sarvaharman_Giuggioli_2023_hexagonal_honeycomb_domains.pdf`** —
  exact spatiotemporal dynamics of lattice walks in hexagonal/honeycomb
  domains. Boundary-geometry reference, especially if we discuss hexagonal
  or triangular domains.
* **`Oropesa_Castro_Lowen_Liarte_2026_adjustable_networks_triangular_rtp.pdf`** —
  active particles doing run-and-tumble motion along links of a triangular
  lattice, but with trail-mediated adjustable networks. Relevant because it
  is very recent and triangular-RTP-adjacent; not the same single-particle
  chiral RTW.
* **`Ouvry_Polychronakos_2019_exclusion_statistics_chiral_triangular_rw.pdf`** —
  "chiral random walk" on a triangular lattice in an exclusion-statistics /
  algebraic-area setting. This is mostly a title trap for us: cite only if
  needed to explain that our chirality is active-matter rotation chirality,
  not the same combinatorial/Hofstadter-style walk.
* **`Ouvry_Polychronakos_2020_lattice_walk_area_combinatorics.pdf`**
  — follow-up title trap. It explicitly discusses triangular lattice chiral
  walks, but in algebraic-area enumeration / Hofstadter combinatorics, not
  active RTW transport. Skim before making novelty claims around
  "triangular chiral random walk."
* **Larralde 1997, `Transport properties of a two-dimensional "chiral"
  persistent random walk`** — historical chiral persistent-walk transport
  paper. Not triangular/six-state, but important citation shield for the
  general statement "chirality modifies persistent-walk transport."
  **PDF not downloaded yet** (APS direct download returned 403).
* **`Batchelor_Henry_2002_triangular_absorbing_boundaries.pdf`** —
  exact triangular absorbing-boundary baseline. Not persistent or active,
  but useful before writing our first-passage section.
* **`Sarmiento_2025_first_passage_asymmetry_biased_RTP.pdf`** — modern FPT
  asymmetry paper for biased RTPs. Relevant to our FPT/search section; not
  triangular or C6.
* **`Gilbert_Sanders_2010_multistep_persistent_walks_lattices.pdf`** —
  multi-step persistent walk diffusion coefficients on lattices. Complements
  the 2009 Gilbert-Sanders triangular Lorentz-gas paper already present.
* **`Delplace_2020_topological_chiral_modes_random_scattering_networks.pdf`**
  — optional topology/edge-mode analogy. Useful later if our OBC spectrum
  discussion needs anomalous Floquet/scattering-network language.

### Direct templates we're extending

Read in this order:

1. **`Mallikarjun_Pal_2023_chiral_run_and_tumble_walker_ref44.pdf`**
   — Ref 44 of TCRW. 2D continuous space and continuous time, with four
   cardinal director states, transport + first-passage + optimal-search
   bias. This is the rate-mechanism and MFPT/search template, but it is
   **not** a square-lattice walker and not a bit-match baseline.
2. **`Hargus_Epstein_Mandadapu_2021_odd_diffusivity_chiral_ref42.pdf`**
   — Ref 42 of TCRW. The "odd diffusivity" antisymmetric tensor: the
   chirality signature in the diffusion tensor. We compute this on
   triangular.
3. **`Sevilla_2016_diffusion_active_chiral_ref43.pdf`** — Ref 43 of
   TCRW. Continuum chiral active particle, 3D. Used as a continuum-limit
   cross-check.

Context paper (already read for the bug-fix work, re-skim for the
edge-mode protocol):

4. **`Osat_Meyberg_Metson_Speck_2026_TCRW.pdf`** — Osat-Speck topological
   chiral random walker. Square-lattice 4-director discrete-time model
   with topological edge modes. Refs 42–44 are inside.

See also **`../LITERATURE_NOVELTY_AUDIT.md`** for the current strict
novelty map and the list of missing citation-shield papers.

## Phase 2 — multi-walker hard-core on triangular (later)

5. **`Whitelam_Klymko_Mandal_2017_lattice_active_MIPS.pdf`** — square
   lattice ABP with hard-core and asymmetric rotation. Phase 2 means
   reproducing this on triangular (with chirality, asking about MIPS).

## Phase 3 — equation of state / pressure with activity (DP2 / thesis)

6. **`Jaleel_Mandal_Thomas_Rajesh_2022_hardcore_triangular_freezing.pdf`**
   — exact-ish work on hard-core lattice gases on triangular (k-NN
   exclusion, freezing). Equilibrium baseline.
7. **`Hancock_Baskaran_2017_self_propelled_hard_spheres.pdf`** — derives
   a Smoluchowski equation for interacting ABPs. Theoretical bridge
   from microscopic to hydrodynamic.
8. **`Solon_etal_2015_pressure_not_state_function_NaturePhys.pdf`** —
   the foundational result that active matter pressure is not a state
   function in general. Motivates the whole Phase 3 angle.
9. **`Panagiotopoulos_2005_lattice_hard_sphere_models.pdf`** —
   lattice hard-sphere thermodynamics and equation-of-state reference from
   Kabir's email thread. Not a triangular active-walker paper, but useful
   background for the later pressure/equation-of-state direction.

## Group / background context

9. **`Jose_Mandal_Barma_Ramola_2022_JMVR_active_walks_1D_2D.pdf`** —
   home reference. The JMVR translation-chirality model in 1D and 2D
   square. The Confinement-2021 draft was the unpublished triangular
   extension; we just finished the bug-fix on its §IV.B.
10. **`Confinement_enhanced_clustering_in_dense_active_systems.pdf`** —
    the unpublished 2021 Mandal-Barma-Ramola draft we just bug-fixed.
11. **`Dolai_single_file_active_2020.pdf`** — single-file active matter.
    Background for the multi-walker direction; not Phase 1.
12. **`Lubensky_2015_Rep._Prog._Phys._78_073901.pdf`** — topological
    mechanics review. Background for the topology side of TCRW.

## Still missing / optional

- Still missing as a local PDF: Larralde 1997,
  `Transport properties of a two-dimensional "chiral" persistent random
  walk`. Keep it on the citation list; direct APS PDF download returned 403.
- Optional if we need more bounded-domain methodology: Giuggioli 2020 PRX,
  `Exact Spatiotemporal Dynamics of Confined Lattice Random Walks...`
  (direct APS/Bristol PDF download returned 403).
- If Kabir sends more hard-hexagon or equation-of-state notes, add them
  under Phase 3 rather than mixing them into the single-walker chiral-RTW
  reading queue.

## Reading priority right now (revised)

Updated reading order based on the priority/novelty audit:

1. **Wójcik-Kalz** (closest lattice prior art, discrete-time coin-step;
   must know exactly what they did)
2. **Mallikarjun-Pal** (rate/MFPT template, but continuum not lattice)
3. **Hargus-Epstein-Mandadapu** (odd-diffusivity tensor formalism)
4. **Gilbert-Sanders + Marris-Giuggioli papers** (quick novelty/FPT audit;
   do not deep-read yet)
5. **Ouvry-Polychronakos 2020 + Sarmiento 2025 + Batchelor-Henry 2002 +
   Gilbert-Sanders 2010** (citation-shield skim; PDFs now present)
6. **Larralde 1997** (citation-shield skim from DOI/abstract until PDF is
   available)
7. **Chiral active in obstacle lattices, Nat Comm 2024** (skim for the
   triangular-is-different evidence)
8. **Sevilla** (continuum-limit sanity check, lower priority)

Everything else (Phase 2, Phase 3, group context) stays in the folder as
background.

## Priority risk note

Wójcik-Kalz published February 2026. If they or any other group write
the triangular extension as their next paper, our priority claim
shrinks. Phase 1 should be a 1–2 month sprint, not a 6-month leisurely
exploration. Preprint within ~3-4 months even if not perfect.
