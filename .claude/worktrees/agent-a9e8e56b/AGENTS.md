This folder contains my TIFR Hyderabad DP1 report and related project context. My real goal is to use this as the starting point for DP2.

My DP1 work is in folder as report and it has all codes i have used and my report has refernces o paper that i extensively used for my DP1, especially jerky active particles, AOUP/ABP-type models, underdamped Langevin dynamics, and related stochastic systems.

Current DP2 motivation:
- I want to continue from jerky active-particle work, especially inspired by Hartmut Lowen’s jerky-particle papers
- I am interested in possible DP2 directions suggested after my report was shared with Stephy Jose(she is working with lowen and they have recent paper Jerky chiral active particles), including:I am sharing the email i got from stephy jose after My PI Kabir Ramola sent my DP1 report to her asking for future work-
  Thank you for sharing the report. It looks nice, and Prashant has done a detailed analysis of different works. Regarding the jerky particles, there are several directions to take forward:

Jerky harmonic oscillator in two dimensions:
One can derive the phase diagram, essentially extending Hartmut’s work to 2D. This is not particularly exciting but would be a relatively straightforward step.
Note: Jerky ABPs and Jerky AOUPs share the same statistics up to second order (same mean and MSD). One could therefore choose whichever model is easier to analyze.

Collective effects in jerky active particles:
This is something I am working on now. A natural question is how jerk influences MIPS. For example, this paper showed how inertia affects MIPS, with large mass destroying MIPS. Similarly, one could study the effect of jerk. I have tried some simulations with jerk, but they are tricky being third order and involving interactions, they are highly unstable. I am still working on this.

Higher-order derivatives and memory effects:
For instance, one could study a “snap-active particle” (fourth derivative) and examine how the MSD exponent depends on the highest-order derivative included in the dynamics.

Adding jerk to Tailleur’s model where exact calculations can be done:
Since it is a lattice model, even incorporating inertia can be different. It may be easy to think about equivalent continuum models with hard-core interactions.

Active harmonic solids:
One can study the effect of jerk on the collective excitations (entropons).

How Codex should help:
- read this report as the high-level summary of my DP1
-understand my complete coding style

- identify the parts of DP1 that naturally lead into DP2
- help me convert broad ideas into a realistic DP2 research plan
- compare possible DP2 topics by feasibility, originality, numerical difficulty, and theoretical value
- help me connect report results with the relevant Fortran (.f90) and gnuplot (.gnu) codes
- help me understand which existing code branches are most useful for extending into DP2
- help me draft proposals, problem statements, literature positioning, and concrete next steps

Important style instructions:
- keep the discussion research-oriented, not generic
- explain physics, equations, and modeling choices clearly and step by step
- if several directions are possible, help me rank them honestly
- distinguish between safe/feasible DP2 directions and ambitious/high-risk ones
- when discussing code, always relate it back to the scientific question

