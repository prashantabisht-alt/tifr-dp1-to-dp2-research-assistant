# Email draft — TCRW authors (Osat & Speck)

**To:** saeedosat13@gmail.com ; thomas.speck@itp4.uni-stuttgart.de
**CC:** (Prof. Kabir Ramola — after his OK)
**Subject:** TCRW (arXiv:2602.12020) — request for simulation code

---

Dear Dr. Osat and Prof. Speck,

I am a first-year PhD student at TIFR Hyderabad, working with Prof. Kabir Ramola on active-matter lattice models. Over the past few weeks I have been reproducing the figures of your paper, "Topological chiral random walker" (arXiv:2602.12020), as part of preparing a DP2 research direction on jerky chiral active particles. My code is public at

https://github.com/prashantabisht-alt/tifr-dp1-to-dp2-research-assistant

Would either of you be willing to share your simulation code — in particular the step-rule implementation and the 4×4 `P(k)` construction used for Figs 2 and 4? I would like to cross-check my conventions (order of operations for the chiral step, handling of blocked translations, and the Bloch-matrix blocks) against your reference before extending the model.

Thank you for the paper.

Best regards,
Prashant Bisht
DP1, TIFR Hyderabad
prashantabisht@gmail.com

---

## What is done (for you, not for them)

- **arXiv ID 2602.12020** — verified against arxiv.org. Title, authors, and v1 submission date (12 Feb 2026) all match. Safe to send.
- **Corner-vs-flat-wall finding** — removed from this short version. It is a real observation and worth raising, but not in the first email. Save it for the reply thread if they write back.
- **GitHub README** — I generated one at `README.md` in this folder. Push it to the repo root before sending the email; otherwise the URL leads to a file dump with no context. See the file for the exact steps.

## What only you can do

1. **Get Kabir's OK before CC-ing him.** One Slack message. Paste this draft.
2. **Push the README to GitHub** (see `README.md` in this folder; 3 git commands).
3. **Send.** Then: if no reply in 10 working days, one polite bump; after that, move on.
