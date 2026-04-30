"""
Shared utilities for Fig 3 cross-check scripts.

Provides:
  - load_summary(path)           : robust .txt loader (repairs Fortran's '1.79+308')
  - tau_relax(L, D_r)            : same formula the Fortran drivers use
  - reconstruct_T_use(...)       : recover per-row T_use to renormalize raw counts
  - rel_error(F, X, mask=...)    : safe relative error (skips |X|<eps)

The exact (authors-Python) physics is in tcrw_fig3_exact.py; this file
just deals with Fortran I/O and per-row arithmetic that's identical
across the six cross-check drivers.
"""
from __future__ import annotations
import numpy as np


# ---------------------------------------------------------------------------
def load_summary(path: str) -> np.ndarray:
    """
    Robust loader for tcrw_fig3*_summary.txt.

    Fortran's `huge(1.0_dp)` writes as '1.79769+308' (no 'E' before the
    sign).  np.loadtxt can't parse that, so we read line by line and
    patch any floating-point token missing the 'E'.

    Returns the data as a 2-D ndarray (rows × cols).
    """
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            fixed = []
            for p in parts:
                if "E" in p or "e" in p:
                    fixed.append(p)
                    continue
                body = p.lstrip("+-")
                sign = p[: len(p) - len(body)]
                ipos = max(body.rfind("+"), body.rfind("-"))
                if ipos > 0 and body[ipos - 1].isdigit():
                    body = body[:ipos] + "E" + body[ipos:]
                fixed.append(sign + body)
            rows.append([float(x) for x in fixed])
    return np.array(rows, dtype=float)


# ---------------------------------------------------------------------------
def tau_relax(L: int, D_r: np.ndarray | float) -> np.ndarray | float:
    """τ_relax = max(L^2, 1/D_r) / D_r  — same as Fortran drivers."""
    return np.maximum(L ** 2, 1.0 / np.asarray(D_r)) / np.asarray(D_r)


def reconstruct_T_use(L: int, D_r, K_meas: float = 100.0,
                      T_floor: float = 1e8) -> np.ndarray:
    """
    Reproduce per-row T_use exactly as the Fortran drivers compute it.
    Used to normalize raw integer current counts to per-step current.
    """
    return np.maximum(T_floor, K_meas * tau_relax(L, D_r))


# ---------------------------------------------------------------------------
def rel_error(F: np.ndarray, X: np.ndarray,
              eps: float = 1e-30,
              keep: np.ndarray | None = None) -> tuple[float, float]:
    """
    Max and mean relative error between Fortran F and exact X, ignoring:
      - rows with |X| < eps   (avoid 0/0)
      - non-finite F or X
      - rows masked by `keep` (False)
    """
    F = np.asarray(F)
    X = np.asarray(X)
    mask = np.isfinite(F) & np.isfinite(X) & (np.abs(X) > eps)
    if keep is not None:
        mask &= keep
    if not mask.any():
        return float("nan"), float("nan")
    rel = np.abs(F[mask] - X[mask]) / np.abs(X[mask])
    return float(np.max(rel)), float(np.mean(rel))


def angle_diff(theta_F: np.ndarray, theta_X: np.ndarray) -> np.ndarray:
    """Wrap angle difference to [-π, π]; useful for θ panels (3e, 3j)."""
    d = np.asarray(theta_F) - np.asarray(theta_X)
    return (d + np.pi) % (2 * np.pi) - np.pi
