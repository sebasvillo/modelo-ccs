"""Design-space search over NaOH concentration and L/G (legacy cells 8, 10 and 12)."""

from collections.abc import Sequence

import numpy as np

from .models import CaseInput, DesignPoint, OptimizationResult, Solvent
from .tea import evaluate_design_point

# LEGACY: multi-objective weights of choose_best_design (cell 8); lower score is better.
LEGACY_WEIGHTS = {
    "LCOC_usd_t": 0.35,
    "capex_total_usd": 0.25,
    "opex_total_usd_y": 0.20,
    "E_total_kWh_t": 0.10,
    "indirect_kgCO2e_t": 0.10,
}


def evaluate_grid(
    case: CaseInput, NaOH_grid_M: Sequence[float], LG_grid_vol: Sequence[float]
) -> tuple[DesignPoint, ...]:
    """Evaluate every (NaOH, L/G) pair, NaOH in the outer loop."""
    return tuple(
        evaluate_design_point(
            case.gas, Solvent(NaOH_M=float(C)), case.absorber, case.cell, case.tea, float(LG)
        )
        for C in NaOH_grid_M
        for LG in LG_grid_vol
    )


def _minmax(values: np.ndarray) -> np.ndarray:
    vmin, vmax = np.nanmin(values), np.nanmax(values)
    if np.isclose(vmax, vmin):
        return np.zeros_like(values)
    return (values - vmin) / (vmax - vmin)


def weighted_scores(objectives: dict[str, np.ndarray], weights: dict[str, float]) -> np.ndarray:
    """Sum of weight × min–max normalised objective (lower is better), in weights order.

    Equations: optimize.weighted_score.
    """
    score = None
    for key, w in weights.items():
        term = w * _minmax(np.asarray(objectives[key], dtype=float))
        score = term if score is None else score + term
    return score


def choose_best(
    points: Sequence[DesignPoint], weights: dict[str, float] | None = None
) -> OptimizationResult:
    """Weighted sum of min–max normalised objectives over the feasible points; lowest wins.

    Only points that reach the capture target are ranked (audit §7, errata E-005). If none
    does, the point with the highest capture is returned with feasible=False and a warning.
    Infeasible points get a NaN score.

    Equations: optimize.weighted_score.
    """
    weights = LEGACY_WEIGHTS if weights is None else weights
    feasible = [i for i, p in enumerate(points) if p.absorber.reached_target]
    scores = np.full(len(points), np.nan)
    if not feasible:
        i_best = max(range(len(points)), key=lambda i: points[i].absorber.capture_achieved)
        top = points[i_best]
        warning = (
            f"no design reaches the capture target; best capture "
            f"{top.absorber.capture_achieved:.1%} at NaOH {top.NaOH_M:.3g} M, L/G {top.LG_vol:.4g}"
        )
        return OptimizationResult(
            best=points[i_best],
            best_score=float("nan"),
            points=tuple(points),
            scores=tuple(float(v) for v in scores),
            feasible=False,
            n_feasible=0,
            warnings=(warning,),
        )
    ranked = [points[i] for i in feasible]
    scores[feasible] = weighted_scores(
        {k: [getattr(p, k) for p in ranked] for k in weights}, weights
    )
    i_best = feasible[int(np.argmin(scores[feasible]))]
    return OptimizationResult(
        best=points[i_best],
        best_score=float(scores[i_best]),
        points=tuple(points),
        scores=tuple(float(v) for v in scores),
        feasible=True,
        n_feasible=len(feasible),
    )


def coarse_grids() -> tuple[np.ndarray, np.ndarray]:
    """Legacy cell 10 grids: 14 NaOH values × 26 L/G values."""
    return np.linspace(0.25, 4.00, 14), np.linspace(0.004, 0.070, 26)


def fine_grids(NaOH_center_M: float, LG_center_vol: float) -> tuple[np.ndarray, np.ndarray]:
    """Legacy cell 12 grids around the coarse optimum: 16 × 20 values."""
    return (
        np.linspace(max(0.10, NaOH_center_M - 0.50), NaOH_center_M + 0.50, 16),
        np.linspace(max(0.001, LG_center_vol - 0.010), LG_center_vol + 0.010, 20),
    )


def optimize_two_stage(case: CaseInput) -> tuple[OptimizationResult, OptimizationResult]:
    """Coarse 2D sweep, then local refinement around its best point (cells 10 and 12)."""
    coarse = choose_best(evaluate_grid(case, *coarse_grids()))
    fine = choose_best(evaluate_grid(case, *fine_grids(coarse.best.NaOH_M, coarse.best.LG_vol)))
    return coarse, fine
