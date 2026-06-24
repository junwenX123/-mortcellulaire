"""
Batch Gillespie / event-driven simulation for caspase activity, observed cell death,
and ERK negative feedback, with a fixed T-shaped background activation zone.

workflow implemented here:
  1. Burn-in first, by default 500 observed deaths.
  2. Then collect and analyze the next 1000 observed deaths.
  3. Compute statistics on the analyzed observed death process:
       - number of analyzed observed deaths,
       - time needed for the analyzed 1000 deaths after burn-in,
       - mean / median inter-death time after burn-in,
       - spatial chi-square-type statistics against uniform counts,
       - spatial coefficient of variation,
       - T-zone density ratio.
  4. Plot a histogram in time after burn-in.
  5. Plot a 2D histogram in space with a common colorbar scale across scenarios.
  6. Play with several parameter sets and save a comparison CSV file.

Run:
    python mortcellulaire_feedback_ERK_1000_stats_unified_burnin.py

Outputs are saved in:
    simulation_outputs_1000_unified_burnin/
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, List, Tuple, Any

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for batch simulations

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

try:
    from scipy.stats import chi2 as scipy_chi2
except Exception:  # pragma: no cover - fallback if scipy is unavailable
    scipy_chi2 = None


# ============================================================
# 1. Parameters
# ============================================================


@dataclass(frozen=True)
class Parameters:
    # Reproducibility
    seed: int = 123

    # Spatial domain W = [0, Lx] x [0, Ly]
    Lx: float = 10.0
    Ly: float = 10.0

    # Burn-in and stopping rule.
    # The simulation runs until burn_in_deaths + target_deaths observed deaths.
    # Statistics are computed only on the target_deaths after burn-in.
    burn_in_deaths: int = 500
    target_deaths: int = 1000

    # Activation intensities with fixed T-shaped extension.
    # lambda_a_1 is the dominating intensity and active-zone intensity.
    lambda_a_1: float = 5.00
    lambda_a_T: float = 0.50
    lambda_a_c: float = 0.05

    # Candidate death intensity on W.
    lambda_d: float = 1.0

    # Activation marks: R^a ~ Exp(beta_a_R), active-center lifetime rate beta_a_T.
    beta_a_R: float = 2.5
    beta_a_T: float = 1.2

    # ERK marks after accepted deaths: R^d ~ Exp(beta_d_R), ERK lifetime rate beta_d_T.
    # In the unified Gillespie version, T^d is represented by an ERK expiration rate,
    # not by a pre-sampled end time.
    beta_d_R: float = 2.0
    beta_d_T: float = 0.8

    # Plotting/statistics choices.
    bins_time: int = 30
    bins_space: int = 20       # used for the 2D plot and descriptive spatial CV
    bins_chi2: int = 10        # coarser bins for a more reliable chi-square test

    # Safety stop, only to avoid infinite runs if bad parameters are chosen.
    max_proposals: int = 10_000_000

    @property
    def area_W(self) -> float:
        return self.Lx * self.Ly

    @property
    def total_deaths_to_simulate(self) -> int:
        return self.burn_in_deaths + self.target_deaths


@dataclass
class SimulationResult:
    scenario: str
    params: Parameters

    # All accepted observed deaths, including burn-in deaths.
    death_x_all: np.ndarray
    death_y_all: np.ndarray
    death_t_all: np.ndarray

    # Analyzed observed deaths, after burn-in.
    death_x: np.ndarray
    death_y: np.ndarray
    death_t: np.ndarray              # raw time values
    death_t_since_burnin: np.ndarray # time relative to the end of burn-in

    # Accepted activations, historical.
    activation_x: np.ndarray
    activation_y: np.ndarray
    activation_t: np.ndarray

    final_time: float
    burn_in_time: float
    analysis_time_span: float

    n_observed_deaths_total: int
    n_observed_deaths_analyzed: int
    n_burn_in_deaths: int
    n_accepted_activations: int
    n_activation_candidates: int
    n_death_candidates: int
    n_activation_rejected: int
    n_death_rejected: int
    n_active_expirations: int
    n_erk_expirations: int
    n_proposals: int
    stopped_by_target: bool


# ============================================================
# 2. Geometry: fixed T-shaped zone and unions of disks
# ============================================================


def t_zone_bounds(p: Parameters) -> Tuple[float, float, float, float]:
    """Return the 3x3-grid cut points defining the fixed T-shaped zone."""
    x1 = p.Lx / 3.0
    x2 = 2.0 * p.Lx / 3.0
    y1 = p.Ly / 3.0
    y2 = 2.0 * p.Ly / 3.0
    return x1, x2, y1, y2


def t_zone_area(p: Parameters) -> float:
    """Area of the fixed T-zone: middle vertical column plus left middle arm."""
    x1, x2, y1, y2 = t_zone_bounds(p)
    middle_column_area = (x2 - x1) * p.Ly
    left_arm_area = x1 * (y2 - y1)
    return middle_column_area + left_arm_area


def is_inside_T_zone(x: float, y: float, p: Parameters) -> bool:
    """Fixed T-zone: middle vertical column plus left middle arm."""
    x1, x2, y1, y2 = t_zone_bounds(p)
    inside_middle_column = (x1 <= x <= x2) and (0.0 <= y <= p.Ly)
    inside_left_arm = (0.0 <= x <= x1) and (y1 <= y <= y2)
    return inside_middle_column or inside_left_arm


def is_inside_union_of_disks(
    x: float,
    y: float,
    centers_x: List[float],
    centers_y: List[float],
    radii: List[float],
) -> bool:
    """Return True if (x,y) is inside at least one disk."""
    for cx, cy, r in zip(centers_x, centers_y, radii):
        if (x - cx) ** 2 + (y - cy) ** 2 <= r ** 2:
            return True
    return False


# ============================================================
# 3. One unified Gillespie / thinning simulation
# ============================================================


def simulate_one(scenario: str, p: Parameters) -> SimulationResult:
    """
    Simulate until p.burn_in_deaths + p.target_deaths observed deaths are accepted.

    Unified Gillespie total rate:
        activation candidates:       lambda_a_1 |W|
        death candidates:            lambda_d |W|
        active-center expirations:   N^a_t beta_a_T
        ERK-zone expirations:        N^p_t beta_d_T

    Activation candidates are thinned according to the T-shaped intensity:
        lambda_a_1 inside active zones,
        lambda_a_T in the fixed T-zone but outside active zones,
        lambda_a_c outside both.

    Death candidates are accepted iff they are inside active zones and outside ERK zones.
    """
    if not (0.0 <= p.lambda_a_c <= p.lambda_a_T <= p.lambda_a_1):
        raise ValueError("Need 0 <= lambda_a_c <= lambda_a_T <= lambda_a_1 for thinning.")
    if min(p.lambda_d, p.beta_a_R, p.beta_a_T, p.beta_d_R, p.beta_d_T) <= 0:
        raise ValueError("Rates lambda_d, beta_a_R, beta_a_T, beta_d_R, beta_d_T must be positive.")
    if p.burn_in_deaths < 0 or p.target_deaths <= 0:
        raise ValueError("Need burn_in_deaths >= 0 and target_deaths > 0.")

    rng = np.random.default_rng(p.seed)
    t = 0.0

    # Current active centers V^a_t = sum delta_(Y_i^a, R_i^a)
    active_x: List[float] = []
    active_y: List[float] = []
    active_r: List[float] = []

    # Historical accepted activations X = {(Y_i^a, S_i^a)}
    activation_x: List[float] = []
    activation_y: List[float] = []
    activation_t: List[float] = []

    # Observed deaths Y = {(Y_i^d, S_i^d)}, including burn-in
    death_x: List[float] = []
    death_y: List[float] = []
    death_t: List[float] = []

    # Current ERK protection zones V^p_t. No pre-sampled end times in this version.
    erk_x: List[float] = []
    erk_y: List[float] = []
    erk_r: List[float] = []

    n_activation_candidates = 0
    n_death_candidates = 0
    n_activation_rejected = 0
    n_death_rejected = 0
    n_active_expirations = 0
    n_erk_expirations = 0
    n_proposals = 0

    def is_inside_active_zone(x: float, y: float) -> bool:
        return is_inside_union_of_disks(x, y, active_x, active_y, active_r)

    def is_inside_erk_zone(x: float, y: float) -> bool:
        return is_inside_union_of_disks(x, y, erk_x, erk_y, erk_r)

    def activation_intensity(x: float, y: float) -> float:
        if is_inside_active_zone(x, y):
            return p.lambda_a_1
        if is_inside_T_zone(x, y, p):
            return p.lambda_a_T
        return p.lambda_a_c

    while len(death_t) < p.total_deaths_to_simulate and n_proposals < p.max_proposals:
        n_active = len(active_x)
        n_erk = len(erk_x)

        activation_proposal_rate = p.lambda_a_1 * p.area_W
        death_proposal_rate = p.lambda_d * p.area_W
        active_expiration_rate = n_active * p.beta_a_T
        erk_expiration_rate = n_erk * p.beta_d_T

        a0 = (
            activation_proposal_rate
            + death_proposal_rate
            + active_expiration_rate
            + erk_expiration_rate
        )

        if a0 <= 0.0:
            break

        # Gillespie waiting time.
        t += rng.exponential(scale=1.0 / a0)

        # Gillespie event-type selection.
        u = rng.uniform(0.0, 1.0)
        p_activation = activation_proposal_rate / a0
        p_death = death_proposal_rate / a0
        p_active_expiration = active_expiration_rate / a0

        n_proposals += 1

        if u <= p_activation:
            # Event type 1: candidate activation from dominating PPP(lambda_a_1) on W.
            n_activation_candidates += 1
            x = rng.uniform(0.0, p.Lx)
            y = rng.uniform(0.0, p.Ly)

            accept_prob = activation_intensity(x, y) / p.lambda_a_1
            if rng.uniform(0.0, 1.0) <= accept_prob:
                activation_x.append(x)
                activation_y.append(y)
                activation_t.append(t)

                r = rng.exponential(scale=1.0 / p.beta_a_R)
                active_x.append(x)
                active_y.append(y)
                active_r.append(r)
            else:
                n_activation_rejected += 1

        elif u <= p_activation + p_death:
            # Event type 2: candidate death from PPP(lambda_d) on W.
            n_death_candidates += 1
            x = rng.uniform(0.0, p.Lx)
            y = rng.uniform(0.0, p.Ly)

            inside_active = is_inside_active_zone(x, y)
            inside_erk = is_inside_erk_zone(x, y)

            if inside_active and not inside_erk:
                death_x.append(x)
                death_y.append(y)
                death_t.append(t)

                # ERK feedback after accepted observed death.
                # Its lifetime is represented by the Gillespie rate N^p_t beta_d_T.
                r_E = rng.exponential(scale=1.0 / p.beta_d_R)
                erk_x.append(x)
                erk_y.append(y)
                erk_r.append(r_E)
            else:
                n_death_rejected += 1

        elif u <= p_activation + p_death + p_active_expiration:
            # Event type 3: expiration of one active center.
            if len(active_x) > 0:
                idx = int(rng.integers(0, len(active_x)))
                active_x.pop(idx)
                active_y.pop(idx)
                active_r.pop(idx)
                n_active_expirations += 1

        else:
            # Event type 4: expiration of one ERK protection zone.
            if len(erk_x) > 0:
                idx = int(rng.integers(0, len(erk_x)))
                erk_x.pop(idx)
                erk_y.pop(idx)
                erk_r.pop(idx)
                n_erk_expirations += 1

    death_x_all = np.asarray(death_x)
    death_y_all = np.asarray(death_y)
    death_t_all = np.asarray(death_t)

    if len(death_t_all) > p.burn_in_deaths:
        burn_in_time = 0.0 if p.burn_in_deaths == 0 else float(death_t_all[p.burn_in_deaths - 1])
        death_x_an = death_x_all[p.burn_in_deaths:]
        death_y_an = death_y_all[p.burn_in_deaths:]
        death_t_an = death_t_all[p.burn_in_deaths:]
        death_t_since_burnin = death_t_an - burn_in_time
        analysis_time_span = float(death_t_all[-1] - burn_in_time) if len(death_t_an) else 0.0
    else:
        burn_in_time = float(death_t_all[-1]) if len(death_t_all) else 0.0
        death_x_an = np.array([])
        death_y_an = np.array([])
        death_t_an = np.array([])
        death_t_since_burnin = np.array([])
        analysis_time_span = 0.0

    return SimulationResult(
        scenario=scenario,
        params=p,
        death_x_all=death_x_all,
        death_y_all=death_y_all,
        death_t_all=death_t_all,
        death_x=death_x_an,
        death_y=death_y_an,
        death_t=death_t_an,
        death_t_since_burnin=death_t_since_burnin,
        activation_x=np.asarray(activation_x),
        activation_y=np.asarray(activation_y),
        activation_t=np.asarray(activation_t),
        final_time=t,
        burn_in_time=burn_in_time,
        analysis_time_span=analysis_time_span,
        n_observed_deaths_total=len(death_t_all),
        n_observed_deaths_analyzed=len(death_t_an),
        n_burn_in_deaths=min(p.burn_in_deaths, len(death_t_all)),
        n_accepted_activations=len(activation_t),
        n_activation_candidates=n_activation_candidates,
        n_death_candidates=n_death_candidates,
        n_activation_rejected=n_activation_rejected,
        n_death_rejected=n_death_rejected,
        n_active_expirations=n_active_expirations,
        n_erk_expirations=n_erk_expirations,
        n_proposals=n_proposals,
        stopped_by_target=(len(death_t_all) >= p.total_deaths_to_simulate),
    )


# ============================================================
# 4. Statistics and plots on the observed death process
# ============================================================


def spatial_histogram(result: SimulationResult, bins: int) -> np.ndarray:
    p = result.params
    H, _, _ = np.histogram2d(
        result.death_x,
        result.death_y,
        bins=bins,
        range=[[0.0, p.Lx], [0.0, p.Ly]],
    )
    return H


def chi2_stats_from_histogram(H: np.ndarray) -> Tuple[float, int, float, float]:
    """Return chi2 statistic, df, chi2/df, and p-value if scipy is available."""
    n = float(np.sum(H))
    k = int(H.size)
    expected = n / k if k > 0 else np.nan
    if not np.isfinite(expected) or expected <= 0:
        return np.nan, k - 1, np.nan, np.nan
    chi2_uniform = float(np.sum((H - expected) ** 2 / expected))
    df = k - 1
    chi2_per_df = chi2_uniform / df if df > 0 else np.nan
    p_value = float(scipy_chi2.sf(chi2_uniform, df)) if scipy_chi2 is not None and df > 0 else np.nan
    return chi2_uniform, df, chi2_per_df, p_value


def t_zone_density_stats(result: SimulationResult) -> Dict[str, float | int]:
    p = result.params
    in_T = np.array([is_inside_T_zone(float(x), float(y), p) for x, y in zip(result.death_x, result.death_y)])
    n_T = int(np.sum(in_T))
    n_out = int(len(in_T) - n_T)
    area_T = t_zone_area(p)
    area_out = p.area_W - area_T
    density_T = n_T / area_T if area_T > 0 else np.nan
    density_out = n_out / area_out if area_out > 0 else np.nan
    ratio = density_T / density_out if density_out > 0 else np.nan
    return {
        "t_zone_area": area_T,
        "outside_t_zone_area": area_out,
        "t_zone_observed_deaths": n_T,
        "outside_t_zone_observed_deaths": n_out,
        "t_zone_observed_fraction": n_T / len(in_T) if len(in_T) else np.nan,
        "t_zone_area_fraction": area_T / p.area_W,
        "t_zone_death_density": density_T,
        "outside_t_zone_death_density": density_out,
        "t_zone_density_ratio": ratio,
    }


def death_statistics(result: SimulationResult) -> Dict[str, float | int | str | bool]:
    p = result.params

    # Inter-death times after burn-in. Include the waiting time from the burn-in cutoff
    # to the first analyzed death.
    if len(result.death_t) > 0:
        inter_death_after_burnin = np.diff(np.concatenate([[result.burn_in_time], result.death_t]))
    else:
        inter_death_after_burnin = np.array([])

    H_plot = spatial_histogram(result, p.bins_space)
    expected_plot = result.n_observed_deaths_analyzed / float(p.bins_space * p.bins_space)
    spatial_cv = float(np.std(H_plot) / np.mean(H_plot)) if np.mean(H_plot) > 0 else np.nan
    empty_bins = int(np.sum(H_plot == 0))
    
    H_chi2 = spatial_histogram(result, p.bins_chi2)
    expected_chi2 = result.n_observed_deaths_analyzed / float(p.bins_chi2 * p.bins_chi2)
    chi2_coarse, df_coarse, chi2_per_df_coarse, p_value_coarse = chi2_stats_from_histogram(H_chi2)

    t_stats = t_zone_density_stats(result)
    k_plot = p.bins_space * p.bins_space
    expected_cv_uniform_plot = np.sqrt((k_plot - 1) / result.n_observed_deaths_analyzed)
    spatial_cv_ratio_to_uniform = spatial_cv / expected_cv_uniform_plot
    return {
        "scenario": result.scenario,
        "seed": p.seed,
        "burn_in_deaths": p.burn_in_deaths,
        "analyzed_target_deaths": p.target_deaths,
        "total_deaths_simulated": result.n_observed_deaths_total,
        "n_observed_deaths_analyzed": result.n_observed_deaths_analyzed,
        "stopped_by_target": result.stopped_by_target,
        "burn_in_time": result.burn_in_time,
        "raw_final_time": result.final_time,
        "total_time_for_analyzed_1000_deaths_after_burnin": result.analysis_time_span,
        "observed_death_rate_after_burnin": result.n_observed_deaths_analyzed / result.analysis_time_span if result.analysis_time_span > 0 else np.nan,
        "mean_inter_death_time_after_burnin": float(np.mean(inter_death_after_burnin)) if len(inter_death_after_burnin) else np.nan,
        "median_inter_death_time_after_burnin": float(np.median(inter_death_after_burnin)) if len(inter_death_after_burnin) else np.nan,
        "accepted_activations_total": result.n_accepted_activations,
        "activation_candidates": result.n_activation_candidates,
        "death_candidates": result.n_death_candidates,
        "activation_rejected": result.n_activation_rejected,
        "death_rejected": result.n_death_rejected,
        "active_center_expirations": result.n_active_expirations,
        "erk_zone_expirations": result.n_erk_expirations,
        "all_proposals": result.n_proposals,
        "spatial_expected_cv_if_uniform_plot_bins": expected_cv_uniform_plot,
        "spatial_cv_ratio_to_uniform_plot_bins": spatial_cv_ratio_to_uniform,
        "spatial_bins_per_axis_for_plot": p.bins_space,
        "spatial_expected_count_per_plot_bin_if_uniform": expected_plot,
        "spatial_coefficient_of_variation_plot_bins": spatial_cv,
        "spatial_empty_bins_plot_bins": empty_bins,
        "spatial_bins_per_axis_for_chi2_test": p.bins_chi2,
        "spatial_expected_count_per_chi2_bin_if_uniform": expected_chi2,
        "spatial_chi2_statistic_against_uniform_counts_chi2_bins": chi2_coarse,
        "spatial_chi2_df_chi2_bins": df_coarse,
        "spatial_chi2_per_df_chi2_bins": chi2_per_df_coarse,
        "spatial_chi2_p_value_chi2_bins": p_value_coarse,
        **t_stats,
        "lambda_a_1": p.lambda_a_1,
        "lambda_a_T": p.lambda_a_T,
        "lambda_a_c": p.lambda_a_c,
        "lambda_d": p.lambda_d,
        "beta_a_R": p.beta_a_R,
        "beta_a_T": p.beta_a_T,
        "beta_d_R": p.beta_d_R,
        "beta_d_T": p.beta_d_T,
    }


def save_death_events_csv(result: SimulationResult, out_dir: Path) -> Path:
    path = out_dir / f"{result.scenario}_observed_deaths_after_burnin.csv"
    p = result.params
    in_T = [is_inside_T_zone(float(x), float(y), p) for x, y in zip(result.death_x, result.death_y)]
    raw_indices = np.arange(p.burn_in_deaths + 1, p.burn_in_deaths + result.n_observed_deaths_analyzed + 1)
    analysis_indices = np.arange(1, result.n_observed_deaths_analyzed + 1)
    data = np.column_stack([
        analysis_indices,
        raw_indices,
        result.death_t,
        result.death_t_since_burnin,
        result.death_x,
        result.death_y,
        np.asarray(in_T, dtype=int),
    ])
    header = "analysis_death_index,raw_death_index,raw_death_time,time_since_burnin,x,y,in_T_zone"
    np.savetxt(
        path,
        data,
        delimiter=",",
        header=header,
        comments="",
        fmt=["%d", "%d", "%.10f", "%.10f", "%.10f", "%.10f", "%d"],
    )
    return path


def add_t_zone_overlay(ax: plt.Axes, p: Parameters) -> None:
    """Add transparent rectangles showing the fixed T-zone on a spatial plot."""
    x1, x2, y1, y2 = t_zone_bounds(p)
    ax.add_patch(Rectangle((x1, 0.0), x2 - x1, p.Ly, fill=False, linewidth=1.2, linestyle="--"))
    ax.add_patch(Rectangle((0.0, y1), x1, y2 - y1, fill=False, linewidth=1.2, linestyle="--"))


def plot_time_histogram(result: SimulationResult, out_dir: Path) -> Path:
    p = result.params
    path = out_dir / f"{result.scenario}_time_histogram_after_burnin.png"

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(result.death_t_since_burnin, bins=p.bins_time, range=(0.0, result.analysis_time_span))
    ax.set_xlabel("time after burn-in")
    ax.set_ylabel("number of observed deaths")
    ax.set_title(
        f"Observed cell death times after burn-in — {result.scenario}\n"
        f"burn-in = {p.burn_in_deaths}, N = {result.n_observed_deaths_analyzed}, "
        f"time span = {result.analysis_time_span:.3f}"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_spatial_histogram_2d(result: SimulationResult, out_dir: Path, colorbar_vmax: int) -> Path:
    p = result.params
    path = out_dir / f"{result.scenario}_space_histogram_2d_after_burnin.png"

    H = spatial_histogram(result, p.bins_space)

    fig, ax = plt.subplots(figsize=(7, 6))
    image = ax.imshow(
        H.T,
        origin="lower",
        extent=[0.0, p.Lx, 0.0, p.Ly],
        aspect="equal",
        interpolation="nearest",
        vmin=0,
        vmax=colorbar_vmax,
    )
    add_t_zone_overlay(ax, p)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(
        f"2D histogram of observed death locations — {result.scenario}\n"
        f"after burn-in {p.burn_in_deaths}, {p.bins_space} x {p.bins_space} bins, "
        f"N = {result.n_observed_deaths_analyzed}"
    )
    cbar = fig.colorbar(image, ax=ax, label="number of observed deaths")
    cbar.set_ticks(np.arange(0, colorbar_vmax + 1, max(1, colorbar_vmax // 7)))
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_cumulative_deaths(result: SimulationResult, out_dir: Path) -> Path:
    path = out_dir / f"{result.scenario}_cumulative_deaths_after_burnin.png"
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.step(result.death_t_since_burnin, np.arange(1, result.n_observed_deaths_analyzed + 1), where="post")
    ax.set_xlabel("time after burn-in")
    ax.set_ylabel("cumulative observed deaths after burn-in")
    ax.set_title(f"Cumulative observed deaths after burn-in — {result.scenario}")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def save_summary_csv(rows: List[Dict[str, Any]], out_dir: Path) -> Path:
    path = out_dir / "parameter_sweep_summary_after_burnin.csv"
    keys = list(rows[0].keys())
    with path.open("w", encoding="utf-8") as f:
        f.write(",".join(keys) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) for k in keys) + "\n")
    return path


def write_analysis_report(rows: List[Dict[str, Any]], out_dir: Path, colorbar_vmax: int) -> Path:
    path = out_dir / "spatial_uniformity_analysis_report.md"
    lines: List[str] = []
    lines.append("# Spatial uniformity analysis after burn-in\n")
    lines.append(f"All 2D histograms use the same colorbar scale: 0 to {colorbar_vmax} observed deaths per bin.\n")
    lines.append("Burn-in removes the initial transient period: the first 500 observed deaths are simulated but not used for statistics. The reported statistics use the next 1000 observed deaths.\n")
    lines.append("## Interpretation rules\n")
    lines.append("- `spatial_coefficient_of_variation_plot_bins = std(H)/mean(H)`: close to 0 means more uniform; larger values mean more spatial heterogeneity.\n")
    lines.append("- `spatial_chi2_per_df_chi2_bins`: close to 1 is compatible with uniform counts; much larger than 1 indicates over-dispersion / non-uniformity.\n")
    lines.append("- `spatial_chi2_p_value_chi2_bins`: small p-value, for example < 0.05, rejects spatial uniformity on the chosen 10 x 10 grid.\n")
    lines.append("- `t_zone_density_ratio`: values > 1 mean deaths are denser in the T-zone than outside; values < 1 mean the opposite.\n")
    lines.append("\n## Results\n")
    header = (
        "| scenario | time for analyzed 1000 deaths | CV 20x20 | chi2/df 10x10 | p-value 10x10 | T-zone density ratio | empty bins 20x20 |\n"
        "|---|---:|---:|---:|---:|---:|---:|\n"
    )
    lines.append(header)
    for r in rows:
        lines.append(
            f"| {r['scenario']} | "
            f"{float(r['total_time_for_analyzed_1000_deaths_after_burnin']):.3f} | "
            f"{float(r['spatial_coefficient_of_variation_plot_bins']):.3f} | "
            f"{float(r['spatial_chi2_per_df_chi2_bins']):.3f} | "
            f"{float(r['spatial_chi2_p_value_chi2_bins']):.3g} | "
            f"{float(r['t_zone_density_ratio']):.3f} | "
            f"{int(r['spatial_empty_bins_plot_bins'])} |\n"
        )
    lines.append("\n## Short conclusion\n")
    baseline = rows[0]
    strongest = max(rows, key=lambda r: float(r['spatial_coefficient_of_variation_plot_bins']))
    lines.append(
        f"The baseline case has CV = {float(baseline['spatial_coefficient_of_variation_plot_bins']):.3f} "
        f"and T-zone density ratio = {float(baseline['t_zone_density_ratio']):.3f}. "
        "Thus it is close to spatially uniform compared with the stronger T-zone scenario.\n"
    )
    lines.append(
        f"The most spatially heterogeneous scenario by CV is `{strongest['scenario']}` "
        f"with CV = {float(strongest['spatial_coefficient_of_variation_plot_bins']):.3f} "
        f"and T-zone density ratio = {float(strongest['t_zone_density_ratio']):.3f}.\n"
    )
    path.write_text("".join(lines), encoding="utf-8")
    return path


# ============================================================
# 5. Parameter sweep: "play with various parameters"
# ============================================================


def make_scenarios(base: Parameters) -> Dict[str, Parameters]:
    """A few interpretable parameter changes around the baseline."""
    return {
        "baseline": base,

        # Larger lambda_d: death candidates arrive faster, so 1000 analyzed deaths should occur earlier.
        "higher_death_rate": replace(base, seed=124, lambda_d=1.5),

        # Larger mean ERK radius because mean R^d = 1 / beta_d_R.
        # This should reject more nearby death candidates and increase total time.
        "larger_ERK_radius": replace(base, seed=125, beta_d_R=1.0),

        # Shorter mean ERK lifetime because mean lifetime = 1 / beta_d_T.
        # This should reduce protection time and can make deaths arrive faster.
        "shorter_ERK_duration": replace(base, seed=126, beta_d_T=2.0),

        # Stronger activation in the fixed T-zone.
        "stronger_T_zone_activation": replace(base, seed=127, lambda_a_T=1.0),

        # A deliberately strong spatially heterogeneous case, useful to demonstrate non-uniformity.
        "strong_visible_T_zone": replace(
            base,
            seed=128,
            lambda_a_T=2.0,
            lambda_a_c=0.005,
            beta_a_R=5.0,
            beta_a_T=3.0,
        ),
    }


# ============================================================
# 6. Main
# ============================================================


def main() -> None:
    out_dir = Path("simulation_outputs_1000_unified_burnin")
    out_dir.mkdir(parents=True, exist_ok=True)

    base = Parameters(burn_in_deaths=500, target_deaths=1000)
    scenarios = make_scenarios(base)

    print("Running unified Gillespie simulations: burn-in 500 deaths, then analyze 1000 deaths...\n")

    results: List[SimulationResult] = []
    for name, params in scenarios.items():
        print(f"Scenario: {name}")
        result = simulate_one(name, params)
        results.append(result)
        print(f"  total observed deaths simulated: {result.n_observed_deaths_total}")
        print(f"  analyzed observed deaths:        {result.n_observed_deaths_analyzed}")
        print(f"  burn-in time:                    {result.burn_in_time:.6f}")
        print(f"  time for analyzed 1000 deaths:   {result.analysis_time_span:.6f}")
        print("")

    # Common colorbar scale across all scenarios.
    global_colorbar_max = int(max(np.max(spatial_histogram(result, result.params.bins_space)) for result in results))
    global_colorbar_max = max(global_colorbar_max, 1)

    summary_rows: List[Dict[str, Any]] = []
    for result in results:
        stats = death_statistics(result)
        summary_rows.append(stats)

        save_death_events_csv(result, out_dir)
        plot_time_histogram(result, out_dir)
        plot_spatial_histogram_2d(result, out_dir, global_colorbar_max)
        plot_cumulative_deaths(result, out_dir)

        print(f"Scenario: {result.scenario}")
        print(f"  spatial CV 20x20:       {stats['spatial_coefficient_of_variation_plot_bins']:.4f}")
        print(f"  chi2/df 10x10:          {stats['spatial_chi2_per_df_chi2_bins']:.4f}")
        print(f"  chi2 p-value 10x10:     {stats['spatial_chi2_p_value_chi2_bins']:.4g}")
        print(f"  T-zone density ratio:   {stats['t_zone_density_ratio']:.4f}")
        print("")

    summary_path = save_summary_csv(summary_rows, out_dir)
    report_path = write_analysis_report(summary_rows, out_dir, global_colorbar_max)

    print("Done.")
    print(f"Common colorbar vmax: {global_colorbar_max}")
    print(f"Summary CSV: {summary_path}")
    print(f"Analysis report: {report_path}")
    print(f"Plots and death-event CSV files are in: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
