"""
Batch rejection-sampling / thinning simulation for the homogeneous-region HMPP model
with active caspase feedback and ERK negative feedback.

This script keeps the original rejection-sampling block structure, but adds the
same analysis workflow as the optimized 1000-death scripts:
  1. Burn-in first, by default 500 observed deaths.
  2. Then collect exactly the next 200 observed deaths and stop.
  3. Compute statistics on the analyzed observed death process:
       - number of analyzed observed events,
       - time needed for the analyzed 200 deaths after burn-in,
       - mean / median inter-death time after burn-in,
       - death candidate rejection rates,
       - spatial coefficient of variation on the plotting grid,
       - chi-square goodness-of-fit test against spatial uniformity on a coarser grid.
  4. Plot a histogram in time after burn-in.
  5. Plot a 2D histogram in space with one common colorbar scale across scenarios.
  6. Save event CSV files, a parameter-sweep summary CSV, and a Markdown analysis report.

Important stopping convention:
  With burn_in_deaths=500 and target_deaths=200, the simulator accepts 700 deaths in total,
  discards the first 500 deaths, analyzes deaths 501--700, and then stops.

Run:
    python homogeneous_region_200_deaths_rejection_burnin_optimized.py

Outputs are saved in:
    simulation_outputs_homogeneous_200_rejection_burnin_optimized/
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Dict, List, Tuple
import csv
import json

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for batch simulations

import matplotlib.pyplot as plt
import numpy as np

try:
    from scipy.stats import chi2 as scipy_chi2
except Exception:  # pragma: no cover - fallback if scipy is unavailable
    scipy_chi2 = None


# ============================================================
# 1. Parameters
# ============================================================


OUTPUT_DIR = Path("simulation_outputs_homogeneous_200_rejection_burnin_optimized")


@dataclass(frozen=True)
class Params:
    name: str = "baseline"
    seed: int = 123

    # Spatial window W = [0, Lx] x [0, Ly]
    Lx: float = 10.0
    Ly: float = 10.0

    # Burn-in and stopping rule.
    # The simulation stops when burn_in_deaths + target_deaths observed deaths
    # have been accepted. Statistics use only the target_deaths after burn-in.
    # For the exercise: burn-in = 500, analyzed deaths = 200.
    burn_in_deaths: int = 500
    target_deaths: int = 200

    # Block simulation parameters.
    # The rejection-sampling structure is kept: candidates are generated block by block.
    block_duration: float = 20.0
    max_blocks: int = 500

    # Activation process V^a_t.
    # lambda_a_1 is the dominating candidate intensity and active-zone intensity;
    # lambda_a_c is the homogeneous background intensity outside active zones.
    lambda_a_1: float = 0.60
    lambda_a_c: float = 0.05

    # Candidate death process Phi^d on W.
    lambda_d: float = 0.25

    # Activation marks: R^a ~ Exp(beta_a_R), T^a ~ Exp(beta_a_T)
    beta_a_R: float = 1.4
    beta_a_T: float = 0.8

    # Protection marks after accepted deaths: R^d ~ Exp(beta_d_R), T^d ~ Exp(beta_d_T)
    beta_d_R: float = 1.2
    beta_d_T: float = 0.6

    # Plotting/statistics choices.
    time_bins: int = 20
    space_bins: int = 10       # used for 2D plots and descriptive spatial CV
    chi2_bins: int = 5         # 200/(5*5)=8 >= 5, better for chi-square approximation

    @property
    def area_W(self) -> float:
        return self.Lx * self.Ly

    @property
    def total_deaths_to_simulate(self) -> int:
        return self.burn_in_deaths + self.target_deaths


@dataclass
class SimulationResult:
    name: str
    params: Params

    # All accepted observed deaths, including burn-in deaths.
    death_x_all: np.ndarray
    death_y_all: np.ndarray
    death_t_all: np.ndarray

    # Analyzed observed deaths, after burn-in.
    death_x: np.ndarray
    death_y: np.ndarray
    death_t: np.ndarray              # raw time values from t=0
    death_t_since_burnin: np.ndarray # time relative to burn-in end time

    # Historical accepted activations.
    activation_x: np.ndarray
    activation_y: np.ndarray
    activation_t: np.ndarray

    burn_in_time: float
    final_time: float
    analysis_time_span: float

    total_deaths_simulated: int
    analyzed_observed_deaths: int
    n_burn_in_deaths: int
    stopped_by_target: bool

    activation_candidates: int
    activation_accepted_total: int
    death_candidates_processed: int
    death_rejected_outside_active: int
    death_rejected_by_erk: int
    death_rejected_total: int
    erk_zones_created_total: int


# ============================================================
# 2. Elementary geometry tests
# ============================================================


def validate_params(p: Params) -> None:
    if not (0.0 <= p.lambda_a_c <= p.lambda_a_1):
        raise ValueError("Need 0 <= lambda_a_c <= lambda_a_1 for rejection sampling.")
    for key in ["beta_a_R", "beta_a_T", "lambda_d", "beta_d_R", "beta_d_T"]:
        if getattr(p, key) <= 0:
            raise ValueError(f"Need {key} > 0.")
    if p.Lx <= 0 or p.Ly <= 0 or p.block_duration <= 0:
        raise ValueError("Need positive domain size and block duration.")
    if p.burn_in_deaths < 0 or p.target_deaths <= 0:
        raise ValueError("Need burn_in_deaths >= 0 and target_deaths > 0.")


def point_in_active_zone_lists(
    x: float,
    y: float,
    s: float,
    X_acc: List[float],
    Y_acc: List[float],
    S_acc: List[float],
    R_acc: List[float],
    Tau_acc: List[float],
) -> bool:
    """Return True if (x, y) belongs to A(V^a_{s-})."""
    for xk, yk, sk, rk, tauk in zip(X_acc, Y_acc, S_acc, R_acc, Tau_acc):
        # State just before s: S_k < s <= S_k + T_k
        if sk < s <= sk + tauk:
            if (x - xk) ** 2 + (y - yk) ** 2 <= rk ** 2:
                return True
    return False


def point_in_active_zone_arrays(
    x: float,
    y: float,
    s: float,
    X_acc: np.ndarray,
    Y_acc: np.ndarray,
    S_acc: np.ndarray,
    R_acc: np.ndarray,
    Tau_acc: np.ndarray,
) -> bool:
    """Return True if (x, y) belongs to A(V^a_{s-}), vectorized over accepted activations."""
    if len(X_acc) == 0:
        return False
    active_before_s = (S_acc < s) & (s <= S_acc + Tau_acc)
    if not np.any(active_before_s):
        return False
    dx = x - X_acc[active_before_s]
    dy = y - Y_acc[active_before_s]
    rr = R_acc[active_before_s]
    return bool(np.any(dx * dx + dy * dy <= rr * rr))


def point_in_protected_zone(
    x: float,
    y: float,
    s: float,
    X_death: List[float],
    Y_death: List[float],
    S_death: List[float],
    R_death: List[float],
    Tau_death: List[float],
) -> bool:
    """Return True if (x, y) belongs to A(V^p_{s-})."""
    for xj, yj, sj, rj, tauj in zip(X_death, Y_death, S_death, R_death, Tau_death):
        # Protection just before s: S_j^d < s <= S_j^d + T_j^d
        if sj < s <= sj + tauj:
            if (x - xj) ** 2 + (y - yj) ** 2 <= rj ** 2:
                return True
    return False


# ============================================================
# 3. Rejection sampling / thinning simulation until burn-in + 200 deaths
# ============================================================


def simulate_one(p: Params) -> SimulationResult:
    """
    Simulate the homogeneous-region model until p.burn_in_deaths + p.target_deaths
    observed deaths are accepted.

    The original rejection-sampling block structure is kept:
      - candidate activations are proposed from a dominating process on W x block;
      - activation candidates are thinned according to active-zone feedback;
      - candidate deaths are proposed on W x block;
      - death candidates are accepted iff they are in active zones and outside ERK zones.
    """
    validate_params(p)
    rng = np.random.default_rng(p.seed)
    area_W = p.area_W

    # Accepted activation centers V^a.
    X_acc: List[float] = []
    Y_acc: List[float] = []
    S_acc: List[float] = []
    R_acc: List[float] = []
    Tau_acc: List[float] = []

    # Accepted observed deaths V^d and their ERK marks V^p.
    X_death: List[float] = []
    Y_death: List[float] = []
    S_death: List[float] = []
    R_death: List[float] = []
    Tau_death: List[float] = []

    activation_candidates_total = 0
    death_candidates_processed = 0
    death_rejected_outside_active = 0
    death_rejected_by_erk = 0

    t0 = 0.0

    for _block in range(p.max_blocks):
        t1 = t0 + p.block_duration

        # ----------------------------------------------------
        # A) Candidate activations in W x [t0, t1]
        # ----------------------------------------------------
        mean_activation_candidates = p.lambda_a_1 * area_W * p.block_duration
        N_act = int(rng.poisson(mean_activation_candidates))
        activation_candidates_total += N_act

        X_star = rng.uniform(0.0, p.Lx, size=N_act)
        Y_star = rng.uniform(0.0, p.Ly, size=N_act)
        S_star = rng.uniform(t0, t1, size=N_act)
        U_star = rng.uniform(0.0, 1.0, size=N_act)
        R_star = rng.exponential(scale=1.0 / p.beta_a_R, size=N_act)
        Tau_star = rng.exponential(scale=1.0 / p.beta_a_T, size=N_act)

        order = np.argsort(S_star)
        for idx in order:
            x = float(X_star[idx])
            y = float(Y_star[idx])
            s = float(S_star[idx])
            u = float(U_star[idx])
            r = float(R_star[idx])
            tau = float(Tau_star[idx])

            inside_active = point_in_active_zone_lists(
                x, y, s, X_acc, Y_acc, S_acc, R_acc, Tau_acc
            )

            # Rejection / thinning rule:
            # p^a_v(x) = 1 if x in A(v), lambda_a_c/lambda_a_1 otherwise.
            if inside_active or (u <= p.lambda_a_c / p.lambda_a_1):
                X_acc.append(x)
                Y_acc.append(y)
                S_acc.append(s)
                R_acc.append(r)
                Tau_acc.append(tau)

        # Convert activations to arrays once per block for fast death tests.
        X_acc_arr = np.asarray(X_acc)
        Y_acc_arr = np.asarray(Y_acc)
        S_acc_arr = np.asarray(S_acc)
        R_acc_arr = np.asarray(R_acc)
        Tau_acc_arr = np.asarray(Tau_acc)

        # ----------------------------------------------------
        # B) Candidate deaths in W x [t0, t1]
        # ----------------------------------------------------
        mean_death_candidates = p.lambda_d * area_W * p.block_duration
        N_death_star = int(rng.poisson(mean_death_candidates))

        X_death_star = rng.uniform(0.0, p.Lx, size=N_death_star)
        Y_death_star = rng.uniform(0.0, p.Ly, size=N_death_star)
        S_death_star = rng.uniform(t0, t1, size=N_death_star)
        R_death_star = rng.exponential(scale=1.0 / p.beta_d_R, size=N_death_star)
        Tau_death_star = rng.exponential(scale=1.0 / p.beta_d_T, size=N_death_star)

        death_order = np.argsort(S_death_star)
        for idx in death_order:
            if len(S_death) >= p.total_deaths_to_simulate:
                break

            death_candidates_processed += 1

            xd = float(X_death_star[idx])
            yd = float(Y_death_star[idx])
            sd = float(S_death_star[idx])
            rd = float(R_death_star[idx])
            taud = float(Tau_death_star[idx])

            inside_active = point_in_active_zone_arrays(
                xd, yd, sd, X_acc_arr, Y_acc_arr, S_acc_arr, R_acc_arr, Tau_acc_arr
            )
            inside_protected = point_in_protected_zone(
                xd, yd, sd, X_death, Y_death, S_death, R_death, Tau_death
            )

            if inside_active and not inside_protected:
                X_death.append(xd)
                Y_death.append(yd)
                S_death.append(sd)
                R_death.append(rd)
                Tau_death.append(taud)
            else:
                if inside_protected:
                    death_rejected_by_erk += 1
                else:
                    death_rejected_outside_active += 1

        if len(S_death) >= p.total_deaths_to_simulate:
            break

        t0 = t1
    else:
        raise RuntimeError(
            f"Only {len(S_death)} deaths after {p.max_blocks} blocks. "
            "Increase block_duration, max_blocks, lambda_d, or active-zone parameters."
        )

    # Keep exactly burn-in + target deaths.
    total = p.total_deaths_to_simulate
    X_death_all = np.asarray(X_death[:total])
    Y_death_all = np.asarray(Y_death[:total])
    S_death_all = np.asarray(S_death[:total])

    if len(S_death_all) > p.burn_in_deaths:
        burn_in_time = 0.0 if p.burn_in_deaths == 0 else float(S_death_all[p.burn_in_deaths - 1])
        X_death_an = X_death_all[p.burn_in_deaths:]
        Y_death_an = Y_death_all[p.burn_in_deaths:]
        S_death_an = S_death_all[p.burn_in_deaths:]
        S_death_since_burnin = S_death_an - burn_in_time
        analysis_time_span = float(S_death_all[-1] - burn_in_time) if len(S_death_an) else 0.0
    else:
        burn_in_time = float(S_death_all[-1]) if len(S_death_all) else 0.0
        X_death_an = np.array([])
        Y_death_an = np.array([])
        S_death_an = np.array([])
        S_death_since_burnin = np.array([])
        analysis_time_span = 0.0

    final_time = float(S_death_all[-1]) if len(S_death_all) else t0

    return SimulationResult(
        name=p.name,
        params=p,
        death_x_all=X_death_all,
        death_y_all=Y_death_all,
        death_t_all=S_death_all,
        death_x=X_death_an,
        death_y=Y_death_an,
        death_t=S_death_an,
        death_t_since_burnin=S_death_since_burnin,
        activation_x=np.asarray(X_acc),
        activation_y=np.asarray(Y_acc),
        activation_t=np.asarray(S_acc),
        burn_in_time=burn_in_time,
        final_time=final_time,
        analysis_time_span=analysis_time_span,
        total_deaths_simulated=len(S_death_all),
        analyzed_observed_deaths=len(S_death_an),
        n_burn_in_deaths=min(p.burn_in_deaths, len(S_death_all)),
        stopped_by_target=(len(S_death_all) >= total),
        activation_candidates=activation_candidates_total,
        activation_accepted_total=len(X_acc),
        death_candidates_processed=death_candidates_processed,
        death_rejected_outside_active=death_rejected_outside_active,
        death_rejected_by_erk=death_rejected_by_erk,
        death_rejected_total=death_rejected_outside_active + death_rejected_by_erk,
        erk_zones_created_total=len(S_death_all),
    )


# ============================================================
# 4. Statistics and plots for the observed death process
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


def death_statistics(result: SimulationResult) -> Dict[str, Any]:
    p = result.params

    if len(result.death_t) > 0:
        inter_death_after_burnin = np.diff(np.concatenate([[result.burn_in_time], result.death_t]))
    else:
        inter_death_after_burnin = np.array([])

    H_plot = spatial_histogram(result, p.space_bins)
    expected_plot = result.analyzed_observed_deaths / float(p.space_bins * p.space_bins)
    spatial_cv = float(np.std(H_plot) / np.mean(H_plot)) if np.mean(H_plot) > 0 else np.nan
    empty_bins = int(np.sum(H_plot == 0))

    k_plot = p.space_bins * p.space_bins
    expected_cv_uniform_plot = (
        float(np.sqrt((k_plot - 1) / result.analyzed_observed_deaths))
        if result.analyzed_observed_deaths > 0
        else np.nan
    )
    spatial_cv_ratio_to_uniform = (
        spatial_cv / expected_cv_uniform_plot
        if expected_cv_uniform_plot > 0
        else np.nan
    )

    H_chi2 = spatial_histogram(result, p.chi2_bins)
    expected_chi2 = result.analyzed_observed_deaths / float(p.chi2_bins * p.chi2_bins)
    chi2_coarse, df_coarse, chi2_per_df_coarse, p_value_coarse = chi2_stats_from_histogram(H_chi2)

    return {
        "scenario": result.name,
        "seed": p.seed,
        "burn_in_deaths": p.burn_in_deaths,
        "analyzed_target_deaths": p.target_deaths,
        "total_deaths_simulated": result.total_deaths_simulated,
        "analyzed_observed_deaths": result.analyzed_observed_deaths,
        "stopped_by_target": result.stopped_by_target,
        "burn_in_time": result.burn_in_time,
        "raw_final_time": result.final_time,
        "total_time_for_analyzed_200_deaths_after_burnin": result.analysis_time_span,
        "observed_death_rate_after_burnin": result.analyzed_observed_deaths / result.analysis_time_span if result.analysis_time_span > 0 else np.nan,
        "mean_inter_death_time_after_burnin": float(np.mean(inter_death_after_burnin)) if len(inter_death_after_burnin) else np.nan,
        "median_inter_death_time_after_burnin": float(np.median(inter_death_after_burnin)) if len(inter_death_after_burnin) else np.nan,
        "activation_candidates": result.activation_candidates,
        "activation_accepted_total": result.activation_accepted_total,
        "death_candidates_processed": result.death_candidates_processed,
        "death_rejected_outside_active": result.death_rejected_outside_active,
        "death_rejected_by_erk": result.death_rejected_by_erk,
        "death_rejected_total": result.death_rejected_total,
        "death_rejection_rate": result.death_rejected_total / result.death_candidates_processed if result.death_candidates_processed > 0 else np.nan,
        "death_acceptance_rate_total": result.total_deaths_simulated / result.death_candidates_processed if result.death_candidates_processed > 0 else np.nan,
        "erk_zones_created_total": result.erk_zones_created_total,
        "spatial_bins_per_axis_for_plot": p.space_bins,
        "spatial_expected_count_per_plot_bin_if_uniform": expected_plot,
        "spatial_coefficient_of_variation_plot_bins": spatial_cv,
        "spatial_expected_cv_if_uniform_plot_bins": expected_cv_uniform_plot,
        "spatial_cv_ratio_to_uniform_plot_bins": spatial_cv_ratio_to_uniform,
        "spatial_empty_bins_plot_bins": empty_bins,
        "spatial_bins_per_axis_for_chi2_test": p.chi2_bins,
        "spatial_expected_count_per_chi2_bin_if_uniform": expected_chi2,
        "spatial_chi2_statistic_against_uniform_counts_chi2_bins": chi2_coarse,
        "spatial_chi2_df_chi2_bins": df_coarse,
        "spatial_chi2_per_df_chi2_bins": chi2_per_df_coarse,
        "spatial_chi2_p_value_chi2_bins": p_value_coarse,
        "lambda_a_1": p.lambda_a_1,
        "lambda_a_c": p.lambda_a_c,
        "lambda_d": p.lambda_d,
        "beta_a_R": p.beta_a_R,
        "beta_a_T": p.beta_a_T,
        "beta_d_R": p.beta_d_R,
        "beta_d_T": p.beta_d_T,
    }


def save_death_events_csv(result: SimulationResult, out_dir: Path) -> Path:
    scenario_dir = out_dir / result.name
    scenario_dir.mkdir(parents=True, exist_ok=True)
    path = scenario_dir / f"{result.name}_observed_deaths_after_burnin.csv"
    p = result.params

    raw_indices = np.arange(p.burn_in_deaths + 1, p.burn_in_deaths + result.analyzed_observed_deaths + 1)
    analysis_indices = np.arange(1, result.analyzed_observed_deaths + 1)
    data = np.column_stack([
        analysis_indices,
        raw_indices,
        result.death_t,
        result.death_t_since_burnin,
        result.death_x,
        result.death_y,
    ])
    header = "analysis_death_index,raw_death_index,raw_death_time,time_since_burnin,x,y"
    np.savetxt(
        path,
        data,
        delimiter=",",
        header=header,
        comments="",
        fmt=["%d", "%d", "%.10f", "%.10f", "%.10f", "%.10f"],
    )
    return path


def save_params_and_summary(result: SimulationResult, stats: Dict[str, Any], out_dir: Path) -> None:
    scenario_dir = out_dir / result.name
    scenario_dir.mkdir(parents=True, exist_ok=True)
    with (scenario_dir / f"{result.name}_parameters.json").open("w", encoding="utf-8") as f:
        json.dump(asdict(result.params), f, indent=2)
    with (scenario_dir / f"{result.name}_summary.json").open("w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)


def plot_time_histogram(result: SimulationResult, out_dir: Path) -> Path:
    p = result.params
    scenario_dir = out_dir / result.name
    scenario_dir.mkdir(parents=True, exist_ok=True)
    path = scenario_dir / f"{result.name}_time_histogram_after_burnin.png"

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(result.death_t_since_burnin, bins=p.time_bins, range=(0.0, result.analysis_time_span))
    ax.set_xlabel("time after burn-in")
    ax.set_ylabel("number of observed deaths")
    ax.set_title(
        f"Observed cell death times after burn-in — {result.name}\n"
        f"burn-in = {p.burn_in_deaths}, N = {result.analyzed_observed_deaths}, "
        f"time span = {result.analysis_time_span:.3f}"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def plot_spatial_histogram_2d(result: SimulationResult, out_dir: Path, colorbar_vmax: int) -> Path:
    p = result.params
    scenario_dir = out_dir / result.name
    scenario_dir.mkdir(parents=True, exist_ok=True)
    path = scenario_dir / f"{result.name}_space_2d_histogram_after_burnin.png"

    H = spatial_histogram(result, p.space_bins)

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
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(
        f"2D histogram of observed death locations — {result.name}\n"
        f"after burn-in {p.burn_in_deaths}, {p.space_bins} x {p.space_bins} bins, "
        f"N = {result.analyzed_observed_deaths}"
    )
    cbar = fig.colorbar(image, ax=ax, label="number of observed deaths")
    cbar.set_ticks(np.arange(0, colorbar_vmax + 1, max(1, colorbar_vmax // 7)))
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def plot_death_locations(result: SimulationResult, out_dir: Path) -> Path:
    p = result.params
    scenario_dir = out_dir / result.name
    scenario_dir.mkdir(parents=True, exist_ok=True)
    path = scenario_dir / f"{result.name}_death_locations_after_burnin.png"

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(result.death_x, result.death_y, s=18)
    ax.set_xlim(0.0, p.Lx)
    ax.set_ylim(0.0, p.Ly)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Observed death locations after burn-in, N = {result.analyzed_observed_deaths} ({result.name})")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def plot_cumulative_deaths(result: SimulationResult, out_dir: Path) -> Path:
    scenario_dir = out_dir / result.name
    scenario_dir.mkdir(parents=True, exist_ok=True)
    path = scenario_dir / f"{result.name}_cumulative_deaths_after_burnin.png"

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.step(result.death_t_since_burnin, np.arange(1, result.analyzed_observed_deaths + 1), where="post")
    ax.set_xlabel("time after burn-in")
    ax.set_ylabel("cumulative observed deaths after burn-in")
    ax.set_title(f"Cumulative observed deaths after burn-in — {result.name}")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def save_summary_csv(rows: List[Dict[str, Any]], out_dir: Path) -> Path:
    path = out_dir / "summary_homogeneous_200_after_burnin.csv"
    keys = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    return path


def plot_comparison_total_time(rows: List[Dict[str, Any]], out_dir: Path) -> Path:
    path = out_dir / "comparison_total_time_after_burnin_T200.png"
    names = [str(r["scenario"]) for r in rows]
    times = [float(r["total_time_for_analyzed_200_deaths_after_burnin"]) for r in rows]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(names, times)
    ax.set_ylabel("time for analyzed 200 observed deaths after burn-in")
    ax.set_title("Effect of parameters on the stopping time after burn-in")
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def write_analysis_report(rows: List[Dict[str, Any]], out_dir: Path, colorbar_vmax: int) -> Path:
    path = out_dir / "spatial_uniformity_analysis_report_homogeneous_200.md"
    lines: List[str] = []

    lines.append("# Homogeneous-region 200-death analysis after burn-in\n\n")
    lines.append(f"All 2D histograms use the same colorbar scale: 0 to {colorbar_vmax} observed deaths per bin.\n\n")
    lines.append(
        "Burn-in convention: the first 500 observed deaths are simulated but not used for statistics. "
        "The reported statistics use the next 200 observed deaths. Therefore each run accepts 700 deaths in total, "
        "then analyzes deaths 501--700.\n\n"
    )

    lines.append("## Statistical interpretation\n\n")
    lines.append(
        "Let the analyzed 200 death locations be binned into a spatial histogram with counts "
        "\\(H_1,\\ldots,H_k\\). For the 10 x 10 plot grid, \\(k=100\\) and \\(n=200\\). "
        "Under spatial uniformity, the reference CV is \\(\\sqrt{(k-1)/n}=\\sqrt{99/200}\\approx 0.704\\). "
        "Thus `spatial_cv_ratio_to_uniform_plot_bins` close to 1 is compatible with uniform random fluctuations, "
        "whereas a value much larger than 1 indicates spatial heterogeneity.\n\n"
    )
    lines.append(
        "For the formal chi-square goodness-of-fit test, we use a coarser 5 x 5 grid. "
        "Then \\(k=25\\), \\(df=24\\), and the expected count per bin under uniformity is \\(200/25=8\\), "
        "which is more reliable than using the 10 x 10 grid where the expected count would be only 2. "
        "The Pearson statistic is \\(X^2=\\sum_i (O_i-E_i)^2/E_i\\), and the p-value is "
        "\\(P(\\chi^2_{df}\\ge X^2_{obs})\\).\n\n"
    )

    lines.append("## Results\n\n")
    lines.append(
        "| scenario | time for analyzed 200 deaths | death rejection rate | CV 10x10 | CV ratio | chi2/df 5x5 | p-value 5x5 | empty bins 10x10 |\n"
        "|---|---:|---:|---:|---:|---:|---:|---:|\n"
    )
    for r in rows:
        lines.append(
            f"| {r['scenario']} | "
            f"{float(r['total_time_for_analyzed_200_deaths_after_burnin']):.3f} | "
            f"{float(r['death_rejection_rate']):.3f} | "
            f"{float(r['spatial_coefficient_of_variation_plot_bins']):.3f} | "
            f"{float(r['spatial_cv_ratio_to_uniform_plot_bins']):.3f} | "
            f"{float(r['spatial_chi2_per_df_chi2_bins']):.3f} | "
            f"{float(r['spatial_chi2_p_value_chi2_bins']):.3g} | "
            f"{int(r['spatial_empty_bins_plot_bins'])} |\n"
        )

    lines.append("\n## Short conclusion\n\n")
    baseline = rows[0]
    most_heterogeneous = max(rows, key=lambda r: float(r["spatial_cv_ratio_to_uniform_plot_bins"]))
    lines.append(
        f"The baseline case has CV ratio = {float(baseline['spatial_cv_ratio_to_uniform_plot_bins']):.3f} "
        f"and chi2/df = {float(baseline['spatial_chi2_per_df_chi2_bins']):.3f}. "
        "These quantities should be interpreted relative to the uniform references, not as isolated numbers.\n\n"
    )
    lines.append(
        f"The most spatially heterogeneous scenario by CV ratio is `{most_heterogeneous['scenario']}`, "
        f"with CV ratio = {float(most_heterogeneous['spatial_cv_ratio_to_uniform_plot_bins']):.3f}.\n"
    )

    path.write_text("".join(lines), encoding="utf-8")
    return path


# ============================================================
# 5. Parameter sweep: play with various parameters
# ============================================================


def make_scenarios(base: Params) -> Dict[str, Params]:
    """Interpretable parameter changes around the homogeneous baseline."""
    return {
        "baseline": base,
        "higher_death_intensity": replace(base, name="higher_death_intensity", seed=124, lambda_d=0.35),
        "stronger_ERK_protection": replace(base, name="stronger_ERK_protection", seed=125, beta_d_R=0.8, beta_d_T=0.4),
        "stronger_activation_feedback": replace(base, name="stronger_activation_feedback", seed=126, lambda_a_1=0.90),
        "weaker_ERK_protection": replace(base, name="weaker_ERK_protection", seed=127, beta_d_R=1.8, beta_d_T=1.2),
    }


# ============================================================
# 6. Main
# ============================================================


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    base = Params(name="baseline", burn_in_deaths=500, target_deaths=200)
    scenarios = make_scenarios(base)

    print("Running homogeneous-region rejection-sampling simulations: burn-in 500 deaths, then analyze 200 deaths...\n")

    results: List[SimulationResult] = []
    for name, params in scenarios.items():
        print(f"Scenario: {name}")
        result = simulate_one(params)
        results.append(result)
        print(f"  total observed deaths simulated: {result.total_deaths_simulated}")
        print(f"  analyzed observed deaths:        {result.analyzed_observed_deaths}")
        print(f"  burn-in time:                    {result.burn_in_time:.6f}")
        print(f"  time for analyzed 200 deaths:    {result.analysis_time_span:.6f}")
        print("")

    # Common colorbar scale across all scenarios.
    global_colorbar_max = int(max(np.max(spatial_histogram(result, result.params.space_bins)) for result in results))
    global_colorbar_max = max(global_colorbar_max, 1)

    summary_rows: List[Dict[str, Any]] = []
    for result in results:
        stats = death_statistics(result)
        summary_rows.append(stats)

        save_death_events_csv(result, OUTPUT_DIR)
        save_params_and_summary(result, stats, OUTPUT_DIR)
        plot_time_histogram(result, OUTPUT_DIR)
        plot_spatial_histogram_2d(result, OUTPUT_DIR, global_colorbar_max)
        plot_death_locations(result, OUTPUT_DIR)
        plot_cumulative_deaths(result, OUTPUT_DIR)

        print(f"Scenario: {result.name}")
        print(f"  death rejection rate:   {stats['death_rejection_rate']:.4f}")
        print(f"  spatial CV 10x10:       {stats['spatial_coefficient_of_variation_plot_bins']:.4f}")
        print(f"  CV ratio to uniform:    {stats['spatial_cv_ratio_to_uniform_plot_bins']:.4f}")
        print(f"  chi2/df 5x5:            {stats['spatial_chi2_per_df_chi2_bins']:.4f}")
        print(f"  chi2 p-value 5x5:       {stats['spatial_chi2_p_value_chi2_bins']:.4g}")
        print("")

    summary_path = save_summary_csv(summary_rows, OUTPUT_DIR)
    comparison_path = plot_comparison_total_time(summary_rows, OUTPUT_DIR)
    report_path = write_analysis_report(summary_rows, OUTPUT_DIR, global_colorbar_max)

    print("Done.")
    print(f"Common colorbar vmax: {global_colorbar_max}")
    print(f"Summary CSV: {summary_path}")
    print(f"Comparison plot: {comparison_path}")
    print(f"Analysis report: {report_path}")
    print(f"Outputs are in: {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
