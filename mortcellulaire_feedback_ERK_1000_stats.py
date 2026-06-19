"""
Gillespie / event-driven simulation for caspase activity, observed cell death,
and ERK negative feedback, with a fixed T-shaped background activation zone.

  1. Stop after 1000 observed cell deaths.
  2. Compute statistics on the observed death process:
       - number of observed deaths,
       - total time needed to reach 1000 observed deaths,
       - mean / median inter-death time,
       - simple spatial non-uniformity statistics.
  3. Plot a histogram of death times.
  4. Plot a 2D histogram of death locations in space.
  5. Play with several parameter sets and record the comparison in a CSV file.

Run:
    python mortcellulaire_feedback_ERK_1000_stats.py

Outputs are saved in:
    simulation_outputs_1000/
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")  # non-interactive backend: faster and safer for batch simulations

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np


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

    # Stopping rule: stop when this number of observed deaths is reached
    target_deaths: int = 1000

    # Activation intensities with fixed T-shaped extension
    # lambda_a_1 is the dominating intensity and active-zone intensity.
    lambda_a_1: float = 5.00
    lambda_a_T: float = 0.50
    lambda_a_c: float = 0.05

    # Candidate death intensity on W
    lambda_d: float = 1.0

    # Activation marks: R^a ~ Exp(beta_a_R), lifetime rate beta_a_T
    beta_a_R: float = 2.5
    beta_a_T: float = 1.2

    # ERK marks after accepted deaths: R^d ~ Exp(beta_d_R), T^d ~ Exp(beta_d_T)
    beta_d_R: float = 2.0
    beta_d_T: float = 0.8

    # Plotting/statistics choices
    bins_time: int = 30
    bins_space: int = 20

    # Safety stop, only to avoid infinite runs if bad parameters are chosen
    max_proposals: int = 5_000_000

    @property
    def area_W(self) -> float:
        return self.Lx * self.Ly


@dataclass
class SimulationResult:
    scenario: str
    params: Parameters
    death_x: np.ndarray
    death_y: np.ndarray
    death_t: np.ndarray
    activation_x: np.ndarray
    activation_y: np.ndarray
    activation_t: np.ndarray
    final_time: float
    n_observed_deaths: int
    n_accepted_activations: int
    n_activation_candidates: int
    n_death_candidates: int
    n_activation_rejected: int
    n_death_rejected: int
    n_expirations: int
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
# 3. One Gillespie / thinning simulation
# ============================================================


def simulate_one(scenario: str, p: Parameters) -> SimulationResult:
    """
    Simulate until p.target_deaths observed death events are accepted.

    Event-driven proposal rates:
        activation candidates: lambda_a_1 |W|
        death candidates:     lambda_d |W|
        active expirations:   N_t beta_a_T

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

    rng = np.random.default_rng(p.seed)
    t = 0.0

    # Active centers V^a_t = sum delta_(Y_i^a, R_i^a)
    active_x: List[float] = []
    active_y: List[float] = []
    active_r: List[float] = []

    # Accepted activations X = {(Y_i^a, S_i^a)}
    activation_x: List[float] = []
    activation_y: List[float] = []
    activation_t: List[float] = []

    # Observed deaths Y = {(Y_i^d, S_i^d)}
    death_x: List[float] = []
    death_y: List[float] = []
    death_t: List[float] = []

    # Active ERK protection zones V^p_t, active on [start, end)
    erk_x: List[float] = []
    erk_y: List[float] = []
    erk_r: List[float] = []
    erk_end_t: List[float] = []

    n_activation_candidates = 0
    n_death_candidates = 0
    n_activation_rejected = 0
    n_death_rejected = 0
    n_expirations = 0
    n_proposals = 0

    def clean_expired_erk(current_t: float) -> None:
        for i in range(len(erk_x) - 1, -1, -1):
            if current_t >= erk_end_t[i]:
                erk_x.pop(i)
                erk_y.pop(i)
                erk_r.pop(i)
                erk_end_t.pop(i)

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

    while len(death_t) < p.target_deaths and n_proposals < p.max_proposals:
        clean_expired_erk(t)

        n_active = len(active_x)
        activation_proposal_rate = p.lambda_a_1 * p.area_W
        death_proposal_rate = p.lambda_d * p.area_W
        expiration_rate = n_active * p.beta_a_T
        a0 = activation_proposal_rate + death_proposal_rate + expiration_rate

        if a0 <= 0.0:
            break

        # Gillespie waiting time
        t += rng.exponential(scale=1.0 / a0)
        clean_expired_erk(t)

        # Gillespie event-type selection
        u = rng.uniform(0.0, 1.0)
        p_activation = activation_proposal_rate / a0
        p_death = death_proposal_rate / a0

        n_proposals += 1

        if u <= p_activation:
            # Candidate activation from dominating PPP(lambda_a_1) on W
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
            # Candidate death from PPP(lambda_d) on W
            n_death_candidates += 1
            x = rng.uniform(0.0, p.Lx)
            y = rng.uniform(0.0, p.Ly)

            inside_active = is_inside_active_zone(x, y)
            inside_erk = is_inside_erk_zone(x, y)

            if inside_active and not inside_erk:
                death_x.append(x)
                death_y.append(y)
                death_t.append(t)

                # ERK feedback after accepted observed death
                r_E = rng.exponential(scale=1.0 / p.beta_d_R)
                T_E = rng.exponential(scale=1.0 / p.beta_d_T)
                erk_x.append(x)
                erk_y.append(y)
                erk_r.append(r_E)
                erk_end_t.append(t + T_E)
            else:
                n_death_rejected += 1

        else:
            # Expiration of one active center
            if len(active_x) > 0:
                idx = int(rng.integers(0, len(active_x)))
                active_x.pop(idx)
                active_y.pop(idx)
                active_r.pop(idx)
                n_expirations += 1

    return SimulationResult(
        scenario=scenario,
        params=p,
        death_x=np.asarray(death_x),
        death_y=np.asarray(death_y),
        death_t=np.asarray(death_t),
        activation_x=np.asarray(activation_x),
        activation_y=np.asarray(activation_y),
        activation_t=np.asarray(activation_t),
        final_time=t,
        n_observed_deaths=len(death_t),
        n_accepted_activations=len(activation_t),
        n_activation_candidates=n_activation_candidates,
        n_death_candidates=n_death_candidates,
        n_activation_rejected=n_activation_rejected,
        n_death_rejected=n_death_rejected,
        n_expirations=n_expirations,
        n_proposals=n_proposals,
        stopped_by_target=(len(death_t) >= p.target_deaths),
    )


# ============================================================
# 4. Statistics and plots on the observed death process
# ============================================================


def death_statistics(result: SimulationResult) -> Dict[str, float | int | str | bool]:
    p = result.params
    death_t = result.death_t
    death_x = result.death_x
    death_y = result.death_y

    inter_death = np.diff(death_t) if len(death_t) >= 2 else np.array([])

    H, _, _ = np.histogram2d(
        death_x,
        death_y,
        bins=p.bins_space,
        range=[[0.0, p.Lx], [0.0, p.Ly]],
    )

    expected = result.n_observed_deaths / float(p.bins_space * p.bins_space)
    chi2_uniform = float(np.sum((H - expected) ** 2 / expected)) if expected > 0 else np.nan
    spatial_cv = float(np.std(H) / np.mean(H)) if np.mean(H) > 0 else np.nan
    empty_bins = int(np.sum(H == 0))

    return {
        "scenario": result.scenario,
        "seed": p.seed,
        "target_deaths": p.target_deaths,
        "n_observed_deaths": result.n_observed_deaths,
        "stopped_by_target": result.stopped_by_target,
        "total_time_for_1000_deaths": result.final_time,
        "observed_death_rate_per_time": result.n_observed_deaths / result.final_time,
        "mean_inter_death_time": float(np.mean(inter_death)) if len(inter_death) else np.nan,
        "median_inter_death_time": float(np.median(inter_death)) if len(inter_death) else np.nan,
        "accepted_activations": result.n_accepted_activations,
        "activation_candidates": result.n_activation_candidates,
        "death_candidates": result.n_death_candidates,
        "activation_rejected": result.n_activation_rejected,
        "death_rejected": result.n_death_rejected,
        "active_center_expirations": result.n_expirations,
        "all_proposals": result.n_proposals,
        "spatial_bins_per_axis": p.bins_space,
        "spatial_expected_count_per_bin_if_uniform": expected,
        "spatial_chi2_statistic_against_uniform_counts": chi2_uniform,
        "spatial_coefficient_of_variation": spatial_cv,
        "spatial_empty_bins": empty_bins,
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
    path = out_dir / f"{result.scenario}_observed_deaths.csv"
    data = np.column_stack([
        np.arange(1, result.n_observed_deaths + 1),
        result.death_t,
        result.death_x,
        result.death_y,
    ])
    header = "death_index,death_time,x,y"
    np.savetxt(path, data, delimiter=",", header=header, comments="", fmt=["%d", "%.10f", "%.10f", "%.10f"])
    return path


def add_t_zone_overlay(ax: plt.Axes, p: Parameters) -> None:
    """Add transparent rectangles showing the fixed T-zone on a spatial plot."""
    x1, x2, y1, y2 = t_zone_bounds(p)
    ax.add_patch(Rectangle((x1, 0.0), x2 - x1, p.Ly, fill=False, linewidth=1.2, linestyle="--"))
    ax.add_patch(Rectangle((0.0, y1), x1, y2 - y1, fill=False, linewidth=1.2, linestyle="--"))


def plot_time_histogram(result: SimulationResult, out_dir: Path) -> Path:
    p = result.params
    path = out_dir / f"{result.scenario}_time_histogram.png"

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(result.death_t, bins=p.bins_time, range=(0.0, result.final_time))
    ax.set_xlabel("time")
    ax.set_ylabel("number of observed deaths")
    ax.set_title(
        f"Observed cell death times — {result.scenario}\n"
        f"N = {result.n_observed_deaths}, total time = {result.final_time:.3f}"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_spatial_histogram_2d(result: SimulationResult, out_dir: Path) -> Path:
    p = result.params
    path = out_dir / f"{result.scenario}_space_histogram_2d.png"

    H, x_edges, y_edges = np.histogram2d(
        result.death_x,
        result.death_y,
        bins=p.bins_space,
        range=[[0.0, p.Lx], [0.0, p.Ly]],
    )

    fig, ax = plt.subplots(figsize=(7, 6))
    image = ax.imshow(
        H.T,
        origin="lower",
        extent=[0.0, p.Lx, 0.0, p.Ly],
        aspect="equal",
        interpolation="nearest",
    )
    add_t_zone_overlay(ax, p)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(
        f"2D histogram of observed death locations — {result.scenario}\n"
        f"{p.bins_space} x {p.bins_space} bins, N = {result.n_observed_deaths}"
    )
    fig.colorbar(image, ax=ax, label="number of observed deaths")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_cumulative_deaths(result: SimulationResult, out_dir: Path) -> Path:
    path = out_dir / f"{result.scenario}_cumulative_deaths.png"
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.step(result.death_t, np.arange(1, result.n_observed_deaths + 1), where="post")
    ax.set_xlabel("time")
    ax.set_ylabel("cumulative observed deaths")
    ax.set_title(f"Cumulative observed deaths — {result.scenario}")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def save_summary_csv(rows: List[Dict[str, float | int | str | bool]], out_dir: Path) -> Path:
    path = out_dir / "parameter_sweep_summary.csv"
    keys = list(rows[0].keys())
    with path.open("w", encoding="utf-8") as f:
        f.write(",".join(keys) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) for k in keys) + "\n")
    return path


# ============================================================
# 5. Parameter sweep: "play with various parameters"
# ============================================================


def make_scenarios(base: Parameters) -> Dict[str, Parameters]:
    """A few interpretable parameter changes around the baseline."""
    return {
        "baseline": base,

        # Larger lambda_d: death candidates arrive faster, so 1000 deaths should occur earlier.
        "higher_death_rate": replace(base, seed=124, lambda_d=1.5),

        # Larger mean ERK radius because mean R^d = 1 / beta_d_R.
        # This should reject more nearby death candidates and increase total time.
        "larger_ERK_radius": replace(base, seed=125, beta_d_R=1.0),

        # Shorter mean ERK duration because mean T^d = 1 / beta_d_T.
        # This should reduce protection time and can make deaths arrive faster.
        "shorter_ERK_duration": replace(base, seed=126, beta_d_T=2.0),

        # Stronger activation in the fixed T-zone.
        # This should create more active zones near the T-shape and may make space less uniform.
        "stronger_T_zone_activation": replace(base, seed=127, lambda_a_T=1.0),
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
    out_dir = Path("simulation_outputs_1000")
    out_dir.mkdir(parents=True, exist_ok=True)

    base = Parameters(target_deaths=1000)
    scenarios = make_scenarios(base)

    summary_rows: List[Dict[str, float | int | str | bool]] = []

    print("Running Gillespie simulations until 1000 observed deaths...\n")

    for name, params in scenarios.items():
        print(f"Scenario: {name}")
        result = simulate_one(name, params)
        stats = death_statistics(result)
        summary_rows.append(stats)

        save_death_events_csv(result, out_dir)
        plot_time_histogram(result, out_dir)
        plot_spatial_histogram_2d(result, out_dir)
        plot_cumulative_deaths(result, out_dir)

        print(f"  observed deaths: {result.n_observed_deaths}")
        print(f"  total time:       {result.final_time:.6f}")
        print(f"  activations:      {result.n_accepted_activations}")
        print(f"  spatial CV:       {stats['spatial_coefficient_of_variation']:.4f}")
        print("")

    summary_path = save_summary_csv(summary_rows, out_dir)

    print("Done.")
    print(f"Summary CSV: {summary_path}")
    print(f"Plots and death-event CSV files are in: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
