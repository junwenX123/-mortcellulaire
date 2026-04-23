import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle


def simulate_caspase_process(
    Lx=10.0,
    Ly=6.0,
    T=12.0,
    lambda_c=0.03,        # intensity of initial active centers in space-time
    lambda_init=0.18,     # offspring intensity of generation-0 centers
    lambda_gain=1.35,     # each caspase-born center gets a higher lambda
    lambda_max=1.20,      # cap to avoid explosion
    beta_R=1.40,          # rate of Exp(beta_R) for R_k
    beta_tau=0.40,        # rate of Exp(beta_tau) for tau_k
    max_generation=4,     # safety cutoff for recursive excitation
    max_centers=4000,     # safety cutoff for recursive excitation
    seed=None,
):
    """
    Recursive self-exciting caspase model.

    1) Initial active centers are born from a Poisson point process on W x [0,T].
    2) Each active center k has a disk radius R_k and lifetime tau_k.
    3) During its active period, it generates caspase events according to a PPP
       with intensity lambda_k inside its own space-time cylinder.
    4) Every generated caspase event becomes a NEW active center.
    5) The new center receives a larger offspring intensity:
           lambda_child = min(lambda_gain * lambda_parent, lambda_max)

    """

    rng = np.random.default_rng(seed)
    area = Lx * Ly

    # ------------------------------------------------------------
    # 1) Initial active centers (generation 0)
    # ------------------------------------------------------------
    n0 = rng.poisson(lambda_c * area * T)

    centers = {
        "x": list(rng.uniform(0.0, Lx, size=n0)),
        "y": list(rng.uniform(0.0, Ly, size=n0)),
        "s": list(rng.uniform(0.0, T, size=n0)),
        "R": list(rng.exponential(scale=1.0 / beta_R, size=n0)),
        "tau": list(rng.exponential(scale=1.0 / beta_tau, size=n0)),
        "lambda": [lambda_init] * n0,
        "generation": [0] * n0,
        "parent": [-1] * n0,
    }

    # all caspase events (every one of these will also become a new center)
    events = {
        "x": [],
        "y": [],
        "t": [],
        "parent_center": [],
        "generation": [],
    }

    # ------------------------------------------------------------
    # 2) Recursive generation of caspase events / new active centers
    # ------------------------------------------------------------
    k = 0
    while k < len(centers["x"]):
        if len(centers["x"]) >= max_centers:
            print(f"[warning] max_centers={max_centers} reached; simulation truncated.")
            break

        x0 = centers["x"][k]
        y0 = centers["y"][k]
        s0 = centers["s"][k]
        r0 = centers["R"][k]
        tau0 = centers["tau"][k]
        lam0 = centers["lambda"][k]
        gen0 = centers["generation"][k]

        # stop recursion after some generation
        if gen0 >= max_generation:
            k += 1
            continue

        mean_k = lam0 * np.pi * r0**2 * tau0
        nk = rng.poisson(mean_k)

        if nk == 0:
            k += 1
            continue

        rr = r0 * np.sqrt(rng.uniform(0.0, 1.0, size=nk))
        theta = rng.uniform(0.0, 2.0 * np.pi, size=nk)
        xk = x0 + rr * np.cos(theta)
        yk = y0 + rr * np.sin(theta)
        tk = s0 + rng.uniform(0.0, tau0, size=nk)

        keep = (
            (0.0 <= xk) & (xk <= Lx) &
            (0.0 <= yk) & (yk <= Ly) &
            (0.0 <= tk) & (tk <= T)
        )

        xk = xk[keep]
        yk = yk[keep]
        tk = tk[keep]

        if xk.size == 0:
            k += 1
            continue

        child_lambda = min(lambda_gain * lam0, lambda_max)

        # register events
        events["x"].extend(xk.tolist())
        events["y"].extend(yk.tolist())
        events["t"].extend(tk.tolist())
        events["parent_center"].extend([k] * xk.size)
        events["generation"].extend([gen0 + 1] * xk.size)

        # every event becomes a new active center
        child_R = rng.exponential(scale=1.0 / beta_R, size=xk.size)
        child_tau = rng.exponential(scale=1.0 / beta_tau, size=xk.size)

        centers["x"].extend(xk.tolist())
        centers["y"].extend(yk.tolist())
        centers["s"].extend(tk.tolist())
        centers["R"].extend(child_R.tolist())
        centers["tau"].extend(child_tau.tolist())
        centers["lambda"].extend([child_lambda] * xk.size)
        centers["generation"].extend([gen0 + 1] * xk.size)
        centers["parent"].extend([k] * xk.size)

        if len(centers["x"]) >= max_centers:
            print(f"[warning] max_centers={max_centers} reached; simulation truncated.")
            break

        k += 1

    # convert to numpy arrays
    centers = {key: np.array(val) for key, val in centers.items()}
    events = {
        "x": np.array(events["x"]),
        "y": np.array(events["y"]),
        "t": np.array(events["t"]),
        "parent_center": np.array(events["parent_center"], dtype=int),
        "generation": np.array(events["generation"], dtype=int),
    }

    return centers, events


def plot_snapshot(
    centers,
    events,
    Lx,
    Ly,
    t_obs=7.0,
    dt_window=1.0,
    show_all_centers=False,
    show_parent_links=False,
):
    """
    Plot one time slice.

    - only active centers are shown as blue '+' and blue circles
    - red points are caspase events near t_obs
    """

    active = (centers["s"] <= t_obs) & (t_obs <= centers["s"] + centers["tau"])
    near_t = np.abs(events["t"] - t_obs) <= dt_window / 2.0

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.add_patch(Rectangle((0, 0), Lx, Ly, fill=False, edgecolor="black", linewidth=2.0))

    if show_all_centers:
        ax.scatter(
            centers["x"], centers["y"],
            marker="x", s=22, color="lightsteelblue", alpha=0.35,
            label="tous les centres actifs (toutes générations)"
        )

    # active centers at t_obs
    first_active = True
    for x0, y0, r in zip(centers["x"][active], centers["y"][active], centers["R"][active]):
        circ = Circle((x0, y0), r, fill=False, edgecolor="steelblue", linewidth=1.6, alpha=0.80)
        ax.add_patch(circ)
        if first_active:
            ax.scatter(x0, y0, marker="+", s=140, linewidths=2.2, color="navy", label="centre actif")
            first_active = False
        else:
            ax.scatter(x0, y0, marker="+", s=140, linewidths=2.2, color="navy")

    # caspase events near t_obs
    if np.any(near_t):
        ax.scatter(
            events["x"][near_t], events["y"][near_t],
            s=28, color="crimson", alpha=0.88, label="caspase"
        )

        if show_parent_links:
            idxs = np.where(near_t)[0]
            for j in idxs:
                p = events["parent_center"][j]
                ax.plot(
                    [centers["x"][p], events["x"][j]],
                    [centers["y"][p], events["y"][j]],
                    color="lightcoral", alpha=0.25, linewidth=0.8,
                )

    ax.set_aspect("equal")
    ax.set_xlim(-0.2, Lx + 0.2)
    ax.set_ylim(-0.2, Ly + 0.2)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(
        rf"Snapshot du processus auto-excitant autour de $t_{{obs}}={t_obs}$"
        + "\n"
        + rf"(points affichés si $|t-t_{{obs}}|\leq {dt_window/2:.2f}$)"
    )
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    Lx, Ly = 10.0, 6.0
    T = 12.0

    centers, events = simulate_caspase_process(
        Lx=Lx,
        Ly=Ly,
        T=T,
        lambda_c=0.02,
        lambda_init=0.10,
        lambda_gain=1.35,
        lambda_max=0.55,
        beta_R=1.40,
        beta_tau=0.40,
        max_generation=4,
        max_centers=3000,
        seed=None,
    )

    print("Nombre total de centres actifs (toutes générations) :", len(centers["x"]))
    print("Nombre total d'événements caspase :", len(events["x"]))
    print("Génération maximale observée :", centers["generation"].max() if len(centers["generation"]) else 0)

    plot_snapshot(
        centers,
        events,
        Lx=Lx,
        Ly=Ly,
        t_obs=7.0,
        dt_window=1.0,
        show_all_centers=True,
        show_parent_links=True,
    )
