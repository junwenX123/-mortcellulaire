import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Circle


# ============================================================
# Offline construction with positive feedback
#
# Goal:
#   Construct V_t by first simulating a dominating PPP of candidates,
#   then accepting/rejecting candidates according to the current active zone.
#
# Positive feedback intensity:
#
#   lambda_nu(y) = lambda_c + (lambda_1 - lambda_c) 1_{y in A(nu)}
#
# where
#
#   A(nu) = W cap union_i B(Y_i, R_i).
#
# Thinning construction:
#   1. Simulate candidate points with maximal intensity lambda_1.
#   2. Process candidates in chronological order.
#   3. If candidate X is inside current active zone: accept.
#   4. If candidate X is outside current active zone: accept with probability lambda_c/lambda_1.
#   5. If accepted, it becomes a caspase event and creates a new active center.
#
# ============================================================


# -----------------------------
# 1. Parameters
# -----------------------------

seed = 123
rng = np.random.default_rng(seed)

# Spatial domain W = [0, Lx] x [0, Ly]
Lx = 10.0
Ly = 10.0
area_W = Lx * Ly

# Time horizon [0, T]
T = 20.0

# Background intensity outside active zones
lambda_c = 0.05

# Higher intensity inside active zones
lambda_1 = 0.60

if lambda_1 <= lambda_c:
    raise ValueError("We need lambda_1 > lambda_c.")

# Marks:
# R_k ~ Exp(beta_R)
# tau_k ~ Exp(beta_tau)
beta_R = 2.5
beta_tau = 1.2

# Animation parameters
n_frames = 250
show_history = True
show_rejected_candidates = True

save_gif = True
gif_name = "vt_positive_feedback_offline.gif"


# -----------------------------
# 2. Helper functions
# -----------------------------

def is_inside_active_zone(x, y, active_x, active_y, active_r):
    """
    Check whether point (x,y) belongs to the current active zone
    
        A(nu) = W cap union_i B(Y_i, R_i).

    Here active_x, active_y, active_r describe the currently active centers.
    """
    if len(active_x) == 0:
        return False

    ax_arr = np.asarray(active_x)
    ay_arr = np.asarray(active_y)
    ar_arr = np.asarray(active_r)

    dist2 = (x - ax_arr) ** 2 + (y - ay_arr) ** 2
    return np.any(dist2 <= ar_arr ** 2)


def remove_expired_centers(current_time, active_x, active_y, active_r, active_death):
    """
    Keep only centers that are still active at current_time.

    A center born at S with lifetime tau is active iff

        S <= current_time < S + tau.

    During chronological construction, all centers in active_* have already been born.
    So we only need to remove those with death_time <= current_time.
    """
    keep_x = []
    keep_y = []
    keep_r = []
    keep_death = []

    for x, y, r, d in zip(active_x, active_y, active_r, active_death):
        if current_time < d:
            keep_x.append(x)
            keep_y.append(y)
            keep_r.append(r)
            keep_death.append(d)

    return keep_x, keep_y, keep_r, keep_death


# -----------------------------
# 3. Offline simulation of dominating PPP
# -----------------------------

# Dominating candidate PPP intensity:
#   lambda_1 dx dt on W x [0,T]
#
# Number of candidate points:
#   N_prop ~ Poisson(lambda_1 |W| T)

mean_number_candidates = lambda_1 * area_W * T
N_prop = rng.poisson(mean_number_candidates)

print(f"Nombre total de candidats simulés : N_prop = {N_prop}")

# Candidate positions and times
cand_x = rng.uniform(0.0, Lx, size=N_prop)
cand_y = rng.uniform(0.0, Ly, size=N_prop)
cand_s = rng.uniform(0.0, T, size=N_prop)

# Process candidates in chronological order
order = np.argsort(cand_s)
cand_x = cand_x[order]
cand_y = cand_y[order]
cand_s = cand_s[order]


# -----------------------------
# 4. Chronological thinning with positive feedback
# -----------------------------

# Accepted caspase events / active centers:
# accepted center k has:
#   position (X_acc[k], Y_acc[k])
#   birth time S_acc[k]
#   radius R_acc[k]
#   lifetime tau_acc[k]
#   death time D_acc[k] = S_acc[k] + tau_acc[k]

X_acc = []
Y_acc = []
S_acc = []
R_acc = []
tau_acc = []
D_acc = []

# Rejected candidates, only for optional display
X_rej = []
Y_rej = []
S_rej = []

# Current active configuration used during thinning
current_active_x = []
current_active_y = []
current_active_r = []
current_active_death = []

for x, y, s in zip(cand_x, cand_y, cand_s):

    # Before deciding whether candidate at time s is accepted,
    # update the current active configuration V_{s-}.
    current_active_x, current_active_y, current_active_r, current_active_death = (
        remove_expired_centers(
            s,
            current_active_x,
            current_active_y,
            current_active_r,
            current_active_death
        )
    )

    # Determine whether candidate lies in current active zone A(V_{s-})
    inside = is_inside_active_zone(
        x,
        y,
        current_active_x,
        current_active_y,
        current_active_r
    )

    # Thinning rule:
    #   if x in A(V_{s-}), accept with probability 1
    #   else accept with probability lambda_c/lambda_1
    if inside:
        accept = True
    else:
        accept = rng.uniform(0.0, 1.0) <= (lambda_c / lambda_1)

    if accept:
        # Draw marks for the accepted center
        r = rng.exponential(scale=1.0 / beta_R)
        tau = rng.exponential(scale=1.0 / beta_tau)
        death_time = s + tau

        # Store accepted center
        X_acc.append(x)
        Y_acc.append(y)
        S_acc.append(s)
        R_acc.append(r)
        tau_acc.append(tau)
        D_acc.append(death_time)

        # Add it immediately to current active configuration
        # since it is born at time s and active on [s, s+tau)
        current_active_x.append(x)
        current_active_y.append(y)
        current_active_r.append(r)
        current_active_death.append(death_time)

    else:
        X_rej.append(x)
        Y_rej.append(y)
        S_rej.append(s)


# Convert lists to arrays
X_acc = np.asarray(X_acc)
Y_acc = np.asarray(Y_acc)
S_acc = np.asarray(S_acc)
R_acc = np.asarray(R_acc)
tau_acc = np.asarray(tau_acc)
D_acc = np.asarray(D_acc)

X_rej = np.asarray(X_rej)
Y_rej = np.asarray(Y_rej)
S_rej = np.asarray(S_rej)

print(f"Nombre de candidats acceptés : {len(X_acc)}")
print(f"Nombre de candidats rejetés  : {len(X_rej)}")


# -----------------------------
# 5. Construct V_t by direct filtering
# -----------------------------

def active_mask(t):
    """
    Active centers at time t:

        S_k <= t < S_k + tau_k.
    """
    if len(S_acc) == 0:
        return np.zeros(0, dtype=bool)
    return (S_acc <= t) & (t < D_acc)


def appeared_mask(t):
    """
    Accepted centers already born before time t:

        S_k <= t.
    """
    if len(S_acc) == 0:
        return np.zeros(0, dtype=bool)
    return S_acc <= t


def rejected_mask(t):
    """
    Rejected candidates that have appeared before time t.
    Only used for optional display.
    """
    if len(S_rej) == 0:
        return np.zeros(0, dtype=bool)
    return S_rej <= t


# -----------------------------
# 6. Figure setup
# -----------------------------

fig, ax = plt.subplots(figsize=(7, 7))

ax.set_xlim(0, Lx)
ax.set_ylim(0, Ly)
ax.set_aspect("equal", adjustable="box")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Offline construction of $V_t$ with positive feedback")

# Rejected candidates, optional
rejected_scatter = ax.scatter(
    [],
    [],
    s=10,
    alpha=0.15,
    label="candidats rejetés"
)

# Accepted but currently inactive centers
history_scatter = ax.scatter(
    [],
    [],
    s=15,
    alpha=0.25,
    label="centres acceptés non actifs"
)

# Currently active centers
active_scatter = ax.scatter(
    [],
    [],
    s=40,
    label="centres actifs $V_t$"
)

time_text = ax.text(
    0.02,
    0.98,
    "",
    transform=ax.transAxes,
    va="top",
    ha="left"
)

count_text = ax.text(
    0.02,
    0.93,
    "",
    transform=ax.transAxes,
    va="top",
    ha="left"
)

rate_text = ax.text(
    0.02,
    0.88,
    "",
    transform=ax.transAxes,
    va="top",
    ha="left"
)

ax.legend(loc="upper right")

# Circles representing active zones
active_circles = []


# -----------------------------
# 7. Animation functions
# -----------------------------

def init():
    rejected_scatter.set_offsets(np.empty((0, 2)))
    history_scatter.set_offsets(np.empty((0, 2)))
    active_scatter.set_offsets(np.empty((0, 2)))
    time_text.set_text("")
    count_text.set_text("")
    rate_text.set_text("")
    return rejected_scatter, history_scatter, active_scatter, time_text, count_text, rate_text


def update(frame):
    global active_circles

    # Current time
    t = frame * T / (n_frames - 1)

    # Remove old active-zone circles
    for circle in active_circles:
        circle.remove()
    active_circles = []

    # Current active centers
    is_active = active_mask(t)

    # Accepted centers already appeared but currently inactive
    has_appeared = appeared_mask(t)
    inactive_appeared = has_appeared & (~is_active)

    # Optional rejected candidates
    if show_rejected_candidates:
        rej = rejected_mask(t)
        rejected_points = np.column_stack([X_rej[rej], Y_rej[rej]])
        rejected_scatter.set_offsets(rejected_points)
    else:
        rejected_scatter.set_offsets(np.empty((0, 2)))

    # History: accepted but inactive
    if show_history:
        history_points = np.column_stack([
            X_acc[inactive_appeared],
            Y_acc[inactive_appeared]
        ])
        history_scatter.set_offsets(history_points)
    else:
        history_scatter.set_offsets(np.empty((0, 2)))

    # Current support of V_t
    active_points = np.column_stack([
        X_acc[is_active],
        Y_acc[is_active]
    ])
    active_scatter.set_offsets(active_points)

    # Draw active balls B(Y_i, R_i)
    for xk, yk, rk in zip(X_acc[is_active], Y_acc[is_active], R_acc[is_active]):
        circle = Circle(
            (xk, yk),
            rk,
            fill=False,
            alpha=0.40,
            linewidth=1.4
        )
        ax.add_patch(circle)
        active_circles.append(circle)

    # Approximate current |A(V_t)| by Monte Carlo for display only
    # This is not used in the simulation.
    n_active = int(np.sum(is_active))
    if n_active > 0:
        n_mc = 1000
        mc_x = rng.uniform(0.0, Lx, size=n_mc)
        mc_y = rng.uniform(0.0, Ly, size=n_mc)

        inside_count = 0
        active_x_now = X_acc[is_active]
        active_y_now = Y_acc[is_active]
        active_r_now = R_acc[is_active]

        for mx, my in zip(mc_x, mc_y):
            if is_inside_active_zone(mx, my, active_x_now, active_y_now, active_r_now):
                inside_count += 1

        area_A_est = area_W * inside_count / n_mc
    else:
        area_A_est = 0.0

    # Real accepted birth rate:
    #   int_W lambda_nu(y) dy
    birth_rate_est = lambda_c * area_W + (lambda_1 - lambda_c) * area_A_est

    time_text.set_text(f"t = {t:.2f}")
    count_text.set_text(
        f"N_t = {n_active} active centers | "
        f"N_acc = {len(X_acc)} accepted"
    )
    rate_text.set_text(
        r"$\lambda_c |W| + (\lambda_1-\lambda_c)|A(V_t)|$"
        f" ≈ {birth_rate_est:.3f}"
    )

    return (
        rejected_scatter,
        history_scatter,
        active_scatter,
        time_text,
        count_text,
        rate_text,
        *active_circles
    )


# Matplotlib animation calls update(frame) repeatedly, one frame at a time.
# This only displays the already-constructed process;
anim = FuncAnimation(
    fig,
    update,
    frames=n_frames,
    init_func=init,
    interval=60,
    blit=False
)

if save_gif:
    writer = PillowWriter(fps=20)
    anim.save(gif_name, writer=writer)
    print(f"Animation sauvegardée dans : {gif_name}")

plt.show()