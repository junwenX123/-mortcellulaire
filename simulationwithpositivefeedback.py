import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


# ============================================================
# Gillespie / event-driven simulation with positive feedback
#
# Inside active zones:      lambda_cas = lambda_1
# Outside active zones:     lambda_cas = lambda_c
#
# Rejection/thinning:
#   Simulate candidate caspase points with rate lambda_1.
#   If candidate is inside active zone: keep it.
#   If candidate is outside active zone: keep it with probability lambda_c/lambda_1.
#
# If a candidate is kept:
#   it becomes a caspase event,
#   and it creates a new active center.
#
# Active centers expire with rate beta_tau.
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

# Final time
T_max = 20.0

# Background intensity outside active zones
lambda_c = 0.05

# Higher intensity inside active zones
lambda_1 = 0.60

if lambda_1 <= lambda_c:
    raise ValueError("We need lambda_1 > lambda_c.")

# Radius distribution: R ~ Exp(beta_R)
# Larger beta_R means smaller average radius.
beta_R = 2.5

# Expiration rate of each active center
# Larger beta_tau means faster expiration.
beta_tau = 1.2

# Pause between visual updates
display_pause = 0.08

# Do not draw rejected candidates
show_rejected_points = False


# -----------------------------
# 2. State variables
# -----------------------------

t = 0.0

# Active centers V_t = sum delta_(Y_i, R_i)
active_x = []
active_y = []
active_r = []

# Accepted caspase events N_cas
cas_x = []
cas_y = []
cas_t = []

# Histories
event_times = [0.0]
active_counts = [0]
cas_counts = [0]


# -----------------------------
# 3. Helper functions
# -----------------------------

def is_inside_active_zone(x, y):
    """
    Return True if point (x,y) is in A_t = union_i B(Y_i,R_i).
    """
    for cx, cy, r in zip(active_x, active_y, active_r):
        if (x - cx) ** 2 + (y - cy) ** 2 <= r ** 2:
            return True
    return False


def add_active_center(x, y):
    """
    Accepted caspase event creates a new active center.
    """
    r = rng.exponential(scale=1.0 / beta_R)
    active_x.append(x)
    active_y.append(y)
    active_r.append(r)


def remove_random_active_center():
    """
    Conditional on an expiration event, every active center has the same rate.
    Therefore choose one uniformly.
    """
    n = len(active_x)

    if n == 0:
        return None

    idx = rng.integers(0, n)

    active_x.pop(idx)
    active_y.pop(idx)
    active_r.pop(idx)

    return idx


# -----------------------------
# 4. Matplotlib real-time setup
# -----------------------------

plt.ion()

fig, ax = plt.subplots(figsize=(7, 7))

ax.set_xlim(0, Lx)
ax.set_ylim(0, Ly)
ax.set_aspect("equal", adjustable="box")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Gillespie simulation with positive feedback")

# Accepted caspase events
cas_scatter = ax.scatter(
    [],
    [],
    s=18,
    label=r"$N_{\mathrm{cas}}$ accepted events"
)

# Active centers
active_scatter = ax.scatter(
    [],
    [],
    s=55,
    label=r"active centers $V_t$"
)

time_text = ax.text(
    0.02,
    0.98,
    "",
    transform=ax.transAxes,
    va="top",
    ha="left"
)

event_text = ax.text(
    0.02,
    0.93,
    "",
    transform=ax.transAxes,
    va="top",
    ha="left"
)

count_text = ax.text(
    0.02,
    0.88,
    "",
    transform=ax.transAxes,
    va="top",
    ha="left"
)

rate_text = ax.text(
    0.02,
    0.83,
    "",
    transform=ax.transAxes,
    va="top",
    ha="left"
)

ax.legend(loc="upper right")

active_circles = []


def redraw(event_name):
    """
    Redraw current state:
      - accepted caspase events
      - active centers
      - active zones
    """
    global active_circles

    # Remove old active-zone circles
    for circle in active_circles:
        circle.remove()
    active_circles = []

    # Draw accepted caspase events
    if len(cas_x) > 0:
        cas_points = np.column_stack([cas_x, cas_y])
        cas_scatter.set_offsets(cas_points)
    else:
        cas_scatter.set_offsets(np.empty((0, 2)))

    # Draw active centers
    if len(active_x) > 0:
        active_points = np.column_stack([active_x, active_y])
        active_scatter.set_offsets(active_points)
    else:
        active_scatter.set_offsets(np.empty((0, 2)))

    # Draw current active zones only
    for x, y, r in zip(active_x, active_y, active_r):
        circle = Circle(
            (x, y),
            r,
            fill=False,
            alpha=0.35,
            linewidth=1.3
        )
        ax.add_patch(circle)
        active_circles.append(circle)

    n_t = len(active_x)
    candidate_rate = lambda_1 * area_W
    expiration_rate = n_t * beta_tau
    a0 = candidate_rate + expiration_rate

    time_text.set_text(f"t = {t:.3f}")
    event_text.set_text(f"dernier événement : {event_name}")
    count_text.set_text(
        f"N_t = {n_t} active centers | "
        f"N_cas = {len(cas_x)} accepted"
    )
    rate_text.set_text(
        f"a0 = λ1|W| + N_tβτ = {a0:.3f}"
    )

    fig.canvas.draw_idle()
    plt.pause(display_pause)


# Initial display
redraw("initialisation")


# -----------------------------
# 5. Gillespie / event-driven loop
# -----------------------------

while t < T_max:

    n_t = len(active_x)

    # Candidate caspase event rate from dominating PPP(lambda_1)
    candidate_rate = lambda_1 * area_W

    # Expiration rate of active centers
    expiration_rate = n_t * beta_tau

    # Total rate
    a0 = candidate_rate + expiration_rate

    if a0 <= 0:
        break

    # Waiting time until next event
    delta_t = rng.exponential(scale=1.0 / a0)

    # Update time
    t += delta_t

    if t > T_max:
        break

    # Choose event type
    u = rng.uniform(0.0, 1.0)

    p_candidate = candidate_rate / a0

    if u <= p_candidate:
        # ----------------------------------------------------
        # Candidate caspase point from PPP(lambda_1)
        # ----------------------------------------------------

        x = rng.uniform(0.0, Lx)
        y = rng.uniform(0.0, Ly)

        inside = is_inside_active_zone(x, y)

        if inside:
            # Inside active zone: keep with probability 1
            accept = True
            event_name = "candidate inside active zone: accepted"

        else:
            # Outside active zone: keep with probability lambda_c/lambda_1
            u_accept = rng.uniform(0.0, 1.0)
            accept = u_accept <= (lambda_c / lambda_1)

            if accept:
                event_name = "candidate outside: accepted by thinning"
            else:
                event_name = "candidate outside: rejected"

        if accept:
            # Add accepted caspase event
            cas_x.append(x)
            cas_y.append(y)
            cas_t.append(t)

            # Positive feedback: accepted event creates new active center
            add_active_center(x, y)

    else:
        # ----------------------------------------------------
        # Expiration of one active center
        # ----------------------------------------------------

        idx = remove_random_active_center()

        if idx is None:
            event_name = "expiration impossible"
        else:
            event_name = f"expiration du centre actif {idx}"

    event_times.append(t)
    active_counts.append(len(active_x))
    cas_counts.append(len(cas_x))

    redraw(event_name)


# -----------------------------
# 6. End
# -----------------------------

redraw("fin de la simulation")

plt.ioff()
plt.show()