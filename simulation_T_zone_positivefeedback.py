import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle


# ============================================================
# Gillespie / event-driven simulation with positive feedback
# and a fixed T-shaped background zone.
#
# Spatial domain W = [0, Lx] x [0, Ly] is divided into thirds
# in both directions.
#
# Intensity rule:
#   1) inside active zones A_t:               lambda_cas = lambda_3
#   2) outside active zones but inside T:     lambda_cas = lambda_2
#   3) outside active zones and outside T:    lambda_cas = lambda_1
#
# with lambda_1 << lambda_2 << lambda_3.
#
# Rejection / thinning:
#   Simulate candidate caspase points on the whole W with rate lambda_3.
#   If candidate is inside active zone: keep it with probability 1.
#   Else if candidate is inside the fixed T zone: keep it with probability lambda_2/lambda_3.
#   Else: keep it with probability lambda_1/lambda_3.
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

# Intensities:
# lambda_1: lowest background intensity outside the fixed T zone and outside active zones
# lambda_2: medium background intensity inside the fixed T zone but outside active zones
# lambda_3: highest intensity inside active zones
lambda_1 = 0.50
lambda_2 = 5.00
lambda_3 = 50.00

if not (lambda_1 < lambda_2 < lambda_3):
    raise ValueError("We need lambda_1 < lambda_2 < lambda_3.")

# Radius distribution: R ~ Exp(beta_R)
# Larger beta_R means smaller average radius.
beta_R = 2.5

# Expiration rate of each active center
# Larger beta_tau means faster expiration.
beta_tau = 1.2

# Pause between visual updates
display_pause = 0.08

# Draw rejected candidates or not
show_rejected_points = False


# -----------------------------
# 2. Fixed T-shaped zone
# -----------------------------

# Divide W into thirds.
x1 = Lx / 3.0
x2 = 2.0 * Lx / 3.0
y1 = Ly / 3.0
y2 = 2.0 * Ly / 3.0

# In this code, the T-shaped zone is the union of:
#   - the middle vertical column:  x in [Lx/3, 2Lx/3], y in [0, Ly]
#   - the middle-left horizontal arm: x in [0, Lx/3], y in [Ly/3, 2Ly/3]
#
# This matches the drawing: the T region has intensity lambda_2,
# while the three outside pieces have intensity lambda_1.

def is_inside_T_zone(x, y):
    """
    Return True if point (x,y) is inside the fixed T-shaped region.
    """
    inside_middle_column = (x1 <= x <= x2) and (0.0 <= y <= Ly)
    inside_middle_left_arm = (0.0 <= x <= x1) and (y1 <= y <= y2)

    return inside_middle_column or inside_middle_left_arm


# -----------------------------
# 3. State variables
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
cas_origin = []

# Rejected candidate points, only used for optional display
rej_x = []
rej_y = []
rej_t = []

# Histories
event_times = [0.0]
active_counts = [0]
cas_counts = [0]


# -----------------------------
# 4. Helper functions
# -----------------------------

def is_inside_active_zone(x, y):
    """
    Return True if point (x,y) is in A_t = union_i B(Y_i,R_i).
    """
    for cx, cy, r in zip(active_x, active_y, active_r):
        if (x - cx) ** 2 + (y - cy) ** 2 <= r ** 2:
            return True
    return False


def local_intensity(x, y):
    """
    Return lambda_nu(x,y), following priority:
      active zone -> lambda_3,
      T zone      -> lambda_2,
      outside     -> lambda_1.
    """
    if is_inside_active_zone(x, y):
        return lambda_3
    if is_inside_T_zone(x, y):
        return lambda_2
    return lambda_1


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
# 5. Matplotlib real-time setup
# -----------------------------

plt.ion()

fig, ax = plt.subplots(figsize=(7, 7))

ax.set_xlim(0, Lx)
ax.set_ylim(0, Ly)
ax.set_aspect("equal", adjustable="box")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Gillespie simulation: T-shaped background + positive feedback")

# Draw the fixed T-shaped lambda_2 background zone.
# Patch 1: middle vertical column.
t_patch_vertical = Rectangle(
    (x1, 0.0),
    x2 - x1,
    Ly,
    fill=True,
    alpha=0.12,
    linewidth=0.0,
    label=r"fixed $T$ zone: $\lambda_2$"
)
ax.add_patch(t_patch_vertical)

# Patch 2: middle-left horizontal arm.
t_patch_arm = Rectangle(
    (0.0, y1),
    x1,
    y2 - y1,
    fill=True,
    alpha=0.12,
    linewidth=0.0
)
ax.add_patch(t_patch_arm)

# Draw the 3 x 3 grid used to define the fixed T zone.
for xx in [x1, x2]:
    ax.axvline(xx, linestyle="--", linewidth=0.8, alpha=0.35)
for yy in [y1, y2]:
    ax.axhline(yy, linestyle="--", linewidth=0.8, alpha=0.35)

# Add text labels for the fixed background intensities.
ax.text(0.5 * x1, 0.5 * (y2 + Ly), r"$\lambda_1$", ha="center", va="center", fontsize=13, alpha=0.75)
ax.text(0.5 * x1, 0.5 * y1, r"$\lambda_1$", ha="center", va="center", fontsize=13, alpha=0.75)
ax.text(0.5 * (x2 + Lx), 0.5 * Ly, r"$\lambda_1$", ha="center", va="center", fontsize=13, alpha=0.75)
ax.text(0.5 * (x1 + x2), 0.5 * Ly, r"$\lambda_2$", ha="center", va="center", fontsize=14, alpha=0.75)

# Accepted caspase events
cas_scatter = ax.scatter(
    [],
    [],
    s=18,
    label=r"$N_{\mathrm{cas}}$ accepted events"
)

# Rejected candidate points, optional
rejected_scatter = ax.scatter(
    [],
    [],
    s=10,
    marker="x",
    alpha=0.30,
    label="rejected candidates"
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

intensity_text = ax.text(
    0.02,
    0.78,
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
      - fixed T-shaped background zone
      - accepted caspase events
      - active centers
      - active zones
    """
    global active_circles

    # Remove old active-zone circles.
    # The fixed T-zone patches are not removed because they do not change.
    for circle in active_circles:
        circle.remove()
    active_circles = []

    # Draw accepted caspase events
    if len(cas_x) > 0:
        cas_points = np.column_stack([cas_x, cas_y])
        cas_scatter.set_offsets(cas_points)
    else:
        cas_scatter.set_offsets(np.empty((0, 2)))

    # Draw rejected candidates if requested
    if show_rejected_points and len(rej_x) > 0:
        rejected_points = np.column_stack([rej_x, rej_y])
        rejected_scatter.set_offsets(rejected_points)
    else:
        rejected_scatter.set_offsets(np.empty((0, 2)))

    # Draw active centers
    if len(active_x) > 0:
        active_points = np.column_stack([active_x, active_y])
        active_scatter.set_offsets(active_points)
    else:
        active_scatter.set_offsets(np.empty((0, 2)))

    # Draw current active zones only.
    # Any point inside one of these circles has intensity lambda_3.
    for x, y, r in zip(active_x, active_y, active_r):
        circle = Circle(
            (x, y),
            r,
            fill=False,
            alpha=0.45,
            linewidth=1.3
        )
        ax.add_patch(circle)
        active_circles.append(circle)

    n_t = len(active_x)

    # Proposal candidate rate from dominating PPP(lambda_3).
    # The real accepted birth intensity is obtained by thinning.
    candidate_rate = lambda_3 * area_W

    # Expiration rate of active centers.
    expiration_rate = n_t * beta_tau

    # Total rate used in the event-driven simulation.
    a0 = candidate_rate + expiration_rate

    time_text.set_text(f"t = {t:.3f}")
    event_text.set_text(f"dernier événement : {event_name}")
    count_text.set_text(
        f"N_t = {n_t} active centers | "
        f"N_cas = {len(cas_x)} accepted"
    )
    rate_text.set_text(
        f"proposal a0 = λ3|W| + N_tβτ = {a0:.3f}"
    )
    intensity_text.set_text(
        r"$\lambda_1 < \lambda_2 < \lambda_3$; "
        r"active zone has priority"
    )

    fig.canvas.draw_idle()
    plt.pause(display_pause)


# Initial display
redraw("initialisation")


# -----------------------------
# 6. Gillespie / event-driven loop
# -----------------------------

while t < T_max:

    n_t = len(active_x)
    
    # Candidate caspase event rate from dominating PPP(lambda_3)
    candidate_rate = lambda_3 * area_W

    # Expiration rate of active centers
    expiration_rate = n_t * beta_tau

    # Total rate for proposal events and expiration events
    a0 = candidate_rate + expiration_rate

    if a0 <= 0:
        break

    # Waiting time until next proposal/expiration event
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
        # Candidate caspase point from dominating PPP(lambda_3)
        # ----------------------------------------------------

        x = rng.uniform(0.0, Lx)
        y = rng.uniform(0.0, Ly)

        inside_active = is_inside_active_zone(x, y)
        inside_T = is_inside_T_zone(x, y)

        if inside_active:
            # Inside active zone: highest intensity lambda_3
            # Since the proposal rate is lambda_3, keep with probability 1.
            accept = True
            origin = "active zone"
            event_name = "candidate inside active zone: accepted with λ3"

        elif inside_T:
            # Outside active zone but inside fixed T zone: intensity lambda_2
            # Keep with probability lambda_2/lambda_3.
            u_accept = rng.uniform(0.0, 1.0)
            accept = u_accept <= (lambda_2 / lambda_3)
            origin = "T zone"

            if accept:
                event_name = "candidate inside T zone: accepted by λ2/λ3 thinning"
            else:
                event_name = "candidate inside T zone: rejected"

        else:
            # Outside active zone and outside fixed T zone: intensity lambda_1
            # Keep with probability lambda_1/lambda_3.
            u_accept = rng.uniform(0.0, 1.0)
            accept = u_accept <= (lambda_1 / lambda_3)
            origin = "outside T"

            if accept:
                event_name = "candidate outside T zone: accepted by λ1/λ3 thinning"
            else:
                event_name = "candidate outside T zone: rejected"

        if accept:
            # Add accepted caspase event
            cas_x.append(x)
            cas_y.append(y)
            cas_t.append(t)
            cas_origin.append(origin)

            # Positive feedback: accepted event creates new active center
            add_active_center(x, y)
        else:
            # Store rejected candidate only for optional display
            rej_x.append(x)
            rej_y.append(y)
            rej_t.append(t)

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
# 7. End
# -----------------------------

redraw("fin de la simulation")

plt.ioff()
plt.show()
