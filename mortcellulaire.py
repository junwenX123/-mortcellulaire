import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle


# ============================================================
# Gillespie / event-driven simulation with:
#   1. fixed T-shaped background zone,
#   2. positive feedback for caspase,
#   3. true cell death point process.
#
# Spatial domain W = [0, Lx] x [0, Ly].
#
# Caspase intensity rule:
#   1) inside active zones A_t:               lambda_cas = lambda_3
#   2) outside active zones but inside T:     lambda_cas = lambda_2
#   3) outside active zones and outside T:    lambda_cas = lambda_1
#
# Cell death intensity:
#   lambda_death(x,t | V_{t-}) = u_d * 1_{x in A(V_{t-})}
#
# Meaning:
#   - cell death happens only inside active zones;
#   - outside active zones, death rate is 0.
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

# Final time
T_max = 20.0

# Caspase intensities
lambda_1 = 0.05   # outside T zone and outside active zones
lambda_2 = 0.50   # inside fixed T zone, outside active zones
lambda_3 = 5.00  # inside active zones

if not (lambda_1 < lambda_2 < lambda_3):
    raise ValueError("We need lambda_1 < lambda_2 < lambda_3.")

# Cell death rate inside active zones
# Outside active zones, death rate is 0.
u_d = 1

# Radius distribution: R ~ Exp(beta_R)
beta_R = 2.5

# Expiration rate of each active center
beta_tau = 1.2

# Pause between visual updates
display_pause = 0.008

# Draw rejected candidates or not
show_rejected_points = False


# -----------------------------
# 2. Fixed T-shaped zone
# -----------------------------

x1 = Lx / 3.0
x2 = 2.0 * Lx / 3.0
y1 = Ly / 3.0
y2 = 2.0 * Ly / 3.0


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

# True cell death events D_t
death_x = []
death_y = []
death_t = []

# Rejected candidate points, only used for optional display
rej_x = []
rej_y = []
rej_t = []
rej_type = []

# Histories
event_times = [0.0]
active_counts = [0]
cas_counts = [0]
death_counts = [0]


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


def local_caspase_intensity(x, y):
    """
    Return lambda_cas(x,y), following priority:
      active zone -> lambda_3,
      T zone      -> lambda_2,
      outside     -> lambda_1.
    """
    if is_inside_active_zone(x, y):
        return lambda_3
    if is_inside_T_zone(x, y):
        return lambda_2
    return lambda_1


def local_death_intensity(x, y):
    """
    Return lambda_death(x,y | V_t).

    Death happens only inside active zones:
        lambda_death = u_d inside A_t,
        lambda_death = 0 outside A_t.
    """
    if is_inside_active_zone(x, y):
        return u_d
    return 0.0


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
ax.set_title("Gillespie simulation: caspase + cell death")

# Draw fixed T-shaped lambda_2 background zone.
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

t_patch_arm = Rectangle(
    (0.0, y1),
    x1,
    y2 - y1,
    fill=True,
    alpha=0.12,
    linewidth=0.0
)
ax.add_patch(t_patch_arm)

# Draw 3 x 3 grid.
for xx in [x1, x2]:
    ax.axvline(xx, linestyle="--", linewidth=0.8, alpha=0.35)
for yy in [y1, y2]:
    ax.axhline(yy, linestyle="--", linewidth=0.8, alpha=0.35)

# Background intensity labels.
ax.text(0.5 * x1, 0.5 * (y2 + Ly), r"$\lambda_1$", ha="center", va="center", fontsize=13, alpha=0.75)
ax.text(0.5 * x1, 0.5 * y1, r"$\lambda_1$", ha="center", va="center", fontsize=13, alpha=0.75)
ax.text(0.5 * (x2 + Lx), 0.5 * Ly, r"$\lambda_1$", ha="center", va="center", fontsize=13, alpha=0.75)
ax.text(0.5 * (x1 + x2), 0.5 * Ly, r"$\lambda_2$", ha="center", va="center", fontsize=14, alpha=0.75)

# Accepted caspase events
cas_scatter = ax.scatter(
    [],
    [],
    s=18,
    label=r"$N_{\mathrm{cas}}$ accepted caspase events"
)

# True cell death events
death_scatter = ax.scatter(
    [],
    [],
    s=55,
    marker="x",
    linewidths=2.0,
    label=r"$D_t$ cell death events"
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

death_text = ax.text(
    0.02,
    0.78,
    "",
    transform=ax.transAxes,
    va="top",
    ha="left"
)

intensity_text = ax.text(
    0.02,
    0.73,
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
      - true cell death events
      - active centers
      - active zones
    """
    global active_circles

    # Remove old active-zone circles.
    for circle in active_circles:
        circle.remove()
    active_circles = []

    # Draw accepted caspase events
    if len(cas_x) > 0:
        cas_points = np.column_stack([cas_x, cas_y])
        cas_scatter.set_offsets(cas_points)
    else:
        cas_scatter.set_offsets(np.empty((0, 2)))

    # Draw true cell death events
    if len(death_x) > 0:
        death_points = np.column_stack([death_x, death_y])
        death_scatter.set_offsets(death_points)
    else:
        death_scatter.set_offsets(np.empty((0, 2)))

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

    # Draw active zones.
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

    # Proposal rates:
    #   caspase candidates are proposed with rate lambda_3 |W|
    #   death candidates are proposed with rate u_d |W|
    #   active centers expire with rate N_t beta_tau
    cas_candidate_rate = lambda_3 * area_W
    death_candidate_rate = u_d * area_W
    expiration_rate = n_t * beta_tau

    a0 = cas_candidate_rate + death_candidate_rate + expiration_rate

    time_text.set_text(f"t = {t:.3f}")
    event_text.set_text(f"dernier événement : {event_name}")
    count_text.set_text(
        f"N_t = {n_t} active centers | "
        f"N_cas = {len(cas_x)}"
    )
    rate_text.set_text(
        f"proposal a0 = λ3|W| + u_d|W| + N_tβτ = {a0:.3f}"
    )
    death_text.set_text(
        f"N_death = {len(death_x)} | death only inside active zones"
    )
    intensity_text.set_text(
        r"$\lambda_d(x,t)=u_d\,\mathbf{1}_{x\in A(V_{t-})}$"
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

    # Proposal rate for caspase candidates from dominating PPP(lambda_3)
    cas_candidate_rate = lambda_3 * area_W

    # Proposal rate for cell death candidates from dominating PPP(u_d)
    # Then we reject outside active zones.
    death_candidate_rate = u_d * area_W

    # Expiration rate of active centers
    expiration_rate = n_t * beta_tau

    # Total proposal rate
    a0 = cas_candidate_rate + death_candidate_rate + expiration_rate

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

    p_cas_candidate = cas_candidate_rate / a0
    p_death_candidate = death_candidate_rate / a0

    if u <= p_cas_candidate:
        # ----------------------------------------------------
        # Event type 1:
        # Candidate caspase point from dominating PPP(lambda_3)
        # ----------------------------------------------------

        x = rng.uniform(0.0, Lx)
        y = rng.uniform(0.0, Ly)

        inside_active = is_inside_active_zone(x, y)
        inside_T = is_inside_T_zone(x, y)

        if inside_active:
            # Inside active zone: intensity lambda_3
            accept = True
            origin = "active zone"
            event_name = "caspase candidate inside active zone: accepted with λ3"

        elif inside_T:
            # Inside T zone but outside active zone: intensity lambda_2
            u_accept = rng.uniform(0.0, 1.0)
            accept = u_accept <= (lambda_2 / lambda_3)
            origin = "T zone"

            if accept:
                event_name = "caspase candidate inside T zone: accepted by λ2/λ3"
            else:
                event_name = "caspase candidate inside T zone: rejected"

        else:
            # Outside T zone and outside active zone: intensity lambda_1
            u_accept = rng.uniform(0.0, 1.0)
            accept = u_accept <= (lambda_1 / lambda_3)
            origin = "outside T"

            if accept:
                event_name = "caspase candidate outside T: accepted by λ1/λ3"
            else:
                event_name = "caspase candidate outside T: rejected"

        if accept:
            # Add accepted caspase event
            cas_x.append(x)
            cas_y.append(y)
            cas_t.append(t)
            cas_origin.append(origin)

            # Positive feedback:
            # accepted caspase event creates new active center
            add_active_center(x, y)

        else:
            # Store rejected candidate only for optional display
            rej_x.append(x)
            rej_y.append(y)
            rej_t.append(t)
            rej_type.append("caspase")

    elif u <= p_cas_candidate + p_death_candidate:
        # ----------------------------------------------------
        # Event type 2:
        # Candidate cell death point from PPP(u_d) on W
        #
        # True target death intensity:
        #   lambda_death(x,t) = u_d 1_{x in A(V_{t-})}
        #
        # Therefore:
        #   propose uniformly in W,
        #   accept only if inside active zone.
        # ----------------------------------------------------

        x = rng.uniform(0.0, Lx)
        y = rng.uniform(0.0, Ly)

        if is_inside_active_zone(x, y):
            # True cell death event
            death_x.append(x)
            death_y.append(y)
            death_t.append(t)

            event_name = "mort cellulaire inside active zone: accepted"

        else:
            # Outside active zone:
            # death rate is 0, so reject.
            rej_x.append(x)
            rej_y.append(y)
            rej_t.append(t)
            rej_type.append("death")

            event_name = "death candidate outside active zone: rejected"

    else:
        # ----------------------------------------------------
        # Event type 3:
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
    death_counts.append(len(death_x))

    redraw(event_name)


# -----------------------------
# 7. End
# -----------------------------

redraw("fin de la simulation")

plt.ioff()
plt.show()