import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle


# ============================================================
# Simulation par REJECTION SAMPLING / THINNING
# Modèle spatio-temporel : activité caspase, mort cellulaire,
# et feedback ERK protecteur.
#
#   V^a_t : processus caché d'activation caspase
#   V^d_t : processus observé de morts cellulaires
#   V^p_t : processus de protection ERK généré après les morts
#
# Il simule d'abord des processus ponctuels de Poisson candidats,
# puis garde/rejette les candidats selon l'état juste avant l'événement.
# ============================================================


# -----------------------------
# 1. Paramètres du modèle
# -----------------------------

seed = 123
rng = np.random.default_rng(seed)

# Domaine spatial W = [0, Lx] x [0, Ly]
Lx = 10.0
Ly = 10.0
area_W = Lx * Ly

# Horizon de simulation.
# En simulation numérique, on prend toujours une fenêtre finie [0, T_max].
T_max = 20.0

# Paramètres de l'activation caspase, comme dans le PDF :
# lambda_a_1 : intensité dominante des candidats d'activation
# lambda_a_c : intensité de fond hors des zones actives
lambda_a_c = 0.05
lambda_a_1 = 0.60

if lambda_a_1 <= lambda_a_c:
    raise ValueError("Il faut lambda_a_1 > lambda_a_c pour la méthode de rejection.")

# Marques des centres d'activation :
# R^a ~ Exp(beta_a_R), T^a ~ Exp(beta_a_T)
beta_a_R = 1.4
beta_a_T = 0.8

if beta_a_R <= 0 or beta_a_T <= 0:
    raise ValueError("Il faut beta_a_R > 0 et beta_a_T > 0.")

# Paramètre du processus candidat de mort :
# lambda_d : intensité des candidats de mort sur W x [0,T_max]
lambda_d = 0.25

if lambda_d <= 0:
    raise ValueError("Il faut lambda_d > 0.")

# Marques associées à une mort acceptée, utilisées pour le feedback ERK :
# R^d ~ Exp(beta_d_R), T^d ~ Exp(beta_d_T)
beta_d_R = 1.2
beta_d_T = 0.6

if beta_d_R <= 0 or beta_d_T <= 0:
    raise ValueError("Il faut beta_d_R > 0 et beta_d_T > 0.")

# Paramètres d'affichage
n_frames = 300
show_history = False
show_rejected_death_candidates = False
show_erk_zones = True
interval_ms = 50


# -----------------------------
# 2. Candidats d'activation caspase : Phi^a
# -----------------------------

# PPP dominant de taux lambda_a_1 sur W x [0,T_max]
# N_star ~ Poisson(lambda_a_1 |W| T_max)
mean_activation_candidates = lambda_a_1 * area_W * T_max
N_star = rng.poisson(mean_activation_candidates)

print(f"Nombre total de candidats caspase simulés : N_star = {N_star}")

# Candidats : (Y_i^a, S_i^a, U_i^a, R_i^a, T_i^a)
X_star = rng.uniform(0.0, Lx, size=N_star)
Y_star = rng.uniform(0.0, Ly, size=N_star)
S_star = rng.uniform(0.0, T_max, size=N_star)
U_star = rng.uniform(0.0, 1.0, size=N_star)
R_star = rng.exponential(scale=1.0 / beta_a_R, size=N_star)
Tau_star = rng.exponential(scale=1.0 / beta_a_T, size=N_star)

# Tri chronologique : l'acceptation dépend de V^a_{s-}
order = np.argsort(S_star)
X_star = X_star[order]
Y_star = Y_star[order]
S_star = S_star[order]
U_star = U_star[order]
R_star = R_star[order]
Tau_star = Tau_star[order]


# -----------------------------
# 3. Rejection / thinning récursif pour V^a_t
# -----------------------------

X_acc = []      # positions x des activations acceptées
Y_acc = []      # positions y des activations acceptées
S_acc = []      # temps S_i^a
R_acc = []      # rayons R_i^a
Tau_acc = []    # durées T_i^a


def point_inside_current_active_zone(x, y, s):
    """
    Teste si (x,y) appartient à A(V^a_{s-}).

    Pour un centre déjà accepté k, actif sur [S_k, S_k + T_k),
    le centre est présent dans l'état à gauche V^a_{s-} si :
        S_k < s <= S_k + T_k.

    La condition stricte S_k < s évite qu'un candidat soit accepté
    à cause de lui-même.
    """
    for xk, yk, sk, rk, tauk in zip(X_acc, Y_acc, S_acc, R_acc, Tau_acc):
        if sk < s <= sk + tauk:
            if (x - xk) ** 2 + (y - yk) ** 2 <= rk ** 2:
                return True
    return False


for x, y, s, u, r, tau in zip(X_star, Y_star, S_star, U_star, R_star, Tau_star):
    inside_active_zone = point_inside_current_active_zone(x, y, s)

    # PDF : p^a_v(x) = 1_{x in A(v)} + (lambda_a_c/lambda_a_1) 1_{x notin A(v)}
    if inside_active_zone:
        accept = True
    else:
        accept = (u <= lambda_a_c / lambda_a_1)

    if accept:
        X_acc.append(x)
        Y_acc.append(y)
        S_acc.append(s)
        R_acc.append(r)
        Tau_acc.append(tau)

X_acc = np.array(X_acc)
Y_acc = np.array(Y_acc)
S_acc = np.array(S_acc)
R_acc = np.array(R_acc)
Tau_acc = np.array(Tau_acc)
activation_end_time = S_acc + Tau_acc

print(f"Nombre de centres de caspase acceptés : N_acc = {len(X_acc)}")


# -----------------------------
# 4. Fonctions associées à V^a_t et A(V^a_{s-})
# -----------------------------


def active_mask(t):
    """Centres de caspase actifs dans V^a_t : S_i^a <= t < S_i^a + T_i^a."""
    return (S_acc <= t) & (t < activation_end_time)


def appeared_mask(t):
    """Centres de caspase déjà apparus avant t."""
    return S_acc <= t


def point_inside_active_zone_from_arrays(x, y, s, X_centers, Y_centers, S_centers, R_centers, Tau_centers):
    """Teste si (x,y) appartient à A(V^a_{s-})."""
    if len(X_centers) == 0:
        return False

    active_before_s = (S_centers < s) & (s <= S_centers + Tau_centers)

    if not np.any(active_before_s):
        return False

    dx = x - X_centers[active_before_s]
    dy = y - Y_centers[active_before_s]
    rr = R_centers[active_before_s]

    return np.any(dx ** 2 + dy ** 2 <= rr ** 2)


# -----------------------------
# 5. Candidats de mort : Phi^d, puis V^d_t et V^p_t
# -----------------------------

# PPP candidat de morts de taux lambda_d sur W x [0,T_max]
# Chaque candidat porte déjà ses marques (R_i^d, T_i^d), comme dans le PDF.
mean_death_candidates = lambda_d * area_W * T_max
N_death_star = rng.poisson(mean_death_candidates)

print(f"Nombre total de candidats de mort simulés : N_death_star = {N_death_star}")

X_death_star = rng.uniform(0.0, Lx, size=N_death_star)
Y_death_star = rng.uniform(0.0, Ly, size=N_death_star)
S_death_star = rng.uniform(0.0, T_max, size=N_death_star)
R_death_star = rng.exponential(scale=1.0 / beta_d_R, size=N_death_star)
Tau_death_star = rng.exponential(scale=1.0 / beta_d_T, size=N_death_star)

# Tri chronologique : une mort acceptée crée une zone ERK qui peut bloquer les suivantes.
death_order = np.argsort(S_death_star)
X_death_star = X_death_star[death_order]
Y_death_star = Y_death_star[death_order]
S_death_star = S_death_star[death_order]
R_death_star = R_death_star[death_order]
Tau_death_star = Tau_death_star[death_order]

# Morts acceptées = processus observé V^d_t
X_death = []
Y_death = []
S_death = []

# Marques de protection associées aux morts acceptées = processus V^p_t
R_death = []
Tau_death = []

# Candidats rejetés, pour diagnostic/visualisation
X_death_rej = []
Y_death_rej = []
S_death_rej = []
reason_death_rej = []


def point_inside_protected_zone(x, y, s):
    """
    Teste si (x,y) appartient à A(V^p_{s-}).

    Une mort acceptée j au temps S_j^d crée une protection active dans le PDF pour :
        S_j^d < t <= S_j^d + T_j^d.
    Donc, pour tester l'état juste avant s, on utilise :
        S_j^d < s <= S_j^d + T_j^d.
    """
    for xj, yj, sj, rj, tauj in zip(X_death, Y_death, S_death, R_death, Tau_death):
        if sj < s <= sj + tauj:
            if (x - xj) ** 2 + (y - yj) ** 2 <= rj ** 2:
                return True
    return False


for xd, yd, sd, rd, taud in zip(
    X_death_star, Y_death_star, S_death_star, R_death_star, Tau_death_star
):
    inside_active_zone = point_inside_active_zone_from_arrays(
        xd, yd, sd, X_acc, Y_acc, S_acc, R_acc, Tau_acc
    )
    inside_protected_zone = point_inside_protected_zone(xd, yd, sd)

    # PDF : V^d garde les candidats qui vérifient :
    #   Y_i^d in A(V^a_{S_i^d-}) et Y_i^d notin A(V^p_{S_i^d-})
    if inside_active_zone and not inside_protected_zone:
        X_death.append(xd)
        Y_death.append(yd)
        S_death.append(sd)
        R_death.append(rd)
        Tau_death.append(taud)
    else:
        X_death_rej.append(xd)
        Y_death_rej.append(yd)
        S_death_rej.append(sd)
        if inside_protected_zone:
            reason_death_rej.append("blocked_by_ERK")
        else:
            reason_death_rej.append("outside_active_zone")

X_death = np.array(X_death)
Y_death = np.array(Y_death)
S_death = np.array(S_death)
R_death = np.array(R_death)
Tau_death = np.array(Tau_death)

X_death_rej = np.array(X_death_rej)
Y_death_rej = np.array(Y_death_rej)
S_death_rej = np.array(S_death_rej)
reason_death_rej = np.array(reason_death_rej, dtype=object)

N_rej_erk = np.sum(reason_death_rej == "blocked_by_ERK")
N_rej_outside = np.sum(reason_death_rej == "outside_active_zone")

print(f"Nombre de vrais événements de mort : N_death = {len(X_death)}")
print(f"Nombre de zones ERK créées : N_ERK = {len(R_death)}")
print(f"Nombre de candidats de mort rejetés : N_death_rej = {len(X_death_rej)}")
print(f"  - rejetés hors zone active : {N_rej_outside}")
print(f"  - rejetés par feedback ERK : {N_rej_erk}")


# -----------------------------
# 6. Masques temporels pour l'animation
# -----------------------------


def death_mask(t):
    """Morts observées jusqu'au temps t : S_i^d <= t."""
    return S_death <= t


def rejected_death_mask(t):
    """Candidats de mort rejetés jusqu'au temps t."""
    return S_death_rej <= t


def protection_active_mask(t):
    """Zones ERK actives dans V^p_t, avec la convention du PDF : S_i^d < t <= S_i^d + T_i^d."""
    return (S_death < t) & (t <= S_death + Tau_death)


# -----------------------------
# 7. Création de la figure
# -----------------------------

fig, ax = plt.subplots(figsize=(9.5, 7))

ax.set_xlim(0, Lx)
ax.set_ylim(0, Ly)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(r"Rejection sampling : $V^a_t$, $V^d_t$ et feedback ERK $V^p_t$")

history_scatter = ax.scatter([], [], s=15, alpha=0.25, label="centres caspase apparus non actifs")
active_scatter = ax.scatter([], [], s=45, label=r"centres actifs $V^a_t$")
death_scatter = ax.scatter([], [], s=60, marker="x", linewidths=1.6, label=r"morts observées $V^d_t$")
rejected_death_scatter = ax.scatter([], [], s=20, marker=".", alpha=0.15, label="candidats mort rejetés")
erk_center_scatter = ax.scatter([], [], s=30, marker="+", linewidths=1.2, label=r"centres ERK actifs $V^p_t$")

time_text = ax.text(0.02, 0.98, "", transform=ax.transAxes, va="top", ha="left")
count_text = ax.text(0.02, 0.93, "", transform=ax.transAxes, va="top", ha="left")
info_text = ax.text(0.02, 0.88, "", transform=ax.transAxes, va="top", ha="left")

ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, framealpha=0.90, borderaxespad=0.0)
fig.subplots_adjust(right=0.72, top=0.90)

active_circles = []
protection_circles = []


# -----------------------------
# 8. Animation
# -----------------------------


def init():
    history_scatter.set_offsets(np.empty((0, 2)))
    active_scatter.set_offsets(np.empty((0, 2)))
    death_scatter.set_offsets(np.empty((0, 2)))
    rejected_death_scatter.set_offsets(np.empty((0, 2)))
    erk_center_scatter.set_offsets(np.empty((0, 2)))
    time_text.set_text("")
    count_text.set_text("")
    info_text.set_text("")
    return (
        history_scatter,
        active_scatter,
        death_scatter,
        rejected_death_scatter,
        erk_center_scatter,
        time_text,
        count_text,
        info_text,
    )


def update(frame):
    global active_circles, protection_circles

    t = frame * T_max / (n_frames - 1)

    for c in active_circles:
        c.remove()
    active_circles = []

    for c in protection_circles:
        c.remove()
    protection_circles = []

    is_active = active_mask(t)
    has_appeared = appeared_mask(t)
    inactive_appeared = has_appeared & (~is_active)

    if show_history:
        history_points = np.column_stack([X_acc[inactive_appeared], Y_acc[inactive_appeared]])
        history_scatter.set_offsets(history_points)
    else:
        history_scatter.set_offsets(np.empty((0, 2)))

    active_points = np.column_stack([X_acc[is_active], Y_acc[is_active]])
    active_scatter.set_offsets(active_points)

    for xk, yk, rk in zip(X_acc[is_active], Y_acc[is_active], R_acc[is_active]):
        circle = Circle((xk, yk), rk, fill=False, alpha=0.45, linewidth=1.5)
        ax.add_patch(circle)
        active_circles.append(circle)

    died_until_t = death_mask(t)
    death_points = np.column_stack([X_death[died_until_t], Y_death[died_until_t]])
    death_scatter.set_offsets(death_points)

    is_protection_active = protection_active_mask(t)
    if show_erk_zones:
        protection_points = np.column_stack([X_death[is_protection_active], Y_death[is_protection_active]])
        erk_center_scatter.set_offsets(protection_points)

        for xj, yj, rj in zip(X_death[is_protection_active], Y_death[is_protection_active], R_death[is_protection_active]):
            circle = Circle((xj, yj), rj, fill=False, alpha=0.35, linewidth=1.4, linestyle="--")
            ax.add_patch(circle)
            protection_circles.append(circle)
    else:
        erk_center_scatter.set_offsets(np.empty((0, 2)))

    if show_rejected_death_candidates:
        rejected_until_t = rejected_death_mask(t)
        rejected_points = np.column_stack([X_death_rej[rejected_until_t], Y_death_rej[rejected_until_t]])
        rejected_death_scatter.set_offsets(rejected_points)
    else:
        rejected_death_scatter.set_offsets(np.empty((0, 2)))

    time_text.set_text(f"t = {t:.2f}")
    count_text.set_text(
        f"N_active = {np.sum(is_active)} | "
        f"N_death(t) = {np.sum(died_until_t)} | "
        f"N_ERK_active = {np.sum(is_protection_active)}"
    )
    info_text.set_text(
        rf"$\lambda_\nu(x)=\lambda^a_1\mathbf{{1}}_{{x\in A(\nu)}}"
        rf"+\lambda^a_c\mathbf{{1}}_{{x\notin A(\nu)}}$, "
        rf"$\lambda^d={lambda_d}$"
    )

    return (
        history_scatter,
        active_scatter,
        death_scatter,
        rejected_death_scatter,
        erk_center_scatter,
        time_text,
        count_text,
        info_text,
        *active_circles,
        *protection_circles,
    )


anim = FuncAnimation(
    fig,
    update,
    frames=n_frames,
    init_func=init,
    interval=interval_ms,
    blit=False,
)

plt.show()
