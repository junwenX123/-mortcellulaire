import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Circle


# ============================================================
# Simulation de V_t construit à partir de Phi_c tilde
# Sans ERK, sans mort cellulaire, sans interaction.
#
# Phi_c_tilde = sum_k delta_(Y_k, S_k, R_k, tau_k)
#
# V_t = sum_{k : S_k <= t < S_k + tau_k} delta_(Y_k, R_k)
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

# Horizon temporel [0, T]
T = 20.0

# Intensité espace-temps constante des centres
# lambda_c = nombre moyen de centres par unité d'aire et par unité de temps
lambda_c = 0.08

# Marques :
# R_k ~ Exp(beta_R)
# tau_k ~ Exp(beta_tau)
beta_R = 1.2
beta_tau = 0.35

# Nombre de frames pour l'animation
n_frames = 250

# Affichage de l'historique :
# True  = montre aussi les centres déjà apparus mais inactifs en gris
# False = montre seulement V_t
show_history = True

# Sauvegarde de l'animation
save_gif = True
gif_name = "vt_simulation.gif"


# -----------------------------
# 2. Simulation de Phi_c tilde
# -----------------------------

# Nombre total de centres dans W x [0,T]
# N ~ Poisson(lambda_c * |W| * T)
mean_number_centers = lambda_c * area_W * T
N = rng.poisson(mean_number_centers)

print(f"Nombre total de centres simulés : N = {N}")

# Positions spatiales Y_k = (X_k, Y_k)
X = rng.uniform(0, Lx, size=N)
Y = rng.uniform(0, Ly, size=N)

# Temps d'apparition S_k
S = rng.uniform(0, T, size=N)

# Rayons R_k
R = rng.exponential(scale=1.0 / beta_R, size=N)

# Durées d'activité tau_k
tau = rng.exponential(scale=1.0 / beta_tau, size=N)

# Temps de disparition S_k + tau_k
death_time = S + tau


# -----------------------------
# 3. Fonction qui construit V_t
# -----------------------------

def active_mask(t):
    """
    Retourne les indices des centres actifs au temps t.

    Condition :
        S_k <= t < S_k + tau_k
    """
    return (S <= t) & (t < death_time)


def appeared_mask(t):
    """
    Centres déjà apparus avant t.

    Condition :
        S_k <= t
    """
    return S <= t


# -----------------------------
# 4. Création de la figure
# -----------------------------

fig, ax = plt.subplots(figsize=(7, 7))

ax.set_xlim(0, Lx)
ax.set_ylim(0, Ly)
ax.set_aspect("equal", adjustable="box")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Simulation dynamique de $V_t$")

# Points historiques gris : centres apparus mais non actifs
history_scatter = ax.scatter([], [], s=15, alpha=0.25, label="centres apparus non actifs")

# Points actifs rouges : support de V_t
active_scatter = ax.scatter([], [], s=35, label="centres actifs $V_t$")

# Texte temps
time_text = ax.text(
    0.02,
    0.98,
    "",
    transform=ax.transAxes,
    va="top",
    ha="left"
)

# Texte nombre actifs
count_text = ax.text(
    0.02,
    0.93,
    "",
    transform=ax.transAxes,
    va="top",
    ha="left"
)

ax.legend(loc="upper right")

# Liste des cercles actifs dessinés à chaque frame
active_circles = []


# -----------------------------
# 5. Fonction d'initialisation
# -----------------------------

def init():
    history_scatter.set_offsets(np.empty((0, 2)))
    active_scatter.set_offsets(np.empty((0, 2)))
    time_text.set_text("")
    count_text.set_text("")
    return history_scatter, active_scatter, time_text, count_text


# -----------------------------
# 6. Fonction de mise à jour
# -----------------------------

def update(frame):
    global active_circles

    # Temps courant
    t = frame * T / (n_frames - 1)

    # Supprimer les anciens cercles actifs
    for c in active_circles:
        c.remove()
    active_circles = []

    # Masques
    is_active = active_mask(t)
    has_appeared = appeared_mask(t)

    # Centres apparus mais non actifs
    inactive_appeared = has_appeared & (~is_active)

    if show_history:
        history_points = np.column_stack([X[inactive_appeared], Y[inactive_appeared]])
        history_scatter.set_offsets(history_points)
    else:
        history_scatter.set_offsets(np.empty((0, 2)))

    # Centres actifs : support de V_t
    active_points = np.column_stack([X[is_active], Y[is_active]])
    active_scatter.set_offsets(active_points)

    # Dessiner les rayons R_k des centres actifs
    for xk, yk, rk in zip(X[is_active], Y[is_active], R[is_active]):
        circle = Circle(
            (xk, yk),
            rk,
            fill=False,
            alpha=0.45,
            linewidth=1.5
        )
        ax.add_patch(circle)
        active_circles.append(circle)

    # Textes
    time_text.set_text(f"t = {t:.2f}")
    count_text.set_text(f"Nombre de centres actifs = {np.sum(is_active)}")

    return history_scatter, active_scatter, time_text, count_text, *active_circles


# -----------------------------
# 7. Animation
# -----------------------------

anim = FuncAnimation(
    fig,
    update,
    frames=n_frames,
    init_func=init,
    interval=60,
    blit=False
)

# Sauvegarde en GIF
if save_gif:
    writer = PillowWriter(fps=20)
    anim.save(gif_name, writer=writer)
    print(f"Animation sauvegardée dans : {gif_name}")

plt.show()