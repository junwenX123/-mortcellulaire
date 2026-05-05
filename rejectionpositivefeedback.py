import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle


# ============================================================
# Simulation de V_t avec positive feedback
# Méthode de rejection / thinning à partir d'un PPP dominant
#
# Intensité cible :
#   lambda_nu(y) = lambda_c + (lambda_1 - lambda_c) 1_{y in A(nu)}
#
# Active zone :
#   A(nu) = W ∩ union_i B(Y_i, R_i)
#
# Construction :
#   1. Simuler un PPP dominant de taux lambda_1 sur W x [0,T].
#   2. Traiter les points candidats par ordre temporel.
#   3. Si le point est dans A(V_{s-}), on l'accepte.
#   4. Sinon, on l'accepte avec probabilité lambda_c / lambda_1.
#   5. Chaque point accepté devient un nouveau centre actif.
#
# Ensuite :
#   V_t = sum_{k : S_k <= t < S_k + tau_k} delta_(Y_k, R_k)
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

# Horizon temporel [0,T]
T = 20.0

# Intensité de fond hors zones actives
lambda_c = 0.05

# Intensité plus forte dans les zones actives
lambda_1 = 0.60

if lambda_1 <= lambda_c:
    raise ValueError("Il faut lambda_1 > lambda_c pour la méthode de rejection.")

# Marques :
# R_k ~ Exp(beta_R)
# tau_k ~ Exp(beta_tau)
beta_R = 1.4
beta_tau = 0.8

# Nombre de frames pour l'animation
n_frames = 300

# Affichage de l'historique :
# True  = montre les centres déjà apparus mais non actifs
# False = montre seulement V_t
show_history = False

# Vitesse de l'animation
interval_ms = 50


# -----------------------------
# 2. Simuler le PPP dominant
# -----------------------------

# PPP dominant de taux lambda_1 sur W x [0,T]
# Nombre de candidats :
# N_star ~ Poisson(lambda_1 |W| T)
mean_candidates = lambda_1 * area_W * T
N_star = rng.poisson(mean_candidates)

print(f"Nombre total de candidats simulés : N_star = {N_star}")

# Positions candidates X_j = (X_star_j, Y_star_j)
X_star = rng.uniform(0.0, Lx, size=N_star)
Y_star = rng.uniform(0.0, Ly, size=N_star)

# Temps candidats S_j
S_star = rng.uniform(0.0, T, size=N_star)

# Uniformes pour rejection / thinning
U_star = rng.uniform(0.0, 1.0, size=N_star)

# Marques candidates
R_star = rng.exponential(scale=1.0 / beta_R, size=N_star)
tau_star = rng.exponential(scale=1.0 / beta_tau, size=N_star)

# On trie les candidats par temps croissant.
# C'est essentiel, car l'acceptation d'un candidat dépend de V_{s-}.
order = np.argsort(S_star)

X_star = X_star[order]
Y_star = Y_star[order]
S_star = S_star[order]
U_star = U_star[order]
R_star = R_star[order]
tau_star = tau_star[order]


# -----------------------------
# 3. Rejection / thinning récursif
# -----------------------------

# Listes des centres acceptés
X_acc = []
Y_acc = []
S_acc = []
R_acc = []
tau_acc = []


def point_inside_current_active_zone(x, y, s):
    """
    Teste si le point spatial (x,y) appartient à la zone active A(V_{s-}).

    V_{s-} contient les centres acceptés qui sont actifs juste avant s :
        S_k < s <= S_k + tau_k.

    On utilise strictement S_k < s pour éviter qu'un candidat soit accepté
    à cause de lui-même.
    """
    for xk, yk, sk, rk, tauk in zip(X_acc, Y_acc, S_acc, R_acc, tau_acc):
        if sk < s <= sk + tauk:
            if (x - xk) ** 2 + (y - yk) ** 2 <= rk ** 2:
                return True
    return False


# Traitement chronologique des candidats
for x, y, s, u, r, tau in zip(X_star, Y_star, S_star, U_star, R_star, tau_star):

    inside = point_inside_current_active_zone(x, y, s)

    if inside:
        # Dans la zone active : intensité cible = lambda_1
        # Donc probabilité d'acceptation = 1.
        accept = True
    else:
        # Hors zone active : intensité cible = lambda_c.
        # Comme le PPP dominant a taux lambda_1,
        # probabilité d'acceptation = lambda_c / lambda_1.
        accept = (u <= lambda_c / lambda_1)

    if accept:
        X_acc.append(x)
        Y_acc.append(y)
        S_acc.append(s)
        R_acc.append(r)
        tau_acc.append(tau)


# Convertir en arrays NumPy pour faciliter les masques
X_acc = np.array(X_acc)
Y_acc = np.array(Y_acc)
S_acc = np.array(S_acc)
R_acc = np.array(R_acc)
tau_acc = np.array(tau_acc)

death_time_acc = S_acc + tau_acc

print(f"Nombre de centres acceptés : N_acc = {len(X_acc)}")


# -----------------------------
# 4. Fonctions pour construire V_t
# -----------------------------

def active_mask(t):
    """
    Centres actifs au temps t.

    Condition :
        S_k <= t < S_k + tau_k
    """
    return (S_acc <= t) & (t < death_time_acc)


def appeared_mask(t):
    """
    Centres déjà apparus avant t.

    Condition :
        S_k <= t
    """
    return S_acc <= t


# -----------------------------
# 5. Création de la figure
# -----------------------------

fig, ax = plt.subplots(figsize=(7, 7))

ax.set_xlim(0, Lx)
ax.set_ylim(0, Ly)
ax.set_aspect("equal", adjustable="box")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(r"Simulation de $V_t$ avec positive feedback")

# Centres déjà apparus mais non actifs
history_scatter = ax.scatter(
    [],
    [],
    s=15,
    alpha=0.25,
    label="centres apparus non actifs"
)

# Centres actifs : support de V_t
active_scatter = ax.scatter(
    [],
    [],
    s=45,
    label=r"centres actifs $V_t$"
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

info_text = ax.text(
    0.02,
    0.88,
    "",
    transform=ax.transAxes,
    va="top",
    ha="left"
)

ax.legend(loc="upper right")

active_circles = []


# -----------------------------
# 6. Initialisation de l'animation
# -----------------------------

def init():
    history_scatter.set_offsets(np.empty((0, 2)))
    active_scatter.set_offsets(np.empty((0, 2)))
    time_text.set_text("")
    count_text.set_text("")
    info_text.set_text("")
    return history_scatter, active_scatter, time_text, count_text, info_text


# -----------------------------
# 7. Fonction de mise à jour
# -----------------------------

def update(frame):
    global active_circles

    # Temps courant
    t = frame * T / (n_frames - 1)

    # Supprimer les anciens cercles
    for c in active_circles:
        c.remove()
    active_circles = []

    # Masques
    is_active = active_mask(t)
    has_appeared = appeared_mask(t)
    inactive_appeared = has_appeared & (~is_active)

    # Historique optionnel
    if show_history:
        history_points = np.column_stack(
            [X_acc[inactive_appeared], Y_acc[inactive_appeared]]
        )
        history_scatter.set_offsets(history_points)
    else:
        history_scatter.set_offsets(np.empty((0, 2)))

    # Centres actifs
    active_points = np.column_stack([X_acc[is_active], Y_acc[is_active]])
    active_scatter.set_offsets(active_points)

    # Cercles actifs : active zone A(V_t)
    for xk, yk, rk in zip(X_acc[is_active], Y_acc[is_active], R_acc[is_active]):
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
    info_text.set_text(
        rf"$\lambda_\nu(y)=\lambda_c+(\lambda_1-\lambda_c)\mathbf{{1}}_{{y\in A(\nu)}}$"
    )

    return (
        history_scatter,
        active_scatter,
        time_text,
        count_text,
        info_text,
        *active_circles
    )


# -----------------------------
# 8. Animation
# -----------------------------

anim = FuncAnimation(
    fig,
    update,
    frames=n_frames,
    init_func=init,
    interval=interval_ms,
    blit=False
)

plt.show()