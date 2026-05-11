import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle


# ============================================================
# Simulation de V_t avec positive feedback
# + vrai processus de mort cellulaire Y_t
#
# 1) Processus des centres actifs de caspase V_t
#
# Intensité cible des naissances de centres de caspase :
#   lambda_nu(y) = lambda_c + (lambda_1 - lambda_c) 1_{y in A(nu)}
#
# Active zone :
#   A(nu) = W ∩ union_i B(Y_i, R_i)
#
# Construction par thinning :
#   - on simule un PPP dominant de taux lambda_1 sur W x [0,T]
#   - on accepte tous les candidats dans A(V_{s-})
#   - hors de A(V_{s-}), on accepte avec probabilité lambda_c/lambda_1
#
# 2) Processus de mort cellulaire Y_t
#
# Intensité conditionnelle de mort :
#   mu_nu(y) = mu_d 1_{y in A(nu)}
#
# Donc :
#   - dans la zone active : death rate = mu_d
#   - hors de la zone active : death rate = 0
#
# Construction :
#   - on simule un PPP candidat de morts de taux mu_d sur W x [0,T]
#   - on conserve seulement les candidats tombant dans A(V_{s-})
#
# Résultat :
#   V_t = sum_{k : S_k <= t < S_k + tau_k} delta_(Y_k, R_k)
#   Y_t = sum_{j : S_j^d <= t, X_j^d in A(V_{S_j^d-})} delta_{X_j^d}
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

# Intensité de fond des centres de caspase hors zones actives
lambda_c = 0.05

# Intensité plus forte des centres de caspase dans les zones actives
lambda_1 = 0.60

if lambda_1 <= lambda_c:
    raise ValueError("Il faut lambda_1 > lambda_c pour la méthode de rejection.")

# Marques des centres de caspase :
# R_k ~ Exp(beta_R)
# tau_k ~ Exp(beta_tau)
beta_R = 1.4
beta_tau = 0.8

# Nouvelle intensité de mort cellulaire dans la zone active
# mu_d = nombre moyen de morts par unité d'aire et par unité de temps
# dans A(V_t). En dehors de A(V_t), l'intensité est 0.
mu_d = 0.25

if mu_d <= 0:
    raise ValueError("Il faut mu_d > 0.")

# Nombre de frames pour l'animation
n_frames = 300

# Affichage de l'historique :
# True  = montre les centres déjà apparus mais non actifs
# False = montre seulement V_t
show_history = False

# Affichage des candidats de mort rejetés ou non
show_rejected_death_candidates = False

# Vitesse de l'animation
interval_ms = 50


# -----------------------------
# 2. Simuler le PPP dominant des centres de caspase
# -----------------------------

# PPP dominant de taux lambda_1 sur W x [0,T]
# Nombre de candidats :
# N_star ~ Poisson(lambda_1 |W| T)
mean_candidates = lambda_1 * area_W * T
N_star = rng.poisson(mean_candidates)

print(f"Nombre total de candidats caspase simulés : N_star = {N_star}")

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
# 3. Rejection / thinning récursif pour les centres de caspase
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


# Traitement chronologique des candidats caspase
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

print(f"Nombre de centres de caspase acceptés : N_acc = {len(X_acc)}")


# -----------------------------
# 4. Fonctions pour construire V_t et A(V_{s-})
# -----------------------------

def active_mask(t):
    """
    Centres de caspase actifs au temps t.

    Condition :
        S_k <= t < S_k + tau_k
    """
    return (S_acc <= t) & (t < death_time_acc)


def appeared_mask(t):
    """
    Centres de caspase déjà apparus avant t.

    Condition :
        S_k <= t
    """
    return S_acc <= t


def point_inside_active_zone_from_arrays(x, y, s, X_centers, Y_centers, S_centers, R_centers, tau_centers):
    """
    Teste si (x,y) est dans A(V_{s-}) pour des centres déjà simulés.

    Pour V_{s-}, la condition d'activité d'un centre k est :
        S_k < s <= S_k + tau_k.

    Cette fonction est utilisée pour filtrer les candidats de mort cellulaire.
    """
    if len(X_centers) == 0:
        return False

    active_before_s = (S_centers < s) & (s <= S_centers + tau_centers)

    if not np.any(active_before_s):
        return False

    dx = x - X_centers[active_before_s]
    dy = y - Y_centers[active_before_s]
    rr = R_centers[active_before_s]

    return np.any(dx ** 2 + dy ** 2 <= rr ** 2)


# -----------------------------
# 5. Vrai processus de mort cellulaire Y_t
# -----------------------------

# On simule un PPP candidat de morts de taux mu_d sur W x [0,T].
# Puis on garde uniquement les candidats qui tombent dans A(V_{s-}).
# Cela correspond à l'intensité réelle :
#     mu_{V_{s-}}(y) = mu_d 1_{y in A(V_{s-})}.

mean_death_candidates = mu_d * area_W * T
N_death_star = rng.poisson(mean_death_candidates)

print(f"Nombre total de candidats de mort simulés : N_death_star = {N_death_star}")

X_death_star = rng.uniform(0.0, Lx, size=N_death_star)
Y_death_star = rng.uniform(0.0, Ly, size=N_death_star)
S_death_star = rng.uniform(0.0, T, size=N_death_star)

# Trier par temps croissant : utile pour lire le processus chronologiquement.
# Ici, les morts ne modifient pas encore V_t, donc ce tri n'est pas mathématiquement
# indispensable, mais il rend le code plus clair et prépare l'ajout futur du feedback ERK.
death_order = np.argsort(S_death_star)
X_death_star = X_death_star[death_order]
Y_death_star = Y_death_star[death_order]
S_death_star = S_death_star[death_order]

# Candidats acceptés = vrais événements de mort cellulaire
X_death = []
Y_death = []
S_death = []

# Candidats rejetés = hors de la zone active au moment où ils apparaissent
X_death_rej = []
Y_death_rej = []
S_death_rej = []

for xd, yd, sd in zip(X_death_star, Y_death_star, S_death_star):
    inside_active_zone = point_inside_active_zone_from_arrays(
        xd,
        yd,
        sd,
        X_acc,
        Y_acc,
        S_acc,
        R_acc,
        tau_acc,
    )

    if inside_active_zone:
        # Vrai événement de mort : on ajoute delta_{X_d} au processus Y_t.
        X_death.append(xd)
        Y_death.append(yd)
        S_death.append(sd)
    else:
        # Hors de la zone active : death rate = 0, donc on rejette.
        X_death_rej.append(xd)
        Y_death_rej.append(yd)
        S_death_rej.append(sd)

X_death = np.array(X_death)
Y_death = np.array(Y_death)
S_death = np.array(S_death)

X_death_rej = np.array(X_death_rej)
Y_death_rej = np.array(Y_death_rej)
S_death_rej = np.array(S_death_rej)

print(f"Nombre de vrais événements de mort : N_death = {len(X_death)}")
print(f"Nombre de candidats de mort rejetés : N_death_rej = {len(X_death_rej)}")


def death_mask(t):
    """
    Morts cellulaires observées jusqu'au temps t.

    Y_t = sum_{j : S_j^d <= t} delta_{X_j^d}
    parmi les candidats qui ont été acceptés car dans A(V_{S_j^d-}).
    """
    return S_death <= t


def rejected_death_mask(t):
    """
    Candidats de mort rejetés jusqu'au temps t.
    Optionnel, seulement pour visualisation.
    """
    return S_death_rej <= t


# -----------------------------
# 6. Création de la figure
# -----------------------------

fig, ax = plt.subplots(figsize=(7, 7))

ax.set_xlim(0, Lx)
ax.set_ylim(0, Ly)
ax.set_aspect("equal", adjustable="box")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(r"Simulation de $V_t$ avec positive feedback et mort cellulaire")

# Centres déjà apparus mais non actifs
history_scatter = ax.scatter(
    [],
    [],
    s=15,
    alpha=0.25,
    label="centres caspase apparus non actifs",
)

# Centres actifs : support de V_t
active_scatter = ax.scatter(
    [],
    [],
    s=45,
    label=r"centres actifs $V_t$",
)

# Vrais événements de mort cellulaire : support de Y_t
death_scatter = ax.scatter(
    [],
    [],
    s=60,
    marker="x",
    linewidths=1.6,
    label=r"morts cellulaires $Y_t$",
)

# Candidats de mort rejetés, optionnel
rejected_death_scatter = ax.scatter(
    [],
    [],
    s=20,
    marker=".",
    alpha=0.15,
    label="candidats mort rejetés",
)

time_text = ax.text(
    0.02,
    0.98,
    "",
    transform=ax.transAxes,
    va="top",
    ha="left",
)

count_text = ax.text(
    0.02,
    0.93,
    "",
    transform=ax.transAxes,
    va="top",
    ha="left",
)

info_text = ax.text(
    0.02,
    0.88,
    "",
    transform=ax.transAxes,
    va="top",
    ha="left",
)

ax.legend(loc="upper right")

active_circles = []


# -----------------------------
# 7. Initialisation de l'animation
# -----------------------------

def init():
    history_scatter.set_offsets(np.empty((0, 2)))
    active_scatter.set_offsets(np.empty((0, 2)))
    death_scatter.set_offsets(np.empty((0, 2)))
    rejected_death_scatter.set_offsets(np.empty((0, 2)))
    time_text.set_text("")
    count_text.set_text("")
    info_text.set_text("")
    return (
        history_scatter,
        active_scatter,
        death_scatter,
        rejected_death_scatter,
        time_text,
        count_text,
        info_text,
    )


# -----------------------------
# 8. Fonction de mise à jour
# -----------------------------

def update(frame):
    global active_circles

    # Temps courant
    t = frame * T / (n_frames - 1)

    # Supprimer les anciens cercles
    for c in active_circles:
        c.remove()
    active_circles = []

    # Centres de caspase actifs
    is_active = active_mask(t)
    has_appeared = appeared_mask(t)
    inactive_appeared = has_appeared & (~is_active)

    # Historique optionnel des centres de caspase apparus mais non actifs
    if show_history:
        history_points = np.column_stack(
            [X_acc[inactive_appeared], Y_acc[inactive_appeared]]
        )
        history_scatter.set_offsets(history_points)
    else:
        history_scatter.set_offsets(np.empty((0, 2)))

    # Centres actifs V_t
    active_points = np.column_stack([X_acc[is_active], Y_acc[is_active]])
    active_scatter.set_offsets(active_points)

    # Cercles actifs : active zone A(V_t)
    for xk, yk, rk in zip(X_acc[is_active], Y_acc[is_active], R_acc[is_active]):
        circle = Circle(
            (xk, yk),
            rk,
            fill=False,
            alpha=0.45,
            linewidth=1.5,
        )
        ax.add_patch(circle)
        active_circles.append(circle)

    # Vrais événements de mort cellulaire Y_t
    died_until_t = death_mask(t)
    death_points = np.column_stack([X_death[died_until_t], Y_death[died_until_t]])
    death_scatter.set_offsets(death_points)

    # Candidats de mort rejetés, optionnel
    if show_rejected_death_candidates:
        death_rej_until_t = rejected_death_mask(t)
        death_rej_points = np.column_stack(
            [X_death_rej[death_rej_until_t], Y_death_rej[death_rej_until_t]]
        )
        rejected_death_scatter.set_offsets(death_rej_points)
    else:
        rejected_death_scatter.set_offsets(np.empty((0, 2)))

    # Textes
    time_text.set_text(f"t = {t:.2f}")
    count_text.set_text(
        f"N_active = {np.sum(is_active)} | "
        f"N_death(t) = {np.sum(died_until_t)}"
    )
    info_text.set_text(
        rf"$\mu_{{\nu}}(y)=\mu_d\mathbf{{1}}_{{y\in A(\nu)}}$, "
        rf"$\mu_d={mu_d}$"
    )

    return (
        history_scatter,
        active_scatter,
        death_scatter,
        rejected_death_scatter,
        time_text,
        count_text,
        info_text,
        *active_circles,
    )


# -----------------------------
# 9. Animation
# -----------------------------

anim = FuncAnimation(
    fig,
    update,
    frames=n_frames,
    init_func=init,
    interval=interval_ms,
    blit=False,
)

plt.show()
