# Analyse de l'uniformité spatiale des morts observées après burn-in

## 1. Objectif de l'expérience

L'objectif de cette série de simulations est d'étudier si les positions des morts observées restent spatialement uniformes ou deviennent concentrées dans la zone en forme de T, appelée **T-zone**.

Pour éviter d'interpréter une phase transitoire du processus, on applique d'abord un **burn-in** de 500 morts observées. Les statistiques sont ensuite calculées uniquement sur les 1000 morts observées suivantes.

Pour chaque jeu de paramètres, plusieurs simulations indépendantes sont réalisées avec différents seeds. Cette procédure est nécessaire car le modèle est stochastique : une seule simulation peut donner une conclusion trompeuse à cause de la variabilité aléatoire.

Les principaux indicateurs utilisés sont :

* `spatial_chi2_p_value_chi2_bins` : p-value du test d'uniformité sur une grille 10 × 10 ;
* `spatial_chi2_per_df_chi2_bins` : statistique (\chi^2) normalisée par le nombre de degrés de liberté ;
* `spatial_coefficient_of_variation_plot_bins` : coefficient de variation spatial sur la grille 20 × 20 ;
* `t_zone_density_ratio` : rapport entre la densité de morts observées dans la T-zone et la densité hors T-zone ;
* `total_time_for_analyzed_1000_deaths_after_burnin` : temps nécessaire pour obtenir les 1000 morts analysées après burn-in.

---

## 2. Résultat de référence : le cas baseline est presque uniforme

Dans le cas baseline, les cinq réplications indépendantes donnent une p-value non significative. La fraction de simulations rejetant l'uniformité au seuil 5% est égale à 0.

Les valeurs moyennes principales sont approximativement :

```text
baseline:
    fraction p-value < 0.05 = 0.0
    mean chi2/df ≈ 1.105
    mean p-value ≈ 0.226
    mean T-zone density ratio ≈ 0.975
    mean time for 1000 deaths ≈ 24.15
```

Ainsi, sous les paramètres de référence, les morts observées après burn-in sont compatibles avec une distribution spatialement uniforme. En particulier, le ratio de densité dans la T-zone est proche de 1, ce qui indique qu'il n'y a pas de concentration nette des morts dans la T-zone.

**Conclusion baseline :**

```text
Under the baseline parameters, the observed death process remains approximately spatially uniform after burn-in.
```

---

## 3. Augmenter seulement le contraste lambda_a_T / lambda_a_c ne suffit pas

On a d'abord testé l'effet du contraste entre l'intensité d'activation dans la T-zone et l'intensité d'activation de fond hors T-zone :

[
\frac{\lambda_{a,T}}{\lambda_{a,c}}.
]

Les valeurs testées incluent des ratios de l'ordre de :

[
10,\ 20,\ 50,\ 100,\ 1000,\ 10000.
]

Cependant, les résultats montrent que l'augmentation de ce ratio ne suffit pas à produire une structure spatiale clairement non uniforme.

Exemples de résultats agrégés :

```text
baseline:
    ratio = 10
    fraction p-value < 0.05 = 0.0
    mean T-zone density ratio ≈ 0.975

T0p5_c0p025_ratio_20:
    ratio = 20
    fraction p-value < 0.05 = 0.0
    mean T-zone density ratio ≈ 0.963

T0p5_c0p01_ratio_50:
    ratio = 50
    fraction p-value < 0.05 = 0.0
    mean T-zone density ratio ≈ 0.968

T0p5_c0p005_ratio_1e2:
    ratio = 100
    fraction p-value < 0.05 = 0.0
    mean T-zone density ratio ≈ 0.993

T0p5_c0p0005_ratio_1e3:
    ratio = 1000
    fraction p-value < 0.05 = 0.4
    mean T-zone density ratio ≈ 0.920

T0p5_c0p00005_ratio_1e4:
    ratio = 10000
    fraction p-value < 0.05 = 0.4
    mean T-zone density ratio ≈ 0.920
```

Même lorsque le ratio atteint (10^3) ou (10^4), le ratio de densité dans la T-zone reste inférieur à 1 en moyenne. Cela signifie que la non-uniformité observée dans certaines réplications ne correspond pas à une concentration des morts dans la T-zone.

**Conclusion :**

```text
Increasing only the contrast lambda_a_T / lambda_a_c is not sufficient to produce a visible T-zone pattern under the baseline active-zone radius and lifetime.
```

Une explication possible est que la T-zone influence seulement la naissance des centres actifs. Les morts observées, elles, résultent ensuite d'un mécanisme indirect : elles doivent apparaître dans une zone active et hors des zones de protection ERK. Si les zones actives sont trop grandes ou trop persistantes, le signal spatial initial de la T-zone peut être diffusé et finalement masqué.

---

## 4. Le scénario larger_ERK_radius rejette l'uniformité, mais ce n'est pas un effet T-zone

Le scénario `larger_ERK_radius` donne un résultat particulier :

```text
larger_ERK_radius:
    fraction p-value < 0.05 = 1.0
    mean p-value ≈ 0.00242
    mean chi2/df ≈ 1.506
    mean T-zone density ratio ≈ 0.894
    mean time for 1000 deaths ≈ 55.50
```

Les cinq réplications rejettent donc l'uniformité spatiale. Cependant, le ratio de densité dans la T-zone est inférieur à 1.

Cela signifie que ce scénario produit une non-uniformité globale, mais pas une concentration des morts dans la T-zone.

L'interprétation la plus naturelle est que l'augmentation du rayon ERK renforce l'effet d'inhibition spatiale autour des morts déjà observées. Comme une mort candidate est acceptée seulement si elle est dans une zone active et hors d'une zone ERK, un rayon ERK plus grand peut créer une structure spatiale non uniforme par exclusion locale, sans révéler la forme de la T-zone.

**Conclusion pour `larger_ERK_radius` :**

```text
The larger_ERK_radius scenario rejects global spatial uniformity, but its T-zone density ratio is below 1. Therefore, this non-uniformity should not be interpreted as a T-zone enrichment effect. It is more likely caused by stronger ERK-mediated spatial inhibition.
```

---

## 5. Le vrai facteur déterminant : la localisation des zones actives

Les résultats les plus importants apparaissent quand on modifie les paramètres des zones actives :

[
\beta_{a,R}
\quad \text{et} \quad
\beta_{a,T}.
]

On a :

[
\mathbb{E}[R^a] = \frac{1}{\beta_{a,R}},
]

et

[
\mathbb{E}[T^a] = \frac{1}{\beta_{a,T}}.
]

Donc :

* lorsque (\beta_{a,R}) augmente, le rayon moyen des zones actives diminue ;
* lorsque (\beta_{a,T}) augmente, la durée de vie moyenne des centres actifs diminue ;
* les zones actives deviennent alors plus locales spatialement et temporellement.

Les scénarios `local_activation` montrent très clairement cet effet :

```text
T0p5_c0p005_ratio_1e2_local_activation:
    beta_a_R = 5.0
    beta_a_T = 3.0
    mean T-zone density ratio ≈ 11.93
    mean chi2/df ≈ 9.78
    fraction p-value < 0.05 = 1.0

T0p5_c0p0005_ratio_1e3_local_activation:
    beta_a_R = 5.0
    beta_a_T = 3.0
    mean T-zone density ratio ≈ 12.93
    mean chi2/df ≈ 10.28
    fraction p-value < 0.05 = 1.0

T1e-1_c1e-4_ratio_1e3_local_activation:
    beta_a_R = 5.0
    beta_a_T = 3.0
    mean T-zone density ratio ≈ 17.69
    mean chi2/df ≈ 10.71
    fraction p-value < 0.05 = 1.0
```

Ces trois scénarios rejettent l'uniformité pour toutes les réplications et montrent une densité de morts beaucoup plus élevée dans la T-zone.

**Conclusion principale :**

```text
The transition from spatially uniform deaths to T-zone concentrated deaths is mainly controlled by the localization of active zones, not by the intensity contrast lambda_a_T / lambda_a_c alone.
```

---

## 6. Analyse de la grille beta pour ratio = 100

Une grille de paramètres a été testée avec :

[
\lambda_{a,T}=0.5,
\qquad
\lambda_{a,c}=0.005,
\qquad
\frac{\lambda_{a,T}}{\lambda_{a,c}}=100.
]

Ce cas est particulièrement intéressant car (\lambda_{a,T}=0.5) reste proche du baseline, tandis que l'on étudie l'effet de la localisation des zones actives.

Résumé qualitatif des résultats :

| (\beta_{a,R}) |                        (\beta_{a,T}=1.2) |          (\beta_{a,T}=2.0) |          (\beta_{a,T}=3.0) |          (\beta_{a,T}=5.0) |
| ------------: | ---------------------------------------: | -------------------------: | -------------------------: | -------------------------: |
|           2.5 |                 T-ratio ≈ 0.99, uniforme |     T-ratio ≈ 1.19, limite |      T-ratio ≈ 2.16, clair |       T-ratio ≈ 4.20, fort |
|           3.5 | T-ratio ≈ 1.29, faible mais significatif |       T-ratio ≈ 4.30, fort |       T-ratio ≈ 6.24, fort |       T-ratio ≈ 8.15, fort |
|           5.0 |                     T-ratio ≈ 8.35, fort |      T-ratio ≈ 11.30, fort |      T-ratio ≈ 11.93, fort |      T-ratio ≈ 13.23, fort |
|           7.5 |               T-ratio ≈ 18.38, très fort | T-ratio ≈ 17.12, très fort | T-ratio ≈ 23.00, très fort | T-ratio ≈ 22.15, très fort |

Le passage important se produit déjà autour de :

[
\beta_{a,R}=2.5,
\qquad
\beta_{a,T}=3.0.
]

À ce point, les cinq réplications rejettent l'uniformité et le ratio de densité dans la T-zone est environ égal à 2.16. Ce scénario peut être vu comme un premier point où la T-zone devient réellement visible.

Un autre point très utile pour les figures est :

[
\beta_{a,R}=3.5,
\qquad
\beta_{a,T}=2.0.
]

Dans ce cas, le ratio de densité dans la T-zone est environ égal à 4.30, avec une structure plus claire mais sans être aussi extrême que les scénarios très localisés.

---

## 7. Analyse de la grille beta pour ratio = 1000

Une seconde grille a été testée avec :

[
\lambda_{a,T}=0.1,
\qquad
\lambda_{a,c}=10^{-4},
\qquad
\frac{\lambda_{a,T}}{\lambda_{a,c}}=1000.
]

Les résultats montrent la même tendance que pour le ratio 100, mais les temps nécessaires pour observer 1000 morts deviennent beaucoup plus grands lorsque les zones actives sont très locales.

Exemples :

```text
beta_a_R = 2.5, beta_a_T = 1.2:
    mean T-zone density ratio ≈ 0.94
    fraction p-value < 0.05 = 0.0
    mean time ≈ 23.70

beta_a_R = 2.5, beta_a_T = 2.0:
    mean T-zone density ratio ≈ 1.05
    fraction p-value < 0.05 = 0.6
    mean time ≈ 33.40

beta_a_R = 2.5, beta_a_T = 3.0:
    mean T-zone density ratio ≈ 2.36
    fraction p-value < 0.05 = 1.0
    mean time ≈ 146.12

beta_a_R = 3.5, beta_a_T = 2.0:
    mean T-zone density ratio ≈ 5.33
    fraction p-value < 0.05 = 1.0
    mean time ≈ 369.19

beta_a_R = 5.0, beta_a_T = 3.0:
    mean T-zone density ratio ≈ 17.69
    fraction p-value < 0.05 = 1.0
    mean time ≈ 2338.41
```

Cela confirme que le ratio (\lambda_{a,T}/\lambda_{a,c}) ne suffit pas à lui seul. La T-zone devient visible seulement si les zones actives sont suffisamment locales.

---

## 8. Point critique proposé

Il y a deux notions possibles de point critique.

### 8.1. Point critique statistique

On peut définir le point critique statistique comme le plus petit jeu de paramètres pour lequel les cinq réplications rejettent l'uniformité spatiale au seuil 5%.

Avec ratio = 100, le scénario :

```text
T0p5_c0p005_ratio_1e2_betaR3p5_betaT1p2
```

rejette déjà l'uniformité dans les cinq réplications. Cependant, le ratio de densité dans la T-zone est seulement environ 1.29. La non-uniformité est donc statistiquement détectable, mais visuellement et biologiquement encore faible.

### 8.2. Point critique visuel/interprétable

Pour obtenir une T-zone clairement visible, il est plus raisonnable de demander à la fois :

```text
fraction p-value < 0.05 = 1.0
```

et

```text
T-zone density ratio > 2.
```

Avec ce critère, un bon point critique est :

```text
T0p5_c0p005_ratio_1e2_betaR2p5_betaT3p0
```

c'est-à-dire :

[
\lambda_{a,T}=0.5,
\qquad
\lambda_{a,c}=0.005,
\qquad
\beta_{a,R}=2.5,
\qquad
\beta_{a,T}=3.0.
]

Ce scénario donne :

```text
mean T-zone density ratio ≈ 2.16
fraction p-value < 0.05 = 1.0
mean time for 1000 deaths ≈ 60.51
```

Un scénario plus clair, mais encore raisonnable, est :

```text
T0p5_c0p005_ratio_1e2_betaR3p5_betaT2p0
```

avec :

```text
mean T-zone density ratio ≈ 4.30
fraction p-value < 0.05 = 1.0
mean time for 1000 deaths ≈ 107.78
```

Un scénario très fort, utile pour illustrer clairement la T-zone, est :

```text
T0p5_c0p005_ratio_1e2_betaR5p0_betaT3p0
```

avec :

```text
mean T-zone density ratio ≈ 11.93
fraction p-value < 0.05 = 1.0
mean time for 1000 deaths ≈ 488.14
```

---

## 9. Conclusion générale

Les simulations montrent que le processus observé de mort cellulaire reste approximativement uniforme sous les paramètres baseline. Augmenter uniquement le contraste d'intensité (\lambda_{a,T}/\lambda_{a,c}) ne suffit pas à créer un motif visible en forme de T lorsque les zones actives restent grandes ou longues.

En revanche, lorsque les zones actives deviennent plus locales, soit par une diminution de leur rayon moyen, soit par une diminution de leur durée de vie moyenne, la structure de la T-zone devient très visible dans les morts observées.

Ainsi, la transition entre un régime spatialement uniforme et un régime concentré dans la T-zone est principalement contrôlée par la localisation des zones actives :

[
\beta_{a,R}
\quad \text{et} \quad
\beta_{a,T},
]

et non uniquement par le ratio :

[
\frac{\lambda_{a,T}}{\lambda_{a,c}}.
]

Une formulation synthétique est :

```text
The simulations show that increasing only the contrast lambda_a_T/lambda_a_c does not reliably produce a visible T-zone pattern under the baseline active-zone radius and lifetime. In contrast, the spatial pattern becomes strongly non-uniform when active zones are made more local by increasing beta_a_R and/or beta_a_T. Therefore, the transition from approximately uniform to T-zone concentrated observed deaths is mainly controlled by active-zone localization, not by the intensity contrast alone.
```

---

## 10. Scénarios recommandés pour le rapport

Pour présenter les résultats clairement, il est préférable de ne pas montrer seulement les cas les plus extrêmes. Une bonne progression est :

| Type de cas | Scénario recommandé                       | Rôle dans l'analyse                                          |
| ----------- | ----------------------------------------- | ------------------------------------------------------------ |
| Baseline    | `baseline`                                | Montre que le modèle de référence est proche de l'uniformité |
| Ratio seul  | `T0p5_c0p005_ratio_1e2`                   | Montre que le ratio seul ne suffit pas                       |
| Transition  | `T0p5_c0p005_ratio_1e2_betaR2p5_betaT3p0` | Premier cas clairement non uniforme                          |
| Cas clair   | `T0p5_c0p005_ratio_1e2_betaR3p5_betaT2p0` | T-zone visible sans être trop extrême                        |
| Cas fort    | `T0p5_c0p005_ratio_1e2_betaR5p0_betaT3p0` | Exemple très net de concentration dans la T-zone             |

Cette progression raconte une histoire cohérente :

```text
baseline uniform
→ ratio alone is not sufficient
→ local activation reveals the T-zone
→ the critical transition depends on active-zone localization.
```
