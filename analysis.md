# Spatial uniformity analysis after burn-in
All 2D histograms use the same colorbar scale: 0 to 12 observed deaths per bin.
Burn-in removes the initial transient period: the first 500 observed deaths are simulated but not used for statistics. The reported statistics use the next 1000 observed deaths.
## Interpretation rules
- `spatial_coefficient_of_variation_plot_bins = std(H)/mean(H)`: close to 0 means more uniform; larger values mean more spatial heterogeneity.
- `spatial_chi2_per_df_chi2_bins`: close to 1 is compatible with uniform counts; much larger than 1 indicates over-dispersion / non-uniformity.
- `spatial_chi2_p_value_chi2_bins`: small p-value, for example < 0.05, rejects spatial uniformity on the chosen 10 x 10 grid.
- `t_zone_density_ratio`: values > 1 mean deaths are denser in the T-zone than outside; values < 1 mean the opposite.

## Results
| scenario | time for analyzed 1000 deaths | CV 20x20 | chi2/df 10x10 | p-value 10x10 | T-zone density ratio | empty bins 20x20 |
|---|---:|---:|---:|---:|---:|---:|
| baseline | 23.329 | 0.623 | 1.063 | 0.316 | 0.943 | 31 |
| higher_death_rate | 22.343 | 0.597 | 1.081 | 0.274 | 0.920 | 25 |
| larger_ERK_radius | 50.491 | 0.639 | 1.186 | 0.1 | 0.974 | 33 |
| shorter_ERK_duration | 16.589 | 0.657 | 0.794 | 0.935 | 0.978 | 35 |
| stronger_T_zone_activation | 25.125 | 0.641 | 1.210 | 0.076 | 0.879 | 31 |
| strong_visible_T_zone | 163.846 | 1.116 | 9.630 | 9.45e-140 | 10.886 | 169 |

## Short conclusion
The baseline case has CV = 0.623 and T-zone density ratio = 0.943. Thus it is close to spatially uniform compared with the stronger T-zone scenario.
The most spatially heterogeneous scenario by CV is `strong_visible_T_zone` with CV = 1.116 and T-zone density ratio = 10.886.
