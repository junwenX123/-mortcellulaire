# Spatial uniformity analysis after burn-in
All 2D histograms use the same colorbar scale: 0 to 9 observed deaths per bin.
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
| shorter_ERK_duration | 16.955 | 0.578 | 0.848 | 0.859 | 1.061 | 28 |
| smaller_ERK_radius | 17.695 | 0.579 | 0.909 | 0.73 | 1.014 | 29 |
| optimized_mild | 16.022 | 0.648 | 1.038 | 0.377 | 1.035 | 38 |
| optimized_strong | 13.337 | 0.643 | 1.111 | 0.211 | 1.031 | 28 |
| moderate_T_zone_activation | 15.910 | 0.606 | 0.853 | 0.852 | 0.990 | 22 |

## Short conclusion
The baseline case has CV = 0.623 and T-zone density ratio = 0.943. Thus it is close to spatially uniform compared with the stronger T-zone scenario.
The most spatially heterogeneous scenario by CV is `optimized_mild` with CV = 0.648 and T-zone density ratio = 1.035.
