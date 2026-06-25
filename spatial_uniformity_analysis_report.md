# Spatial uniformity analysis after burn-in
All 2D histograms use the same colorbar scale: 0 to 14 observed deaths per bin.
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
| T0p5_c0p005_ratio_1e2 | 25.613 | 0.631 | 1.327 | 0.0163 | 1.052 | 36 |
| T0p5_c0p0005_ratio_1e3 | 24.180 | 0.622 | 1.083 | 0.269 | 0.970 | 35 |
| T0p5_c0p00005_ratio_1e4 | 25.445 | 0.618 | 1.087 | 0.261 | 1.002 | 29 |
| T0p5_c0p025_ratio_20 | 25.190 | 0.623 | 0.939 | 0.651 | 0.943 | 23 |
| T0p5_c0p01_ratio_50 | 26.139 | 0.627 | 0.956 | 0.606 | 0.994 | 31 |
| T1e-1_c1e-2_ratio_1e1 | 24.682 | 0.618 | 1.069 | 0.302 | 1.061 | 32 |
| T1e-1_c1e-3_ratio_1e2 | 25.335 | 0.602 | 0.994 | 0.498 | 0.935 | 20 |
| T1e-1_c1e-4_ratio_1e3 | 26.282 | 0.622 | 1.257 | 0.043 | 0.869 | 24 |
| T1e-2_c1e-3_ratio_1e1 | 24.174 | 0.613 | 0.921 | 0.699 | 1.108 | 36 |
| T1e-2_c1e-4_ratio_1e2 | 23.776 | 0.612 | 0.937 | 0.656 | 1.019 | 29 |
| T1e-3_c1e-4_ratio_1e1 | 24.385 | 0.642 | 1.154 | 0.141 | 1.182 | 35 |
| T0p5_c0p005_ratio_1e2_local_activation | 516.270 | 1.162 | 10.378 | 3e-154 | 13.810 | 175 |
| T0p5_c0p0005_ratio_1e3_local_activation | 514.041 | 1.163 | 10.188 | 1.48e-150 | 14.984 | 185 |
| T1e-1_c1e-4_ratio_1e3_local_activation | 2342.542 | 1.178 | 10.067 | 3.36e-148 | 12.639 | 178 |

## Short conclusion
The baseline case has CV = 0.623 and T-zone density ratio = 0.943. Thus it is close to spatially uniform compared with the stronger T-zone scenario.
The most spatially heterogeneous scenario by CV is `T1e-1_c1e-4_ratio_1e3_local_activation` with CV = 1.178 and T-zone density ratio = 12.639.
