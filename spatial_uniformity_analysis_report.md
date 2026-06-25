# Spatial uniformity analysis after burn-in
All 2D histograms use the same colorbar scale: 0 to 17 observed deaths per bin.
Burn-in removes the initial transient period: the first 500 observed deaths are simulated but not used for statistics. The reported statistics use the next 1000 observed deaths.
## Interpretation rules
- `spatial_coefficient_of_variation_plot_bins = std(H)/mean(H)`: close to 0 means more uniform; larger values mean more spatial heterogeneity.
- `spatial_chi2_per_df_chi2_bins`: close to 1 is compatible with uniform counts; much larger than 1 indicates over-dispersion / non-uniformity.
- `spatial_chi2_p_value_chi2_bins`: small p-value, for example < 0.05, rejects spatial uniformity on the chosen 10 x 10 grid.
- `t_zone_density_ratio`: values > 1 mean deaths are denser in the T-zone than outside; values < 1 mean the opposite.

## Results
| scenario | time for analyzed 1000 deaths | CV 20x20 | chi2/df 10x10 | p-value 10x10 | T-zone density ratio | empty bins 20x20 |
|---|---:|---:|---:|---:|---:|---:|
| baseline__rep01 | 24.347 | 0.620 | 1.079 | 0.278 | 1.131 | 32 |
| baseline__rep02 | 22.789 | 0.631 | 1.107 | 0.219 | 0.932 | 30 |
| baseline__rep03 | 25.512 | 0.623 | 1.133 | 0.172 | 0.966 | 27 |
| baseline__rep04 | 23.451 | 0.639 | 1.099 | 0.235 | 0.890 | 39 |
| baseline__rep05 | 24.660 | 0.597 | 1.105 | 0.223 | 0.955 | 24 |
| higher_death_rate__rep01 | 18.756 | 0.602 | 1.055 | 0.336 | 1.048 | 26 |
| higher_death_rate__rep02 | 19.977 | 0.610 | 1.040 | 0.372 | 0.935 | 31 |
| higher_death_rate__rep03 | 21.913 | 0.593 | 0.804 | 0.924 | 0.932 | 36 |
| higher_death_rate__rep04 | 21.669 | 0.635 | 1.244 | 0.0502 | 1.117 | 25 |
| higher_death_rate__rep05 | 19.577 | 0.647 | 1.301 | 0.0237 | 0.966 | 26 |
| larger_ERK_radius__rep01 | 56.894 | 0.683 | 1.446 | 0.00245 | 0.913 | 38 |
| larger_ERK_radius__rep02 | 52.666 | 0.667 | 1.552 | 0.000363 | 0.840 | 37 |
| larger_ERK_radius__rep03 | 52.519 | 0.688 | 1.370 | 0.00861 | 0.876 | 47 |
| larger_ERK_radius__rep04 | 50.316 | 0.664 | 1.525 | 0.000597 | 0.932 | 39 |
| larger_ERK_radius__rep05 | 65.104 | 0.661 | 1.636 | 6.69e-05 | 0.909 | 40 |
| shorter_ERK_duration__rep01 | 16.900 | 0.578 | 0.857 | 0.845 | 0.970 | 23 |
| shorter_ERK_duration__rep02 | 16.088 | 0.611 | 0.887 | 0.782 | 0.943 | 34 |
| shorter_ERK_duration__rep03 | 17.354 | 0.620 | 0.844 | 0.866 | 1.086 | 31 |
| shorter_ERK_duration__rep04 | 16.520 | 0.613 | 0.895 | 0.764 | 0.935 | 39 |
| shorter_ERK_duration__rep05 | 16.543 | 0.611 | 0.840 | 0.873 | 1.002 | 30 |
| stronger_T_zone_activation__rep01 | 22.775 | 0.622 | 1.085 | 0.265 | 0.909 | 33 |
| stronger_T_zone_activation__rep02 | 26.479 | 0.652 | 1.032 | 0.393 | 0.731 | 40 |
| stronger_T_zone_activation__rep03 | 27.158 | 0.609 | 1.083 | 0.269 | 0.990 | 33 |
| stronger_T_zone_activation__rep04 | 24.111 | 0.611 | 0.923 | 0.694 | 0.943 | 31 |
| stronger_T_zone_activation__rep05 | 23.688 | 0.615 | 0.905 | 0.74 | 1.002 | 31 |
| strong_visible_T_zone__rep01 | 166.164 | 1.218 | 10.475 | 3.88e-156 | 15.417 | 179 |
| strong_visible_T_zone__rep02 | 169.827 | 1.168 | 10.182 | 1.95e-150 | 12.955 | 174 |
| strong_visible_T_zone__rep03 | 167.982 | 1.104 | 9.210 | 1.18e-131 | 11.376 | 166 |
| strong_visible_T_zone__rep04 | 158.214 | 1.131 | 9.624 | 1.24e-139 | 10.011 | 164 |
| strong_visible_T_zone__rep05 | 165.274 | 1.126 | 9.655 | 3.21e-140 | 11.005 | 175 |
| T0p5_c0p005_ratio_1e2__rep01 | 24.380 | 0.621 | 1.166 | 0.124 | 0.994 | 31 |
| T0p5_c0p005_ratio_1e2__rep02 | 25.268 | 0.633 | 1.117 | 0.2 | 0.861 | 36 |
| T0p5_c0p005_ratio_1e2__rep03 | 24.605 | 0.611 | 1.024 | 0.414 | 1.108 | 32 |
| T0p5_c0p005_ratio_1e2__rep04 | 23.052 | 0.607 | 0.974 | 0.555 | 0.982 | 28 |
| T0p5_c0p005_ratio_1e2__rep05 | 24.713 | 0.651 | 0.970 | 0.567 | 1.019 | 37 |
| T0p5_c0p0005_ratio_1e3__rep01 | 27.621 | 0.639 | 1.317 | 0.0189 | 0.958 | 35 |
| T0p5_c0p0005_ratio_1e3__rep02 | 24.418 | 0.621 | 0.943 | 0.64 | 0.982 | 33 |
| T0p5_c0p0005_ratio_1e3__rep03 | 27.135 | 0.594 | 0.804 | 0.924 | 0.883 | 28 |
| T0p5_c0p0005_ratio_1e3__rep04 | 25.182 | 0.650 | 1.554 | 0.000349 | 0.854 | 37 |
| T0p5_c0p0005_ratio_1e3__rep05 | 25.168 | 0.613 | 0.824 | 0.898 | 0.924 | 34 |
| T0p5_c0p00005_ratio_1e4__rep01 | 27.621 | 0.639 | 1.317 | 0.0189 | 0.958 | 35 |
| T0p5_c0p00005_ratio_1e4__rep02 | 24.418 | 0.621 | 0.943 | 0.64 | 0.982 | 33 |
| T0p5_c0p00005_ratio_1e4__rep03 | 27.135 | 0.594 | 0.804 | 0.924 | 0.883 | 28 |
| T0p5_c0p00005_ratio_1e4__rep04 | 25.182 | 0.650 | 1.554 | 0.000349 | 0.854 | 37 |
| T0p5_c0p00005_ratio_1e4__rep05 | 25.168 | 0.613 | 0.824 | 0.898 | 0.924 | 34 |
| T0p5_c0p025_ratio_20__rep01 | 23.172 | 0.595 | 0.651 | 0.997 | 0.982 | 33 |
| T0p5_c0p025_ratio_20__rep02 | 23.638 | 0.620 | 0.881 | 0.796 | 0.894 | 31 |
| T0p5_c0p025_ratio_20__rep03 | 26.476 | 0.668 | 1.186 | 0.1 | 1.100 | 30 |
| T0p5_c0p025_ratio_20__rep04 | 24.879 | 0.628 | 0.923 | 0.694 | 0.898 | 32 |
| T0p5_c0p025_ratio_20__rep05 | 25.149 | 0.650 | 1.137 | 0.165 | 0.939 | 38 |
| T0p5_c0p01_ratio_50__rep01 | 25.199 | 0.593 | 0.768 | 0.958 | 1.069 | 24 |
| T0p5_c0p01_ratio_50__rep02 | 25.268 | 0.633 | 1.117 | 0.2 | 0.861 | 36 |
| T0p5_c0p01_ratio_50__rep03 | 23.651 | 0.633 | 1.055 | 0.336 | 0.932 | 42 |
| T0p5_c0p01_ratio_50__rep04 | 23.052 | 0.607 | 0.974 | 0.555 | 0.982 | 28 |
| T0p5_c0p01_ratio_50__rep05 | 22.087 | 0.607 | 1.006 | 0.464 | 0.994 | 31 |
| T1e-1_c1e-2_ratio_1e1__rep01 | 23.138 | 0.582 | 1.103 | 0.227 | 0.962 | 21 |
| T1e-1_c1e-2_ratio_1e1__rep02 | 23.618 | 0.653 | 1.174 | 0.114 | 0.861 | 36 |
| T1e-1_c1e-2_ratio_1e1__rep03 | 27.377 | 0.623 | 1.002 | 0.475 | 1.078 | 31 |
| T1e-1_c1e-2_ratio_1e1__rep04 | 25.296 | 0.642 | 1.143 | 0.156 | 0.861 | 39 |
| T1e-1_c1e-2_ratio_1e1__rep05 | 23.175 | 0.629 | 0.974 | 0.555 | 0.909 | 38 |
| T1e-1_c1e-3_ratio_1e2__rep01 | 23.138 | 0.582 | 1.103 | 0.227 | 0.962 | 21 |
| T1e-1_c1e-3_ratio_1e2__rep02 | 23.618 | 0.653 | 1.174 | 0.114 | 0.861 | 36 |
| T1e-1_c1e-3_ratio_1e2__rep03 | 24.724 | 0.646 | 1.172 | 0.117 | 0.998 | 36 |
| T1e-1_c1e-3_ratio_1e2__rep04 | 24.047 | 0.625 | 0.867 | 0.825 | 0.872 | 25 |
| T1e-1_c1e-3_ratio_1e2__rep05 | 22.971 | 0.636 | 0.901 | 0.75 | 0.998 | 29 |
| T1e-1_c1e-4_ratio_1e3__rep01 | 23.138 | 0.582 | 1.103 | 0.227 | 0.962 | 21 |
| T1e-1_c1e-4_ratio_1e3__rep02 | 23.618 | 0.653 | 1.174 | 0.114 | 0.861 | 36 |
| T1e-1_c1e-4_ratio_1e3__rep03 | 24.724 | 0.646 | 1.172 | 0.117 | 0.998 | 36 |
| T1e-1_c1e-4_ratio_1e3__rep04 | 24.047 | 0.625 | 0.867 | 0.825 | 0.872 | 25 |
| T1e-1_c1e-4_ratio_1e3__rep05 | 22.971 | 0.636 | 0.901 | 0.75 | 0.998 | 29 |
| T1e-2_c1e-3_ratio_1e1__rep01 | 27.243 | 0.604 | 1.210 | 0.076 | 0.876 | 24 |
| T1e-2_c1e-3_ratio_1e1__rep02 | 24.226 | 0.615 | 0.875 | 0.809 | 0.858 | 36 |
| T1e-2_c1e-3_ratio_1e1__rep03 | 28.304 | 0.616 | 1.358 | 0.0104 | 1.078 | 28 |
| T1e-2_c1e-3_ratio_1e1__rep04 | 26.765 | 0.643 | 1.119 | 0.196 | 0.920 | 32 |
| T1e-2_c1e-3_ratio_1e1__rep05 | 25.947 | 0.620 | 1.079 | 0.278 | 1.010 | 32 |
| T1e-2_c1e-4_ratio_1e2__rep01 | 27.243 | 0.604 | 1.210 | 0.076 | 0.876 | 24 |
| T1e-2_c1e-4_ratio_1e2__rep02 | 24.226 | 0.615 | 0.875 | 0.809 | 0.858 | 36 |
| T1e-2_c1e-4_ratio_1e2__rep03 | 25.285 | 0.639 | 1.103 | 0.227 | 0.939 | 29 |
| T1e-2_c1e-4_ratio_1e2__rep04 | 26.765 | 0.643 | 1.119 | 0.196 | 0.920 | 32 |
| T1e-2_c1e-4_ratio_1e2__rep05 | 25.947 | 0.620 | 1.079 | 0.278 | 1.010 | 32 |
| T1e-3_c1e-4_ratio_1e1__rep01 | 27.378 | 0.632 | 0.974 | 0.555 | 0.958 | 41 |
| T1e-3_c1e-4_ratio_1e1__rep02 | 24.837 | 0.671 | 1.386 | 0.00668 | 0.898 | 34 |
| T1e-3_c1e-4_ratio_1e1__rep03 | 25.369 | 0.663 | 1.295 | 0.0258 | 0.876 | 34 |
| T1e-3_c1e-4_ratio_1e1__rep04 | 25.689 | 0.608 | 0.956 | 0.606 | 0.924 | 37 |
| T1e-3_c1e-4_ratio_1e1__rep05 | 21.549 | 0.611 | 1.059 | 0.326 | 0.920 | 27 |
| T0p5_c0p005_ratio_1e2_local_activation__rep01 | 500.963 | 1.138 | 10.204 | 7.2e-151 | 12.191 | 173 |
| T0p5_c0p005_ratio_1e2_local_activation__rep02 | 495.324 | 1.123 | 9.673 | 1.43e-140 | 12.795 | 171 |
| T0p5_c0p005_ratio_1e2_local_activation__rep03 | 471.502 | 1.145 | 9.459 | 1.94e-136 | 11.505 | 167 |
| T0p5_c0p005_ratio_1e2_local_activation__rep04 | 488.157 | 1.167 | 9.804 | 4.13e-143 | 11.908 | 166 |
| T0p5_c0p005_ratio_1e2_local_activation__rep05 | 484.775 | 1.130 | 9.756 | 3.58e-142 | 11.250 | 167 |
| T0p5_c0p0005_ratio_1e3_local_activation__rep01 | 477.195 | 1.138 | 10.026 | 2.04e-147 | 11.771 | 168 |
| T0p5_c0p0005_ratio_1e3_local_activation__rep02 | 516.548 | 1.183 | 10.487 | 2.25e-156 | 16.111 | 178 |
| T0p5_c0p0005_ratio_1e3_local_activation__rep03 | 513.689 | 1.187 | 10.044 | 9.06e-148 | 12.955 | 176 |
| T0p5_c0p0005_ratio_1e3_local_activation__rep04 | 482.786 | 1.164 | 10.651 | 1.44e-159 | 14.182 | 177 |
| T0p5_c0p0005_ratio_1e3_local_activation__rep05 | 481.140 | 1.134 | 10.198 | 9.45e-151 | 9.620 | 175 |
| T1e-1_c1e-4_ratio_1e3_local_activation__rep01 | 2194.206 | 1.180 | 10.521 | 4.81e-157 | 13.631 | 180 |
| T1e-1_c1e-4_ratio_1e3_local_activation__rep02 | 2348.940 | 1.205 | 11.034 | 4.49e-167 | 19.242 | 185 |
| T1e-1_c1e-4_ratio_1e3_local_activation__rep03 | 2441.805 | 1.198 | 10.851 | 1.78e-163 | 21.898 | 190 |
| T1e-1_c1e-4_ratio_1e3_local_activation__rep04 | 2538.260 | 1.202 | 10.929 | 5.12e-165 | 18.911 | 183 |
| T1e-1_c1e-4_ratio_1e3_local_activation__rep05 | 2168.858 | 1.171 | 10.220 | 3.49e-151 | 14.776 | 182 |

## Short conclusion
The baseline case has CV = 0.620 and T-zone density ratio = 1.131. Thus it is close to spatially uniform compared with the stronger T-zone scenario.
The most spatially heterogeneous scenario by CV is `strong_visible_T_zone__rep01` with CV = 1.218 and T-zone density ratio = 15.417.
