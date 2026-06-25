# T-shaped rejection-sampling analysis after burn-in

All 2D histograms use the same colorbar scale: 0 to 13 observed deaths per bin.

The simulation keeps the block rejection sampling / thinning structure. The first 500 observed deaths are burn-in. The reported statistics use the next 200 observed deaths.

## T-shaped activation intensity

For activation candidates, the acceptance probability is `lambda_a(x | V_a) / lambda_a_1`, where

$$
\lambda^a(x\mid V^a_{s-}) =
\begin{cases}
\lambda^a_1, & x\in A(V^a_{s-}),\\
\lambda^a_T, & x\notin A(V^a_{s-})\text{ and }x\in\mathcal{T},\\
\lambda^a_c, & \text{otherwise}.
\end{cases}
$$

The fixed T-zone is the middle vertical column plus the left middle arm of a 3 x 3 grid.

## Spatial coefficient of variation

Let `H_i` be the observed death count in spatial bin `i`. The spatial coefficient of variation is

$$ CV_{\mathrm{space}} = \frac{\mathrm{sd}(H_1,\ldots,H_k)}{\mathrm{mean}(H_1,\ldots,H_k)} $$

Under spatial uniformity, a useful finite-sample reference is

$$ CV_{\mathrm{uniform}} \approx \sqrt{\frac{k-1}{n}}. $$

For the 10 x 10 plot grid, `k = 100` and `n = 200`, so `CV_uniform ≈ sqrt((100-1)/200) = 0.704`.

## Chi-square uniformity test

The 10 x 10 grid has expected count `200/100 = 2.00` per bin. The formal chi-square test uses the coarser 5 x 5 grid, where expected count is `200/25 = 8.00` per bin.

The Pearson statistic is

$$ X^2 = \sum_{i=1}^{k} \frac{(O_i-E_i)^2}{E_i} $$

Under spatial uniformity, `X^2` is approximately chi-square with `df = 24`. Thus `X^2/df` should fluctuate around 1 with standard deviation approximately `sqrt(2/df) = 0.289`.

## Results

| scenario | lambda_a_T | lambda_a_c | time for 200 deaths | death rejection | outside-active | ERK | CV ratio | chi2/df | p-value | T density ratio | T deaths | empty bins |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | 0.1 | 0.01 | 10 | 32.713 | 0.723 | 0.345 | 0.378 | 1.159 | 2.719 | 1.12e-05 | 1.382 | 0.0134 | 105 | 15 |
| higher_death_rate | 0.1 | 0.01 | 10 | 28.000 | 0.747 | 0.266 | 0.481 | 1.142 | 1.448 | 0.0722 | 1.225 | 0.086 | 99 | 20 |
| larger_ERK_radius | 0.1 | 0.01 | 10 | 33.165 | 0.770 | 0.320 | 0.449 | 0.916 | 1.250 | 0.185 | 1.086 | 0.303 | 93 | 10 |
| shorter_ERK_duration | 0.1 | 0.01 | 10 | 14.590 | 0.564 | 0.372 | 0.192 | 0.835 | 0.646 | 0.905 | 1.497 | 0.00271 | 109 | 8 |
| stronger_T_zone_activation | 0.3 | 0.01 | 30 | 24.293 | 0.644 | 0.258 | 0.386 | 1.092 | 1.667 | 0.0214 | 1.410 | 0.00923 | 106 | 14 |
| strong_visible_T_zone | 0.3 | 0.0001 | 3000 | 84.347 | 0.912 | 0.777 | 0.135 | 1.741 | 7.260 | 9.13e-25 | 5.893 | 1.53e-28 | 165 | 38 |
| ratio_T_over_c_1e1 | 0.1 | 0.01 | 10 | 20.345 | 0.667 | 0.308 | 0.359 | 1.054 | 1.438 | 0.0762 | 1.250 | 0.0658 | 100 | 16 |
| ratio_T_over_c_1e2 | 0.1 | 0.001 | 100 | 25.459 | 0.713 | 0.335 | 0.378 | 0.953 | 0.948 | 0.535 | 1.108 | 0.255 | 94 | 14 |
| ratio_T_over_c_1e3 | 0.1 | 0.0001 | 1000 | 24.775 | 0.672 | 0.303 | 0.369 | 1.115 | 1.958 | 0.00335 | 1.131 | 0.212 | 95 | 14 |
| ratio_T_over_c_1e4 | 0.1 | 1e-05 | 1e+04 | 30.332 | 0.737 | 0.357 | 0.380 | 1.064 | 1.562 | 0.039 | 1.657 | 0.000242 | 114 | 15 |
| high_T_ratio_1e3 | 0.3 | 0.0003 | 1000 | 20.632 | 0.647 | 0.271 | 0.376 | 1.082 | 1.667 | 0.0214 | 1.225 | 0.086 | 99 | 14 |

## Short conclusion

The baseline T-shaped case has T-zone density ratio 1.382. A ratio above 1 means death locations are denser in the T-zone than outside.

The fastest scenario is `shorter_ERK_duration`, with 14.590 time units for the analyzed 200 deaths.

The strongest T-zone enrichment is `strong_visible_T_zone`, with T-zone density ratio 5.893.

The largest spatial CV ratio is `strong_visible_T_zone`, with CV ratio 1.741.
