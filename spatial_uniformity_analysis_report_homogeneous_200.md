# Homogeneous-region observed death analysis after burn-in

All 2D histograms use the same colorbar scale: 0 to 8 observed deaths per bin.

The simulation keeps the rejection sampling / thinning structure. The first 500 observed deaths are used as burn-in and are not used for statistics. The reported statistics use the next 200 observed deaths.

## Mathematical/statistical interpretation

### Spatial coefficient of variation

Let the spatial window be divided into equal-area bins and let `H_i` be the observed death count in bin `i`. The spatial coefficient of variation is

$$
CV_{\mathrm{space}}
=
\frac{\mathrm{sd}(H_1,\ldots,H_k)}
{\mathrm{mean}(H_1,\ldots,H_k)}.
$$
Under spatial uniformity, the bin counts are approximately multinomial with probability `1/k` per bin. Thus a useful finite-sample reference is

$$
CV_{\mathrm{uniform}}
\approx
\sqrt{\frac{k-1}{n}}.
$$

For the 10 x 10 plot grid, `k = 100` and `n = 200`, so `CV_uniform ≈ sqrt((100-1)/200) = 0.704`. The ratio `CV_observed / CV_uniform` is therefore easier to interpret than the raw CV alone.

### Chi-square uniformity test

The 10 x 10 grid is useful for visualization, but its expected count under uniformity is `200/100 = 2.00` per bin. For the formal chi-square test, the script uses a coarser 5 x 5 grid, where the expected count is `200/25 = 8.00` per bin.

The Pearson statistic is

$$
X^2
=
\sum_{i=1}^{k}
\frac{(O_i-E_i)^2}{E_i}.
$$

Under spatial uniformity, `X^2` is approximately chi-square with `df = k-1 = 24` degrees of freedom. Therefore `X^2/df` should fluctuate around 1, with standard deviation approximately `sqrt(2/df) = 0.289`. The p-value is the right-tail probability `P(chi2_df >= X_obs^2)`.

## Results

| scenario | time for analyzed 200 deaths | death rejection rate | outside-active rejection | ERK rejection | CV ratio to uniform | chi2/df | p-value | empty bins |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | 40.746 | 0.735 | 0.323 | 0.412 | 1.181 | 1.677 | 0.0201 | 18 |
| higher_death_intensity | 18.916 | 0.734 | 0.279 | 0.455 | 0.974 | 1.479 | 0.0613 | 17 |
| more_active_coverage | 13.989 | 0.513 | 0.010 | 0.503 | 0.959 | 1.156 | 0.271 | 12 |
| weaker_ERK | 14.259 | 0.462 | 0.385 | 0.077 | 0.969 | 1.323 | 0.133 | 11 |
| balanced_optimized | 9.052 | 0.146 | 0.062 | 0.084 | 0.835 | 0.604 | 0.935 | 9 |
| low_rejection_test | 9.083 | 0.101 | 0.013 | 0.089 | 1.115 | 1.042 | 0.406 | 17 |

## Short conclusion

The baseline analyzes exactly 200 deaths after burn-in. Its CV ratio is 1.181; values close to 1 are compatible with uniform multinomial fluctuations on the plotting grid.

The fastest scenario is `balanced_optimized`, with 9.052 time units for the analyzed 200 deaths.

The largest spatial CV ratio is obtained by `baseline`, with CV ratio 1.181.
