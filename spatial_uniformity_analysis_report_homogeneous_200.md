# Homogeneous-region 200-death analysis after burn-in

All 2D histograms use the same colorbar scale: 0 to 8 observed deaths per bin.

Burn-in convention: the first 500 observed deaths are simulated but not used for statistics. The reported statistics use the next 200 observed deaths. Therefore each run accepts 700 deaths in total, then analyzes deaths 501--700.

## Statistical interpretation

Let the analyzed 200 death locations be binned into a spatial histogram with counts \(H_1,\ldots,H_k\). For the 10 x 10 plot grid, \(k=100\) and \(n=200\). Under spatial uniformity, the reference CV is \(\sqrt{(k-1)/n}=\sqrt{99/200}\approx 0.704\). Thus `spatial_cv_ratio_to_uniform_plot_bins` close to 1 is compatible with uniform random fluctuations, whereas a value much larger than 1 indicates spatial heterogeneity.

For the formal chi-square goodness-of-fit test, we use a coarser 5 x 5 grid. Then \(k=25\), \(df=24\), and the expected count per bin under uniformity is \(200/25=8\), which is more reliable than using the 10 x 10 grid where the expected count would be only 2. The Pearson statistic is \(X^2=\sum_i (O_i-E_i)^2/E_i\), and the p-value is \(P(\chi^2_{df}\ge X^2_{obs})\).

## Results

| scenario | time for analyzed 200 deaths | death rejection rate | CV 10x10 | CV ratio | chi2/df 5x5 | p-value 5x5 | empty bins 10x10 |
|---|---:|---:|---:|---:|---:|---:|---:|
| baseline | 20.845 | 0.672 | 0.714 | 1.015 | 1.375 | 0.104 | 12 |
| higher_death_intensity | 25.073 | 0.740 | 0.700 | 0.995 | 1.385 | 0.0989 | 13 |
| stronger_ERK_protection | 32.566 | 0.785 | 0.624 | 0.888 | 0.885 | 0.624 | 11 |
| stronger_activation_feedback | 19.553 | 0.597 | 0.718 | 1.020 | 1.292 | 0.154 | 10 |
| weaker_ERK_protection | 19.703 | 0.579 | 0.735 | 1.044 | 0.865 | 0.653 | 12 |

## Short conclusion

The baseline case has CV ratio = 1.015 and chi2/df = 1.375. These quantities should be interpreted relative to the uniform references, not as isolated numbers.

The most spatially heterogeneous scenario by CV ratio is `weaker_ERK_protection`, with CV ratio = 1.044.
