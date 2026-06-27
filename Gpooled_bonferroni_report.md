# Pooled T-zone binomial test with Bonferroni correction

## Method

For each parameter setting, the five independent replicates were pooled. If a parameter setting has five replicates of 1000 analyzed observed deaths, then the pooled sample size is N = 5000. The pooled number of T-zone deaths was tested against the null probability p0 = |T|/|W| using a one-sided exact binomial test.

The hypotheses are:

$$H_0:\pi=p_0,\qquad H_A:\pi>p_0.$$

Since 52 parameter settings were tested, Bonferroni correction was applied. With global alpha = 0.05, the Bonferroni threshold is:

$$\alpha_\mathrm{Bonferroni} = \frac{0.05}{52} = 0.000961538.$$

A parameter setting is considered significant only if:

$$p_\mathrm{pooled} < 0.000961538.$$

## Results

| parameter setting | replicates | pooled k/n | pooled fraction | pooled p-value | Bonferroni adjusted p-value | significant after Bonferroni |
|---|---:|---:|---:|---:|---:|---:|
| strong_visible_T_zone | 5 | 4525/5000 | 0.9050 | 0 | 0 | True |
| T0p5_c0p005_ratio_1e2_local_activation | 5 | 4525/5000 | 0.9050 | 0 | 0 | True |
| T0p5_c0p0005_ratio_1e3_local_activation | 5 | 4548/5000 | 0.9096 | 0 | 0 | True |
| T1e-1_c1e-4_ratio_1e3_local_activation | 5 | 4661/5000 | 0.9322 | 0 | 0 | True |
| T0p5_c0p005_ratio_1e2_betaR2p5_betaT5p0 | 5 | 3848/5000 | 0.7696 | 0 | 0 | True |
| T0p5_c0p005_ratio_1e2_betaR3p5_betaT2p0 | 5 | 3865/5000 | 0.7730 | 0 | 0 | True |
| T0p5_c0p005_ratio_1e2_betaR3p5_betaT3p0 | 5 | 4161/5000 | 0.8322 | 0 | 0 | True |
| T0p5_c0p005_ratio_1e2_betaR3p5_betaT5p0 | 5 | 4328/5000 | 0.8656 | 0 | 0 | True |
| T0p5_c0p005_ratio_1e2_betaR5p0_betaT1p2 | 5 | 4340/5000 | 0.8680 | 0 | 0 | True |
| T0p5_c0p005_ratio_1e2_betaR5p0_betaT2p0 | 5 | 4498/5000 | 0.8996 | 0 | 0 | True |
| T0p5_c0p005_ratio_1e2_betaR5p0_betaT3p0 | 5 | 4525/5000 | 0.9050 | 0 | 0 | True |
| T0p5_c0p005_ratio_1e2_betaR5p0_betaT5p0 | 5 | 4566/5000 | 0.9132 | 0 | 0 | True |
| T0p5_c0p005_ratio_1e2_betaR7p5_betaT1p2 | 5 | 4662/5000 | 0.9324 | 0 | 0 | True |
| T0p5_c0p005_ratio_1e2_betaR7p5_betaT2p0 | 5 | 4658/5000 | 0.9316 | 0 | 0 | True |
| T0p5_c0p005_ratio_1e2_betaR7p5_betaT3p0 | 5 | 4738/5000 | 0.9476 | 0 | 0 | True |
| T0p5_c0p005_ratio_1e2_betaR7p5_betaT5p0 | 5 | 4725/5000 | 0.9450 | 0 | 0 | True |
| T1e-1_c1e-4_ratio_1e3_betaR2p5_betaT5p0 | 5 | 3908/5000 | 0.7816 | 0 | 0 | True |
| T1e-1_c1e-4_ratio_1e3_betaR3p5_betaT2p0 | 5 | 4043/5000 | 0.8086 | 0 | 0 | True |
| T1e-1_c1e-4_ratio_1e3_betaR3p5_betaT3p0 | 5 | 4297/5000 | 0.8594 | 0 | 0 | True |
| T1e-1_c1e-4_ratio_1e3_betaR3p5_betaT5p0 | 5 | 4429/5000 | 0.8858 | 0 | 0 | True |
| T1e-1_c1e-4_ratio_1e3_betaR5p0_betaT1p2 | 5 | 4398/5000 | 0.8796 | 0 | 0 | True |
| T1e-1_c1e-4_ratio_1e3_betaR5p0_betaT2p0 | 5 | 4518/5000 | 0.9036 | 0 | 0 | True |
| T1e-1_c1e-4_ratio_1e3_betaR5p0_betaT3p0 | 5 | 4661/5000 | 0.9322 | 0 | 0 | True |
| T1e-1_c1e-4_ratio_1e3_betaR5p0_betaT5p0 | 5 | 4627/5000 | 0.9254 | 0 | 0 | True |
| T1e-1_c1e-4_ratio_1e3_betaR7p5_betaT1p2 | 5 | 4725/5000 | 0.9450 | 0 | 0 | True |
| T1e-1_c1e-4_ratio_1e3_betaR7p5_betaT2p0 | 5 | 4720/5000 | 0.9440 | 0 | 0 | True |
| T1e-1_c1e-4_ratio_1e3_betaR7p5_betaT3p0 | 5 | 4743/5000 | 0.9486 | 0 | 0 | True |
| T1e-1_c1e-4_ratio_1e3_betaR7p5_betaT5p0 | 5 | 4749/5000 | 0.9498 | 0 | 0 | True |
| T1e-1_c1e-4_ratio_1e3_betaR2p5_betaT3p0 | 5 | 3252/5000 | 0.6504 | 1.5e-188 | 7.82e-187 | True |
| T0p5_c0p005_ratio_1e2_betaR2p5_betaT3p0 | 5 | 3165/5000 | 0.6330 | 4.48e-158 | 2.33e-156 | True |
| T1e-1_c1e-4_ratio_1e3_betaR3p5_betaT1p2 | 5 | 2541/5000 | 0.5082 | 8.74e-20 | 4.54e-18 | True |
| T0p5_c0p005_ratio_1e2_betaR3p5_betaT1p2 | 5 | 2532/5000 | 0.5064 | 8.67e-19 | 4.51e-17 | True |
| T0p5_c0p005_ratio_1e2_betaR2p5_betaT2p0 | 5 | 2435/5000 | 0.4870 | 8.38e-10 | 4.36e-08 | True |
| T1e-1_c1e-4_ratio_1e3_betaR2p5_betaT2p0 | 5 | 2286/5000 | 0.4572 | 0.0359 | 1 | False |
| higher_death_rate | 5 | 2219/5000 | 0.4438 | 0.542 | 1 | False |
| T0p5_c0p005_ratio_1e2 | 5 | 2210/5000 | 0.4420 | 0.641 | 1 | False |
| T0p5_c0p005_ratio_1e2_betaR2p5_betaT1p2 | 5 | 2210/5000 | 0.4420 | 0.641 | 1 | False |
| shorter_ERK_duration | 5 | 2205/5000 | 0.4410 | 0.693 | 1 | False |
| baseline | 5 | 2187/5000 | 0.4374 | 0.845 | 1 | False |
| T0p5_c0p01_ratio_50 | 5 | 2179/5000 | 0.4358 | 0.893 | 1 | False |
| T0p5_c0p025_ratio_20 | 5 | 2172/5000 | 0.4344 | 0.926 | 1 | False |
| T1e-2_c1e-3_ratio_1e1 | 5 | 2153/5000 | 0.4306 | 0.976 | 1 | False |
| T1e-1_c1e-3_ratio_1e2 | 5 | 2142/5000 | 0.4284 | 0.989 | 1 | False |
| T1e-1_c1e-4_ratio_1e3 | 5 | 2142/5000 | 0.4284 | 0.989 | 1 | False |
| T1e-1_c1e-4_ratio_1e3_betaR2p5_betaT1p2 | 5 | 2142/5000 | 0.4284 | 0.989 | 1 | False |
| T1e-1_c1e-2_ratio_1e1 | 5 | 2135/5000 | 0.4270 | 0.994 | 1 | False |
| T0p5_c0p0005_ratio_1e3 | 5 | 2119/5000 | 0.4238 | 0.998 | 1 | False |
| T0p5_c0p00005_ratio_1e4 | 5 | 2119/5000 | 0.4238 | 0.998 | 1 | False |
| T1e-2_c1e-4_ratio_1e2 | 5 | 2119/5000 | 0.4238 | 0.998 | 1 | False |
| T1e-3_c1e-4_ratio_1e1 | 5 | 2113/5000 | 0.4226 | 0.999 | 1 | False |
| stronger_T_zone_activation | 5 | 2107/5000 | 0.4214 | 1 | 1 | False |
| larger_ERK_radius | 5 | 2084/5000 | 0.4168 | 1 | 1 | False |
