import pandas as pd
import numpy as np
from scipy.stats import binomtest

# ============
# Input file
# ============
input_csv = "parameter_sweep_summary_after_burnin.csv"

# ============
# Output files
# ============
output_csv = "pooled_bonferroni_summary.csv"
output_md = "pooled_bonferroni_report.md"

alpha = 0.05

df = pd.read_csv(input_csv)

required_cols = [
    "base_scenario",
    "t_zone_binomial_test_k",
    "t_zone_binomial_test_n",
    "t_zone_binomial_null_probability_area_fraction",
]

missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns in CSV: {missing}")

rows = []

for base_scenario, g in df.groupby("base_scenario", sort=False):
    pooled_k = int(g["t_zone_binomial_test_k"].sum())
    pooled_n = int(g["t_zone_binomial_test_n"].sum())
    p0 = float(g["t_zone_binomial_null_probability_area_fraction"].iloc[0])

    pooled_fraction = pooled_k / pooled_n

    pooled_p_value = binomtest(
        k=pooled_k,
        n=pooled_n,
        p=p0,
        alternative="greater",
    ).pvalue

    # pooled density ratio:
    # density in T-zone / density outside T-zone
    n_out = pooled_n - pooled_k
    if n_out > 0 and 0 < p0 < 1:
        pooled_density_ratio = (pooled_k / p0) / (n_out / (1 - p0))
    else:
        pooled_density_ratio = np.nan

    rows.append({
        "base_scenario": base_scenario,
        "n_replicates": int(len(g)),
        "pooled_t_zone_deaths_k": pooled_k,
        "pooled_total_deaths_n": pooled_n,
        "t_zone_area_fraction_p0": p0,
        "pooled_t_zone_fraction": pooled_fraction,
        "pooled_t_zone_density_ratio": pooled_density_ratio,
        "pooled_binomial_p_value_greater": pooled_p_value,
    })

pooled = pd.DataFrame(rows)

# Number of tested base scenarios
m = len(pooled)

# Bonferroni correction
bonferroni_alpha = alpha / m
pooled["bonferroni_alpha"] = bonferroni_alpha
pooled["bonferroni_adjusted_p_value"] = np.minimum(
    pooled["pooled_binomial_p_value_greater"] * m,
    1.0,
)
pooled["significant_after_bonferroni_0p05"] = (
    pooled["pooled_binomial_p_value_greater"] < bonferroni_alpha
)

# Sort: strongest evidence first
pooled = pooled.sort_values(
    by=["significant_after_bonferroni_0p05", "pooled_binomial_p_value_greater"],
    ascending=[False, True],
)

pooled.to_csv(output_csv, index=False)

# ============
# Markdown report
# ============
with open(output_md, "w", encoding="utf-8") as f:
    f.write("# Pooled T-zone binomial test with Bonferroni correction\n\n")

    f.write("## Method\n\n")
    f.write(
        "For each parameter setting, the five independent replicates were pooled. "
        "If a parameter setting has five replicates of 1000 analyzed observed deaths, "
        "then the pooled sample size is N = 5000. The pooled number of T-zone deaths "
        "was tested against the null probability p0 = |T|/|W| using a one-sided exact "
        "binomial test.\n\n"
    )

    f.write("The hypotheses are:\n\n")
    f.write(r"\[" + "\n")
    f.write(r"H_0:\pi=p_0,\qquad H_A:\pi>p_0." + "\n")
    f.write(r"\]" + "\n\n")

    f.write(
        f"Since {m} parameter settings were tested, Bonferroni correction was applied. "
        f"With global alpha = {alpha}, the Bonferroni threshold is:\n\n"
    )
    f.write(r"\[" + "\n")
    f.write(rf"\alpha_\mathrm{{Bonferroni}} = \frac{{0.05}}{{{m}}} = {bonferroni_alpha:.6g}." + "\n")
    f.write(r"\]" + "\n\n")

    f.write("A parameter setting is considered significant only if:\n\n")
    f.write(r"\[" + "\n")
    f.write(rf"p_\mathrm{{pooled}} < {bonferroni_alpha:.6g}." + "\n")
    f.write(r"\]" + "\n\n")

    f.write("## Results\n\n")
    f.write(
        "| parameter setting | replicates | pooled k/n | pooled fraction | "
        "pooled p-value | Bonferroni adjusted p-value | significant after Bonferroni |\n"
    )
    f.write("|---|---:|---:|---:|---:|---:|---:|\n")

    for _, r in pooled.iterrows():
        f.write(
            f"| {r['base_scenario']} | "
            f"{int(r['n_replicates'])} | "
            f"{int(r['pooled_t_zone_deaths_k'])}/{int(r['pooled_total_deaths_n'])} | "
            f"{float(r['pooled_t_zone_fraction']):.4f} | "
            f"{float(r['pooled_binomial_p_value_greater']):.3g} | "
            f"{float(r['bonferroni_adjusted_p_value']):.3g} | "
            f"{bool(r['significant_after_bonferroni_0p05'])} |\n"
        )

print(f"Saved: {output_csv}")
print(f"Saved: {output_md}")
print(f"Number of parameter settings: {m}")
print(f"Bonferroni alpha: {bonferroni_alpha:.6g}")
print(
    "Significant after Bonferroni:",
    int(pooled["significant_after_bonferroni_0p05"].sum()),
)