"""P3 detectability table — checked-in computation (replaces in-writing arithmetic).

Fisher-z power for IC detection under Bonferroni(1450) + 80% power, and
Hanley-McNeil AUC detectability at n1 events vs n2 non-events.
"""
import numpy as np
from scipy.stats import norm

N_TESTS = 1450
ALPHA, POWER = 0.05, 0.80
SPEARMAN = 1.06     # the search's ICs are rank-residualized Spearman; Fisher se ~ 1.06/sqrt(n-3)
za = norm.ppf(1 - ALPHA / N_TESTS / 2)      # two-sided Bonferroni
zp = norm.ppf(POWER)
zn = norm.ppf(1 - ALPHA / 2)                # nominal two-sided

print(f"z_alpha(Bonf) = {za:.3f}   z_power = {zp:.3f}   z_nominal = {zn:.3f}")
print(f"\nIC detectability (Fisher-Spearman, se = {SPEARMAN}/sqrt(n_eff - 3)):")
for h, n in [(10, 1168 // 10), (5, 1168 // 5)]:
    se = SPEARMAN / np.sqrt(n - 3)
    print(f"  h={h:>2}  n_eff={n:>4}  |IC|min Bonf+80% = {np.tanh((za+zp)*se):.3f}"
          f"   nominal = {np.tanh((zn+zp)*se):.3f}")

# Episode-clustered IC branch: if tail information is concentrated in ~n_ep episode-level
# draws (the abstract's "honest unit"), the Fisher se is SPEARMAN/sqrt(n_ep - 3).
print(f"\nIC detectability, episode-clustered (n_eff = distinct loss episodes):")
for n_ep in (15, 30):
    se = SPEARMAN / np.sqrt(n_ep - 3)
    print(f"  n_ep={n_ep:>3}  |IC|min Bonf+80% = {np.tanh((za+zp)*se):.3f}"
          f"   nominal = {np.tanh((zn+zp)*se):.3f}")

def hm_se(auc, n1, n2):
    q1, q2 = auc / (2 - auc), 2 * auc**2 / (1 + auc)
    return np.sqrt((auc*(1-auc) + (n1-1)*(q1-auc**2) + (n2-1)*(q2-auc**2)) / (n1*n2))

n1, n2 = 15, 1150
# CORRECTED 2026-08-29 (fresh-context referee pass): the rejection threshold must use the
# NULL SE (A = 0.5); only the power condition uses the alternative SE. The earlier version
# used SE(A) for both, giving 0.828 -- understating the certification floor.
se0 = hm_se(0.5, n1, n2)
crit = 0.5 + za * se0
print(f"\nAUC detectability (Hanley-McNeil, n1={n1} events, n2={n2}):")
print(f"  null SE {se0:.4f}  critical value (Bonf) {crit:.4f}")
for auc in (0.70, 0.75, 0.80, 0.85, 0.90):
    se = hm_se(auc, n1, n2)
    power = norm.cdf((auc - crit) / se)
    print(f"  AUC {auc:.2f}: se {se:.4f}  power at Bonf critical value {power:.3f}")
# solve: smallest A with 80% power against the null-SE critical value
lo, hi = 0.51, 0.999
for _ in range(60):
    mid = (lo + hi) / 2
    if mid - zp * hm_se(mid, n1, n2) > crit: hi = mid
    else: lo = mid
print(f"  threshold AUC (Bonf + 80% power): {hi:.3f}")
