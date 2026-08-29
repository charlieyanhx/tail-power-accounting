"""P3 detectability table — checked-in computation (replaces in-writing arithmetic).

Fisher-z power for IC detection under Bonferroni(1450) + 80% power, and
Hanley-McNeil AUC detectability at n1 events vs n2 non-events.
"""
import numpy as np
from scipy.stats import norm

N_TESTS = 1450
ALPHA, POWER = 0.05, 0.80
za = norm.ppf(1 - ALPHA / N_TESTS / 2)      # two-sided Bonferroni
zp = norm.ppf(POWER)
zn = norm.ppf(1 - ALPHA / 2)                # nominal two-sided

print(f"z_alpha(Bonf) = {za:.3f}   z_power = {zp:.3f}   z_nominal = {zn:.3f}")
print(f"\nIC detectability (Fisher, se = 1/sqrt(n_eff - 3)):")
for h, n in [(10, 1168 // 10), (5, 1168 // 5)]:
    se = 1 / np.sqrt(n - 3)
    print(f"  h={h:>2}  n_eff={n:>4}  |IC|min Bonf+80% = {np.tanh((za+zp)*se):.3f}"
          f"   nominal = {np.tanh((zn+zp)*se):.3f}")

def hm_se(auc, n1, n2):
    q1, q2 = auc / (2 - auc), 2 * auc**2 / (1 + auc)
    return np.sqrt((auc*(1-auc) + (n1-1)*(q1-auc**2) + (n2-1)*(q2-auc**2)) / (n1*n2))

n1, n2 = 15, 1150
print(f"\nAUC detectability (Hanley-McNeil, n1={n1} events, n2={n2}):")
for auc in (0.70, 0.75, 0.80, 0.85, 0.90):
    se = hm_se(auc, n1, n2)
    z = (auc - 0.5) / se
    print(f"  AUC {auc:.2f}: se {se:.4f}  z {z:5.2f}  "
          f"{'detectable' if z > za + zp else 'NOT detectable'} at Bonf+80%")
# solve threshold
lo, hi = 0.51, 0.999
for _ in range(60):
    mid = (lo + hi) / 2
    if (mid - 0.5) / hm_se(mid, n1, n2) > za + zp: hi = mid
    else: lo = mid
print(f"  threshold AUC (Bonf + 80% power): {hi:.3f}")
