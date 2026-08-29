# Tail Power Accounting

Code and data-availability statement for:

> **How Much Tail Prediction Could We Have Detected? A Power Accounting for a 1,450-Test Search**
> Charlie Yan, 2026. [`paper/p3_tail_power.pdf`](paper/p3_tail_power.pdf)

## What the paper claims

A null from a large signal search is uninterpretable unless the searcher reports what the search could have detected. This supplies that accounting for a 1,450-test tail-predictor search, and proposes a tail-IC/mean-IC ≥ 2 pre-registration gate to remove premium-confound false positives.

**Headline result.** At n_eff ≈ 116 (ten-day overlap) and Bonferroni over 1,450 tests, 80% power requires **|IC| ≥ 0.44**; event classification requires **AUC ≥ 0.828** on ~15 episodes. The best candidates sat at 0.21–0.23.

## Reproducibility

**FULLY REPRODUCIBLE — no data needed.** `python code/p3_power_table.py` regenerates the paper's detectability table analytically (Fisher-z for IC, Hanley–McNeil for AUC). `python code/make_p3_figure.py` regenerates Exhibit 1.

## What is here

`code/p3_power_table.py` — detectability thresholds.
`code/make_p3_figure.py` — the frontier figure.
`paper/` — paper and exhibit.

## Evidence conventions used throughout

Every performance figure in the paper carries its accounting basis inline. Unless labelled
otherwise: **line 3** = full cross-spread fills (buy at ask, sell at bid, every leg both ways),
ex-commission, marked to market daily, padded to the full business calendar. Figures labelled
**screen** are descriptive or information-coefficient statistics and are never annualised into a
Sharpe ratio. Numbers marked **invalid** appear only as invalidated examples, with the corrected
figure alongside.

All tests reported in the paper were pre-registered — horizons, controls, nulls and decision bars
fixed before execution — and deviations are recorded rather than edited away. Where pre-registration
documents exist in this repository they are included verbatim.

## Citation

```bibtex
@techreport{yan2026tailpoweraccounting,
  title  = {How Much Tail Prediction Could We Have Detected? A Power Accounting for a 1,450-Test Search},
  author = {Yan, Charlie},
  year   = {2026},
  type   = {Working paper}
}
```

## License

Code MIT (see `LICENSE`). The paper PDF is © 2026 Charlie Yan, all rights reserved.
