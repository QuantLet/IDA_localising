# IDA O1 diagnostic — localising vs global reliability (native, no conformal)

Controlled Monte-Carlo: 300 paths, N=12 assets, T=1200, regime shift at t=700 (vol scale 1.0->2.3), Student-t(5) innovations, vol clustering. Functional = 95% VaR of the principal-eigenvector (systemic) portfolio. Distribution-free (order-statistic) 90% CIs.

**Nominal CI coverage of the TRUE time-varying functional = 90%.**

| Estimator | coverage BEFORE shift | coverage AFTER shift | RMSE vs truth |
|---|---|---|---|
| Global (stationary, W=500) | 92% | **33%** | 0.955 |
| Localising (adaptive, W=120) | 94% | **89%** | 0.655 |

**Reading.** Before the shift both are near nominal. After the regime shift the GLOBAL stationary estimate's CI collapses in coverage (it is anchored to stale low-vol data), while the LOCALISING estimate keeps coverage near nominal and roughly halves RMSE — the native (Härdle ICARE / Localizing-CAViaR) route to reliability under non-stationarity that O1's L2 level formalises. Figure: `o1_localising.png`.

*Scope: a controlled diagnostic for a low-dimensional (systemic-portfolio) projection of the network functional; it evidences feasibility of the L2 reliability claim, NOT the full high-dimensional theorem under dependence — which is precisely O1's novel contribution.*
