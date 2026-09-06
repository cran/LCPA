# LCPA 1.0.4

-   We gratefully thank Sungbo Sim (`pposam@naver.com`) for using LCPA and for the careful, detailed issue reports and suggestions that helped improve this release.
-   Unified `LCPA()` and `LTA()` around `type.analysis = "XZ"` and `"ZY"`. The measurement model is selected by `type.model`, its estimator by `method.model`, and its detailed settings by `control.model`. 
-   Added ML and BCH three-step analyses for covariate effects on latent groups (XZ) and latent-group effects on continuous or categorical dependent variables (ZY). Regression estimation and standard errors can be analytic, numerical, or bootstrap-based; continuous dependent variables must be standardized before analysis.
-   Extended LTA with state-based and complete-path ZY analyses, time-specific or time-invariant effects, optional pooled Step 1 estimation, participant-level bootstrap, and compiled forward-backward evaluation.
-   Expanded ZY results to report conditional means, variances, category probabilities, standard errors, confidence intervals, and omnibus Wald tests. Step 3 sandwich standard errors now account for Step 2 CEP estimation uncertainty where applicable, which can change standard errors relative to earlier releases.
-   Reimplemented `LRT.test.VLMR()` using the robust-sandwich weighted chi-square reference distribution for Mplus TECH11, evaluated by Imhof's method while retaining negative eigenvalue weights. The adjusted LMR uses the general parameter-difference correction from Lo, Mendell, and Rubin (2001).
-   Improved the bootstrap likelihood-ratio test so refits inherit the original fitting controls, with sequential stopping and diagnostics for negative bootstrap LRT statistics.
-   Corrected the LCA parameter count, two-sided p-values and confidence-interval rounding, CEP orientation and normalization, LCA EM synchronization, polytomous smoothing, and several likelihood, posterior-probability, model-selection, and simulation-label calculations.
-   Improved analytic-gradient optimization, boundary and singular-matrix diagnostics, Louis-information and sandwich standard errors, and standardized fitted-object, simulation, progress, and S3 output. Optional SEM backends through `flexmix`, `Rmixmod`, and `RMixtComp` and the LCPA/LTA documentation and examples were also expanded.

# LCPA 1.0.3

-   Improved NNE performance and corrected the LPA EM algorithm.
-   Corrected documentation errors.

# LCPA 1.0.2

-   Added `adjust.model()` for aligning LCA or LPA solutions and `plot()` methods for both model types.
-   Improved Mplus variable handling and removed unnecessary Mplus statements.
-   Improved Python dependency selection and GPU-based neural-network estimation.
-   Improved bootstrap likelihood-ratio testing and corrected documentation errors.

# LCPA 1.0.1

-   Added `use.attention` control for neural-network estimation.
-   Corrected `control.NNE` configuration and LCA parameter names.

# LCPA 1.0.0

-   Initial release
