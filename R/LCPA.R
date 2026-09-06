#' Three-step latent class/profile analysis
#'
#' Fits one of two independent auxiliary-variable paths. `type.analysis = "XZ"`
#' estimates the effect of observed covariates on latent class/profile membership.
#' `type.analysis = "ZY"` estimates class/profile-specific distributions of
#' external observed dependent variables; these variables are never treated as
#' indicators.
#'
#' @param response An \eqn{N \times I} numeric matrix or data frame containing
#'   the Step 1 latent class/profile indicators. For `type.model = "LCA"`, each
#'   column is a categorical indicator and follows the same category-mapping
#'   requirements as the `response` argument of \code{\link[LCPA]{LCA}()}. For
#'   `type.model = "LPA"`, all columns are continuous, missing values are not
#'   allowed, and the indicators should be standardized with `scale()` or
#'   \code{\link[LCPA]{normalize}()} before analysis, exactly as required by \code{\link[LCPA]{LPA}()}.
#' @param L Integer number of latent classes/profiles in the Step 1 measurement
#'   model (default: 2; must be at least 2). It has the same meaning as `L` in
#'   \code{\link[LCPA]{LCA}()} or \code{\link[LCPA]{LPA}()}.
#' @param type.analysis Character string selecting the independent Step 3 path:
#'   \itemize{
#'     \item `"XZ"`: estimate the effect of observed covariates
#'       \eqn{\boldsymbol{\zeta}} on latent class/profile membership \eqn{Z};
#'       `covariates` is used and
#'       `dependent.variables` is ignored.
#'     \item `"ZY"`: estimate class/profile-specific distributions of external
#'       observed dependent variables \eqn{\mathbf{Y}}; `dependent.variables`
#'       and `family` are used and
#'       `covariates` and `ref.class` are ignored.
#'   }
#'   The default is `"XZ"`.
#' @param type.model Character string selecting the Step 1 measurement model:
#'   `"LCA"` for categorical indicators or `"LPA"` for continuous indicators
#'   (default: `"LCA"`). Its meaning, response-data requirements, and fitted
#'   measurement parameters are the same as in \code{\link[LCPA]{LCA}()} and \code{\link[LCPA]{LPA}()},
#'   respectively.
#' @param covariates Used only when `type.analysis = "XZ"`. An
#'   \eqn{N\times(U+1)} numeric matrix/data frame representing
#'   \eqn{\boldsymbol{\zeta}_n=(1,\zeta_{n1},\ldots,\zeta_{nU})^\top}.
#'   Its first column is an unstandardized all-ones intercept; standardize the
#'   \eqn{U} observed covariates before analysis and construct interactions
#'   from the standardized variables. If `NULL` (default), the
#'   function automatically creates an \eqn{N \times 1} all-ones design, so an
#'   intercept-only class-membership model remains estimable. The argument is
#'   ignored for `"ZY"`; that path automatically uses latent-class-specific
#'   intercepts as the dependent-variable predictors.
#' @param ref.class Integer from 1 to `L` selecting the reference category in
#'   the `X -> Z` multinomial regression (default: `L`). Its coefficient vector
#'   is fixed to zero. If `control.model$is.sort = TRUE`, the value refers to the class/profile
#'   position after Step 1 sorting. This argument is not used for `"ZY"` because
#'   that path estimates a separate dependent-variable distribution for every class.
#' @param dependent.variables Used only when `type.analysis = "ZY"`. An observed
#'   dependent-variable vector, or an \eqn{N\times V} matrix/data frame
#'   representing
#'   \eqn{\mathbf{Y}_n=(Y_{n1},\ldots,Y_{nV})^\top}, with one dependent variable
#'   per column and one participant per row. These are external dependent
#'   variables (for example, depression or anxiety), not Step 1 latent
#'   class/profile indicators. Missing dependent-variable values are allowed and are
#'   omitted separately for each dependent variable; the indicator data in `response`
#'   remain subject to the requirements of \code{\link[LCPA]{LCA}()} or \code{\link[LCPA]{LPA}()}.
#'   Before calling `LCPA()`, standardize every continuous dependent-variable
#'   column assigned `family = "gaussian"` over the \eqn{N} participants,
#'   preferably with \code{\link[base]{scale}()}, so that its observed sample mean
#'   is 0 and sample standard deviation is 1. Do not standardize columns assigned
#'   `family = "categorical"`; retain their original category values.
#' @param family Used only when `type.analysis = "ZY"`. Either one character
#'   string applied to every dependent variable or a character vector with one value per
#'   dependent-variable column. `"gaussian"` (default) estimates class-specific means for
#'   a standardized numeric continuous dependent variable; `"categorical"` estimates
#'   class-specific probabilities for the observed categories. Missing values are
#'   excluded separately for each dependent variable.
#' @param method.model Character string selecting the Step 1 parameter estimator.
#'   It has exactly the same meaning and available values as the `method`
#'   argument of \code{\link[LCPA]{LCA}()} or \code{\link[LCPA]{LPA}()}: `"EM"` (default), `"NNE"`, `"Mplus"`,
#'   `"flexmix"`, `"Rmixmod"`, or `"RMixtComp"`. The corresponding
#'   `control.*` argument is passed to the selected Step 1 function. It is not
#'   used to estimate the Step 3 regression/dependent-variable model.
#' @param control.model Optional named list of common Step 1 measurement-model
#'   settings. If `NULL` (default), all settings below use their defaults. Supply
#'   only the elements to override:
#'   \describe{
#'     \item{`params`}{Optional fixed Step 1 parameter list (default: `NULL`).
#'       When supplied, Step 1 fitting is skipped. For LCA it must contain `par`,
#'       `P.Z`, and `category.levels`; for LPA it must contain `means`, `covs`,
#'       and `P.Z`, with the definitions and dimensions returned by
#'       \code{\link[LCPA]{LCA}()} and \code{\link[LCPA]{LPA}()}.}
#'     \item{`par.ini`}{Initialization used when `params = NULL` (default:
#'       `"random"`). It accepts `"random"`, `"kmeans"`, or the model-specific
#'       parameter-list forms documented for \code{\link[LCPA]{LCA}()} and
#'       \code{\link[LCPA]{LPA}()}.}
#'     \item{`constraint`}{LPA covariance structure (default: `"VV"`). It accepts
#'       `"UE"`, `"UV"`, `"E0"`, `"V0"`, `"EE"`, `"VV"`, `"VE"`, `"EV"`, or
#'       a custom equality list as documented for \code{\link[LCPA]{LPA}()}.
#'       It is ignored when `type.model = "LCA"`.}
#'     \item{`is.sort`}{Logical (default: `TRUE`). Order classes/profiles by
#'       decreasing Step 1 prior probability and consistently permute Step 1,
#'       CEP, and Step 3 class-specific results.}
#'     \item{`starts`}{Positive integer number of independently initialized Step 1
#'       warm-up analyses (default: 100).}
#'     \item{`maxiter.warmup`}{Positive integer maximum number of Step 1 iterations
#'       per warm-up analysis (default: 20).}
#'     \item{`nrep`}{Positive integer not exceeding `starts` (default: 20). The
#'       best warm-up states are continued to final Step 1 fits.}
#'   }
#'   Element names must be unique; unknown or unnamed elements are rejected.
#' @param control.EM Optional Step 1 control list for the EM estimator in
#'   \code{\link[LCPA]{LCA}()} or \code{\link[LCPA]{LPA}()}. `maxiter` sets the
#'   maximum number of EM iterations, `tol` sets the convergence tolerance, and
#'   the LPA covariance-floor element bounds small covariance eigenvalues away
#'   from zero. See `control.EM` in those functions for the complete list.
#' @param control.Mplus Optional list passed to the selected \code{\link[LCPA]{LCA}()}/\code{\link[LCPA]{LPA}()}
#'   Mplus backend. It has the same `maxiter`, `tol`, `files.path`,
#'   `files.clean`, and LPA covariance-floor semantics documented there.
#' @param control.NNE Optional list passed to the selected \code{\link[LCPA]{LCA}()}/\code{\link[LCPA]{LPA}()} NNE
#'   backend. Network architecture, attention, annealing, optimizer, plotting,
#'   and device fields have exactly the meanings documented for `control.NNE`
#'   in those functions.
#' @param control.flexmix Optional list passed to the selected \code{\link[LCPA]{LCA}()}/\code{\link[LCPA]{LPA}()}
#'   flexmix SEM backend. `maxiter`, `minprior`, `tol`, and the LPA covariance
#'   floor have the same meanings as in \code{\link[LCPA]{LCA}()}/\code{\link[LCPA]{LPA}()}.
#' @param control.Rmixmod Optional list passed to the selected \code{\link[LCPA]{LCA}()}/\code{\link[LCPA]{LPA}()}
#'   Rmixmod backend. The `path`, native strategy, algorithm, initialization,
#'   iteration, tolerance, and replication fields retain the meanings and
#'   restrictions documented in those functions.
#' @param control.RMixtComp Optional list passed to the selected \code{\link[LCPA]{LCA}()}/\code{\link[LCPA]{LPA}()}
#'   RMixtComp SEM backend. Burn-in, SEM/Gibbs iterations, stability,
#'   initialization, criterion, replication, and core controls have the same
#'   meanings and restrictions documented in those functions.
#' @param method.3step Character string selecting the Step 2--3 correction:
#'   \itemize{
#'     \item `"ML"`: maximize a likelihood that treats the modal Step 1 class
#'       as measured with error through the CEP matrix.
#'     \item `"BCH"`: invert the CEP matrix and estimate Step 3 with BCH
#'       pseudo-weights, leaving the Step 1 class/profile definition fixed.
#'   }
#'   If `NULL` (default), `"ML"` is selected for `"XZ"` and `"BCH"` for
#'   `"ZY"`. Both methods are available for both paths. ML/CEP estimates a
#'   classification-error-corrected likelihood, whereas BCH estimates
#'   inverse-CEP-weighted score equations. BCH is generally preferred for
#'   distal dependent variables, and ML is generally preferred for
#'   class-membership regression. Vermunt (2010) develops both BCH and ML
#'   corrections for covariates predicting class membership, including the
#'   BCH-XZ specification. Bakk, Tekle, and Vermunt (2013) develop the
#'   bias-adjusted ML formulation for class membership predicting a distal
#'   dependent variable; Nylund-Gibson, Grimm, and Masyn (2019) provide a
#'   worked manual ML three-step distal-outcome analysis.
#' @param CEP.error Logical. If `TRUE` (default and recommended), Step 2
#'   estimates the classification-error probability matrix
#'   \eqn{\mathrm{CEP}(l,k)=P(\widehat{Z}=k\mid Z=l)} using
#'   \code{\link[LCPA]{get.CEP}()}, and Step 3 applies the selected correction.
#'   If `FALSE`, an identity CEP matrix is used, reducing the analysis to naive
#'   modal-class assignment without classification-error correction.
#' @param method.regression Character string controlling Step 3 point estimation
#'   (default: `"Analytic"`). For ML-XZ, BCH-XZ, and ML-ZY, `"Analytic"`
#'   supplies the exact analytic gradient to a numerical optimizer; these models
#'   do not generally have closed-form coefficient estimates. For BCH-ZY,
#'   `"Analytic"` directly evaluates the weighted closed-form class means or
#'   categorical probabilities. `"Numeric"` optimizes the same likelihood
#'   without the supplied gradient or numerically solves the same BCH estimating
#'   equations. This choice is independent of `method.SE`.
#' @param maxiter Positive integer giving the maximum number of Step 3 `"XZ"`
#'   optimization iterations (default: 5000). It is distinct from Step 1
#'   `control.*$maxiter` and `control.model$maxiter.warmup`, and is ignored for `"ZY"`.
#' @param tol Positive finite convergence tolerance for the Step 3 `"XZ"`
#'   optimization (default: \eqn{10^{-4}}). It is separate from the Step 1
#'   tolerance in `control.EM`, `control.Mplus`, or another backend control and
#'   is ignored for `"ZY"`, whose dependent-variable estimators use their own equations.
#' @param lower,upper Finite lower and upper bounds applied to every `"XZ"`
#'   multinomial regression coefficient (defaults: -10 and 10). They are ignored
#'   for `"ZY"`. Inspect the returned bound diagnostics when an estimate reaches
#'   a limit.
#' @param method.SE Character string selecting Step 3 uncertainty estimation:
#'   \itemize{
#'     \item `"Analytic"`: for ML-XZ, use the analytic observed-information
#'       bread and an empirical sandwich meat that includes the influence of
#'       estimating the CEP matrix; for ML-ZY, use Louis observed information;
#'       for BCH, use the analytic estimating-equation bread and empirical
#'       sandwich meat.
#'     \item `"Numeric"`: replace the analytic information/bread with a
#'       numerical Hessian or Jacobian while retaining the corresponding
#'       empirical meat and Step 3 estimand.
#'     \item `"Bootstrap"` (default): resample individuals, keep the Step 1
#'       measurement parameters fixed, recompute posterior assignments and CEP
#'       or BCH weights, and re-estimate Step 3.
#'   }
#' @param nrep.bootstrap Integer number of nonparametric bootstrap replications
#'   used only when `method.SE = "Bootstrap"` (default: 100; minimum: 2). Only
#'   successful Step 3 replications contribute to the empirical covariance;
#'   larger values such as 500--1000 are advisable for final publication-level
#'   inference when computationally feasible.
#' @param vis Logical (default: `TRUE`). If enabled, display the Step 1
#'   measurement-model progress, Step 2 posterior/CEP preparation, and Step 3
#'   regression or distal-dependent-variable estimation progress. Each Step 2
#'   and Step 3 heading identifies the selected `X -> Z` or `Z -> Y` path. For
#'   `Z -> Y`, the output also gives the number and family of dependent-variable
#'   models, convergence and iteration information, and the selected
#'   standard-error or bootstrap progress.
#'
#' @details
#' The notation distinguishes the Step 1 indicators from the Step 3 auxiliary
#' variables. Write the indicator matrix as
#' \eqn{\mathbf{X}=(X_{ni})_{N\times I}}, where
#' \eqn{n=1,2,\ldots,N} indexes participants and
#' \eqn{i=1,2,\ldots,I} indexes observed indicators. Participant \eqn{n}'s
#' indicator vector is
#' \eqn{\mathbf{X}_n=(X_{n1},\ldots,X_{nI})^\top}, and
#' \eqn{Z_n\in\{1,2,\ldots,L\}}, with
#' \eqn{l=1,2,\ldots,L} indexing latent classes/profiles.
#'
#' The covariate vector is
#' \eqn{\boldsymbol{\zeta}_n=(1,\zeta_{n1},\ldots,\zeta_{nU})^\top}, where
#' \eqn{u=1,2,\ldots,U} indexes the \eqn{U} observed covariates and the leading
#' 1 is the intercept. The dependent-variable vector is
#' \eqn{\mathbf{Y}_n=(Y_{n1},\ldots,Y_{nV})^\top}, where
#' \eqn{v=1,2,\ldots,V} indexes the \eqn{V} external observed dependent
#' variables. Neither \eqn{\boldsymbol{\zeta}_n} nor \eqn{\mathbf{Y}_n} is part
#' of the indicator vector \eqn{\mathbf{X}_n}.
#' `"XZ"` is the function-interface label for the covariate-to-latent path;
#' the formulas use \eqn{\boldsymbol{\zeta}} for its covariates because
#' \eqn{\mathbf{X}} is reserved for the LCA/LPA indicator data.
#' `type.analysis = "XZ"` uses \eqn{\boldsymbol{\zeta}_n} but not
#' \eqn{\mathbf{Y}_n}; `type.analysis = "ZY"` uses \eqn{\mathbf{Y}_n} but not
#' \eqn{\boldsymbol{\zeta}_n}. These are separate Step 3 analyses rather than a
#' jointly estimated mediation model. Run both analyses when both the
#' covariate-to-class and class-to-dependent-variable associations are required.
#' For a ZY analysis, standardization applies to the continuous dependent variables
#' in \eqn{\mathbf{Y}_n}, not to the Step 1 indicator matrix \eqn{\mathbf{X}}.
#' Each Gaussian dependent variable must be transformed before model fitting to
#' have observed sample mean 0 and sample standard deviation 1. Consequently,
#' its class/profile-specific estimates and standard errors are expressed in
#' observed-standard-deviation units. Categorical dependent variables retain
#' their original category values.
#'
#' @section Methodology overview:
#' The cross-sectional three-step analysis proceeds as follows.
#'
#' Step 1 -- Unconditional measurement model. Fit an unconditional
#' \code{\link[LCPA]{LCA}()} or \code{\link[LCPA]{LPA}()} to `response`. Let
#' \eqn{\pi_l=P(Z_n=l)}. For LCA, the Step 1 observed-data log-likelihood is
#' \deqn{\log\mathcal{L}_{\mathrm{LCA}}=
#' \sum_{n=1}^N\log\left\{\sum_{l=1}^L\pi_l
#' \prod_{i=1}^I P(X_{ni}=x_{ni}\mid Z_n=l)\right\}.}
#' For LPA, it is
#' \deqn{\log\mathcal{L}_{\mathrm{LPA}}=
#' \sum_{n=1}^N\log\left\{\sum_{l=1}^L\pi_l
#' \mathcal{N}(\mathbf{X}_n\mid\boldsymbol{\mu}_l,
#' \boldsymbol{\Sigma}_l)\right\}.}
#' These are the likelihoods defined in
#' \code{\link[LCPA]{get.Log.Lik.LCA}()} and
#' \code{\link[LCPA]{get.Log.Lik.LPA}()}. Bayes' theorem gives
#' \deqn{\tau_{nl}=P(Z_n=l\mid\mathbf{X}_n)=
#' \frac{\pi_l\prod_{i=1}^I P(X_{ni}=x_{ni}\mid Z_n=l)}
#' {\sum_{h=1}^L\pi_h\prod_{i=1}^I
#' P(X_{ni}=x_{ni}\mid Z_n=h)}}
#' for LCA and
#' \deqn{\tau_{nl}=P(Z_n=l\mid\mathbf{X}_n)=
#' \frac{\pi_l\mathcal{N}(\mathbf{X}_n\mid\boldsymbol{\mu}_l,
#' \boldsymbol{\Sigma}_l)}
#' {\sum_{h=1}^L\pi_h\mathcal{N}(\mathbf{X}_n\mid
#' \boldsymbol{\mu}_h,\boldsymbol{\Sigma}_h)}}
#' for LPA.
#' The modal assignment is
#' \eqn{\widehat{Z}_n=\arg\max_l\tau_{nl}}. `method.model` selects the estimator of
#' this measurement model, and `control.model` supplies its initialization,
#' covariance-constraint, sorting, and replication settings. If
#' `control.model$params` is supplied, those fixed measurement parameters are
#' used to calculate \eqn{\tau_{nl}}.
#'
#' Step 2 -- Classification-error probabilities. The \eqn{L\times L} CEP matrix
#' has rows indexed by latent class/profile \eqn{l} and columns indexed by modal
#' assignment \eqn{k}, so
#' \eqn{\mathrm{CEP}(l,k)=P(\widehat{Z}_n=k\mid Z_n=l)}. The modal assignment, posterior-weight
#' estimator, matrix orientation, and pooling rules are defined in
#' \code{\link[LCPA]{get.CEP}()}. With `CEP.error = FALSE`, \eqn{\mathrm{CEP}}
#' is replaced
#' by the identity matrix and Step 3 becomes an uncorrected modal-assignment
#' analysis. The BCH and ML corrections for this classification error follow
#' Bolck, Croon, and Hagenaars (2004) and Vermunt (2010).
#'
#' Step 3A -- Covariates predicting latent membership (XZ). With reference
#' class `ref.class` denoted by \eqn{l_0}, the multinomial-logit model is
#' \deqn{P(Z_n=l\mid\boldsymbol{\zeta}_n)=
#' \frac{\exp(\boldsymbol{\zeta}_n^\top\boldsymbol{\beta}_l)}
#' {1+\sum_{h\ne l_0}\exp(\boldsymbol{\zeta}_n^\top
#' \boldsymbol{\beta}_h)},\quad l\ne l_0,}
#' with \eqn{\boldsymbol{\beta}_{l_0}=0}. If `covariates = NULL`,
#' \eqn{\boldsymbol{\zeta}_n=1} and the model
#' contains class-specific intercepts only. Vermunt's (2010) ML/CEP estimator
#' maximizes
#' \deqn{\ell_{\mathrm{ML}}(\boldsymbol{\beta})=\sum_{n=1}^N
#' \log\left\{\sum_{l=1}^L\mathrm{CEP}(l,\widehat{Z}_n)
#' P(Z_n=l\mid\boldsymbol{\zeta}_n)\right\}.}
#' For each non-reference class/profile \eqn{l\ne l_0}, the first derivative of
#' the observed-data log-likelihood with respect to the coefficient vector
#' \eqn{\boldsymbol{\beta}_l} is
#' \deqn{\frac{\partial\ell_{\mathrm{ML}}(\boldsymbol{\beta})}
#' {\partial\boldsymbol{\beta}_l}=
#' \sum_{n=1}^N\boldsymbol{\zeta}_n
#' P(Z_n=l\mid\boldsymbol{\zeta}_n)
#' \left\{\frac{\mathrm{CEP}(l,\widehat{Z}_n)}
#' {\sum_{h=1}^L\mathrm{CEP}(h,\widehat{Z}_n)
#' P(Z_n=h\mid\boldsymbol{\zeta}_n)}-1\right\}.}
#' This derivative is an \eqn{(U+1)\times 1} vector: its entries correspond to
#' the intercept and the \eqn{U} covariate coefficients in
#' \eqn{\boldsymbol{\beta}_l}. At an interior maximum, the ML estimates jointly
#' satisfy \eqn{\partial\ell_{\mathrm{ML}}/
#' \partial\boldsymbol{\beta}_l=\mathbf{0}} for every \eqn{l\ne l_0}; `lower`
#' and `upper` define the permitted coefficient range.
#'
#' Vermunt's (2010) BCH-XZ estimator instead solves, for class
#' \eqn{l\ne l_0},
#' \deqn{\sum_{n=1}^N\boldsymbol{\zeta}_n
#' \left\{(\mathrm{CEP}^{-1})_{\widehat{Z}_n,l}
#' -P(Z_n=l\mid\boldsymbol{\zeta}_n)
#' \sum_{h=1}^L(\mathrm{CEP}^{-1})_{\widehat{Z}_n,h}
#' \right\}=\mathbf{0}.}
#' This is a BCH estimating equation, not the derivative of the ML corrected
#' likelihood above. Its \eqn{U+1} equations correspond to the intercept and
#' covariate coefficients in \eqn{\boldsymbol{\beta}_l}.
#' Thus ML and BCH estimate the same multinomial-logit parameters but use
#' different corrections for modal-classification error.
#'
#' Step 3B -- Latent membership predicting dependent variables (ZY). The
#' bias-adjusted ML three-step formulation follows Bakk, Tekle, and Vermunt
#' (2013) and Nylund-Gibson, Grimm, and Masyn (2019), with the BCH
#' secondary-model formulation described by Asparouhov and Muthén (2014b). The model
#' estimates the conditional distribution of \eqn{Y_{nv}} given \eqn{Z_n=l}
#' separately for \eqn{v=1,\ldots,V}. For `family = "gaussian"`,
#' \eqn{Y_{nv}\mid Z_n=l\sim N(\mu_{lv},\sigma_{lv}^2)}. For
#' `family = "categorical"`,
#' \eqn{P(Y_{nv}=q\mid Z_n=l)=p_{lvq}}, where \eqn{q} indexes the observed
#' categories of dependent variable \eqn{v}. No design matrix is required because the
#' model contains a separate intercept for every class/profile.
#'
#' ML/CEP maximizes
#' \deqn{\ell_{\mathrm{ML},v}=
#' \begin{cases}
#' \sum_{n=1}^N\log\left\{\sum_{l=1}^L
#' \mathrm{CEP}(l,\widehat{Z}_n)\pi_l
#' \mathcal{N}(Y_{nv}\mid\mu_{lv},\sigma_{lv}^2)\right\},
#' & \text{for a Gaussian dependent variable},\\
#' \sum_{n=1}^N\log\left\{\sum_{l=1}^L
#' \mathrm{CEP}(l,\widehat{Z}_n)\pi_l
#' \prod_q p_{lvq}^{\mathbb{1}(Y_{nv}=q)}\right\},
#' & \text{for a categorical dependent variable},
#' \end{cases}}
#' The reported
#' class-shift rate compares the Step 3 modal class with
#' \eqn{\widehat{Z}_n} from Step 1.
#'
#' BCH Gaussian means solve
#' \deqn{\sum_{n=1}^N
#' (\mathrm{CEP}^{-1})_{\widehat{Z}_n,l}
#' (Y_{nv}-\mu_{lv})=0,}
#' and the corresponding Gaussian variances solve
#' \deqn{\sum_{n=1}^N
#' (\mathrm{CEP}^{-1})_{\widehat{Z}_n,l}
#' \{(Y_{nv}-\mu_{lv})^2-\sigma_{lv}^2\}=0.}
#' and categorical probabilities solve
#' \deqn{\sum_{n=1}^N
#' (\mathrm{CEP}^{-1})_{\widehat{Z}_n,l}
#' \{\mathbb{1}(Y_{nv}=q)-p_{lvq}\}=0.}
#' @section Parameter estimation and uncertainty:
#' With `method.regression = "Analytic"`, ML-XZ, BCH-XZ, and ML-ZY use their
#' exact scores or gradients within numerical optimization; BCH-ZY evaluates
#' the closed-form weighted estimates shown above. `"Numeric"` evaluates the
#' same likelihoods without supplied gradients or minimizes the squared BCH-ZY
#' estimating equations numerically.
#'
#' For ML-XZ and BCH, the covariance has sandwich form
#' \eqn{A^{-1}BA^{-\top}}. In ML-XZ, \eqn{A} is the observed information and
#' \eqn{B} is formed from individual likelihood scores plus the influence
#' function of the estimated CEP matrix. In BCH, \eqn{A} is the
#' estimating-equation Jacobian and \eqn{B} is the empirical covariance of
#' individual estimating-function contributions. For BCH-ZY Gaussian models,
#' the mean and variance equations are stacked so their sandwich covariance
#' includes the covariance between \eqn{\widehat{\mu}_{lv}} and
#' \eqn{\widehat{\sigma}_{lv}^2}. ML-ZY uses the inverse Louis
#' observed-information matrix; the variance standard error follows by applying
#' the delta method to the fitted log-standard-deviation parameter.
#' `method.SE = "Numeric"` evaluates the required
#' Hessian or Jacobian numerically. `"Bootstrap"` resamples individuals, recalculates posterior
#' assignments and CEP/BCH weights, and re-estimates Step 3 while holding the
#' Step 1 measurement parameters fixed. For a Gaussian dependent variable,
#' `estimate` and `se` report \eqn{\mu_{lv}} and its standard error, whereas
#' `variance` and `variance.se` report \eqn{\sigma_{lv}^2} and its standard
#' error. Separate omnibus Wald tests assess equality of the conditional means
#' and equality of the conditional variances across classes/profiles. For a
#' categorical dependent variable, `estimate` and `se` report every
#' class/profile-specific category probability and its standard error; the
#' omnibus Wald test assesses equality of the complete conditional category
#' distributions.
#'
#' @section Method selection:
#' ML is the default for XZ because its likelihood directly represents the
#' error-prone modal assignment through \eqn{\mathrm{CEP}}. BCH is the default for ZY
#' because its weights are calculated without using \eqn{\mathbf{Y}}, so the Step 1
#' class/profile definition is not changed by the dependent variable. ML-ZY
#' provides a corrected-likelihood sensitivity analysis and reports class
#' shifts. Use a dedicated DCAT procedure when the DCAT estimand is required.
#'
#' @return An object of class `"LCPA"`. The selected result is available from
#'   `analysis$XZ` or `analysis$ZY`; posterior probabilities, modal assignments,
#'   and CEP matrices use the same list structure as \code{\link[LCPA]{LTA}()}.
#'   For `type.analysis = "ZY"`, `dependent.variables$t1` contains one fitted
#'   model per observed dependent variable. A Gaussian model reports
#'   class/profile-specific `estimate`, `se`, `variance`, `variance.se`, their
#'   covariance matrices, omnibus Wald tests for both means and variances,
#'   group weight masses, observations, omitted values, iterations, and
#'   convergence. A categorical model reports class/profile-by-category
#'   `estimate` and `se` matrices, their covariance matrix, an omnibus Wald
#'   test of equality of the conditional category distributions, group weight
#'   masses, observations, omitted values, iterations, and convergence.
#'
#' @examples
#' \donttest{
#' library(LCPA)
#'
#' set.seed(1245)
#' N <- 2000
#' L <- 3
#' I <- 6
#'
#' # Two observed covariates plus the required intercept
#' covariates <- cbind(
#'   Intercept = 1,
#'   Zeta.1 = as.numeric(scale(rnorm(N))),
#'   Zeta.2 = rbinom(N, 1, 0.5)
#' )
#' beta <- matrix(c(
#'    0.70,  0.30, 0,
#'    0.40, -0.20, 0,
#'   -0.30,  0.30, 0
#' ), ncol = L, byrow = TRUE)
#' rownames(beta) <- colnames(covariates)
#'
#' data.LCPA <- sim.LTA(
#'   N = N, I = I, L = L, times = 1, type = "LPA",
#'   constraint = "VE", mean.range = c(-3, 3),
#'   covs.range = c(0.4, 0.8),
#'   covariates = list(covariates), ref.class = 3,
#'   beta = beta, is.sort = TRUE
#' )
#' control.model <- list(
#'   constraint = "VE", is.sort = TRUE,
#'   starts = 10, maxiter.warmup = 10, nrep = 3
#' )
#'
#' # Covariates predicting latent profiles: XZ analysis
#' fit.LCPA.XZ <- LCPA(
#'   response = data.LCPA$responses[[1]], L = L,
#'   type.analysis = "XZ", type.model = "LPA",
#'   covariates = covariates, ref.class = 3,
#'   method.model = "EM", control.model = control.model,
#'   method.3step = "ML", method.regression = "Analytic",
#'   method.SE = "Analytic", maxiter = 500, vis = TRUE
#' )
#' round(cbind(
#'   "True Class 1" = beta[, 1],
#'   "Estimate Class 1" = fit.LCPA.XZ$beta[, 1],
#'   "True Class 2" = beta[, 2],
#'   "Estimate Class 2" = fit.LCPA.XZ$beta[, 2]
#' ), 3)
#'
#' # Latent profiles predicting two dependent variables: ZY analysis
#' true.mean <- rbind(
#'   "Class 1" = c(Depression = 8, Anxiety = 12),
#'   "Class 2" = c(Depression = 10, Anxiety = 10),
#'   "Class 3" = c(Depression = 13, Anxiety = 8)
#' )
#' dependent.variables <- scale(
#'   true.mean[data.LCPA$Zs[[1]], ] +
#'     matrix(rnorm(N * 2, sd = 1.5), N, 2)
#' )
#' true.mean.standardized <- sweep(
#'   sweep(true.mean, 2, attr(dependent.variables, "scaled:center"), "-"),
#'   2, attr(dependent.variables, "scaled:scale"), "/"
#' )
#' true.variance.standardized <-
#'   (1.5 / attr(dependent.variables, "scaled:scale"))^2
#' dependent.variables <- as.data.frame(dependent.variables)
#' fit.LCPA.ZY <- LCPA(
#'   response = data.LCPA$responses[[1]], L = L,
#'   type.analysis = "ZY", type.model = "LPA",
#'   dependent.variables = dependent.variables,
#'   family = "gaussian",
#'   method.model = "EM", control.model = control.model,
#'   method.3step = "BCH", method.regression = "Analytic",
#'   method.SE = "Analytic", vis = TRUE
#' )
#' round(cbind(
#'   True.Depression = true.mean.standardized[, "Depression"],
#'   Estimate.Depression =
#'     fit.LCPA.ZY$dependent.variables$t1$Depression$estimate,
#'   True.Anxiety = true.mean.standardized[, "Anxiety"],
#'   Estimate.Anxiety =
#'     fit.LCPA.ZY$dependent.variables$t1$Anxiety$estimate
#' ), 3)
#' round(cbind(
#'   True.Variance.Depression = rep(
#'     true.variance.standardized["Depression"], L
#'   ),
#'   Estimate.Variance.Depression =
#'     fit.LCPA.ZY$dependent.variables$t1$Depression$variance,
#'   True.Variance.Anxiety = rep(
#'     true.variance.standardized["Anxiety"], L
#'   ),
#'   Estimate.Variance.Anxiety =
#'     fit.LCPA.ZY$dependent.variables$t1$Anxiety$variance
#' ), 3)
#' }
#'
#' @references
#' Asparouhov, T., & Muthén, B. (2014a). Auxiliary variables in mixture
#' modeling: Three-step approaches using Mplus. *Structural Equation Modeling:
#' A Multidisciplinary Journal, 21*(3), 329--341.
#' \doi{10.1080/10705511.2014.915181}
#'
#' Asparouhov, T., & Muthén, B. (2014b). *Auxiliary variables in mixture
#' modeling: Using the BCH method in Mplus to estimate a distal outcome model
#' and an arbitrary secondary model* (Mplus Web Note No. 21, Version 2).
#' \url{https://www.statmodel.com/examples/webnotes/webnote21.pdf}
#'
#' Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the association
#' between latent class membership and external variables using bias-adjusted
#' three-step approaches. *Sociological Methodology, 43*(1), 272--311.
#' \doi{10.1177/0081175012470644}
#'
#' Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent structure
#' models with categorical variables: One-step versus three-step estimators.
#' *Political Analysis, 12*(1), 3--27.
#' \doi{10.1093/pan/mph001}
#'
#' Nylund-Gibson, K., Grimm, R. P., & Masyn, K. E. (2019). Prediction from
#' latent classes: A demonstration of different approaches to include distal
#' outcomes in mixture models. *Structural Equation Modeling: A
#' Multidisciplinary Journal, 26*(6), 967--985.
#' \doi{10.1080/10705511.2019.1590146}
#'
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved
#' three-step approaches. *Political Analysis, 18*(4), 450--469.
#' \doi{10.1093/pan/mpq025}
#'
#' @seealso \code{\link[LCPA]{LCA}()}, \code{\link[LCPA]{LPA}()}, \code{\link[LCPA]{LTA}()}, \code{\link[LCPA]{get.CEP}()}
#'
#' @export
LCPA <- function(response, L = 2,
                 type.analysis = c("XZ", "ZY"),
                 type.model = c("LCA", "LPA"),
                 covariates = NULL,
                 ref.class = L,
                 dependent.variables = NULL,
                 family = "gaussian",
                 method.model = "EM",
                 control.model = NULL,
                 control.EM = NULL,
                 control.Mplus = NULL,
                 control.NNE = NULL,
                 control.flexmix = NULL,
                 control.Rmixmod = NULL,
                 control.RMixtComp = NULL,
                 method.3step = NULL,
                 CEP.error = TRUE,
                 method.regression = "Analytic",
                 maxiter = 5000, tol = 1e-4,
                 lower = -10, upper = 10,
                 method.SE = "Bootstrap", nrep.bootstrap = 100,
                 vis = TRUE) {

  call <- match.call()
  type.analysis <- match.arg(type.analysis)
  type.model <- match.arg(type.model)
  control.model <- .three.step.control.model(control.model)
  params <- control.model$params
  constraint <- control.model$constraint
  is.sort <- control.model$is.sort
  par.ini <- control.model$par.ini
  starts <- control.model$starts
  maxiter.warmup <- control.model$maxiter.warmup
  nrep <- control.model$nrep
  if(is.null(method.3step)){
    method.3step <- if(type.analysis == "XZ") "ML" else "BCH"
  }else{
    method.3step <- match.arg(method.3step, c("ML", "BCH"))
  }
  method.regression <- match.arg(method.regression, c("Analytic", "Numeric"))
  method.SE <- match.arg(method.SE, c("Analytic", "Numeric", "Bootstrap"))

  if(type.analysis == "XZ"){
    res <- XZ.LCPA(
      response = response, L = L,
      ref.class = ref.class, type.model = type.model,
      covariates = covariates, CEP.error = CEP.error,
      par.ini = par.ini, params = params, is.sort = is.sort,
      constraint = constraint,
      method.model = method.model, tol = tol,
      method.regression = method.regression,
      lower = lower, upper = upper,
      method.SE = method.SE, nrep.bootstrap = nrep.bootstrap,
      maxiter = maxiter, starts = starts,
      maxiter.warmup = maxiter.warmup, nrep = nrep,
      vis = vis,
      control.EM = control.EM, control.Mplus = control.Mplus,
      control.NNE = control.NNE, control.flexmix = control.flexmix,
      control.Rmixmod = control.Rmixmod,
      control.RMixtComp = control.RMixtComp,
      method.3step = method.3step
    )
    res$P.Z.Xns <- list(res$P.Z.Xn)
    res$P.Zs <- list(res$P.Z)
    res$Zs <- list(res$Z)
    if(is.null(res$CEP)){
      res$CEP <- if(CEP.error){
        get.CEP(res$P.Z.Xns, CEP.time.cross = FALSE)
      }else list(diag(L))
    }
    res$XZ <- list(
      beta = res$beta,
      beta.se = res$beta.se,
      beta.Z.sta = res$beta.Z.sta,
      beta.p.value.tail1 = res$beta.p.value.tail1,
      beta.p.value.tail2 = res$beta.p.value.tail2,
      vcov = res$vcov,
      information = res$information,
      SE.diagnostics = res$SE.diagnostics
    )
    res$ZY <- NULL
  }else{
    if(is.null(dependent.variables)){
      stop("dependent.variables must be supplied when type.analysis = 'ZY'")
    }
    res <- ZY.LCPA(
      response = response, dependent.variables = dependent.variables, L = L,
      type.model = type.model, family = family,
      CEP.error = CEP.error,
      par.ini = par.ini, params = params, is.sort = is.sort,
      constraint = constraint, method.model = method.model,
      method.regression = method.regression,
      method.SE = method.SE, nrep.bootstrap = nrep.bootstrap,
      starts = starts, maxiter.warmup = maxiter.warmup, nrep = nrep,
      vis = vis,
      control.EM = control.EM, control.Mplus = control.Mplus,
      control.NNE = control.NNE, control.flexmix = control.flexmix,
      control.Rmixmod = control.Rmixmod,
      control.RMixtComp = control.RMixtComp,
      method.3step = method.3step
    )
    res$XZ <- NULL
    res$ZY <- res$dependent.variables
  }

  res$type.analysis <- type.analysis
  res$type.model <- type.model
  res$method.3step <- method.3step
  res$analysis <- list(type = type.analysis, XZ = res$XZ, ZY = res$ZY)
  dependent.variables.stored <- if(type.analysis == "ZY"){
    res$arguments$dependent.variables
  }else dependent.variables
  covariates.stored <- if(type.analysis == "XZ"){
    res$arguments$covariates
  }else covariates
  family.stored <- if(type.analysis == "ZY") res$arguments$family else family
  res$arguments <- list(
    response = response,
    L = L,
    type.analysis = type.analysis,
    type.model = type.model,
    covariates = covariates.stored,
    ref.class = ref.class,
    dependent.variables = dependent.variables.stored,
    family = family.stored,
    method.model = method.model,
    control.model = control.model,
    control.EM = control.EM,
    control.Mplus = control.Mplus,
    control.NNE = control.NNE,
    control.flexmix = control.flexmix,
    control.Rmixmod = control.Rmixmod,
    control.RMixtComp = control.RMixtComp,
    method.3step = method.3step,
    CEP.error = CEP.error,
    method.regression = method.regression,
    maxiter = maxiter,
    tol = tol,
    lower = lower,
    upper = upper,
    method.SE = method.SE,
    nrep.bootstrap = nrep.bootstrap,
    vis = vis
  )
  res$call <- call
  class(res) <- "LCPA"
  res
}
