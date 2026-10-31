# Burr XII inference with doubly Type-II censored fuzzy observations

R code for a Monte Carlo simulation study of classical and Bayesian estimation of the two shape parameters of the Burr XII distribution under doubly Type-II censoring and triangular fuzzy observations.

**Author:** Dr. Riyadh Rustam Al-Mosawi, Department of Mathematics, University of Thi-Qar, Nasiriyah, Iraq.

## Repository files

| File | Purpose |
| --- | --- |
| `BFuz_DT2C_Sim_R1.R` | Main simulation script: selects a setting, generates samples, runs estimation methods, summarizes performance, and saves results. |
| `Library.r` | User-defined functions for the model, fuzzy likelihood, estimation, simulation settings, and reporting. |
| `README.md` | Setup instructions, function reference, and output documentation. |

Keep the two R files in the same directory. Filenames are case-sensitive on some systems; the main script uses `source("Library.r")`.

## Statistical model

For x > 0 and alpha, beta > 0, the model is

```math
f(x;\alpha,\beta)=\alpha\beta x^{\alpha-1}(1+x^\alpha)^{-(\beta+1)},
\qquad F(x;\alpha,\beta)=1-(1+x^\alpha)^{-\beta}.
```

The censoring convention used by the code is:

- `n`: original sample size.
- `r`: number of left-censored observations.
- `s`: index of the last retained order statistic.
- `s - r`: number of retained observations, stored as the sorted vector `E`.
- `n - s`: number of right-censored observations.

Each retained value `g` is assigned a triangular membership function with support `[(1-f1)*g, (1+f2)*g]` and peak at `g`. The defaults `f1 = f2 = 0.05` give 5% uncertainty on each side. Use positive fuzzification widths; setting them to zero does not produce an exact-data implementation because the membership functions divide by these widths.

`LogLike()` integrates the Burr XII density against each membership function and adds the censoring contributions at the first and last retained centers. Up to a parameter-independent factor, its likelihood is

```math
L(\alpha,\beta)\propto F(E_1)^r
\prod_{j=1}^{m}\left[\int f(w;\alpha,\beta)\mu_j(w)\,dw\right]
[1-F(E_m)]^{n-r-m},\qquad m=s-r.
```

## Estimation methods

| Output label | Method |
| --- | --- |
| `MLE` / `NR_MLE` | Direct likelihood maximization with `optim(method = "L-BFGS-B")`; the historical function name is `NR.fun`. |
| `EM` / `EM_MLE` | EM with numerical integration in the E-step and numerical optimization in the M-step. |
| `MCEM` / `SM_MLE` | Monte Carlo EM with membership-weighted draws for fuzzy observations and conditional draws for censored observations. Retains the iterate with the largest observed log-likelihood. |
| `LN.NON`, `LN.INF` | Lindley approximations to posterior means. |
| `TK.NON`, `TK.INF` | Tierney–Kadane approximations to posterior means. |
| `MH.NON`, `MH.INF` | Metropolis–Hastings posterior means and 95% HPD intervals. |

Independent Gamma priors use the **shape–rate** parameterization: alpha ~ Gamma(a,b), beta ~ Gamma(c,d). `NON` denotes the proper Gamma(0.001, 0.001) priors used by the code; it is a naming convention, not an improper flat prior. `INF` uses the setting-specific hyperparameters below.

The MH sampler proposes a symmetric normal random walk in `(log(alpha), log(beta))`. Its log target is `LogPostDist(exp(z), ...) + sum(z)`, including the transformation Jacobian. It converts the original-scale MLE covariance to the log scale using the delta method. Each call runs one chain for each prior specification; these are not multiple chains for the same posterior.

## Requirements

Install R and, optionally, RStudio. Required packages are:

```r
install.packages(c("numDeriv", "coda", "knitr"),
                 repos = "https://cloud.r-project.org")
```

- `numDeriv`: numerical derivatives for information matrices and Lindley approximation.
- `coda`: MCMC objects and HPD intervals.
- `knitr`: loaded by the main script; used by the earlier, overridden export function. The active LaTeX exporter writes the table directly.

The main script also installs missing packages automatically. First-time installation requires internet access and a writable R package library. No external data file is required.

## How to run the main script

1. Download or clone the repository and keep both R files together.
2. Open `BFuz_DT2C_Sim_R1.R`. Replace its existing computer-specific `setwd(...)` line with the absolute path to your local repository, for example:

   ```r
   setwd("C:/Users/YourName/Documents/BurrXII-Fuzzy-DT2C")
   source("Library.r")
   ```

   Alternatively, remove the `setwd(...)` line, open an RStudio project in the repository directory, and retain `source("Library.r")`. Merely changing the console working directory will not override the existing `setwd(...)` line inside the script.

3. Set `index` and the simulation controls near the beginning of the main file.
4. Run the entire main script in RStudio or execute:

   ```r
   source("BFuz_DT2C_Sim_R1.R")
   ```

   From a terminal in the repository directory, after correcting or removing the hard-coded path:

   ```bash
   Rscript BFuz_DT2C_Sim_R1.R
   ```

The script runs **one setting per execution**. Change `index` and rerun for other settings. It prints replication progress, usable estimates, summaries, and elapsed time. Numerical integration, numerical third derivatives, and two MCMC chains per successful replication can make the full simulation time-consuming.

### Default controls

| Variable | Default | Meaning |
| --- | --- | --- |
| `index` | `7` | Selected simulation setting. |
| `alp_ini`, `bet_ini` | `0.5`, `0.5` | Initial values for numerical estimation. |
| `f1`, `f2` | `0.05`, `0.05` | Relative left and right membership widths. |
| `Min_Dif` | `0.0001` | EM stopping tolerance for the Euclidean parameter change. |
| `Max_Iter` | `30` | Maximum EM iterations. |
| `N` | `5100` | MH iterations per prior per replication. |
| `P` | `100` | Discarded MH burn-in iterations; 5,000 draws remain by default. |
| `M` | `100` | MCEM iterations; distinct from the MH control `N`. |
| `K` | `5` | Monte Carlo draws per fuzzy observation or censored conditional component per MCEM iteration. |
| `iterations` | `500` | Number of simulation replications. |
| `set.seed(...)` | `2021` | Random seed set before the simulation loop. |

For an initial execution check, change `iterations = 2`, `N = 300`, `P = 100`, and `M = 5` **inside the main file**, before running it. These reduced settings are only for checking execution, not for reported inference. Values assigned in the console before sourcing are overwritten by the script.

### Available settings

| Index, first parameter pair | Index, second parameter pair | n | r | s | Retained | Right-censored |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 10 | 25 | 0 | 15 | 15 | 10 |
| 2 | 11 | 25 | 5 | 15 | 10 | 10 |
| 3 | 12 | 25 | 10 | 15 | 5 | 10 |
| 4 | 13 | 50 | 0 | 30 | 30 | 20 |
| 5 | 14 | 50 | 10 | 30 | 20 | 20 |
| 6 | 15 | 50 | 20 | 30 | 10 | 20 |
| 7 | 16 | 100 | 0 | 60 | 60 | 40 |
| 8 | 17 | 100 | 20 | 60 | 40 | 40 |
| 9 | 18 | 100 | 40 | 60 | 20 | 40 |

| Indices | True alpha | True beta | a | b | c | d |
| --- | --- | --- | --- | --- | --- | --- |
| 1–9 | 0.5 | 1.5 | 1 | 2 | 3 | 2 |
| 10–18 | 1.5 | 1 | 3 | 2 | 3 | 3 |

The informative prior means equal the corresponding true parameter values. The default setting 7 uses `(n,r,s) = (100,0,60)` and `(alpha,beta) = (0.5,1.5)`.

## Function reference: Library.r

### Model, sampling, and likelihood

| Function | Purpose |
| --- | --- |
| `pdf(x, alp, bet)` | Evaluates the Burr XII density. |
| `cdf(x, alp, bet)` | Evaluates the Burr XII cumulative distribution function. |
| `qnt(p, alp, bet)` | Evaluates the inverse CDF. |
| `membership_fun(g)` | Returns the two triangular membership arms and their integration limits; uses global `f1` and `f2`. Its nested `memb.fun1()` and `memb.fun2()` evaluate the increasing and decreasing arms on their respective supports. |
| `Gen.Data(n, r, s, alp, bet)` | Generates the retained doubly Type-II censored sample through a beta/uniform order-statistic construction; returns a data frame with column `E`. |
| `LogLike(y, E)` | Evaluates the fuzzy observed-data log-likelihood at `y = c(alpha,beta)`, including left/right censoring. |
| `Like(y, E)` | Exponentiates `LogLike()` to obtain the likelihood. |
| `LogCompLike(y, E)` | Computes the expected complete-data log-likelihood for EM; `y` contains proposed alpha/beta followed by current alpha/beta. Local `logpdf()` and `oldpdf()` evaluate the proposed log-density and current density. |
| `CompLike(y, E)` | Exponentiates `LogCompLike()`. |

### Priors, posterior, and derivatives

| Function | Purpose |
| --- | --- |
| `Prior1(x, a, b)` | Gamma prior density for alpha. |
| `Prior2(x, c, d)` | Gamma prior density for beta. |
| `Rho(y, a, b, c, d)` | Joint log-prior density. |
| `LogPostDist(x, E, a, b, c, d)` | Log-likelihood plus independent Gamma log-priors, up to a parameter-independent constant. |
| `PostDist(x, E, a, b, c, d)` | Exponentiates the log-posterior. |
| `H(y, E, a, b, c, d)` | Log-posterior divided by `n`, used in Tierney–Kadane approximation. |
| `Ha(y, E, a, b, c, d)` | Adds `log(alpha)/n` to `H()` for the posterior mean of alpha. |
| `Hb(y, E, a, b, c, d)` | Adds `log(beta)/n` to `H()` for the posterior mean of beta. |
| `dLogLike(y, E)` | Numerical gradient of the observed-data log-likelihood. |
| `ddLogLike(y, E)` | Numerical Hessian of the observed-data log-likelihood. |

### Estimation

| Function | Purpose and returned information |
| --- | --- |
| `NR.fun(alp, bet, E)` | Bounded direct likelihood optimization; returns MLEs, observed-information covariance, standard errors, log-transformed 95% confidence intervals, coverage indicators, and `Err_NR`. |
| `EM.fun(alp, bet, E, max_iter, tol)` | Numerical EM; returns converged estimates, standard errors, 95% confidence intervals, coverage, and `Err_EM`. |
| `SM.fun(alp, bet, E)` | MCEM using global `M` and `K`; returns estimates, standard errors, 95% confidence intervals, coverage, and `Err.SM`. Nested `Q()` evaluates the fixed Monte Carlo EM objective and `logpdf()` evaluates its proposed complete-data log-density. |
| `TK.fun(alp, bet, E, a, b, c, d)` | Tierney–Kadane posterior mean approximations and `Err_TK`; no intervals are computed. Nested `fit_TK()`, `valid_fit()`, and `logdet()` optimize objectives, validate fits, and compute positive-definite log-determinants. |
| `LN.fun(alp, bet, E, a, b, c, d)` | Lindley posterior mean approximations and `Err_LN`; inputs must be the observed-data MLEs. Local `loglike()` and `Hess()` support numerical second and third derivatives. No intervals are computed. |
| `MH.fun(alp, bet, E, varcov, N, P, a, b, c, d, verbose)` | Runs the two prior-specific log-scale MH chains; returns `MH.Res_non` and `MH.Res_inf`, each containing posterior means, posterior SDs, 95% HPD intervals, coverage, acceptance rate, and retained `MH_sample`. Nested `run_chain()` runs one prior-specific chain and `log_target()` includes the Jacobian. |

For estimation functions with error flags, `FALSE` indicates successful completion of the required checks. A point estimate can sometimes be present even when covariance estimation fails and the flag remains `TRUE`.

### Settings and reporting

| Function | Purpose |
| --- | --- |
| `Select.Set(index)` | Returns the sample/censoring setting, true parameters, and informative prior hyperparameters for an index from 1 to 18. |
| `safe_mean(x)` | Computes a mean after excluding nonfinite values; returns `NA` if no usable values remain. |
| `parameter_summary(x, truth, parameter)` | Summarizes one parameter for one method, including successful counts, bias, absolute estimated bias, MSE, RMSE, interval endpoints, interval length, and coverage. |
| `Summary.Results(...)` | Combines alpha/beta summaries for the nine method/prior combinations. |
| `Pres.fun(...)` | Prints a compact summary; displays coverage as a percentage. |
| `Save.Results(...)` | The final definition saves raw matrices, updates the cumulative summary for the current sample size/true parameter pair, and writes a formatted LaTeX table. Nested `number()` and `cells()` format numeric entries and parameter-specific table cells. |

## Helper functions in the main script

| Function | Purpose |
| --- | --- |
| `load_or_install(pkg)` | Installs a missing package from CRAN, then loads it. |
| `empty_results(mh = FALSE)` | Allocates an `NA`-filled replication-by-statistic matrix; adds an acceptance column for MH. |
| `safe_run(expr)` | Catches errors, prints the message, and returns `NULL` so the simulation can continue. |
| `successful(result, flag)` | Checks that an estimator returned a list with its error flag exactly `FALSE`. |
| `extract_results(result, prefix, bayes = FALSE)` | Extracts estimates, SDs, interval endpoints, and coverage into the common ten-column format. |
| `extract_points(result, prefix)` | Stores Lindley/Tierney–Kadane point estimates and leaves interval-related columns as `NA`. |

## Saved results and interpretation

The **last** `Save.Results()` definition in `Library.r` overrides the earlier definition. The active outputs, written to the working directory, are:

| File pattern | Contents |
| --- | --- |
| `Results_n{n}_alpha{alpha}_beta{beta}.rds` | Cumulative full-precision summary data frame for available censoring schemes with this sample size and true parameter pair. |
| `Results_n{n}_alpha{alpha}_beta{beta}_r{r}_s{s}.RData` | Nine raw replication-level result matrices for the current censoring scheme. |
| `Results_n{n}_alpha{alpha}_beta{beta}.tex` | Formatted LaTeX table combining the saved censoring schemes. |

For the default setting, the prefix is `Results_n100_alpha0.5_beta1.5` and the raw-results suffix is `_r0_s60.RData`. Rerunning a scheme replaces its previous summary rows and overwrites its raw-results file. Use a separate output directory when varying controls or prior settings that are not encoded in the filename.

The earlier `Simulation_*.csv` exporter is overridden and does not run when the complete function library is sourced. The main script calls the active save function twice; the second call updates the same outputs.

Read exported data with:

```r
summary <- readRDS("Results_n100_alpha0.5_beta1.5.rds")
print(summary)

# Load raw matrices into a separate environment.
raw <- new.env()
load("Results_n100_alpha0.5_beta1.5_r0_s60.RData", envir = raw)
head(raw$NR_MLE)
```

To include a generated table in a LaTeX manuscript:

```latex
% In the preamble:
\usepackage{graphicx}
\usepackage{multirow}

% In the document body:
\input{Results_n100_alpha0.5_beta1.5.tex}
```

Each raw matrix has ten columns, in this order:

```text
alp, sd.alp, lower.alp, upper.alp, coverage.alp,
bet, sd.bet, lower.bet, upper.bet, coverage.bet
```

The MH matrices have an additional `acceptance` column. Classical `sd` entries are estimated standard errors; MH `sd` entries are posterior standard deviations. Lindley and Tierney–Kadane matrices contain point estimates only.

Summary definitions:

- `N_est`, `N_interval`: numbers of usable point estimates and intervals, respectively.
- `Bias`: mean estimate minus the true parameter.
- `AEB`: absolute value of the estimated bias, **not** the mean absolute error.
- `MSE`, `RMSE`: mean squared error and its square root.
- `Mean_lower`, `Mean_upper`: separately averaged valid interval endpoints; together they summarize intervals across replications rather than forming an interval from pooled data.
- `Mean_length` / `LEN`: average valid interval length.
- `ECP`: empirical coverage among valid intervals; stored as a proportion, though `Pres.fun()` prints a percentage.

Failures and skipped methods remain `NA`. Summary denominators use available finite results rather than automatically using all requested replications. Review the usable counts when interpreting comparisons. The raw MH chains are returned by `MH.fun()` but are not included in the nine saved simulation matrices; only their scalar summaries are stored there.

## Implementation notes and troubleshooting

- These scripts use global settings. Sourcing `Library.r` defines functions but does not initialize a full simulation. The main script sets `n`, `r`, `s`, `alp_true`, `bet_true`, `f1`, `f2`, the initial values, and EM/MCEM controls before estimation.
- The estimators use `alp_true` and `bet_true` to calculate simulation coverage. Adapting them to real observations requires handling these simulation-only fields.
- If `Library.r` cannot be found, check the working directory and filename capitalization.
- If EM fails to converge within `Max_Iter`, that replication is recorded as unsuccessful. Inspect numerical messages and consider a larger iteration limit or different initial values.
- MH is skipped when the direct MLE fit does not pass its covariance checks. Lindley can run when positive MLE point estimates are available even if the covariance check failed.
- The main simulation saves acceptance rates but does not automatically report ESS, Monte Carlo standard errors, multiple-chain R-hat, or diagnostic plots. Evaluate convergence separately before using MCMC summaries in final reporting.
- The saved files do not capture every simulation control or package version. For reproducibility, retain the exact scripts, chosen index, controls, and `sessionInfo()` output alongside the results.

## Author and citation

Developed by Dr. Riyadh Rustam Al-Mosawi. If you use the code in research, acknowledge the author and cite the associated Burr XII study when its final bibliographic details are available.

A reuse license has not been specified in these files. Add a `LICENSE` file reflecting the author's chosen terms before distributing the project under an open-source license.
