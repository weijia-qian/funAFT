
<!-- README.md is generated from README.Rmd. Please edit README.Rmd -->

# funAFT

An R package for **functional accelerated failure time (AFT) models**:

- **Linear functional AFT (lfAFT)** with a functional predictor entering
  as  
  $$
    \int X(s)\,\beta(s)\,ds
  $$ and estimated via penalized splines.
- **Additive functional AFT (afAFT)** where the functional predictor
  enters as  
  $$
    \int F\{X(s), s\}\,ds
  $$ with an unknown bivariate function $F(\cdot,\cdot)$, fitted using
  tensor-product smooths in **mgcv**.

The package currently supports:

- lfAFT: lognormal and loglogistic AFT models.  
- afAFT: **lognormal AFT family only**.

## Installation

You can install the development version of funAFT from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("weijia-qian/funAFT")
```

## Data structure

Both lfAFT and afAFT assume:

- `y`: positive survival times
- `delta`: censoring indicator (1 = observed event, 0 = right-censored)
- `X`: **wide functional covariate matrix** (one row per subject, one
  column per grid point)
- (optional) `z`: scalar covariates

For convenience, you typically store everything in one `data.frame`:

``` r
head(dat)
#>   time status age  X1   X2   X3  ...  X50
#> 1 ...   ...   ...  ...  ...  ...
```

You refer to functional columns either by explicit names
(e.g. `"X1", "X2", ...`) or via a regex pattern (e.g. `"^X"` with
`x_as_regex = TRUE`).

------------------------------------------------------------------------

## 1. Linear functional AFT model: `fit_lfAFT()`

The linear functional AFT model assumes $$
  \log T_i = \beta_0 + Z_i^\top \beta_Z + \int X_i(s)\,\beta(s)\,ds + \sigma \varepsilon_i,
$$ with lognormal or loglogistic errors. The coefficient function
$\beta(s)$ is represented using spline basis functions, and a roughness
penalty is applied to control smoothness.

### Automatic selection of the smoothing parameter

Instead of specifying a single `lambda`, you supply a **grid**
`lambda_grid`. Internally:

- `optimize_lambda()` evaluates GCV for each candidate `lambda`
- the `lambda` minimizing GCV is chosen
- the final model is refit at that optimal `lambda`

### Example: simulated lognormal lfAFT

``` r
set.seed(1)

n  <- 100        # subjects
nS <- 50         # grid points
s  <- seq(0, 1, length.out = nS)

# true coefficient function
beta_true <- sin(2 * pi * s)

# functional predictor
X <- matrix(rnorm(n * nS), nrow = n, ncol = nS)

# scalar covariate
age <- rnorm(n, 60, 8)

# linear predictor and lognormal survival times
mu <- 0.5 + 0.01 * age + as.vector(X %*% beta_true / nS)
Y  <- exp(mu + rnorm(n, sd = 0.5))
Delta <- rbinom(n, size = 1, prob = 0.8)

dat <- data.frame(
  time   = Y,
  status = Delta,
  age    = age
)
colnames(X) <- paste0("X", seq_len(nS))
dat <- cbind(dat, X)

# fit lfAFT with automatic lambda selection
lambda_grid <- exp(seq(log(1000), log(10000), length.out = 50))

fit_lf <- fit_lfAFT(
  data        = dat,
  y           = "time",
  delta       = "status",
  x           = "^X",          # select X1,...,X50 by regex
  x_as_regex  = TRUE,
  z           = "age",
  family      = "lognormal",
  k           = 15,
  lambda_grid = lambda_grid,
  se          = TRUE
)

fit_lf$lambda_opt
#> NULL
```

### Inspecting $\hat\beta(s)$

``` r
plot(fit_lf$s_grid, fit_lf$betaX_hat, type = "l",
     xlab = "s", ylab = "Estimated beta(s)")
lines(s, beta_true, col = "grey", lty = 2)
legend("topright", legend = c("Estimated", "True"),
       col = c("black", "grey"), lty = c(1, 2), bty = "n")
```

<img src="man/figures/README-lfAFT-beta-plot-1.png" width="100%" />

### Prediction with `predict_lfAFT()`

``` r
newdata <- dat[1:5, ]

# linear predictor mu
mu_new <- predict_lfAFT(fit_lf, newdata, type = "link")

# median survival time exp(mu)
tmed_new <- predict_lfAFT(fit_lf, newdata, type = "response")

# survival probabilities at t = 1, 2, 5
S_new <- predict_lfAFT(fit_lf, newdata,
                       type  = "survival",
                       times = c(1, 2, 5))
```

------------------------------------------------------------------------

## 2. Additive functional AFT model: `fit_afAFT()`

The additive functional AFT model allows the functional predictor to
enter through a more flexible bivariate function: $$
  \log T_i = \alpha_0 + Z_i^\top \alpha_Z +
    \int F\{X_i(s), s\}\,ds + \sigma \varepsilon_i,
$$ where $F(\cdot,\cdot)$ is unknown and estimated using a
tensor-product smooth via **mgcv**.

In this package version, `fit_afAFT()` supports **only the lognormal AFT
family**, using `mgcv::cnorm()`.

### Implementation sketch

- Build matrices
  - `X`: functional covariate (n × nS)  
  - `S`: replicated grid (n × nS)  
  - `L`: quadrature weights (n × nS) via the trapezoid rule
- Fit a tensor-product smooth $$
    \texttt{ti(S, X, by = L, bs = basis, k = k, mc = (FALSE, TRUE))}
  $$ where `mc = c(FALSE, TRUE)` imposes marginal identifiability
  constraints in the functional covariate direction.

### Example: lognormal afAFT

``` r
library(mgcv)

set.seed(2)
n  <- 100
nS <- 50
s  <- seq(0, 1, length.out = nS)

Xfun <- matrix(rnorm(n * nS), nrow = n, ncol = nS)
colnames(Xfun) <- paste0("X", seq_len(nS))

age <- rnorm(n, 60, 8)

beta_true <- cos(2 * pi * s)
mu_af     <- 0.5 + 0.02 * age + as.vector(Xfun %*% beta_true / nS)
Y_af      <- exp(mu_af + rnorm(n, sd = 0.5))
Delta_af  <- rbinom(n, 1, 0.8)

dat_af <- data.frame(
  time   = Y_af,
  status = Delta_af,
  age    = age
)
dat_af <- cbind(dat_af, Xfun)

fit_af <- fit_afAFT(
  data   = dat_af,
  y      = "time",
  delta  = "status",
  x      = grep("^X", names(dat_af), value = TRUE),
  z      = "age",
  family = "lognormal",
  k      = c(10, 10),
  basis  = c("cr", "cr")
)

summary(fit_af)
#> 
#> Family: cnorm(0.485) 
#> Link function: identity 
#> 
#> Formula:
#> logY ~ age + ti(S, X, by = L, bs = basis, k = k, mc = c(FALSE, 
#>     TRUE))
#> 
#> Parametric coefficients:
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept) 0.501188   0.444832   1.127   0.2599   
#> age         0.021347   0.007217   2.958   0.0031 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Approximate significance of smooth terms:
#>             edf Ref.df Chi.sq p-value
#> ti(S,X):L 3.495  4.087  7.595   0.109
#> 
#> R-sq.(adj) =  0.144   Deviance explained = 13.6%
#> -REML = -74.889  Scale est. = 1         n = 100
```

### Prediction with `predict_afAFT()`

``` r
new_af <- dat_af[1:5, ]

eta_af <- predict_afAFT(fit_af, new_af, type = "link")
tmed_af <- predict_afAFT(fit_af, new_af, type = "response")
S_af <- predict_afAFT(fit_af, new_af, type = "survival", times = c(1, 2, 5))
```

------------------------------------------------------------------------

## References

- Qian, W., Cui, E., Brooks-Russell, A., & Wrobel, J. (2025).  
  *Functional Accelerated Failure Time Models for Predicting Time Since
  Cannabis Use.*  
  arXiv preprint arXiv:2510.22343.  
  <https://doi.org/10.48550/arXiv.2510.22343>

- Wood S N (2017). *Generalized Additive Models: An Introduction with
  R.* 2nd ed. Chapman & Hall/CRC.

- `mgcv::cnorm` documentation:  
  <https://rdrr.io/cran/mgcv/man/cnorm.html>
