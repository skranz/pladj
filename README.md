# pladj


Partial leverage adjusted DOF for HC1 and HC2 standard errors, as proposed by Kranz (2024)

Main function is dof_pladj that can take the returned object from a call to `fixest:feols` or `lm` as argument.

Example:

```r
  library(pladj)
  library(fixest)
  # --- Setup --------------------------------------------------------------
  set.seed(123)

  # --- Simulate panel with two-way FEs and weights -----------------------
  N = 100
  G_firm = 5
  G_year = 10

  firm = sample.int(G_firm, N, replace = TRUE)
  year = sample.int(G_year, N, replace = TRUE)

  # Fixed effects (drawn once per level)
  alpha = rnorm(G_firm, sd = 1.0)[firm]
  gamma = rnorm(G_year, sd = 0.5)[year]


  # Non-FE regressors
  x1 = rnorm(N) + alpha + gamma
  x2 = rt(N,df = 1)


  # x2 has much wider tails
  # we will expect much more concentrated partial leverages and much
  # smaller partial leverag adjusted dof.
  sd(x1);sd(x2)
  hist(x2)
  # Heteroskedastic error term
  sigma = 2*abs(gamma)^3
  u = rnorm(N, sd = sigma)

  # for the moment no weights
  #w = runif(N, min = 0.2, max = 2.0)
  w = rep(1, N)

  # Outcome
  beta1 = 1.5
  beta2 = -0.7
  y = beta1 * x1 + beta2 * x2 + alpha + gamma + u

  # Build data
  dat = data.frame(y = y, x1 = x1, x2 = x2, firm = firm, year = year, w = w)
  rownames(dat) = sprintf("i%05d", seq_len(N))

  # --- Fit the three models ----------------------------------------------
  # 1) fixest: absorb FEs
  fit_fe = feols(y ~ x1 + x2 | firm + year, data = dat, weights = ~ w)
  # degrees of freedom
  dof_pladj(fit_fe)

  # 2) lm: LSDV (factor dummies)
  fit_lsdv = lm(y ~ x1 + x2 + factor(firm) + factor(year), data = dat, weights = w)
  dof_pladj(fit_lsdv)[2:3]



  # 3) lm: within/demeaned data (weighted demeaning)
  dm = fixest::demean(cbind(y = dat$y, x1 = dat$x1, x2 = dat$x2),
                      f = list(dat$firm, dat$year),
                      weights = dat$w, as.matrix = TRUE)
  dat_dm = data.frame(y = dm[, "y"], x1 = dm[, "x1"], x2 = dm[, "x2"])
  rownames(dat_dm) = rownames(dat)  # keep identical ordering and names
  fit_dm = lm(y ~ 0 + x1 + x2, data = dat_dm, weights = dat$w)
  dof_pladj(fit_dm)


  # Test partial leverages
  pl_fe = partial_leverages.feols(fit_fe, method = "auto")  # N x 2 (x1, x2)
  colSums(pl_fe)
  pl_dm = partial_leverages.lm(fit_dm, method = "auto")     # N x 2
  colSums(pl_dm)
  pl_lsdv = partial_leverages.lm(fit_lsdv, method = "auto", select = c("x1", "x2"))
  colSums(pl_lsdv)
  pl_lsdv = pl_lsdv[, colnames(pl_fe), drop = FALSE]

  # my own simple implementation
  X = model.matrix(~x1 + x2 + factor(firm) + factor(year),dat)
  pl_my = .partial_leverages_my_test(X)[,colnames(pl_fe)]
  colSums(pl_my)

  cat("Max abs diff (feols vs demeaned-lm): ", format(max(abs(pl_fe - pl_dm)), digits = 6), "\n")
  cat("Max abs diff (feols vs my implementation)   : ", format(max(abs(pl_fe - pl_my)), digits = 6), "\n")
  cat("Mean abs diff (feols vs LSDV-lm)   : ", format(max(abs(pl_fe - pl_lsdv)), digits = 6), "\n")

  # Numeric problems or substantial difference between LSDV and demeaned?
  col="x1"
  plot(pl_lsdv[,col], pl_fe[,col]); abline(a=0,b=1, col="blue")
  plot(pl_my[,col], pl_fe[,col]); abline(a=0,b=1, col="blue")
  col="x2"
  plot(pl_lsdv[,col], pl_fe[,col]); abline(a=0,b=1, col="blue")
```

See also the `example` function in `pl.R`.
