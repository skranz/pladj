example = function() {
  library(pladj)
  # --- Setup --------------------------------------------------------------
  set.seed(123)
  suppressPackageStartupMessages({
    library(fixest)
  })

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


  #cat("All equal within numerical tolerance. ✅\n")
}


# =========================
# Core: partial_leverages
# =========================
# Fast N x m partial leverage matrix for columns of X (all columns or a subset).
# Works for OLS and WLS. If weights are provided, we internally use sqrt(w) * X.
# method: "auto" (default), "chol", or "qr".
# select: integer or character vector of columns to compute; default NULL = all.
# chol_G: optional K x K Cholesky R of crossprod(sqrt(w)*X) to reuse.
# qr_X:   optional qr() object of sqrt(w)*X to reuse (e.g., lm$qr).
# block_size: optional number of rows to process at a time to reduce temporaries.
# return_extra: if TRUE, attaches attributes Sdiag (vector of S_jj for returned cols).
partial_leverages = function(
  X,
  weights = NULL,
  method = c("auto", "chol", "qr"),
  select = NULL,
  chol_G = NULL,
  qr_X = NULL,
  block_size = NULL,
  allow_zero_weights = FALSE,
  return_extra = FALSE
) {
  method = match.arg(method)
  X = as.matrix(X)
  N = nrow(X); K = ncol(X)
  if (K == 0L) stop("X must have at least one column.")
  if (N < K) stop("Require N >= K and full column rank.")

  # Column selection
  if (is.null(select)) {
    sel_idx = seq_len(K)
  } else {
    sel_idx = if (is.character(select)) match(select, colnames(X)) else as.integer(select)
    if (anyNA(sel_idx)) stop("select contains names/indices not found in X.")
    if (length(unique(sel_idx)) != length(sel_idx)) stop("select has duplicates.")
  }
  m = length(sel_idx)

  # Weights: build sqrt(w) and optionally drop zero-weight rows
  if (is.null(weights)) {
    sqrtw = rep(1.0, N)
    keep = rep(TRUE, N)
  } else {
    if (length(weights) != N) stop("weights must have length nrow(X).")
    w = as.numeric(weights)
    if (!allow_zero_weights && any(w <= 0, na.rm = TRUE)) {
      stop("All weights must be strictly positive unless allow_zero_weights = TRUE.")
    }
    keep = !is.na(w) & (allow_zero_weights | (w > 0))
    if (!all(keep)) {
      # We drop zero/NA weight rows from computations, will pad back zeros later
      X_full = X
      X = X_full[keep, , drop = FALSE]
      w = w[keep]
      N_keep = nrow(X)
    } else {
      X_full = NULL
      N_keep = N
    }
    sqrtw = sqrt(pmax(w, 0))
  }

  if (nrow(X) < ncol(X)) stop("After filtering rows, N < K; cannot proceed.")

  # Preweighted design
  Xw = X * sqrtw

  # Helper: process in blocks of rows to limit temporaries
  row_blocks = function(n, bs) {
    if (is.null(bs) || !is.numeric(bs) || bs <= 0 || bs >= n) return(list(seq_len(n)))
    starts = seq.int(1L, n, by = bs)
    ends = pmin(starts + bs - 1L, n)
    mapply(seq.int, starts, ends, SIMPLIFY = FALSE)
  }
  blocks = row_blocks(nrow(Xw), block_size)

  # Allocate output (for kept rows only, selected columns only)
  out = matrix(NA_real_, nrow = nrow(Xw), ncol = m,
               dimnames = list(rownames(Xw), colnames(X)[sel_idx]))

  # ---------- Cholesky path (preferred) ----------
  chol_path = function() {
    # Use provided chol_G or compute it
    R = if (!is.null(chol_G)) {
      if (!is.matrix(chol_G) || nrow(chol_G) != K || ncol(chol_G) != K)
        stop("chol_G must be K x K upper triangular.")
      chol_G
    } else {
      chol(crossprod(Xw))  # may fail -> caller will fallback to QR
    }

    # If computing all columns: precompute Sdiag via R^{-1}, then stream blocks
    if (identical(sel_idx, seq_len(K))) {
      Rinv = backsolve(R, diag(K))
      Sdiag = rowSums(Rinv * Rinv)  # diag(S)

      for (idx in blocks) {
        XB = Xw[idx, , drop = FALSE]
        Y = forwardsolve(t(R), t(XB))
        ZS_t = backsolve(R, Y)
        ZS = t(ZS_t)                      # N_block x K
        out[idx, ] = sweep(ZS * ZS, 2L, Sdiag, "/")
      }

      attr(out, "Sdiag") = Sdiag
      return(out)
    }

    # Selected columns only: solve S E_sel and multiply by Xw
    I_sel = diag(K)[, sel_idx, drop = FALSE]
    Ysel = forwardsolve(t(R), I_sel)
    U = backsolve(R, Ysel)                # U = S %*% E_sel, K x m

    # Sdiag for selected columns = take row j of column j
    Sdiag_sel = U[cbind(sel_idx, seq_len(m))]

    for (idx in blocks) {
      XB = Xw[idx, , drop = FALSE]
      XS_sel = XB %*% U                   # N_block x m
      out[idx, ] = sweep(XS_sel * XS_sel, 2L, Sdiag_sel, "/")
    }

    attr(out, "Sdiag") = Sdiag_sel
    out
  }

  # ---------- QR path (pivoted) ----------
  qr_path = function(QR = NULL) {
    QR = if (!is.null(QR)) QR else qr(Xw)  # pivoted QR
    rnk = QR$rank
    if (rnk < K) stop("Rank deficiency detected in X (after weighting).")
    Q = qr.Q(QR, complete = FALSE)
    R = qr.R(QR, complete = FALSE)
    piv = QR$pivot

    # helper to map original column indices to pivot positions
    pos_in_piv = integer(K); pos_in_piv[piv] = seq_len(K)

    if (identical(sel_idx, seq_len(K))) {
      Rt_inv = backsolve(R, diag(K), transpose = TRUE)  # K x K (pivoted order)
      Zp = Q %*% Rt_inv                                 # N x K (pivoted columns)
      Sd_p = colSums(Rt_inv * Rt_inv)
      # place back to original order
      out_full = matrix(NA_real_, nrow = nrow(Xw), ncol = K,
                        dimnames = list(rownames(Xw), colnames(X)))
      out_full[, piv] = sweep(Zp * Zp, 2L, Sd_p, "/")
      attr(out_full, "Sdiag") = { Sd = numeric(K); Sd[piv] = Sd_p; Sd }
      return(out_full[, sel_idx, drop = FALSE])
    }

    # Selected columns only
    ppos = pos_in_piv[sel_idx]            # pivot positions for selected
    I_sel_p = diag(K)[, ppos, drop = FALSE]
    Rt_inv_sel = backsolve(R, I_sel_p, transpose = TRUE)  # K x m
    Zp_sel = Q %*% Rt_inv_sel                              # N x m
    Sd_sel_p = colSums(Rt_inv_sel * Rt_inv_sel)

    out_sel = sweep(Zp_sel * Zp_sel, 2L, Sd_sel_p, "/")
    attr(out_sel, "Sdiag") = Sd_sel_p
    colnames(out_sel) = colnames(X)[sel_idx]
    rownames(out_sel) = rownames(Xw)
    out_sel
  }

  # Choose method
  result_kept =
    if (method == "chol") {
      chol_path()
    } else if (method == "qr") {
      qr_path(qr_X)
    } else {
      # auto: try chol then fallback to QR
      tryCatch(chol_path(),
               error = function(e) qr_path(qr_X))
    }

  # If we dropped some rows due to zero/NA weights, pad them back with zeros
  if (!is.null(X_full)) {
    res = matrix(0.0, nrow = nrow(X_full), ncol = ncol(result_kept),
                 dimnames = list(rownames(X_full), colnames(result_kept)))
    res[keep, ] = result_kept
    if (!is.null(attr(result_kept, "Sdiag"))) attr(res, "Sdiag") = attr(result_kept, "Sdiag")
    result_kept = res
  }

  if (!return_extra) attr(result_kept, "Sdiag") = NULL
  result_kept
}


# =========================
# Wrapper: feols
# =========================
# Computes partial leverages for the non-FE regressors in a feols model.
# By default we demean via fixest::demean() (or reuse X_demeaned if available),
# use the estimation sample obs(), and pass weights(fit) subset to that sample.
partial_leverages.feols = function(
  fit,
  method = c("auto", "chol", "qr"),
  select = NULL,
  block_size = NULL,
  allow_zero_weights = FALSE,
  use_stored_demeaned = TRUE
) {
  method = match.arg(method)
  if (!inherits(fit, "fixest")) stop("fit must be a fixest object from feols().")
  library(fixest)

  # Which RHS? If IV, use the second-stage RHS
  rhs_type = if (!is.null(fit$fml_all) && !is.null(fit$fml_all$iv)) "iv.rhs2" else "rhs"

  # Estimation sample indices and RHS in that sample
  idx = fixest::obs(fit)                                        # integer vector of rows used
  Z = stats::model.matrix(fit, type = rhs_type, sample = "estimation", as.matrix = TRUE)
  N = nrow(Z)
  if (N == 0L) stop("No rows in the estimation sample.")
  K = ncol(Z)
  if (K == 0L) stop("No non-FE regressors found.")

  # Weights aligned to estimation sample
  wf = tryCatch(stats::weights(fit), error = function(e) NULL)  # returns length original, with NA for dropped
  w = if (is.null(wf)) rep(1.0, length(idx)) else as.numeric(wf[idx])

  # Demeaned RHS:
  # 1) reuse X_demeaned if available (and align columns), else
  # 2) use fixest::demean(fit, sample="estimation") and align, else
  # 3) if no FE, Z is already fine.
  Z_dm = NULL
  if (use_stored_demeaned && !is.null(fit$X_demeaned)) {
    Xd = fit$X_demeaned
    if (!is.null(colnames(Xd))) {
      mi = match(colnames(Z), colnames(Xd))
      if (!anyNA(mi)) Z_dm = Xd[, mi, drop = FALSE]
    }
  }
  if (is.null(Z_dm)) {
    dm_all = tryCatch(fixest::demean(fit, sample = "estimation", as.matrix = TRUE),
                      error = function(e) NULL)
    if (!is.null(dm_all) && !is.null(colnames(dm_all))) {
      mi = match(colnames(Z), colnames(dm_all))
      if (!anyNA(mi)) Z_dm = dm_all[, mi, drop = FALSE]
    }
  }
  if (is.null(Z_dm)) {
    # no FE or demeaning not needed
    Z_dm = Z
  }

  # Column selection by name is natural here (e.g., only your non-FE subset)
  if (is.character(select)) {
    if (is.null(colnames(Z_dm))) stop("select is character but Z_dm lacks column names.")
    miss = setdiff(select, colnames(Z_dm))
    if (length(miss)) stop("select contains unknown columns: ", paste(miss, collapse = ", "))
  }

  # Compute partial leverages with the core
  L = partial_leverages(
    X = Z_dm,
    weights = w,
    method = method,
    select = select,
    block_size = block_size,
    allow_zero_weights = allow_zero_weights,
    return_extra = FALSE
  )
  L
}


# =========================
# Wrapper: lm
# =========================
# Uses the model matrix of the fitted lm and reuses fit$qr for speed when possible.
partial_leverages.lm = function(
  fit,
  method = c("auto", "chol", "qr"),
  select = NULL,
  block_size = NULL,
  allow_zero_weights = FALSE
) {
  method = match.arg(method)
  if (!inherits(fit, "lm")) stop("fit must be an 'lm' object.")
  X = stats::model.matrix(fit)                  # estimation-sample design
  w = tryCatch(stats::weights(fit), error = function(e) NULL)

  # Reuse lm's QR (of the weighted design) when QR is chosen/needed
  qr_obj = if (!is.null(fit$qr)) fit$qr else NULL

  L = partial_leverages(
    X = X,
    weights = w,
    method = method,
    select = select,
    qr_X = if (!is.null(qr_obj) && method != "chol") qr_obj else NULL,
    block_size = block_size,
    allow_zero_weights = allow_zero_weights,
    return_extra = FALSE
  )
  L
}


# My own naive implementation: only works without weights
.partial_leverages_my_test <- function(X) {
  if (!is.matrix(X)) stop("X must be a matrix.")
  N <- nrow(X)
  K <- ncol(X)
  # Initialize the results matrix
  pl <- matrix(0.0, nrow = N, ncol = K)
  colnames(pl) <- colnames(X)

  k = 1
  # Loop over each regressor k
  for (k in 1:K) {
    x_k <- X[, k, drop = FALSE]
    X_not_k <- X[, -k, drop = FALSE]
    res = lm.fit(x=X_not_k, y=x_k)
    x_tilde = res$residuals
    pl[,k] = x_tilde^2 / sum(x_tilde^2)
  }

  return(pl)
}
