# Older implementations that can be used for testing newer versions


#' Compute Partial Leverages using the Naive FWL Method
#'
#' This function calculates the N x K matrix of partial leverages by iterating
#' through each regressor and performing a separate FWL-style regression.
#' This approach is computationally expensive, with complexity O(N*K^3).
#'
#' @param X An N x K numeric matrix of regressors. Should include an intercept.
#' @param use_qr Logical. If TRUE, use QR decomposition for the intermediate
#'   regressions, which is more numerically stable. Defaults to FALSE.
#' @param qr_X A pre-computed QR decomposition of X. This argument is ignored
#'   by the naive function but included for signature consistency.
#'
#' @return An N x K matrix where the (i, k)-th element is the partial
#'   leverage of observation i for regressor k.
#' @export
partial_leverages_naive <- function(X, use_qr = FALSE, qr_X = NULL) {
  if (!is.matrix(X)) stop("X must be a matrix.")

  N <- nrow(X)
  K <- ncol(X)

  # Initialize the results matrix
  P <- matrix(0.0, nrow = N, ncol = K)
  colnames(P) <- colnames(X)

  # Loop over each regressor k
  for (k in 1:K) {
    X_k <- X[, k, drop = FALSE]
    X_not_k <- X[, -k, drop = FALSE]

    # Check for the case where there are no other regressors (K=1)
    if (ncol(X_not_k) == 0) {
      residuals <- X_k - mean(X_k) # Residuals from intercept-only model
    } else {
      # Get the residuals of X_k regressed on X_not_k
      if (use_qr) {
        # Numerically stable approach using QR
        qr_not_k <- qr(X_not_k)
        residuals <- qr.resid(qr_not_k, X_k)
      } else {
        # Standard approach using solve(crossprod(X))
        beta_hat <- solve(crossprod(X_not_k), crossprod(X_not_k, X_k))
        residuals <- X_k - X_not_k %*% beta_hat
      }
    }

    # Calculate sum of squared residuals
    ss_res <- sum(residuals^2)

    # Avoid division by zero in case of perfect collinearity
    if (ss_res > 1e-12) {
      P[, k] <- residuals^2 / ss_res
    } else {
      # If X_k is perfectly explained by others, leverages are undefined/zero
      P[, k] <- 0
    }
  }

  return(P)
}


#' Compute Partial Leverages using the Efficient Matrix Method
#'
#' This function calculates the N x K matrix of partial leverages by first
#' computing the inverse of (X'X) and then using its elements to construct
#' all necessary residuals. This approach is much more efficient, with
#' complexity O(N*K^2 + K^3).
#'
#' @param X An N x K numeric matrix of regressors. Should include an intercept.
#' @param use_qr Logical. If TRUE, use QR decomposition to find the inverse
#'   of (X'X), which is more numerically stable. Defaults to FALSE.
#' @param qr_X Optional. A pre-computed QR decomposition of X. If provided,
#'   this will be used instead of re-computing it, saving time.
#'
#' @return An N x K matrix where the (i, k)-th element is the partial
#'   leverage of observation i for regressor k.
#' @export
partial_leverage_g25_1 <- function(X, use_qr = FALSE, qr_X = NULL) {
  if (!is.matrix(X)) stop("X must be a matrix.")

  N <- nrow(X)
  K <- ncol(X)

  # Initialize the results matrix
  P <- matrix(0.0, nrow = N, ncol = K)
  colnames(P) <- colnames(X)

  # Step 1: Compute C = (X'X)^-1, the key matrix
  if (use_qr) {
    # If qr_X is not provided, compute it
    if (is.null(qr_X)) {
      qr_X <- qr(X)
    }
    # Check for rank deficiency
    if(qr_X$rank < K) stop("X is rank-deficient, cannot compute inverse.")
    C <- solve(qr.R(qr_X)) %*% t(solve(qr.R(qr_X)))
  } else {
    C <- solve(crossprod(X))
  }

  # Step 2: Loop over each regressor k
  for (k in 1:K) {
    C_kk <- C[k, k]

    # Calculate the residual vector X_k^* using the elements of C.
    # The regression coefficients of X_k on X_~k are -C_[-k,k] / C_kk
    # So, X_k^* = X_k - X_~k %*% beta = X_k + X_~k %*% (C_[-k,k] / C_kk)
    # Note: Handle K=1 case where there are no other regressors
    if (K > 1) {
        X_k_star <- X[, k] + X[, -k, drop=FALSE] %*% (C[-k, k] / C_kk)
    } else {
        X_k_star <- X[, k] - mean(X[,k])
    }

    # The partial leverages are (X_k^*)^2 divided by sum((X_k^*)^2).
    # The denominator sum((X_k^*)^2) is equal to 1/C_kk.
    # So, P_ik = (X_ik^*)^2 * C_kk
    P[, k] <- X_k_star^2 * C_kk
  }

  return(P)
}

# Compute the N x K matrix of partial leverages (fast: one small factorization).
# You may optionally provide:
# - qr_X: a stats::qr(X) object (can be pivoted). Used when method = "qr".
# - chol_G: an unpivoted Cholesky R of G = t(X) %*% X, i.e. t(R) %*% R = G. Used when method = "chol".
# method = "auto" tries Cholesky first (fast) and falls back to QR if needed

partial_leverages_gpt_1 = function(X, method = c("auto", "chol", "qr"), qr_X = NULL, chol_G = NULL) {
  method = match.arg(method)
  X = as.matrix(X)
  N = nrow(X); K = ncol(X)
  if (K < 1L || N < K) stop("X must be N x K with N >= K and full column rank.")

  # Decide path if auto
  if (method == "auto") {
    # Prefer Cholesky on G for speed if it works cleanly
    if (!is.null(chol_G)) {
      method = "chol"
    } else {
      G = crossprod(X)
      chol_ok = TRUE
      R = tryCatch(chol(G), error = function(e) { chol_ok <<- FALSE; NULL })
      if (chol_ok) {
        chol_G = R
        method = "chol"
      } else {
        method = "qr"
      }
    }
  }

  if (method == "chol") {
    # Use provided chol_G if available; otherwise compute it
    R = if (!is.null(chol_G)) {
      if (!is.matrix(chol_G) || nrow(chol_G) != K || ncol(chol_G) != K)
        stop("chol_G must be a K x K upper-triangular Cholesky factor of t(X) %*% X.")
      chol_G
    } else {
      chol(crossprod(X))
    }
    # Solve G Z^T = X^T using two triangular solves: R^T Y = X^T, then R Z^T = Y
    Y = forwardsolve(t(R), t(X))
    Zt = backsolve(R, Y)
    Z = t(Zt)  # Z = X %*% S (but without forming S)

    # diag(S) from R^{-1}: diag(S) = rowSums((R^{-1})^2)
    Rinv = backsolve(R, diag(K))
    Sdiag = rowSums(Rinv * Rinv)

    # L = (Z^2) scaled columnwise by 1 / Sdiag
    Lmat = sweep(Z * Z, 2L, Sdiag, "/")
    dimnames(Lmat) = list(rownames(X), colnames(X))
    return(Lmat)
  }

  # method == "qr"
  QR = if (!is.null(qr_X)) {
    qr_X
  } else {
    qr(X)  # pivoted by default; stable
  }
  if (is.null(QR$rank) || QR$rank < K) stop("X appears rank-deficient under QR.")
  # Extract thin Q (N x K) and R (K x K)
  Q = qr.Q(QR, complete = FALSE)
  R = qr.R(QR, complete = FALSE)
  piv = QR$pivot
  if (is.null(piv)) piv = seq_len(K)

  # Compute Rt_inv = R^{-T} via a single backsolve with multiple RHS
  Rt_inv = backsolve(R, diag(K), transpose = TRUE)  # K x K

  # Zp corresponds to pivoted column order; then unpermute to original
  Zp = Q %*% Rt_inv                       # N x K (pivoted order)
  Z = matrix(0.0, nrow = N, ncol = K)
  Z[, piv] = Zp

  # diag(S), in pivoted order: colSums(Rt_inv^2). Unpermute to original order.
  Sd_p = colSums(Rt_inv * Rt_inv)
  Sdiag = numeric(K); Sdiag[piv] = Sd_p

  Lmat = sweep(Z * Z, 2L, Sdiag, "/")
  dimnames(Lmat) = list(rownames(X), colnames(X))
  Lmat
}

