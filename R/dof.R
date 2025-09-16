#' Degrees of freedom using the partial leverages adjustment proposed by Kranz (2024)
#'
#' Use these degree of freedoms for the t-test of HC1 or HC2 standard errors
#' Only makes sense for these hetereoksedasticity robust standard errors.
#' Does not work for cluster robust standard errors
#'
#' @param x A regression result from lm or feols. Alternatively, a matrix of partial leverage generated with partial_leverages.
dof_pladj = function(x) {

  if (is_pl_mat(x)) {
   pl_mat = x
  } else if (is(x,"fixest")) {
   # TO DO: check if from feols
   pl_mat = partial_leverages.feols(x)
  } else if (is(x,"lm")) {
   pl_mat = partial_leverages.lm(x)
  } else {
   stop("x must be the result of lm, feols or a partial leverages matrix")
  }
  # Herfindahl-Hirschman Index
  HHI = colSums(pl_mat^2)
  # Partial leverage adjusted sample sizes
  n_adj = 1/HHI
  # Degrees of freedom
  dof = n_adj-1
  dof

}

is_pl_mat = function(x) {
  # maybe test more strictly
  # problem X matrix has same dimensions than pl matrix
  # we could test: sum up to 1 all positive, but need tolerances
  # and takes a bit time
  is.matrix(x)
}
