#' Log-Likelihood Contributions for the rGLO Model
#'
#' Computes the observation-wise log-likelihood contributions for the
#' r-largest generalized logistic distribution (rGLO) model.
#'
#' @param data A numeric vector, matrix, or data frame of observations.
#'   If a vector is supplied, it is treated as a one-column matrix.
#'   If a matrix or data frame is supplied, each row is treated as one
#'   observation and columns represent decreasing order statistics.
#' @param par A numeric vector of length 3 giving the location, scale,
#'   and shape parameters, respectively.
#'
#' @return A numeric vector of log-likelihood contributions, one for each row
#'   of \code{data}. If the parameter combination is invalid, the function
#'   returns \code{Inf}.
#'
#' @details
#' This function is mainly intended for internal likelihood evaluation.
#' Invalid parameter combinations return \code{Inf}, which is often more
#' robust than stopping with an error when used inside iterative procedures.
#'
#'@references
#' Ahmad, M. I., Sinclair, C. D., and Werritty, A. (1988).
#' Log-logistic flood frequency analysis.
#' \emph{Journal of Hydrology}.
#' \doi{10.1016/0022-1694(88)90015-7}
#'
#' Coles, S. (2001).
#' An Introduction to Statistical Modeling of Extreme Values.
#' Springer.
#'
#' Shin, Y., & Park, J-S. (2024).
#' Generalized logistic model for r-largest order statistics with
#' hydrological application.
#' \emph{Stochastic Environmental Research and Risk Assessment}.
#' \doi{10.1007/s00477-023-02642-7}
#'
#' @export
#'
#' @examples
#' x <- rglor(n=50, r=3, loc = 10, scale = 2, shape = 0.1)
#' fit <- rglo.fit(x$rmat, num_inits = 5)
#' rgloLh(data=fit$data,par=fit$mle)
rgloLh <- function(data, par) {

  tol <- .Machine$double.eps^0.5

  if (is.vector(data)) {
    data <- matrix(data, ncol = 1)
  } else {
    data <- as.matrix(data)
  }

  nr <- nrow(data)

  if (!is.numeric(data)) {
    return(rep(1e6, nr))
  }

  if (!is.numeric(par) || length(par) != 3 || any(!is.finite(par))) {
    return(rep(1e6, nr))
  }

  R <- ncol(data)

  mu <- par[1]
  sc <- par[2]
  xi <- par[3]

  if (!is.finite(sc) || sc <= 0) {
    return(rep(1e6, nr))
  }

  if (!is.finite(xi) || abs(xi) < tol) {
    return(rep(1e6, nr))
  }

  ri <- R - seq_len(R)
  cr <- 1 + ri

  if (any(cr <= 0, na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  y <- 1 - xi * ((data - mu) / sc)
  yr <- 1 - xi * ((data[, R] - mu) / sc)

  if (any(y <= 0, na.rm = TRUE) || any(yr <= 0, na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  f <- 1 + yr^(1 / xi)

  if (any(!is.finite(f), na.rm = TRUE) || any(f <= 0, na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  log.den <- -R * log(sc) +
    sum(log(cr)) +
    (1 + R) * log(1 / (1 + yr^(1 / xi))) +
    rowSums(((1 / xi) - 1) * log(y), na.rm = TRUE)

  if (any(!is.finite(log.den), na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  log.den
}
