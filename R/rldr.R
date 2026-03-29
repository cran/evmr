#' Random Generation from the Logistic Distribution for r-Largest Order Statistics
#'
#' Generates random samples from the logistic distribution for
#' \eqn{r}-largest order statistics.
#'
#' @param n A positive integer specifying the number of observations.
#' @param r A positive integer specifying the number of order statistics
#'   for each observation.
#' @param loc A numeric value specifying the location parameter.
#' @param scale A positive numeric value specifying the scale parameter.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{umat}: an \code{n x r} matrix of independent uniform random numbers
#'   \item \code{wmat}: an \code{n x r} matrix of transformed uniform variables
#'   \item \code{rmat}: an \code{n x r} matrix of simulated \eqn{r}-largest order statistics
#' }
#'
#' @details
#' The function first generates independent uniform random variables and then
#' constructs decreasing transformed variables recursively. These are
#' transformed by the logistic quantile function \code{\link{qld}}.
#'
#' @references
#'
#' Coles, S. (2001).
#' An Introduction to Statistical Modeling of Extreme Values.
#' Springer.
#'
#' Bader, B., Yan, J., & Zhang, X. (2017).
#' Automated selection of \eqn{r} for the r-largest order statistics approach.
#' \emph{Statistics and Computing}.
#' \doi{10.1007/s11222-016-9697-3}
#'
#' Shin, Y., & Park, J-S. (2024).
#' Generalized logistic model for r-largest order statistics with
#' hydrological application.
#' \emph{Stochastic Environmental Research and Risk Assessment}.
#' \doi{10.1007/s00477-023-02642-7}
#'
#'
#' @export
#'
#' @examples
#' x <- rldr(n=50, r=3, loc = 0, scale = 1)
#' x$rmat
#'
rldr <- function(n, r, loc = 0, scale = 1) {

  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0 || n != as.integer(n)) {
    stop("'n' must be a positive integer.")
  }

  if (!is.numeric(r) || length(r) != 1 || !is.finite(r) || r <= 0 || r != as.integer(r)) {
    stop("'r' must be a positive integer.")
  }

  if (!is.numeric(loc) || length(loc) != 1 || !is.finite(loc)) {
    stop("'loc' must be a finite numeric scalar.")
  }

  if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale)) {
    stop("'scale' must be a finite numeric scalar.")
  }

  if (scale <= 0) {
    stop("'scale' must be positive.")
  }

  z <- list()

  umat <- matrix(0, nrow = n, ncol = r)
  wmat <- matrix(0, nrow = n, ncol = r)

  colnames(umat) <- paste0("u", seq_len(r))
  colnames(wmat) <- paste0("w", seq_len(r))

  i <- 1
  while (i <= n) {
    u <- stats::runif(r)

    w <- numeric(r)
    w[1] <- u[1]

    if (r >= 2) {
      for (j in 2:r) {
        w[j] <- w[j - 1] * (u[j]^(1 / j))
      }
    }

    if (all(diff(w) < 0)) {
      umat[i, ] <- u
      wmat[i, ] <- w
      i <- i + 1
    }
  }

  z$umat <- umat
  z$wmat <- wmat

  z$rmat <- matrix(
    qld(wmat, loc = loc, scale = scale),
    nrow = nrow(wmat),
    ncol = ncol(wmat)
  )

  colnames(z$rmat) <- paste0("r", seq_len(r))

  invisible(z)
}
