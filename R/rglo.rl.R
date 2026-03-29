#' Return Levels for the Generalized Logistic Distribution
#'
#' Computes return levels and their standard errors for a stationary
#' generalized logistic model fitted by \code{\link{rglo.fit}}.
#'
#' @param z An object returned by \code{\link{rglo.fit}}. The fitted model
#'   should represent a stationary model.
#' @param year A numeric vector of return periods for which return levels
#'   are to be computed.
#' @param show Logical. If \code{TRUE}, the estimated return levels and
#'   their standard errors are printed.
#'
#' @return The input object \code{z} with two additional components:
#' \item{rl}{A numeric vector of estimated return levels.}
#' \item{rlse}{A numeric vector of standard errors of the estimated return levels.}
#'
#' @details
#' For a return period \eqn{T}, the return level is defined as the quantile
#' exceeded with probability \eqn{1/T}. Under the generalized logistic
#' distribution, the return level is
#' \deqn{x_T = \mu + \frac{\sigma}{\xi} \left[1 - \left(\frac{1 - 1/T}{1/T}\right)^{-\xi}\right],}
#' which is equivalently written in the implementation as
#' \deqn{x_T = \mu + \frac{\sigma}{\xi} - \frac{\sigma}{\xi}
#' \left(\frac{1/T}{1 - 1/T}\right)^{\xi}.}
#' Standard errors are obtained using the delta method.
#'
#' @references
#'
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
#' @seealso \code{\link{rglo.fit}}, \code{\link{rglo.prof}}
#' @export
#'
#' @examples
#' x <- rglor(n = 50, r = 3, loc = 10, scale = 2, shape = 0.1)
#' fit <- rglo.fit(x$rmat)
#' out <- rglo.rl(fit, year = c(20, 50, 100, 200))
rglo.rl <- function(z, year = c(20, 50, 100, 200), show=TRUE) {

  if (!inherits(z, "rglo.fit")) {
    warning("'z' does not inherit from class 'rglo.fit'.")
  }

  if (isTRUE(z$trans)) {
    stop("Return levels are only implemented for stationary models.")
  }

  if (!is.numeric(year) || any(!is.finite(year)) || any(year <= 1)) {
    stop("'year' must be a numeric vector with values greater than 1.")
  }

  if (is.null(z$mle) || length(z$mle) < 3) {
    stop("'z$mle' must contain location, scale, and shape estimates.")
  }

  if (is.null(z$cov) || all(is.na(z$cov))) {
    stop("Covariance matrix is not available in 'z$cov'.")
  }

  mu <- z$mle[1]
  sig <- z$mle[2]
  xi <- z$mle[3]

  if (!is.finite(sig) || sig <= 0) {
    stop("The fitted scale parameter must be positive.")
  }

  if (!is.finite(xi) || abs(xi) < .Machine$double.eps^0.5) {
    stop("The fitted shape parameter is too close to 0 for 'rglo.rl()'.")
  }

  del <- matrix(nrow = 3, ncol = length(year))

  f <- 1 - (1 / year)
  f_inv <- 1 / f
  g <- exp(xi * log(f_inv - 1))

  del[1, ] <- 1
  del[2, ] <- (1 - g) / xi
  del[3, ] <- -(sig / xi^2) * (-g * (log(f_inv - 1) * xi - 1) - 1)

  del.t <- t(del)

  z$rl <- mu + sig / xi - sig / xi * g
  z$rlse <- sqrt(diag(del.t %*% z$cov %*% del))

  names(z$rl) <- paste0(year, "y")
  names(z$rlse) <- paste0(year, "y")

  if(show) {
    print(z[c("rl", "rlse")])}

  z
}
