#' Compute the CDF of a normal mixture distribution
#'
#' @param q Numeric value(s) where the CDF is evaluated.
#' @param mean_1 Mean of the first normal distribution.
#' @param sd_1 Standard deviation of the first normal distribution.
#' @param mean_2 Mean of the second normal distribution.
#' @param sd_2 Standard deviation of the second normal distribution.
#' @param mixprob Probability weight for the first normal distribution (between 0 and 1).
#'
#' @return A numeric value representing the CDF at `q`.
#' @export
pnormmix <- function(q, mean_1, sd_1, mean_2, sd_2, mixprob) {
  p_1 <- pnorm(q, mean = mean_1, sd = sd_1)
  p_2 <- pnorm(q, mean = mean_2, sd = sd_2)
  return(mixprob * p_1 + (1 - mixprob) * p_2)
}

#' Compute the PDF of a normal mixture distribution
#'
#' @param x Numeric value(s) where the PDF is evaluated.
#' @inheritParams pnormmix
#' @return A numeric value representing the PDF at `x`.
#' @export
dnormmix <- function(x, mean_1, sd_1, mean_2, sd_2, mixprob) {
  d_1 <- dnorm(x, mean = mean_1, sd = sd_1)
  d_2 <- dnorm(x, mean = mean_2, sd = sd_2)
  return(mixprob * d_1 + (1 - mixprob) * d_2)
}

#' Generate random numbers from a normal mixture distribution
#'
#' @param n Number of random samples to generate.
#' @inheritParams pnormmix
#' @return A vector of `n` random values from the mixture distribution.
#' @export
rnormmix <- function(n, mean_1, sd_1, mean_2, sd_2, mixprob) {
  u <- runif(n)
  x <- ifelse(u < mixprob, rnorm(n, mean = mean_1, sd = sd_1),
              rnorm(n, mean = mean_2, sd = sd_2))
  return(x)
}

#' Compute the quantile function (inverse CDF) of a normal mixture distribution
#'
#' @param p Probability value(s) for quantile calculation.
#' @inheritParams pnormmix
#' @return A numeric value representing the quantile at probability `p`.
#' @export
qnormmix <- function(p, mean_1, sd_1, mean_2, sd_2, mixprob) {
  quantiles <- numeric(length(p))
  for (i in seq_along(p)) {
    if (p[i] == 0) {
      quantiles[i] <- -Inf
    } else if (p[i] == 1) {
      quantiles[i] <- Inf
    } else {
      qf <- function(x) pnormmix(x, mean_1, sd_1, mean_2, sd_2, mixprob) - p[i]
      quantiles[i] <- uniroot(qf, interval = c(-10, 10), extendInt = "yes")$root
    }
  }
  return(quantiles)
}
