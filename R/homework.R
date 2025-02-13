#'  CDF of a normal mixture distribution
#'
#' @param q 
#' @param mean_1 
#' @param sd_1 
#' @param mean_2 
#' @param sd_2 
#' @param mixprob 
#'
#' @return 
#' @export
pnormmix <- function(q, mean_1, sd_1, mean_2, sd_2, mixprob) {
  p_1 <- pnorm(q, mean = mean_1, sd = sd_1)
  p_2 <- pnorm(q, mean = mean_2, sd = sd_2)
  return(mixprob * p_1 + (1 - mixprob) * p_2)
}

#'PDF of a normal mixture distribution
#'
#' @param x 
#' @inheritParams pnormmix
#' @return
#' @export
dnormmix <- function(x, mean_1, sd_1, mean_2, sd_2, mixprob) {
  d_1 <- dnorm(x, mean = mean_1, sd = sd_1)
  d_2 <- dnorm(x, mean = mean_2, sd = sd_2)
  return(mixprob * d_1 + (1 - mixprob) * d_2)
}

#' random numbers from a normal mixture distribution
#'
#' @param n 
#' @inheritParams pnormmix
#' @return
#' @export
rnormmix <- function(n, mean_1, sd_1, mean_2, sd_2, mixprob) {
  u <- runif(n)
  x <- ifelse(u < mixprob, rnorm(n, mean = mean_1, sd = sd_1),
              rnorm(n, mean = mean_2, sd = sd_2))
  return(x)
}

#' quantile function (inverse CDF) of a normal mixture distribution
#'
#' @param p 
#' @inheritParams pnormmix
#' @return
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
