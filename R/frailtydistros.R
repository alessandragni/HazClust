#' @title No Frailty Distribution
#' @description Computes the Laplace transform for a model with no frailty.
#' @param s Argument of the Laplace transform.
#' @param what Character; either `"logLT"` for \eqn{\log[ (-1)^k L^{{(k)}}(s) ]} or `"tau"` the Kendall's Tau.
#' @return Numeric value. `-s` if `what = "logLT"`, otherwise `NULL`.
#' @keywords internal
fr.none <- function(s, what = "logLT") {
  if (what == "logLT")
    return(-s)
  else if (what == "tau")
    return(NULL)
}


#' @title Gamma Frailty Distribution
#' @description Computes the log-Laplace transform or Kendall’s tau for the Gamma frailty model.
#' @param k Integer; derivative order.
#' @param s Numeric; argument of the Laplace transform.
#' @param theta Numeric > 0; heterogeneity parameter.
#' @param what Character; either `"logLT"` for \eqn{\log[ (-1)^k L^{{(k)}}(s) ]} or `"tau"` the Kendall's Tau.
#' @return Numeric value of log-Laplace transform or tau.
#' @keywords internal
fr.gamma <- function(k, s, theta, what = "logLT") {
  if (what == "logLT") {
    res <- ifelse(k == 0,
                  -1 / theta * log(1 + theta * s),
                  - (k + 1 / theta) * log(1 + theta * s) +
                    sum(log(1 + (seq(from = 0, to = k - 1, by = 1) * theta))))
    return(res)
  } else if (what == "tau") {
    return(theta / (theta + 2))
  }
}


#' @title Inverse Gaussian Frailty Distribution
#' @description Computes the log-Laplace transform or Kendall’s tau for the Inverse Gaussian frailty.
#' @param k Integer; derivative order.
#' @param s Numeric > 0; Laplace argument.
#' @param theta Numeric > 0; heterogeneity parameter.
#' @param what Character; either `"logLT"` for \eqn{\log[ (-1)^k L^{{(k)}}(s) ]} or `"tau"` the Kendall's Tau.
#' @return Numeric value or message if tau is unstable.
#' @keywords internal
fr.ingau <- function(k, s, theta, what = "logLT") {
  if (what == "logLT") {
    z <- theta^(-0.5) * sqrt(2 * s + theta^(-1))
    res <- ifelse(k == 0,
                  1 / theta * (1 - sqrt(1 + 2 * theta * s)),
                  -k / 2 * log(2 * theta * s + 1) +
                    log(besselK(z, k - 0.5)) -
                    (log(pi / (2 * z)) / 2 - z) +
                    1 / theta * (1 - sqrt(1 + 2 * theta * s)))
    return(res)
  } else if (what == "tau") {
    integrand <- function(u) exp(-u) / u
    int <- integrate(integrand, lower = 2 / theta, upper = Inf)$value
    tau <- 0.5 - 1 / theta + (2 * theta^(-2) * exp(2 / theta) * int)
    if (is.nan(tau) || tau < 0)
      tau <- "The value of 'theta' is too small for computing Kendall's Tau numerically!"
    return(tau)
  }
}


#' @title Polynomial Omega Sum for Positive Stable
#' @description Computes helper sum over Omega matrix for Positive Stable frailty.
#' @param k Order of derivative.
#' @param s Argument for Laplace.
#' @param nu Stability parameter.
#' @param Omega Matrix of omega coefficients.
#' @param correct Correction factor for large clusters.
#' @return Numeric value of polynomial sum.
#' @keywords internal
J <- function(k, s, nu, Omega, correct) {
  if (k == 0) {
    sum <- 10^-correct
  } else {
    sum <- 0
    for (m in 0:(k - 1)) {
      sum <- sum + (Omega[k, m + 1] * s^(-m * (1 - nu)))
    }
  }
  return(sum)
}


#' @title Positive Stable Frailty Distribution
#' @description Computes log-Laplace transform or tau for Positive Stable frailty.
#' @param k Derivative order.
#' @param s Laplace argument.
#' @param nu Stability parameter in (0,1).
#' @param Omega Omega matrix of coefficients.
#' @param what Character; either `"logLT"` for \eqn{\log[ (-1)^k L^{{(k)}}(s) ]} or `"tau"` the Kendall's Tau.
#' @param correct Correction factor for large clusters.
#' @return Numeric value.
#' @keywords internal
fr.possta <- function(k, s, nu, Omega, what = "logLT", correct) {
  if (what == "logLT") {
    res <- k * (log(1 - nu) - nu * log(s)) - s^(1 - nu) +
      log(J(k, s, nu, Omega, correct))
    return(res)
  } else if (what == "tau") {
    return(nu)
  }
}


#' @title Lognormal Frailty Distribution
#' @description Computes the log-Laplace transform or Kendall's tau for the Lognormal frailty.
#' @param k Order of derivative.
#' @param s Laplace argument.
#' @param sigma2 Variance parameter (positive).
#' @param what Character; either `"logLT"` for \eqn{\log[ (-1)^k L^{{(k)}}(s) ]} or `"tau"` the Kendall's Tau.
#' @return Numeric approximation of log-LT or Kendall's tau.
#' @keywords internal
fr.lognormal <- function(k, s, sigma2, what = "logLT") {
  if (what == "logLT") {
    WARN <- getOption("warn")
    options(warn = -1)
    wTilde <- nlm(f = g, p = 0, k = k, s = s, sigma2 = sigma2)$estimate
    options(warn = WARN)
    res <- -g(w = wTilde, k = k, s = s, sigma2 = sigma2) -
      log(sigma2 * g2(w = wTilde, k = k, s = s, sigma2 = sigma2)) / 2
    return(res)
  } else if (what == "tau") {
    intTau <- Vectorize(function(x, intTau.sigma2 = sigma2) {
      x * Lapl(s = x, k = 0, sigma2 = intTau.sigma2) *
        Lapl(s = x, k = 2, sigma2 = intTau.sigma2)
    }, "x")
    
    tauRes <- 4 * integrate(f = intTau, lower = 0, upper = Inf,
                            intTau.sigma2 = sigma2)$value - 1
    return(tauRes)
  }
}

# Internal helper functions for lognormal Laplace approximation
g <- function(w, k, s, sigma2) {
  -k * w + exp(w) * s + w^2 / (2 * sigma2)
}

g1 <- function(w, k, s, sigma2) {
  -k + exp(w) * s + w / sigma2
}

g2 <- function(w, k, s, sigma2) {
  exp(w) * s + 1 / sigma2
}

#' @keywords internal
Lapl <- Vectorize(function(s, k, sigma2) {
  WARN <- getOption("warn")
  options(warn = -1)
  wTilde <- optimize(f = g, c(-1e10, 1e10), maximum = FALSE,
                     k = k, s = s, sigma2 = sigma2)$minimum
  options(warn = WARN)
  
  res <- (-1)^k *
    exp(-g(w = wTilde, k = k, s = s, sigma2 = sigma2)) /
    sqrt(sigma2 * g2(w = wTilde, k = k, s = s, sigma2 = sigma2))
  return(res)
}, "s")
