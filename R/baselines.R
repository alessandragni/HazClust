#' Weibull Baseline Hazard Function
#'
#' Computes the cumulative hazard or log-hazard for a Weibull baseline hazard.
#'
#' @param pars Numeric vector of parameters: \code{pars[1] = rho > 0}, \code{pars[2] = lambda > 0}.
#' @param t Numeric scalar. Time point at which the hazard is evaluated.
#' @param what Character. Either \code{"H"} for cumulative hazard or \code{"lh"} for log-hazard.
#'
#' @return A numeric value of the requested quantity.
#' @keywords internal
weibull <- function(pars, t, what) {
  if (what == "H")
    return(pars[2] * t ^ (pars[1]))
  else if (what == "lh")
    return(log(pars[1]) + log(pars[2]) + ((pars[1] - 1) * log(t)))
}


#' Fréchet / Inverse Weibull Baseline Hazard Function
#'
#' Computes the cumulative hazard or log-hazard for a Fréchet (inverse Weibull) baseline hazard.
#'
#' @param pars Numeric vector of parameters: \code{pars[1] = rho > 0}, \code{pars[2] = lambda > 0}.
#' @param t Numeric scalar. Time point.
#' @param what Character. Either \code{"H"} or \code{"lh"}.
#'
#' @return A numeric value of the requested quantity.
#' @keywords internal
inweibull <- frechet <- function(pars, t, what) {
  if (what == "H")
    return(-log(1 - exp(-pars[2] * (t ^ -pars[1]))))
  else if (what == "lh")
    return(sum(log(pars[1:2])) - log(t) * (pars[1] + 1) -
             log(exp(pars[2] * (t ^ -pars[1])) - 1))
}

#' Exponential Baseline Hazard Function
#'
#' Computes the cumulative hazard or log-hazard for an exponential baseline hazard.
#'
#' @param pars Numeric scalar: \code{lambda > 0}.
#' @param t Numeric scalar. Time point.
#' @param what Character. Either \code{"H"} or \code{"lh"}.
#'
#' @return A numeric value of the requested quantity.
#' @keywords internal
exponential <- function(pars, t, what) {
  if (what == "H")
    return(pars * t)
  else if (what == "lh")
    return(log(pars))
}

#' Gompertz Baseline Hazard Function
#'
#' Computes the cumulative hazard or log-hazard for a Gompertz baseline hazard.
#'
#' @param pars Numeric vector of parameters: \code{pars[1] = gamma > 0}, \code{pars[2] = lambda > 0}.
#' @param t Numeric scalar. Time point.
#' @param what Character. Either \code{"H"} or \code{"lh"}.
#'
#' @return A numeric value of the requested quantity.
#' @keywords internal
gompertz <- function(pars, t, what) {
  if (what == "H")
    return(pars[2] / pars[1] * (exp(pars[1] * t) - 1))
  else if (what == "lh")
    return(log(pars[2]) + pars[1] * t)
}

#' Lognormal Baseline Hazard Function
#'
#' Computes the cumulative hazard or log-hazard for a lognormal baseline hazard.
#'
#' @param pars Numeric vector: \code{pars[1] = mu}, \code{pars[2] = sigma > 0}.
#' @param t Numeric scalar. Time point.
#' @param what Character. Either \code{"H"} or \code{"lh"}.
#'
#' @return A numeric value of the requested quantity.
#' @keywords internal
lognormal <- function(pars, t, what) {
  if (what == "H")
    return(-log(1 - plnorm(t, meanlog = pars[1], sdlog = pars[2])))
  else if (what == "lh")
    return(dlnorm(t, meanlog = pars[1], sdlog = pars[2], log = TRUE) -
             log(1 - plnorm(t, meanlog = pars[1], sdlog = pars[2])))
}

#' Loglogistic Baseline Hazard Function
#'
#' Computes the cumulative hazard or log-hazard for a loglogistic baseline hazard.
#'
#' @param pars Numeric vector: \code{pars[1] = alpha}, \code{pars[2] = kappa > 0}.
#' @param t Numeric scalar. Time point.
#' @param what Character. Either \code{"H"} or \code{"lh"}.
#'
#' @return A numeric value of the requested quantity.
#' @keywords internal
loglogistic <- function(pars, t, what) {
  if (what == "H")
    return(log(1 + exp(pars[1]) * t ^ (pars[2])))
  else if (what == "lh")
    return(pars[1] + log(pars[2]) + (pars[2] - 1) * log(t) -
             log(1 + exp(pars[1]) * t ^ pars[2]))
}

#' Log-Skew Normal Baseline Hazard Function
#'
#' Computes the cumulative hazard or log-hazard for a log-skew normal distribution.
#'
#' @param pars Numeric vector of parameters:
#'   \code{pars[1] = xi}, location;
#'   \code{pars[2] = omega > 0}, scale;
#'   \code{pars[3] = alpha}, shape.
#' @param t Numeric scalar. Time point.
#' @param what Character. Either \code{"H"} or \code{"lh"}.
#'
#' @details
#' The distribution is based on the skew-normal distribution by Azzalini (1985).
#' Requires the \pkg{sn} package for density and CDF calculations.
#'
#' @return A numeric value of the requested quantity.
#' @references Azzalini, A. (1985). A class of distributions which includes the normal ones. Scandinavian Journal of Statistics, 12:171–178.
#' @keywords internal
logskewnormal <- function(pars, t, what) {
  # library(sn)
  dlsn <- Vectorize(function(t, pars) {
    dsn(log(t), xi = pars[1], omega = pars[2], alpha = pars[3]) / t
  }, 't')
  plsn <- Vectorize(function(t, pars) {
    psn(log(t), xi = pars[1], omega = pars[2], alpha = pars[3])
  }, 't')
  
  if (what == "H")
    return(-log(1 - plsn(t, pars = pars)))
  else if (what == "lh")
    return(log(dlsn(t, pars = pars)) - log(1 - plsn(t, pars = pars)))
}
