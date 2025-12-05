#' Compute Distance Matrix for Silhouette
#'
#' Fits parametric frailty models for clustered survival data by constructing a modified 
#' distance matrix used in the penalized optimization problem (as in \code{Mloglikelihoodl1SIM}).
#' This matrix incorporates covariate information, frailty effects, and an optional penalty matrix.
#'
#' @param formula A \code{\link{formula}} specifying the survival model. The left-hand side 
#'   should be a survival object created by \code{\link[survival]{Surv}}; the right-hand side 
#'   contains covariates.
#' @param cluster Optional character string naming the clustering variable in \code{data}.
#' @param strata Optional character string naming the stratification variable in \code{data}.
#' @param data A \code{data.frame} containing all variables referenced in \code{formula}.
#' @param dist Character string specifying the baseline hazard distribution. One of 
#'   \code{"weibull"}, \code{"inweibull"}, \code{"frechet"}, \code{"exponential"}, 
#'   \code{"gompertz"}, \code{"loglogistic"}, \code{"lognormal"}, or \code{"logskewnormal"}.
#' @param frailty Character string specifying the frailty distribution. One of 
#'   \code{"none"}, \code{"gamma"}, \code{"ingau"}, \code{"possta"}, \code{"lognormal"}, or \code{"loglogistic"}.
#' @param p Numeric vector of model parameters, including frailty, baseline hazard, and regression parameters.
#' @param transform Logical indicating whether parameter transformations (e.g., exponentiation) 
#'   should be applied to improve optimization stability. Default is \code{TRUE}.
#' @param correct Optional numeric correction factor applied only when \code{frailty = "possta"} 
#'   to adjust likelihood calculations for clusters with many events.
#'
#' @return A numeric matrix representing the modified squared distance matrix used in penalized likelihood optimization.
#'
#' @details
#' Frailty distributions supported include Gamma, Inverse Gaussian, Positive Stable, Lognormal, and Loglogistic, 
#' each with specific parameter transformations applied if \code{transform} is \code{TRUE}.
#'
#' The baseline hazard parameters are parsed according to the specified distribution \code{dist}, 
#' with support for multiple strata.
#'
#' @seealso \code{\link{Mloglikelihoodl1SIM}}, \code{\link[survival]{Surv}}, frailty distribution functions (\code{fr.gamma}, \code{fr.ingau}, etc.)
#'
#' @importFrom stats dist aggregate
#' @keywords internal
computeDistforSilhouette <- function(formula, cluster, strata, data, dist,
                                     frailty, p,
                                     transform = TRUE,
                                     correct){
  
  #----- create obsdata ------------------------------------------#
  
  obs = obsdata_creator(formula = formula, data = data, cluster = cluster, 
                        strata = strata, frailty = frailty, dist = dist)
  nRpar <- obs$nRpar
  nBpar <- obs$nBpar
  nFpar <- obs$nFpar
  
  
  # ---- Assign the number of frailty parameters 'obs$nFpar' --------------- #
  # ---- and compute Sigma for the Positive Stable frailty ----------------- #
  
  if (frailty %in% c("gamma", "ingau")) {
    theta <- ifelse(transform, exp(p[1]), p[1])
  } else if (frailty == "lognormal") {
    sigma2 <- ifelse(transform, exp(p[1]), p[1])
  } else if (frailty == "possta") {
    nu <- ifelse(transform, exp(-exp(p[1])), p[1])
    D <- max(obs$dqi)
    Omega <- Omega(D, correct = correct, nu = nu)
  }
  
  
  # ---- Baseline hazard --------------------------------------------------- #
  if (frailty == 'none') obs$nFpar <- 0
  
  # baseline parameters
  if (dist %in% c("weibull", "inweibull", "frechet")) {
    if (transform) {
      pars <- cbind(rho    = exp(p[obs$nFpar + 1:obs$nstr]),
                    lambda = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]))
    } else {
      pars <- cbind(rho    = p[obs$nFpar + 1:obs$nstr],
                    lambda = p[obs$nFpar + obs$nstr + 1:obs$nstr])
    }
    beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
  } else if (dist == "exponential") {
    if (transform) {
      pars <- cbind(lambda = exp(p[obs$nFpar + 1:obs$nstr]))
    } else {
      pars <- cbind(lambda = p[obs$nFpar + 1:obs$nstr])
    }
    beta <- p[-(1:(obs$nFpar + obs$nstr))]
  } else if (dist == "gompertz") {
    if (transform) {
      pars <- cbind(gamma  = exp(p[obs$nFpar + 1:obs$nstr]),
                    lambda = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]))
    } else {
      pars <- cbind(gamma  = p[obs$nFpar + 1:obs$nstr],
                    lambda = p[obs$nFpar + obs$nstr + 1:obs$nstr])
    }  
    beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
  } else if (dist == "lognormal") {
    if (transform) {
      pars <- cbind(mu    = p[obs$nFpar + 1:obs$nstr],
                    sigma = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]))
    } else {
      pars <- cbind(mu    = p[obs$nFpar + 1:obs$nstr],
                    sigma = p[obs$nFpar + obs$nstr + 1:obs$nstr])
    }
    beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
  } else if (dist == "loglogistic") {
    if (transform) {
      pars <- cbind(alpha = p[obs$nFpar + 1:obs$nstr],
                    kappa = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]))
    } else  {
      pars <- cbind(alpha = p[obs$nFpar + 1:obs$nstr],
                    kappa = p[obs$nFpar + obs$nstr + 1:obs$nstr])
    }
    beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
  } else if (dist == "logskewnormal") {
    if (transform) {
      pars <- cbind(mu    = p[obs$nFpar + 1:obs$nstr],
                    sigma = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]),
                    alpha = exp(p[obs$nFpar + 2 * obs$nstr + 1:obs$nstr]))
    } else {
      pars <- cbind(mu    = p[obs$nFpar + 1:obs$nstr],
                    sigma = p[obs$nFpar + obs$nstr + 1:obs$nstr],
                    alpha = p[obs$nFpar + 2 * obs$nstr + 1:obs$nstr])
    }
    beta <- p[-(1:(obs$nFpar + 3 * obs$nstr))]
  }
  rownames(pars) <- levels(as.factor(obs$strata))
  
  # baseline: from string to the associated function
  dist <- eval(parse(text = dist))
  
  
  # ---- Cumulative Hazard by cluster and by strata ------------------------- #
  
  cumhaz <- NULL
  cumhaz <- matrix(unlist(
    sapply(levels(as.factor(obs$strata)),
           function(x) {t(
             cbind(dist(pars[x, ], obs$time[obs$strata == x], what = "H"
             ) * exp(as.matrix(obs$x)[
               obs$strata == x, -1, drop = FALSE] %*% as.matrix(beta)),
             obs$cluster[obs$strata == x]))
           })), ncol = 2, byrow = TRUE)
  # cumhaz_i <- cumhaz[,1]
  cumhaz <- aggregate(cumhaz[, 1], by = list(cumhaz[, 2]), 
                      FUN = sum)[, 2, drop = FALSE]
  ### NO FRAILTY
  if (frailty == "none") cumhaz <- sum(cumhaz)
  
  
  
  # ---- log[ (-1)^di L^(di)(cumhaz) ]-------------------------------------- #
  logSurv <- NULL
  if (frailty == "gamma") {
    logSurv <- mapply(fr.gamma, 
                      k = obs$di, s = as.numeric(cumhaz[[1]]), 
                      theta = rep(theta, obs$ncl), 
                      what = "logLT") 
  } else if (frailty == "ingau") {
    logSurv <- mapply(fr.ingau, 
                      k = obs$di, s = as.numeric(cumhaz[[1]]), 
                      theta = rep(theta, obs$ncl), 
                      what = "logLT") 
  } else if (frailty == "possta") {
    logSurv <- sapply(1:obs$ncl, 
                      function(x) fr.possta(k = obs$di[x], 
                                            s = as.numeric(cumhaz[[1]])[x], 
                                            nu = nu, Omega = Omega, 
                                            what = "logLT",
                                            correct = correct))
  } else if (frailty == "lognormal") {
    logSurv <- mapply(fr.lognormal, 
                      k = obs$di, s = as.numeric(cumhaz[[1]]), 
                      sigma2 = rep(sigma2, obs$ncl), 
                      what = "logLT")
  } else if (frailty == "none") {
    logSurv <- mapply(fr.none, s = cumhaz, what = "logLT")
  }
  
  
  # ---- log(\hat u_g) = log[ (-1)^(di+1) L^(di+1)(cumhaz) ] - log[ (-1)^(di) L^(di)(cumhaz) ] --- #
  DIFFlogSurv <- NULL
  if (frailty == "gamma") {
    DIFFlogSurv <- mapply(fr.gamma, 
                          k = obs$di+1, s = as.numeric(cumhaz[[1]]), 
                          theta = rep(theta, obs$ncl), 
                          what = "logLT") - logSurv
  } else if (frailty == "ingau") {
    DIFFlogSurv <- mapply(fr.ingau, 
                          k = obs$di+1, s = as.numeric(cumhaz[[1]]), 
                          theta = rep(theta, obs$ncl), 
                          what = "logLT") - logSurv
  } else if (frailty == "possta") {
    DIFFlogSurv <- sapply(1:obs$ncl, 
                          function(x) fr.possta(k = obs$di[x]+1, 
                                                s = as.numeric(cumhaz[[1]])[x], 
                                                nu = nu, Omega = Omega, 
                                                what = "logLT",
                                                correct = correct)) - logSurv
  } else if (frailty == "lognormal") {
    DIFFlogSurv <- mapply(fr.lognormal, 
                          k = obs$di+1, s = as.numeric(cumhaz[[1]]), 
                          sigma2 = rep(sigma2, obs$ncl), 
                          what = "logLT") - logSurv
  } else if (frailty == "none") {
    DIFFlogSurv <- mapply(fr.none, s = cumhaz, what = "logLT") - logSurv
  }
  
  
  # ---- X'_gi * beta and Xonly_i by cluster --------------------------------------------- #
  Xbeta <- Xbeta_i <- NULL
  if (frailty != "none")  {
    Xbeta <- matrix(unlist(
      sapply(levels(as.factor(obs$strata)),
             function(x) {
               t(cbind( as.matrix(obs$x)[
                 obs$strata == x, -1, drop = FALSE] %*% 
                   as.matrix(beta),
                 obs$cluster[obs$strata == x]))
             })), ncol = 2, byrow = TRUE)
    Xbeta_i <- Xbeta[,1] 
  } else {
    Xbeta_i <- unlist(lapply(1:nrow(pars), function(i) {
      stratum <- rownames(pars)[i]
      x_mat <- as.matrix(obs$x[obs$strata == stratum, -1, drop = FALSE])
      drop(x_mat %*% as.matrix(beta))
    }))
  }
  
  dist2 = Xbeta_i + 
    + as.vector(DIFFlogSurv[match(obs[["cluster"]], unique(obs[["cluster"]]))]) # making ug as long as gi
  
  D2 <- (outer(dist2, dist2, "-"))^2
  
  return(D2)
}