#' Minus Log-Likelihood and other terms for Clustered Survival Models
#'
#' Internal function to compute the negative log-likelihood for frailty survival models
#' with l1 penalty and similarity regularization.
#'
#' @param p Numeric vector of parameters (frailty, baseline, regression).
#' @param obs A list of observed data (see \code{\link{obsdata_creator}}).
#' @param dist Character; name of the baseline distribution.
#' @param frailty Character; name of the frailty distribution.
#' @param correct Numeric; correction factor for positive stable frailty.
#' @param transform Logical; whether to transform parameters to valid space (default: \code{TRUE}).
#' @param gammapar Numeric; similarity penalty weight.
#' @param S Numeric matrix; similarity matrix for regularization.
#'
#' @return A numeric scalar (negative log-likelihood), with attributes \code{cumhaz}, \code{loghaz}, and \code{logSurv}.
#'
#' @keywords internal

Mloglikelihoodl1SIM <- function(p,
                                obs,
                                dist,
                                frailty,
                                correct,
                                transform = TRUE,
                                gammapar,
                                S) { 
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
  cumhaz <- aggregate(cumhaz[, 1], by = list(cumhaz[, 2]), 
                      FUN = sum)[, 2, drop = FALSE]
  ### NO FRAILTY
  if (frailty == "none") cumhaz <- sum(cumhaz)
  
  
  # ---- log-hazard by cluster --------------------------------------------- #
  loghaz <- NULL
  if (frailty != "none")  {
    loghaz <- matrix(unlist(
      sapply(levels(as.factor(obs$strata)),
             function(x) {
               t(cbind(obs$event[obs$strata == x] * (
                 dist(pars[x, ], obs$time[obs$strata == x],
                      what = "lh") + 
                   as.matrix(obs$x)[
                     obs$strata == x, -1, drop = FALSE] %*% 
                   as.matrix(beta)),
                 obs$cluster[obs$strata == x]))
             })), ncol = 2, byrow = TRUE)
    # loghaz_i <- loghaz[,1] 
    loghaz <- aggregate(loghaz[, 1], by = list(loghaz[, 2]), FUN = sum)[
      , 2, drop = FALSE]
  } else {
    loghaz <- sum(apply(cbind(rownames(pars), pars), 1,
                        function(x) {
                          sum(obs$event[obs$strata == x[1]] * (
                            dist(as.numeric(x[-1]), 
                                 obs$time[obs$strata == x[1]],
                                 what = "lh") + 
                              as.matrix(obs$x[
                                obs$strata == x[1], -1, drop = FALSE]
                              ) %*% as.matrix(beta)))
                        }))
    # loghaz_i <- c(apply(cbind(rownames(pars), pars), 1,
    #                     function(x) {
    #                       obs$event[obs$strata == x[1]] * (
    #                         dist(as.numeric(x[-1]), 
    #                              obs$time[obs$strata == x[1]],
    #                              what = "lh") + 
    #                           as.matrix(obs$x[
    #                             obs$strata == x[1], -1, drop = FALSE]
    #                           ) %*% as.matrix(beta))
    #                     }))
  }
  
  
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
    DIFFlogSurv <- rep(mapply(fr.none, s = cumhaz, what = "logLT") - logSurv, obs$ncl)
  }
  
  
  # ---- X'_gi * beta by cluster --------------------------------------------- #
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
    #Xbeta <- aggregate(Xbeta[, 1], by = list(Xbeta[, 2]), FUN = sum)[
    #  , 2, drop = FALSE]
  } else {
    #Xbeta <- sum(apply(cbind(rownames(pars), pars), 1,
    #                    function(x) {
    #                      sum(
    #                        as.matrix(obs$x[
    #                            obs$strata == x[1], -1, drop = FALSE]
    #                          ) %*% as.matrix(beta))
    #                    }))
    #Xbeta_i <- c(apply(cbind(rownames(pars), pars), 1,
    #                    function(x) {
    #                          as.matrix(obs$x[
    #                            obs$strata == x[1], -1, drop = FALSE]
    #                          ) %*% as.matrix(beta)
    #                    }))
    Xbeta_i <- unlist(lapply(1:nrow(pars), function(i) {
      stratum <- rownames(pars)[i]
      x_sub <- as.matrix(obs$x[obs$strata == stratum, -1, drop = FALSE])
      drop(x_sub %*% as.matrix(beta))  # drop() converts to vector
    }))
    
  }
  
      
  distan = Xbeta_i + 
    as.vector(DIFFlogSurv[match(obs[["cluster"]], unique(obs[["cluster"]]))]) # making ug long as gi
  
  
  # ---- Minus the log likelihood and other terms ----------------------------- #
  Mloglik <- - sum(as.numeric(loghaz[[1]]) + logSurv) + # Minus the log likelihood 
             # + eta*sum(abs(beta))  + # L1-penalization on beta
             + gammapar*sum( (outer(distan, distan, "-"))^2 * S ) # distance/similarity term
    
  attr(Mloglik, "cumhaz") <- as.numeric(cumhaz[[1]])
  attr(Mloglik, "loghaz") <- as.numeric(loghaz[[1]])
  attr(Mloglik, "logSurv") <- logSurv
  
  return(Mloglik)
}


#' Minus Log-Likelihood for Optimization
#'
#' Wrapper for `Mloglikelihoodl1SIM()` to remove attributes and return a pure numeric
#' value, suitable for use in optimizers like `optimx()`.
#'
#' @inheritParams Mloglikelihoodl1SIM
#'
#' @return The numeric value of the minus log-likelihood.
#' @keywords internal
optMloglikelihoodl1SIM <- function(p, obs, dist, frailty, correct, gammapar, S) {
  res <- Mloglikelihoodl1SIM(p = p, obs = obs, dist = dist, 
                             frailty = frailty, correct = correct,
                             # eta = eta, 
                             gammapar = gammapar, S = S)
  as.numeric(res)}


