#' Fit Parametric Frailty Models with Modified Optimization
#'
#' This function fits parametric frailty models for clustered survival data using
#' a modified optimization problem as implemented in \code{Mloglikelihoodl1SIM}.
#' It supports multiple baseline hazard distributions and frailty types, with
#' options for stratification and correction for many events per cluster.
#'
#' @param formula A formula with a survival object created by \code{Surv()} on the left-hand side,
#'   and covariates on the right-hand side.
#' @param cluster Optional. Character string specifying the name of the clustering variable in \code{data}.
#' @param strata Optional. Character string specifying the name of the stratification variable in \code{data}.
#' @param data A data frame containing the variables used in the model.
#' @param inip Optional. Numeric vector of initial values for baseline hazard and regression parameters.
#' @param iniFpar Optional. Initial value(s) for the frailty distribution parameter(s).
#' @param dist Character string specifying the baseline hazard distribution. One of:
#'   \code{"weibull"}, \code{"inweibull"}, \code{"frechet"}, \code{"exponential"},
#'   \code{"gompertz"}, \code{"loglogistic"}, \code{"lognormal"}, \code{"logskewnormal"}.
#' @param frailty Character string specifying the frailty distribution. One of:
#'   \code{"none"}, \code{"gamma"}, \code{"ingau"}, \code{"possta"}, \code{"lognormal"}, \code{"loglogistic"}.
#' @param method Optimization method passed to \code{optimx}, e.g., \code{"nlminb"}.
#' @param maxit Maximum number of iterations for the optimizer.
#' @param Fparscale Scaling factor for the frailty parameter in optimization; the algorithm
#'   optimizes \code{Fpar / Fparscale}.
#' @param showtime Logical; if \code{TRUE}, prints execution time.
#' @param correct Numeric scalar used only when \code{frailty = "possta"}. A correction factor
#'   for the likelihood calculation in cases with many events per cluster.
#' @param gammapar Numeric; additional parameter used in the modified likelihood function.
#' @param S Numeric matrix; additional data-dependent matrix used in the modified likelihood.
#'
#' @return A \code{parfm} class object (a named matrix) with columns:
#'   \item{ESTIMATE}{Parameter estimates for frailty, baseline hazard, and regression coefficients.}
#'   \item{SE}{Standard errors of the estimates.}
#'   \item{p-val}{p-values for regression coefficients (if any).}
#'
#'   The object also contains attributes including the call, convergence status, 
#'   number of iterations, execution time (if requested), log-likelihood at optimum,
#'   baseline cumulative hazard, cluster information, frailty type, and Fisher information matrix.
#'
#' @details
#' The model is fit by maximizing the modified log-likelihood function provided by \code{Mloglikelihoodl1SIM}.
#' If \code{frailty = "none"}, no random effects are included.
#' Stratification allows for different baseline hazards in different strata.
#' The function supports multiple parametric baseline hazards and frailty distributions.
#'
#' @seealso
#' \code{\link{optimx}}, \code{\link{Surv}}, \code{\link{Mloglikelihoodl1SIM}}
#'
#' @examples
#' \dontrun{
#' MYparfm(Surv(time, status) ~ age + sex, cluster = "id", data = kidney,
#'        dist = "weibull", frailty = "gamma",
#'        gammapar = 1, 
#'        S = matrix(1/nrow(kidney), nrow(kidney), nrow(kidney)))
#' }
#'
#' @importFrom survival Surv
#' @importFrom stats model.matrix optimHess pnorm dlnorm integrate nlm plnorm terms
#' @importFrom utils data
#' @importFrom optimx optimx
#' @importFrom sn dsn psn
#' @keywords internal
MYparfm <- function(formula,
                    cluster,
                    strata,
                    data,
                    inip,
                    iniFpar,
                    dist,
                    frailty,
                    method,
                    maxit,
                    Fparscale,
                    showtime,
                    correct,
                    gammapar, 
                    S){
  
  #----- create obsdata ------------------------------------------#
  obsdata = obsdata_creator(formula = formula, data = data, cluster = cluster, strata = strata, frailty = frailty, dist = dist)
  nRpar <- obsdata$nRpar
  nBpar <- obsdata$nBpar
  nFpar <- obsdata$nFpar
  
  #----- Initial parameters -------------------------------------------------#
  if (!is.null(inip)) {
    #if they are specified, then 
    # (1) check the dimension,
    # (2) check whether they lie in their parameter space, and
    # (3) reparametrise them so that they take on values on the whole real line
    
    # (1)
    if (length(inip) != nBpar * obsdata$nstr + nRpar) {
      stop(paste("number of initial parameters 'inip' must be", 
                 nBpar * obsdata$nstr + nRpar))
    }
    p.init <- inip
    
    # (2)-(3)
    if (dist %in% c("exponential", "weibull", "inweibull", "frechet", "gompertz")) {
      # 1st initial par: log(lambda), log(rho), or log(gamma)
      if (any(p.init[1:obsdata$nstr] <= 0)) {
        stop(paste("with that baseline, the 1st parameter has to be > 0"))
      }
      p.init[1:obsdata$nstr] <- log(p.init[1:obsdata$nstr]) 
    }
    if (dist %in% c("weibull", "inweibull", "frechet", "gompertz", 
                    "lognormal", "loglogistic", "logskewnormal")) {
      #2nd initial par: log(lambda), log(lambda), 
      #                 log(sigma), log(kappa), or log(omega)
      if (any(p.init[obsdata$nstr + 1:obsdata$nstr] <= 0)) {
        stop(paste("with that baseline, the 2nd parameter has to be > 0"))
      }
      p.init[obsdata$nstr + 1:obsdata$nstr] <- 
        log(p.init[obsdata$nstr + 1:obsdata$nstr]) 
    }
  } else {
    inires <- optimx(
      par = rep(0, nRpar + nBpar),
      fn = optMloglikelihoodl1SIM, 
      method = method,
      obs = obsdata, dist = dist, frailty = 'none',
      correct = correct,
      gammapar = gammapar, S = S,
      hessian = FALSE,
      control = list(maxit = maxit,
                     starttests = FALSE,
                     trace = 0,
                     dowarn = FALSE)
    )
    p.init <- inires[1:(nRpar + nBpar)]
    rm(inires)
  }
  
  # --- frailty parameters initialisation --- #
  if (frailty == "none") {
    pars <- NULL
  } else if (frailty %in% c("gamma", "ingau", "lognormal")) {
    if (is.null(iniFpar)) {
      iniFpar <- 1
    } else if (iniFpar <= 0) {
      stop("initial heterogeneity parameter (theta) has to be > 0")
    }
    pars <- log(iniFpar)
  } else if (frailty == "possta") {
    if (is.null(iniFpar)) {
      iniFpar <- 0.5
    } else if (iniFpar <= 0 || iniFpar >= 1) {
      stop("initial heterogeneity parameter (nu) must lie in (0, 1)")
    }
    pars <- log(-log(iniFpar))
  }
  
  pars <- c(pars, unlist(p.init))
  res <- NULL
  
  #--------------------------------------------------------------------------#
  #----- Minimise Mloglikelihoodl1SIM() ------------------------------------------#
  #--------------------------------------------------------------------------#
  todo <- expression({
    res <- optimx(
      par = pars, 
      fn = optMloglikelihoodl1SIM, 
      method = method,
      obs = obsdata, dist = dist, frailty = frailty,
      correct = correct,
      gammapar = gammapar, S = S,
      hessian = FALSE,
      control = list(maxit = maxit,
                     starttests = FALSE,
                     trace = 0,
                     dowarn = FALSE)
    )
  })
  if (showtime) {
    extime <- system.time(eval(todo))[1]
  } else {
    eval(todo)
    extime <- NULL
  }
  #--------------------------------------------------------------------------#
  
  if (res$convcode > 0) {
    warning("optimisation procedure did not converge,
              conv = ", bquote(.(res$convergence)), ": see ?optimx for details")
  }
  it <- res$niter   #number of iterations
  lL <- -res$value     # value of the marginal loglikelihood
  if (frailty == "possta") {
    lL <- lL + correct * log(10) * obsdata$ncl
  }
  
  
  #--------------------------------------------------------------------------#
  #----- Recover the estimates ----------------------------------------------#
  #--------------------------------------------------------------------------#
  estim_par <- as.numeric(res[1:(nFpar + nBpar  * obsdata$nstr + nRpar)])
  
  #heterogeneity parameter
  if (frailty %in% c("gamma", "ingau")) {
    theta <- exp(estim_par[1:nFpar])
    sigma2 <- NULL
    nu <- NULL
  } else if (frailty == "lognormal") {
    theta <- NULL
    sigma2 <- exp(estim_par[1:nFpar])
    nu <- NULL
  } else if (frailty == "possta") {
    theta <- NULL
    sigma2 <- NULL
    nu <- exp(-exp(estim_par[1:nFpar]))
  } else if (frailty == "none") {
    theta <- NULL
    sigma2 <- NULL
    nu <- NULL
  }
  
  #baseline hazard parameter(s)
  if (dist == "exponential") {
    lambda <- exp(estim_par[nFpar + 1:obsdata$nstr])
    ESTIMATE <- c(lambda = lambda)
  } else if (dist %in% c("weibull", "inweibull", "frechet")) {
    rho <- exp(estim_par[nFpar + 1:obsdata$nstr])
    lambda <- exp(estim_par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(rho = rho, lambda = lambda)
  } else if (dist == "gompertz") {
    gamma <- exp(estim_par[nFpar + 1:obsdata$nstr])
    lambda <- exp(estim_par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(gamma = gamma, lambda = lambda)
  } else if (dist == "lognormal") {
    mu <- estim_par[nFpar + 1:obsdata$nstr]
    sigma <- exp(estim_par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(mu = mu, sigma = sigma)
  } else if (dist == "loglogistic") {
    alpha <- estim_par[nFpar + 1:obsdata$nstr]
    kappa <- exp(estim_par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(alpha = alpha, kappa = kappa)
  } else if (dist == "logskewnormal") {
    xi <- estim_par[nFpar + 1:obsdata$nstr]
    omega <- exp(estim_par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    alpha <- estim_par[nFpar + 2 * obsdata$nstr + 1:obsdata$nstr]
    ESTIMATE <- c(xi = xi, omega = omega, alpha = alpha)
  }
  
  #regression parameter(s)
  if (nRpar == 0) {
    beta <- NULL
  } else {
    beta <- estim_par[-(1:(nFpar + nBpar * obsdata$nstr))]
    names(beta) <- paste("beta", names(obsdata$x), sep=".")[-1]
  }
  
  #all together
  ESTIMATE <- c(theta = theta,
                sigma2 = sigma2,
                nu = nu,
                ESTIMATE,
                beta = beta)
  #--------------------------------------------------------------------------#
  
  #. print('ok up to here')
  #. #--------------------------------------------------------------------------#
  #. #----- Recover the standard errors ----------------------------------------#
  #. #--------------------------------------------------------------------------#
  #. resHessian <- # attr(res, 'details')[1, 'nhatend'][[1]]
  #.   suppressWarnings(
  #.     optimHess(par = ESTIMATE, fn = Mloglikelihoodl1SIM, 
  #.               obs = obsdata, dist = dist, frailty = frailty,
  #.               correct = correct, transform = FALSE,
  #.               gammapar = gammapar, S = S)
  #.   )    
  #. 
  #. var <- try(diag(solve(resHessian)), silent=TRUE)
  #. if (inherits(var, "try-error") | any(is.nan(var))) {
  #.   warning(var[1])
  #.   STDERR <- rep(NA, nFpar + nBpar * obsdata$nstr + nRpar)
  #.   PVAL <- rep(NA, nFpar + nBpar * obsdata$nstr + nRpar)
  #. } else {
  #.   if (any(var <= 0)) {
  #.     warning(paste("negative variances have been replaced by NAs\n",
  #.                   "Please, try other initial values",
  #.                   "or another optimisation method"))
  #.   }
  #.   
  #.   # heterogeneity (frailty distribution) parameter(s)
  #.   if (frailty %in% c("gamma", "ingau")) {
  #.     seTheta <- sapply(1:nFpar, function(x){
  #.       ifelse(var[x] > 0, sqrt(var[x]), NA)
  #.     })
  #.     seSigma2 <- seNu <- NULL
  #.   } else if (frailty == "lognormal") {
  #.     seSigma2 <- sapply(1:nFpar, function(x){
  #.       ifelse(var[x] > 0, sqrt(var[x]), NA)
  #.     })
  #.     seTheta <- seNu <- NULL
  #.   } else if (frailty == "possta") {
  #.     seNu <- sapply(1:nFpar, function(x){
  #.       ifelse(var[x] > 0, sqrt(var[x]), NA)
  #.     })
  #.     seTheta <- seSigma2 <- NULL
  #.   }
  #.   
  #.   # baseline hazard parameter(s)
  #.   if (dist == "exponential") {
  #.     seLambda <- sapply(1:obsdata$nstr, function(x){
  #.       ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), NA)
  #.     })
  #.     STDERR <- c(seLambda = seLambda)
  #.   } else if (dist %in% c("weibull", "inweibull", "frechet")) {
  #.     seRho <- sapply(1:obsdata$nstr, function(x){
  #.       ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), NA)
  #.     })
  #.     seLambda <- sapply(1:obsdata$nstr, function(x){
  #.       ifelse(var[nFpar + obsdata$nstr + x] > 0, 
  #.              sqrt(var[nFpar + obsdata$nstr + x]), NA)
  #.     })
  #.     STDERR <- c(seRho = seRho, seLambda = seLambda)
  #.   } else if (dist == "gompertz") {
  #.     seGamma <- sapply(1:obsdata$nstr, function(x){
  #.       ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), NA)
  #.     })
  #.     seLambda <- sapply(1:obsdata$nstr, function(x){
  #.       ifelse(var[nFpar + obsdata$nstr + x] > 0,
  #.              sqrt(var[nFpar + obsdata$nstr + x]), NA)
  #.     })
  #.     STDERR <- c(seGamma = seGamma, seLambda = seLambda)
  #.   } else if (dist == "lognormal") {
  #.     seMu <- sapply(1:obsdata$nstr, function(x){
  #.       ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), NA)
  #.     })
  #.     seSigma <- sapply(1:obsdata$nstr, function(x){
  #.       ifelse(var[nFpar + obsdata$nstr + x] > 0,
  #.              sqrt(var[nFpar + obsdata$nstr + x]), NA)
  #.     })
  #.     STDERR <- c(seMu=seMu, seSigma=seSigma)
  #.   } else if (dist == "loglogistic") {
  #.     seAlpha <- sapply(1:obsdata$nstr, function(x){
  #.       ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), NA)
  #.     })
  #.     seKappa <- sapply(1:obsdata$nstr, function(x){
  #.       ifelse(var[nFpar + obsdata$nstr + x] > 0,
  #.              sqrt(var[nFpar + obsdata$nstr + x]), NA)
  #.     })
  #.     STDERR <- c(seAlpha=seAlpha, seKappa=seKappa)
  #.   } else if (dist == "logskewnormal") {
  #.     seXi <- sapply(1:obsdata$nstr, function(x){
  #.       ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), NA)
  #.     })
  #.     seOmega <- sapply(1:obsdata$nstr, function(x){
  #.       ifelse(var[nFpar + obsdata$nstr + x] > 0,
  #.              sqrt(var[nFpar + obsdata$nstr + x]), NA)
  #.     })
  #.     seAlpha <- sapply(1:obsdata$nstr, function(x){
  #.       ifelse(var[nFpar + 2 * obsdata$nstr + x] > 0,
  #.              sqrt(var[nFpar + 2 * obsdata$nstr + x]), NA)
  #.     })
  #.     STDERR <- c(seXi    = seXi, 
  #.                 seOmega = seOmega, 
  #.                 seAlpha = seAlpha)
  #.   }
  #.   
  #.   #regression parameter(s)
  #.   if (nRpar == 0) {
  #.     seBeta <- NULL
  #.   } else {
  #.     seBeta <- numeric(nRpar)
  #.     varBeta <- var[-(1:(nFpar + nBpar * obsdata$nstr))]
  #.     for (i in 1:nRpar) {
  #.       seBeta[i] <- ifelse(varBeta[i] > 0, sqrt(varBeta[i]), NA)
  #.     }
  #.     PVAL <- c(rep(NA, nFpar + nBpar * obsdata$nstr), 
  #.               2 * pnorm(q = -abs(beta / seBeta)))
  #.   }
  #.   
  #.   #all together
  #.   STDERR <- c(STDERR, se.beta = seBeta)
  #.   if (frailty != "none") {
  #.     STDERR <- c(se.theta  = seTheta,
  #.                 se.sigma2 = seSigma2,
  #.                 se.nu     = seNu,
  #.                 STDERR)
  #.   }
  #. }
  
  
  
  #--------------------------------------------------------------------------#
  #----- Output -------------------------------------------------------------#
  #--------------------------------------------------------------------------#
  resmodel <- cbind(ESTIMATE = ESTIMATE) #, 
                    # SE       = STDERR)
  rownames(resmodel) <- gsub("beta.","", rownames(resmodel))
  
  # if (nRpar > 0)
  #   resmodel <- cbind(resmodel, "p-val" = PVAL)
  
  class(resmodel) <- c("MYparfm", class(resmodel))
  
  ### - Checks - #############################################################
  Call <- match.call()
  if (!match("formula", names(Call), nomatch = 0))
    stop("A formula argument is required")
  
  Terms <- terms(formula, data = data)
  ###################################################### - End of Checks - ###
  
  attributes(resmodel) <- c(attributes(resmodel), list(
    call        = Call,
    convergence = res$convergence,
    it          = it,
    extime      = extime,
    nobs        = nrow(data),
    shared      = (nrow(data) > obsdata$ncl),
    loglik      = lL,
    dist        = dist,
    cumhaz      = attributes(Mloglikelihoodl1SIM(p       = estim_par,
                                            obs     = obsdata, 
                                            dist    = dist, 
                                            frailty = frailty,
                                            correct = correct,
                                            gammapar = gammapar, S = S))$cumhaz,
    di          = obsdata$di,
    dq          = obsdata$dq,
    dqi         = obsdata$dqi,
    frailty     = frailty,
    clustname   = cluster,
    stratname   = strata,
    correct     = correct,
    formula     = as.character(Call[match("formula", names(Call), 
                                          nomatch = 0)]),
    terms       = attr(Terms, "term.labels") #,
    #FisherI     = resHessian
  ))
  if (frailty != "none") {
    names(attr(resmodel, "cumhaz")) <-
      names(attr(resmodel, "di")) <- 
      unique(obsdata$cluster)
  }
  if (showtime)
    cat("\nExecution time:", extime, "second(s) \n")
  
  return(resmodel)
}
