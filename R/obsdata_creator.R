#' Prepare Data for Parametric Frailty Models
#'
#' Internal utility to extract and structure survival data with optional frailty and stratification.
#'
#' @param formula A \code{Surv()}-based formula of the form \code{Surv(time, event)} or \code{Surv(entry, time, event)}.
#' @param data A \code{data.frame} containing variables in the formula and optionally the \code{cluster} and \code{strata}.
#' @param cluster Optional character string naming the clustering variable.
#' @param strata Optional character string naming the stratification variable.
#' @param frailty A character string specifying the frailty distribution: \code{"none"}, \code{"gamma"}, \code{"ingau"}, \code{"possta"}, or \code{"lognormal"}.
#' @param dist A character string for the baseline distribution: \code{"exponential"}, \code{"weibull"}, \code{"inweibull"}, \code{"frechet"}, \code{"gompertz"}, \code{"lognormal"}, \code{"loglogistic"}, or \code{"logskewnormal"}.
#'
#' @return A list with time/event/covariate matrices, cluster/strata indicators, event counts, and parameter dimensions.
#'
#' @keywords internal
#' @importFrom stats model.matrix aggregate as.formula xtabs

obsdata_creator <- function(formula, data, cluster, strata, frailty, dist){
  obsdata <- NULL
  #time
  if (length(formula[[2]]) == 3) {# --> without left truncation
    obsdata$time <- eval(
      formula[[2]][[2]], 
      envir = data
    )
    obsdata$event <- eval(
      formula[[2]][[3]],
      envir = data
    )
  } else if (length(formula[[2]]) == 4) {# --> with left truncation
    obsdata$trunc <- eval(
      formula[[2]][[2]],
      envir = data
    )    
    obsdata$time <- eval(
      formula[[2]][[3]] ,
      envir = data
    )
    obsdata$event <- eval(
      formula[[2]][[4]], 
      envir = data
    )
  }
  if (!all(levels(as.factor(obsdata$event)) %in% 0:1)) {
    stop(paste("The status indicator 'event' in the Surv object",
               "in the left-hand side of the formula object",
               "must be either 0 (no event) or 1 (event)."))
  }
  
  #covariates (an intercept is automatically added)
  obsdata$x <- as.data.frame(model.matrix(formula, data=data))
  
  #cluster
  if (is.null(cluster)) {
    if (frailty != "none") {
      stop(paste("if you specify a frailty distribution,\n",
                 "then you have to specify the cluster variable as well"))
    } else {
      obsdata$cluster <- rep(1, nrow(data))
    }
    #number of clusters
    obsdata$ncl <- 1
    #number of events in each cluster
    obsdata$di <- sum(obsdata$event)
  } else {
    if (! cluster %in% names(data)) {
      stop(paste0("object '", cluster, "' not found"))
    }
    obsdata$cluster <- eval(#parse(text=paste0("data$", 
      as.name(cluster),
      envir = data #))
    )
    #number of clusters
    obsdata$ncl <- length(levels(as.factor(obsdata$cluster)))
    #number of events in each cluster
    obsdata$di <- aggregate(obsdata$event,
                            by=list(obsdata$cluster), 
                            FUN=sum)[,, drop=FALSE]
    cnames <- obsdata$di[,1]
    obsdata$di <- as.vector(obsdata$di[,2])
    names(obsdata$di) <- cnames
  }
  
  #strata
  if (is.null(strata)) {
    obsdata$strata <- rep(1, length(obsdata$time))
    #number of strata
    obsdata$nstr <- 1
    #number of events in each stratum
    obsdata$dq <- sum(obsdata$event)
  } else {
    if (!strata %in% names(data)) {
      stop(paste0("object '", strata, "' not found"))
    }
    obsdata$strata <- eval(
      as.name(strata), 
      envir = data
    )
    #number of strata
    obsdata$nstr <- length(levels(as.factor(obsdata$strata)))
    #number of events in each stratum
    obsdata$dq <- aggregate(obsdata$event, 
                            by = list(obsdata$strata), 
                            FUN = sum)[, , drop = FALSE]
    snames <- obsdata$dq[,1]
    obsdata$dq <- as.vector(obsdata$dq[,2])
    names(obsdata$dq) <- snames
  }
  
  #cluster+strata
  if (!is.null(cluster) && !is.null(strata)) {
    #number of events in each cluster for each stratum
    obsdata$dqi <- xtabs(x~Group.1+Group.2, data = aggregate(
      obsdata$event, 
      by = list(obsdata$cluster, obsdata$strata), 
      FUN = sum))
    dimnames(obsdata$dqi) <- list(cluster = dimnames(obsdata$dqi)[[1]], 
                                  strata  = dimnames(obsdata$dqi)[[2]])
  } else if (!is.null(cluster)) {
    obsdata$dqi <- obsdata$di
  } else if (!is.null(strata)) {
    obsdata$dqi <- obsdata$dq
  } else {
    obsdata$dqi <- sum(obsdata$event)
  }
  
  
  
  #----- Dimensions ---------------------------------------------------------#
  #nFpar: number of heterogeneity parameters
  if (frailty == "none") {
    nFpar <- 0
  } else if (frailty %in% c("gamma", "ingau", "possta", "lognormal")) {
    nFpar <- 1
  }
  obsdata$nFpar <- nFpar
  
  #nBpar: number of parameters in the baseline hazard
  if (dist == "exponential") {
    nBpar <- 1
  } else if (dist %in% c("weibull", "inweibull", "frechet", "gompertz",
                         "lognormal", "loglogistic")) {
    nBpar <- 2
  } else if (dist %in% c("logskewnormal")) {
    nBpar <- 3
  }
  obsdata$nBpar <- nBpar
  
  #nRpar: number of regression parameters
  nRpar <- ncol(obsdata$x) - 1
  obsdata$nRpar <- nRpar  
  obsdata
}

