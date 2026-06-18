# Two-Step Baseline: Parametric Frailty Model + Clustering
# Load packages
library(survival, lib.loc="/u/ragni/Rlibs")
library(cluster,  lib.loc="/u/ragni/Rlibs")
# library(parfm,    lib.loc="/u/ragni/Rlibs")
library(igraph, lib.loc="/u/ragni/Rlibs")
library(kernlab,  lib.loc="/u/ragni/Rlibs")
library(optimx, lib.loc="/u/ragni/Rlibs")

# Set working directory
setwd("/u/ragni/HazClust")
r_files <- list.files("R", full.names = TRUE, pattern = "\\.[rR]$")
sapply(r_files, source)
source("/u/ragni/HazClust/simulations/enhanced_simulate_data.R")

# Parse input parameters
args     <- commandArgs(trailingOnly = TRUE)
seed <- as.numeric(args[1])
gammapar <- 0
k <- 20
c <- as.numeric(args[2])
censor <- args[3]
scenario <- args[4]  # NEW
tabfig <- args[5]  # NEW
frailty_intensity <- as.numeric(args[6])  # NEW

# ----- Define Scenarios -----
scenario_configs <- list(
  
  # Baseline: correctly specified
  baseline = list(
    dgp_frailty = "gamma",
    dgp_baseline = "weibull",
    fit_frailty = "gamma",
    fit_baseline = "weibull",
    separation = "high",
    balance = "balanced"
  ),
  
  # Misspecification: Frailty distribution
  misspec_frailty = list(
    dgp_frailty = "lognormal",
    dgp_baseline = "weibull",
    fit_frailty = "gamma",
    fit_baseline = "weibull",
    separation = "high",
    balance = "balanced"
  ),
  
  # Misspecification: Baseline hazard
  misspec_baseline = list(
    dgp_frailty = "gamma",
    dgp_baseline = "gompertz",
    fit_frailty = "gamma",
    fit_baseline = "weibull",
    separation = "high",
    balance = "balanced"
  ),
  
  # Misspecification: Both
  misspec_both = list(
    dgp_frailty = "lognormal",
    dgp_baseline = "loglogistic",
    fit_frailty = "gamma",
    fit_baseline = "weibull",
    separation = "high",
    balance = "balanced"
  ),
  
  # Weak separation
  weak_separation = list(
    dgp_frailty = "gamma",
    dgp_baseline = "weibull",
    fit_frailty = "gamma",
    fit_baseline = "weibull",
    separation = "low",
    balance = "balanced"
  ),
  
  # Medium separation
  medium_separation = list(
    dgp_frailty = "gamma",
    dgp_baseline = "weibull",
    fit_frailty = "gamma",
    fit_baseline = "weibull",
    separation = "medium",
    balance = "balanced"
  ),
  
  # Imbalanced (moderate)
  imbalanced_moderate = list(
    dgp_frailty = "gamma",
    dgp_baseline = "weibull",
    fit_frailty = "gamma",
    fit_baseline = "weibull",
    separation = "high",
    balance = "imbalanced_moderate"
  ),
  
  # Imbalanced (severe)
  imbalanced_severe = list(
    dgp_frailty = "gamma",
    dgp_baseline = "weibull",
    fit_frailty = "gamma",
    fit_baseline = "weibull",
    separation = "high",
    balance = "imbalanced_severe"
  ),
  
  # Worst case
  worst_case = list(
    dgp_frailty = "lognormal",
    dgp_baseline = "gompertz",
    fit_frailty = "gamma",
    fit_baseline = "weibull",
    separation = "medium",
    balance = "imbalanced_severe"
  )
)

# Get configuration for current scenario
config <- scenario_configs[[scenario]]

# Set fitting parameters based on scenario
dist <- config$fit_baseline
frailty <- config$fit_frailty

#### DATA SIMULATION ####
sim_data <- paper_simulate_data_enhanced(
  seed = seed, 
  verbose = FALSE, 
  dgp_frailty_type = config$dgp_frailty,
  dgp_baseline_type = config$dgp_baseline,
  censoring_mode = censor,
  cluster_separation = config$separation,
  cluster_balance = config$balance,
  frailty_intensity = frailty_intensity,
  mismatch = TRUE
)

# -------------------------------------------------------------------------
#### STEP 1: FIT STANDARD PARAMETRIC FRAILTY MODEL (I still use HazClust but penal to zero) ####
# -------------------------------------------------------------------------
cluster = "group"
data = sim_data
formula = Surv(time, status) ~ X1
strata = NULL
transform = TRUE


frailty_model <- tryCatch({
  HazClust(
    formula,
    cluster   = cluster,
    strata    = NULL,
    data      = data,
    inip      = NULL,
    iniFpar   = NULL,
    dist      = dist,
    frailty   = frailty,
    method    = "Nelder-Mead",
    maxitparfm = 500,
    maxit     = 200,
    tolS      = 1e-2,
    tolll     = 1,
    Fparscale = 1,
    showtime  = FALSE,
    correct   = 0,
    transform = TRUE,
    c = c,
    k = k,
    gammapar = gammapar,
    lambda0 = NULL,
    S0 = NULL,
    trace = TRUE
  )
}, error = function(e) {
  message(paste("Error at seed =", seed, "c =", c, "k =", k, "gammapar =", gammapar, ":", e$message))
  return(NULL)
}
)



computeetahat <- function(formula, cluster, strata, data, dist,
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
  
  return(dist2)
}


# -------------------------------------------------------------------------
#### STEP 2: EXTRACT RISK SCORES ####
# -------------------------------------------------------------------------
if (!is.null(frailty_model)) {
  
  # eta_hat = x'*beta + log(u_hat): same metric as eq. (2.5) in the paper
  sim_data$eta_hat <- computeetahat(formula = formula, cluster = cluster, 
                                    strata = strata, data = data, 
                                    dist = dist, frailty = frailty, 
                                    p = frailty_model$estimate,
                                    transform = transform)

  # Pairwise squared distance matrix on eta_hat
  D_mat  <- as.matrix(dist(sim_data$eta_hat))^2   # N x N
  
  # -----------------------------------------------------------------------
  #### STEP 3a: K-MEANS CLUSTERING ####
  # -----------------------------------------------------------------------
  # set.seed(seed)
  # km_fit <- kmeans(sim_data$eta_hat, centers = c, nstart = 20, iter.max = 300)
  # sim_data$cluster_kmeans <- km_fit$cluster
  
  # ward.D2
  # set.seed(seed)
  # d <- dist(sim_data$eta_hat)
  # hc <- hclust(d, method = "ward.D2")
  # sim_data$cluster_hc_ward <- cutree(hc, k = c)

  # -----------------------------------------------------------------------
  #### STEP 3b: HIERARCHICAL CLUSTERING ####
  # -----------------------------------------------------------------------
  set.seed(seed)
  d <- dist(sim_data$eta_hat)
  hc <- hclust(d, method = "single")
  sim_data$cluster_hc <- cutree(hc, k = c)
  
  # -----------------------------------------------------------------------
  #### STEP 4: SAVE RESULTS ####
  # -----------------------------------------------------------------------
  out_dir <- paste0("/u/ragni/HazClust/simulations/output/tabS62/", scenario)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  name <- paste0(out_dir, "/twostep_seed", seed,
                 "_c", c,
                 "_frailty", frailty, "_censor", censor, "_frailty_intensity", frailty_intensity, ".Rds")
  
  saveRDS(list(
    # identifiers
    seed          = seed,
    scenario      = scenario,
    c             = c,
    # parameter estimates
    estimate      = coef(frailty_model),
    true = sim_data$cluster,
    # k-means
    #clusters_kmeans    = sim_data$cluster_kmeans,
    # hierarchical
    clusters_hc_sin        = sim_data$cluster_hc,
    #clusters_hc_ward        = sim_data$cluster_hc_ward,
    # scenario metadata
    dgp_frailty    = config$dgp_frailty,
    dgp_baseline   = config$dgp_baseline,
    fit_frailty    = config$fit_frailty,
    fit_baseline   = config$fit_baseline,
    separation     = config$separation,
    balance        = config$balance,
    censoring_rate = mean(data$status == 0)), name)
}