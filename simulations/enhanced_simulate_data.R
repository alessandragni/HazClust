# Enhanced Simulation: Adding misspecification scenarios for reviewer response

paper_simulate_data_enhanced <- function(
    seed,
    verbose = FALSE,
    
    # --- Data generation parameters ---
    dgp_frailty_type,      # "gamma" or "lognormal" or "ingau" or "none"
    dgp_baseline_type,     # "weibull" or "gompertz" or "loglogistic"
    censoring_mode,        # "administrative" or "normal"
    
    # --- Misspecification scenarios ---
    cluster_separation = "high",  # "high", "medium", "low"
    cluster_balance = "balanced", # "balanced", "imbalanced_moderate", "imbalanced_severe"
    
    # --- Optional: model fitting parameters (what you assume when fitting) ---
    fit_frailty_type = NULL,    # If NULL, assumes same as dgp_frailty_type
    fit_baseline_type = NULL,    # If NULL, assumes same as dgp_baseline_type
    
    frailty_intensity = 0.5,
    mismatch = FALSE
){
  set.seed(seed)
  
  # --- Basic parameters ---
  N_clusters   <- 3
  N_groups     <- 10
  n_per_group  <- 50
  N            <- N_groups * n_per_group
  beta         <- log(2)   # scalar
  eta_sd       <- 1
  theta        <- frailty_intensity # default frailty variance
  
  # --- Cluster separation scenarios ---
  if (cluster_separation == "high") {
    cluster_centers <- c(-8, 0, 8)
  } else if (cluster_separation == "medium") {
    cluster_centers <- c(-6, 0, 6)
  } else if (cluster_separation == "low") {
    cluster_centers <- c(-4, 0, 4)
  }
  
  # --- Cluster balance scenarios ---
  if (cluster_balance == "balanced") {
    cluster_probs <- c(1/3, 1/3, 1/3)
    cluster_assignments <- sample(1:N_clusters, N, replace = TRUE) #, prob = cluster_probs)
  } else if (cluster_balance == "imbalanced_moderate") {
    cluster_probs <- c(0.5, 0.3, 0.2)
    cluster_assignments <- sample(1:N_clusters, N, replace = TRUE, prob = cluster_probs)
  } else if (cluster_balance == "imbalanced_severe") {
    cluster_probs <- c(0.7, 0.2, 0.1)
    cluster_assignments <- sample(1:N_clusters, N, replace = TRUE, prob = cluster_probs)
  }
  
  # ----- 1. Cluster assignments for latent eta -----
  #cluster_assignments <- sample(1:N_clusters, N, replace = TRUE, prob = cluster_probs)
  
  # ----- 2. Simulate latent predictor (eta) -----
  eta <- rnorm(N, mean = cluster_centers[cluster_assignments], sd = eta_sd)
  
  # ----- 3. Simulate group-level frailties (DGP) -----
  group_ids <- rep(1:N_groups, each = n_per_group)
  
  if (dgp_frailty_type == "gamma") {
    # Gamma frailty with mean 1 and var = theta
    shape <- 1/theta
    scale <- theta
    u_groups <- rgamma(N_groups, shape = shape, scale = scale)
    log_u_groups <- log(u_groups)
  } else if (dgp_frailty_type == "lognormal") {
    # Log-normal frailty: log(u) ~ N(-theta/2, theta) to have E[u]=1
    log_u_groups <- rnorm(N_groups, mean = -theta/2, sd = sqrt(theta))
    u_groups <- exp(log_u_groups)
  } else if (dgp_frailty_type == "ingau") {
    # Inverse Gaussian frailty (less common, for misspecification)
    # Using statmod::rinvgauss if available, or approximation
    if (requireNamespace("statmod", quietly = TRUE)) {
      # mu = 1, lambda controls variance
      lambda <- 1/theta
      u_groups <- statmod::rinvgauss(N_groups, mean = 1, shape = lambda)
    } else {
      # Fallback to gamma if statmod not available
      warning("statmod package not available, using gamma instead of inverse Gaussian")
      shape <- 1/theta
      scale <- theta
      u_groups <- rgamma(N_groups, shape = shape, scale = scale)
    }
    log_u_groups <- log(u_groups)
  } else { # none
    log_u_groups <- rep(0, N_groups)
    u_groups <- rep(1, N_groups)
  }
  
  log_u <- log_u_groups[group_ids]
  u     <- u_groups[group_ids]
  
  # ----- 4. Covariate from latent predictors -----
  if (mismatch == TRUE){
    eps_x <- rnorm(N, mean = 0, sd = 1.5) 
    x <- (eta - log_u) / beta + eps_x
  } else {
    x <- (eta - log_u) / beta
  }
  # x <- (eta - log_u) / beta
  X <- matrix(x, ncol = 1)
  colnames(X) <- "Z1"
  
  # ----- 5. Baseline hazard (DGP) via inverse cumulative hazard -----
  U <- runif(N)
  # if (mismatch == TRUE){
  #   beta1 <- beta/100   
  #   lp <- beta1 * x + log_u
  # } else {
  #   lp <- beta * x + log_u
  # }
  lp <- beta * x + log_u
  s <- -log(U) / exp(lp)
  
  if (dgp_baseline_type == "weibull") {
    # pars[1]=rho, pars[2]=lambda
    # rho <- 2.5; lambda <- 0.01
    # t_event <- (s / lambda)^(1 / rho)
    xi <- 0.01; rho <- 2.5
    t_event <- (s^(1/rho)) / xi
    
  } else if (dgp_baseline_type == "gompertz") {
    # pars[1]=gamma (alpha), pars[2]=lambda (gamma)
    alpha_param <- 0.05; gamma_param <- 0.01
    t_event <- (1/alpha_param) * log(1 + alpha_param * s / gamma_param)
    
  } else if (dgp_baseline_type == "loglogistic") {
    # Coherent with your loglogistic() function: H(t) = log(1 + exp(pars[1]) * t^pars[2])
    # Inverse: t = [ (exp(s) - 1) / exp(pars[1]) ]^(1/pars[2])
    pars1 <- -8  # log(scale)
    pars2 <- 2   # shape
    t_event <- ( (exp(s) - 1) / exp(pars1) )^(1 / pars2)
  }
  
  # ----- 6. Censoring -----
  lambda_cens <- 0.003
  t_censor_exp <- rexp(N, rate = lambda_cens)
  admin_censor_time <- 100
  censor_mu <- 130
  censor_sd <- 15
  t_censor_norm <- rnorm(N, mean = censor_mu, sd = censor_sd)
  censor_unif_lower <- 80
  censor_unif_upper <- 160
  t_censor_unif <- runif(N, min = censor_unif_lower, max = censor_unif_upper)
  
  if (censoring_mode == "exponential") {
    t_censor <- t_censor_exp
  } else if (censoring_mode == "administrative") {
    t_censor <- rep(admin_censor_time, N)
  } else if (censoring_mode == "normal") {
    t_censor <- t_censor_norm 
  } else if (censoring_mode == "uniform") {
    t_censor <- t_censor_unif
  } else {
    t_censor <- rep(Inf, N)
  }
  
  time_obs <- pmin(t_event, t_censor)
  event    <- as.integer(t_event <= t_censor)
  
  # ----- 7. Assemble output -----
  simdat <- data.frame(
    id = 1:N,
    cluster = cluster_assignments,
    group = group_ids,
    X1 = x,
    log_u = log_u,
    u = u,
    eta = eta,
    true_time = t_event,
    censor_time = t_censor,
    time = time_obs,
    status = event
  )
  
  # Store metadata about DGP
  attr(simdat, "dgp_frailty") <- dgp_frailty_type
  attr(simdat, "dgp_baseline") <- dgp_baseline_type
  attr(simdat, "cluster_separation") <- cluster_separation
  attr(simdat, "cluster_balance") <- cluster_balance
  attr(simdat, "cluster_centers") <- cluster_centers
  attr(simdat, "cluster_probs") <- cluster_probs
  
  if(verbose){
    if(requireNamespace("cowplot", quietly = TRUE) && 
       requireNamespace("survminer", quietly = TRUE)){
      library(cowplot) 
      cat("N =", N, "(M =", N_groups, "ng =", n_per_group, ")\n")
      cat("DGP: Baseline =", dgp_baseline_type, ", Frailty =", dgp_frailty_type, "\n")
      cat("Separation =", cluster_separation, ", Balance =", cluster_balance, "\n")
      cat("Censoring proportion:", mean(simdat$status == 0), "\n")
      cat("Cluster sizes:", table(simdat$cluster), "\n")
      
      km_clust <- survfit(Surv(time, status) ~ cluster, data = simdat)
      line_types <- c("dashed", "dotdash", "solid")
      
      surv_plot <- ggsurvplot(
        km_clust, 
        data = simdat, 
        linetype = line_types,
        size = 1,
        risk.table = TRUE,
        palette = "Dark2"
      )
      print(surv_plot)
    }
  }
  
  return(simdat)
}
