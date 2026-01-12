# Simulation: only age and sex; latent clusters as log-hazard shifts; compare models

# --- packages ----------------------------------------------------------------
# library(survival)
# library(survminer)

# library(survival)
# library(frailtySurv)
# library(igraph)
# library(survminer)
# library(mvtnorm)



paper_simulate_data <- function(
    seed,
    verbose = FALSE,
    frailty_type,   # "gamma" or "lognormal" or "none"
    censoring_mode  # "exponential" or "administrative" or "normal" or "none"
){
  set.seed(seed)
  
  N_clusters   <- 3
  N_groups     <- 10
  n_per_group  <- 50
  N            <- N_groups * n_per_group
  beta         <- log(2)   # scalar
  theta        <- 0.5      # frailty variance
  eta_sd       <- 1
  cluster_centers <- c(-8, 0, 8)
  
  # ----- 1. Cluster assignments for latent eta -----
  cluster_assignments <- sample(1:N_clusters, N, replace = TRUE)
  
  # ----- 2. Simulate latent predictor (eta) -----
  eta <- rnorm(N, mean = cluster_centers[cluster_assignments], sd = eta_sd)
  
  # ----- 3. Simulate group-level frailties -----
  group_ids <- rep(1:N_groups, each = n_per_group)
  
  if (frailty_type == "gamma") {
    # Gamma frailty with mean 1 and var = theta
    shape <- 1/theta
    scale <- theta
    u_groups <- rgamma(N_groups, shape = shape, scale = scale)  # mean=1
    log_u_groups <- log(u_groups)
  } else if (frailty_type == "lognormal") {
    # Log-normal frailty: log(u) ~ N(0, theta)
    log_u_groups <- rnorm(N_groups, mean = 0, sd = sqrt(theta))
    u_groups <- exp(log_u_groups)
  } else { # none
    log_u_groups <- rep(0, N_groups)
    u_groups <- rep(1, N_groups)
  }
  
  log_u <- log_u_groups[group_ids]
  u     <- u_groups[group_ids]
  
  # ----- 4. Covariate from latent predictors -----
  # We want eta = beta * x + log_u  =>  x = (eta - log_u) / beta
  x <- (eta - log_u) / beta
  
  # Single-covariate matrix (if needed)
  X <- matrix(x, ncol = 1)
  colnames(X) <- "Z1"
  
  # ----- 5. Weibull baseline via inverse cumulative hazard -----
  # Your function: Lambda_0_inv(s) = (s^(1/rho)) / c
  # So if s = -log(U) / (u * exp(beta * x)), then t = Lambda_0_inv(s)
  c_param <- 0.01
  rho     <- 2.5
  
  Lambda0_inv <- function(s, c = c_param, rho = rho) {
    (s^(1 / rho)) / c
  }
  
  # ----- 6. Simulate event times by inversion -----
  U <- runif(N)
  # linear predictor including frailty: lp = beta * x + log(u)
  lp <- beta * x + log_u
  # note: u * exp(beta*x) = exp(lp)
  s <- -log(U) / exp(lp)
  t_event <- Lambda0_inv(s, c = c_param, rho = rho)
  
  # ----- 7. Censoring (choose method) -----
  # Option A: exponential censoring with rate lambda_cens
  lambda_cens <- 0.003
  t_censor_exp <- rexp(N, rate = lambda_cens)
  
  # Option B: administrative (fixed) censoring time
  admin_censor_time <- 100  # change as desired
  
  # Option C: normal censoring, similar to genfrail()
  censor_mu <- 130     # choose mean
  censor_sd <- 15      # choose SD
  t_censor_norm <- rnorm(N, mean = censor_mu, sd = censor_sd)
  
  # Option D: uniform censoring, similar to genfrail()
  censor_unif_lower <- 80
  censor_unif_upper <- 160
  t_censor_unif <- runif(N, min = censor_unif_lower, max = censor_unif_upper)
  
  # Choose which censoring to use:
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
  
  # ----- 8. Assemble output -----
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
  
  
  
  #. ##### --- 6. Control parameters --- #####
  #. genfrail.control <- function(...) list(
  #.   crowther.subdivisions = 100L,
  #.   crowther.reltol = 1e-6,
  #.   censor.subdivisions = 100L,
  #.   censor.reltol = 1e-6
  #. )
  #. control <- genfrail.control()
  #. 
  #. ##### --- 7. Simulate survival data --- #####
  #. simdat <- genfrail(
  #.   N = N_groups,
  #.   K = n_per_group,
  #.   beta = beta,
  #.   frailty = "gamma", # "lognormal", # "none", 
  #.   theta = theta,
  #.   covar.matrix = X,
  #.   Lambda_0_inv = Lambda_0_inv,
  #.   control = control
  #. )
  #. 
  #. ##### --- 8. Add latent variables and cluster info --- #####
  #. simdat$cluster <- cluster_assignments
  #. simdat$log_u <- log_u
  #. simdat$eta <- eta
  #. simdat$X1 <- x  # covariate for later use
  #. simdat$group <- group_ids
  
  
  if(verbose){
    library(cowplot) 
    cat("N =", N, "(M =", N_groups, "ng =", n_per_group, ")\n")
    cat("Censoring proportion:", mean(simdat$status == 0), "\n")
    print(table(simdat$cluster))
    km_clust <- survfit(Surv(time, status) ~ cluster, data = simdat)
    
    line_types <- c("dashed", "dotdash", "solid")  # distinguishable in B&W
    
    # Generate the plot object
    surv_plot <- ggsurvplot(
      km_clust, 
      data = simdat, 
      linetype = line_types,
      size = 1,
      risk.table = TRUE,
      palette = "Dark2"
    )
    
    # Print to console
    print(surv_plot)
    
    # print(ggsurvplot(km_clust, data = simdat, risk.table = TRUE))
    
    combined_plot <- plot_grid(surv_plot$plot, surv_plot$table, ncol = 1, rel_heights = c(3, 1))
    ggsave(
      filename = paste0("simulations/Plots/km_clusters_cens_", censoring_mode, ".pdf"),
      plot = combined_plot,
      width = 6,
      height = 5.5
    )
    
  }
  
  
  return(simdat)
}
