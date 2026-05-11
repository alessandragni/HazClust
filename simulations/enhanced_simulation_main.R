# Enhanced Simulation - Simple version following original structure
# Load packages
library(survival, lib.loc="/u/ragni/Rlibs")
library(igraph, lib.loc="/u/ragni/Rlibs")
library(cluster, lib.loc="/u/ragni/Rlibs")
library(optimx, lib.loc="/u/ragni/Rlibs")

# Set working directory
setwd("/u/ragni/HazClust")
r_files <- list.files("R", full.names = TRUE, pattern = "\\.[rR]$")
sapply(r_files, source)
source("/u/ragni/HazClust/simulations/enhanced_simulate_data.R")

# Parse input parameters
args <- commandArgs(trailingOnly = TRUE)
seed <- as.numeric(args[1])
gammapar <- as.numeric(args[2])
k <- as.numeric(args[3])
c <- as.numeric(args[4])
censor <- args[5]
scenario <- args[6]  # NEW
tabfig <- args[7]  # NEW
frailty_intensity <- as.numeric(args[8])  # NEW


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
    balance = "imbalanced_moderate"
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
  mismatch = FALSE #TRUE
)

#### SURVCLUST ####
cluster = "group"
data = sim_data
formula = Surv(time, status) ~ X1
strata = NULL
transform = TRUE

xx <- tryCatch({
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


# Create directories
scenario_dir <- paste0("/u/ragni/HazClust/simulations/output/", tabfig, "/", scenario)
dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)

if (!is.null(xx)) {
  
  if (gammapar == 0) {
    
    scenario_dir <- paste0("/u/ragni/HazClust/simulations/output/", tabfig, "/", scenario, "/notpenalized")
    dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
    
    name = paste0(scenario_dir, "/sim_seed", seed, 
                  "_c", c, "_k", k, "_gammapar", gammapar, 
                  "_frailty", frailty, "_censor", censor, "_frailty_intensity", frailty_intensity, ".Rds")
    
    saveRDS(list(seed = seed,
                 scenario = scenario,
                 c = c,
                 k = k,
                 gammapar = gammapar,
                 estimate = xx$estimate, 
                 loglik = xx$loglik,
                 iterations = xx$iterations,
                 lambdaoptim = xx$lambda,
                 # Scenario info
                 dgp_frailty = config$dgp_frailty,
                 dgp_baseline = config$dgp_baseline,
                 fit_frailty = config$fit_frailty,
                 fit_baseline = config$fit_baseline,
                 separation = config$separation,
                 balance = config$balance), name)
  } else {
    
    S = as.matrix(xx$S)
    S_sym <- (S + t(S)) / 2
    d <- rowSums(S_sym)
    D <- diag(d)
    L <- D - S_sym
    n_components <- ncol(L) - qr(L)$rank
    
    name = paste0(scenario_dir, "/sim_seed", seed, 
                  "_c", c, "_k", k, "_gammapar", gammapar, 
                  "_frailty", frailty, "_censor", censor, "_frailty_intensity", frailty_intensity, ".Rds")
    
    g <- graph_from_adjacency_matrix(S, mode = "undirected", weighted = TRUE, diag = FALSE)
    data$cluster_after = components(g)$membership
    
    DistforSil = sqrt(computeDistforSilhouette(formula = formula, cluster = cluster, 
                                               strata = strata, data = data, 
                                               dist = dist, frailty = frailty, 
                                               p = xx$estimate,
                                               transform = transform))
    
    # Compute silhouette
    sil <- summary(silhouette(data$cluster_after, as.dist(DistforSil)))
    vec = unlist(sil)[4]
    
    saveRDS(list(seed = seed,
                 scenario = scenario,
                 c = c,
                 k = k,
                 gammapar = gammapar,
                 estimate = xx$estimate, 
                 n_components = n_components,
                 clusters = cbind(data$cluster, data$cluster_after),
                 S = as.matrix(xx$S),
                 loglik = xx$loglik,
                 convergence = xx$convergence,
                 iterations = xx$iterations,
                 lambdaoptim = xx$lambda,
                 silhouette = vec,
                 # Scenario info
                 dgp_frailty = config$dgp_frailty,
                 dgp_baseline = config$dgp_baseline,
                 fit_frailty = config$fit_frailty,
                 fit_baseline = config$fit_baseline,
                 separation = config$separation,
                 balance = config$balance,
                 censoring_rate = mean(data$status == 0)), name)
  }
  
}