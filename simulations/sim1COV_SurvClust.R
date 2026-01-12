
# Load packages
library(survival, lib.loc="/u/ragni/Rlibs")
library(igraph, lib.loc="/u/ragni/Rlibs")
library(cluster, lib.loc="/u/ragni/Rlibs")
library(optimx, lib.loc="/u/ragni/Rlibs")

# Set working directory
setwd("/u/ragni/SurvClust")
r_files <- list.files("R", full.names = TRUE, pattern = "\\.[rR]$")
sapply(r_files, source)
source("/u/ragni/SurvClust/simulations/1COVsimulate_data.R")


# for simulation
baseline_type = "weibull"    # also "lognormal" or "weibull"
frailty_type = "gamma"   # also "lognormal" or "gamma"

dist = baseline_type 

# Parse input parameters
args <- commandArgs(trailingOnly = TRUE)
seed <- as.numeric(args[1])
gammapar <- as.numeric(args[2])
frailty <- args[3]
k <- as.numeric(args[4])
c <- as.numeric(args[5])
censor <- args[6]

if(FALSE){ # if I am working on my laptop
  library(survival)
  library(frailtySurv)
  library(igraph)
  library(survminer)
  library(roxygen2)
  library(cluster)
  
  setwd("~/Documents/DATA/POLITECNICO/PHD/CODE_REPO/SurvClust")
  source("~/Documents/DATA/POLITECNICO/PHD/CODE_REPO/SurvClust/simulations/1COVsimulate_data.R")
  
  seed = 2
  gammapar = 0.00001 #10
  # for SurvClust
  # for simulation
  baseline_type = "weibull"    # also "lognormal" or "weibull"
  frailty_type = "gamma"   # also "lognormal" or "gamma"
  dist = baseline_type 
  frailty = frailty_type # to try also "none"
  k = 20
  c = 3
  censor = "normal" # to try also "exponential" or "administrative" or "normal" or "uniform"
}



#### DATA SIMULATION ####
sim_data <- paper_simulate_data(seed = seed, verbose = FALSE, 
                                frailty_type = frailty,   # "gamma" or "lognormal" or "none"
                                censoring_mode = censor)  # "exponential" or "administrative" or "none")

#### SURVCLUST ####
cluster = "group"
data = sim_data
formula = Surv(time, status) ~ X1
strata = NULL
transform = TRUE

xx <- tryCatch({
        SurvClust(
          formula,
          cluster   = cluster,
          strata    = NULL,
          data      = data,
          inip      = NULL,
          iniFpar   = NULL,
          dist      = dist,
          frailty   = frailty,
          method    = "Nelder-Mead", #"nlminb",
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
    return(NULL)  # return NULL on error
    }
  )


# mygraph <- graph_from_adjacency_matrix(xx$S, mode = "undirected", weighted = TRUE, diag = FALSE)
# components(mygraph)$membership


if(!dir.exists("/u/ragni/SurvClust/simulations/output")) dir.create("/u/ragni/SurvClust/simulations/output")
if(!dir.exists("/u/ragni/SurvClust/simulations/output/errors")) dir.create("/u/ragni/SurvClust/simulations/output/errors")
if(!dir.exists("/u/ragni/SurvClust/simulations/output/notpenalized")) dir.create("/u/ragni/SurvClust/simulations/output/notpenalized")
if(!dir.exists("/u/ragni/SurvClust/simulations/output/penalized")) dir.create("/u/ragni/SurvClust/simulations/output/penalized")

if (!is.null(xx)) {
  
  if (gammapar == 0) {
    
    name = paste0("/u/ragni/SurvClust/simulations/output/notpenalized/sim_ExpGamma_seed", seed, 
                  "_c", c, "_k", k, "_gammapar", gammapar, 
                  "_frailty", frailty, "_censor", censor, ".Rdata")
    
    saveRDS(list(seed = seed,
                 c = c,
                 k = k,
                 gammapar = gammapar,
                 estimate = xx$estimate, 
                 # se = xx$se, 
                 # pval = xx$pval,
                 loglik = xx$loglik,
                 iterations = xx$iterations,
                 lambdaoptim = xx$lambda), name)
  } else {
    
    S = as.matrix(xx$S)
    S_sym <- (S + t(S)) / 2
    d <- rowSums(S_sym) # degree matrix
    D <- diag(d)
    L <- D - S_sym
    n_components <- ncol(L) - qr(L)$rank
    
    name = paste0("/u/ragni/SurvClust/simulations/output/penalized/sim_ExpGamma_seed", seed, 
                  "_c", c, "_k", k, "_gammapar", gammapar, 
                  "_frailty", frailty, "_censor", censor, ".Rdata")
    
    g <- graph_from_adjacency_matrix(S, mode = "undirected", weighted = TRUE, diag = FALSE)
    data$cluster_after = components(g)$membership
    
    DistforSil = sqrt(computeDistforSilhouette(formula = formula, cluster = cluster, 
                                               strata = strata, data = data, 
                                               dist = dist, frailty = frailty, 
                                               p = xx$estimate,
                                               transform = transform))
    
    # Compute silhouette
    library(cluster)
    sil <- summary(silhouette(data$cluster_after, as.dist(DistforSil)))
    #vec <- mean(unlist(sil[2]$clus.avg.widths))
    vec = unlist(sil)[4]
  
    
    saveRDS(list(seed = seed,
                 c = c,
                 k = k,
                 gammapar = gammapar,
                 estimate = xx$estimate, 
                 # se = xx$se, 
                 # pval = xx$pval,
                 n_components = n_components,
                 clusters = cbind(data$cluster, data$cluster_after),
                 S = as.matrix(xx$S),
                 loglik = xx$loglik,
                 convergence = xx$convergence,
                 iterations = xx$iterations,
                 lambdaoptim = xx$lambda,
                 silhouette = vec), name)
    }
  
}
      
