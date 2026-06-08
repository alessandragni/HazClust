
# Load packages
library(survival)
library(igraph)
library(cluster)
library(optimx)


# Parse input parameters
gammapar <- 0.01
frailty <- "ingau"
k <- 150
c <-3
base <- "weibull"


sim_data = readRDS("casestudy/data/sim_data.Rdata")

sim_data$RESP = as.numeric(sim_data$BRONCHITE) + 
  (as.numeric(sim_data$INSUF_RESPIRATORIA)-1) + 
  (as.numeric(sim_data$COPD)-1) + 
  (as.numeric(sim_data$POLMONITE)-1)

sim_data$RESPbin = as.factor(ifelse(sim_data$RESP>0, 1, 0))


#### HAZCLUST ####
cluster = "COD_OSPEDALE"
data = sim_data
formula = Surv(time, cens) ~ ETA_AL_RICOVERO + SESSO + MCS + RESP
strata = NULL
transform = TRUE

# parfm::parfm(formula = formula, cluster = cluster, data = sim_data, dist = base, frailty = frailty)

xx <- tryCatch({
  HazClust(
    formula,
    cluster   = cluster,
    strata    = NULL,
    data      = data,
    inip      = NULL,
    iniFpar   = NULL,
    dist      = base,
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
  message(paste("Error at base =", base, "frailty =", frailty, "c =", c, "k =", k, "gammapar =", gammapar, ":", e$message))
  return(NULL)  # return NULL on error
}
)




    S = as.matrix(xx$S)
    S_sym <- (S + t(S)) / 2
    d <- rowSums(S_sym) # degree matrix
    D <- diag(d)
    L <- D - S_sym
    n_components <- ncol(L) - qr(L)$rank
    
    name = paste0("casestudy/result/output_gammapar", gammapar, 
                  "_frailty", frailty, "_base", base, "_k", k, "_c", c, "SEs.Rdata")
    
    g <- graph_from_adjacency_matrix(S, mode = "undirected", weighted = TRUE, diag = FALSE)
    data$cluster_after = components(g)$membership
    
    DistforSil = sqrt(computeDistforSilhouette(formula = formula, cluster = cluster, 
                                               strata = strata, data = data, 
                                               dist = base, frailty = frailty, 
                                               p = xx$estimate,
                                               transform = transform))
    
    # Compute silhouette
    library(cluster)
    sil <- summary(silhouette(data$cluster_after, as.dist(DistforSil)))
    vec = sil$avg.width
    
    
    saveRDS(list(base = base,
                 frailty = frailty,
                 c = c,
                 k = k,
                 gammapar = gammapar,
                 estimate = xx$estimate, 
                 # se = xx$se, 
                 # pval = xx$pval,
                 n_components = n_components,
                 clusters = data$cluster_after,
                 S = as.matrix(xx$S),
                 loglik = xx$loglik,
                 convergence = xx$convergence,
                 iterations = xx$iterations,
                 lambdaoptim = xx$lambda,
                 silhouette = vec), name)
  
