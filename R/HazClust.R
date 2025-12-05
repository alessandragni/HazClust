#' HazClust algorithm
#'
#' @param formula A formula with a survival object created by \code{Surv()} on the left and 
#' covariates on the right.
#' @param cluster Optional. The name (character string) of the variable in \code{data} indicating cluster IDs.
#' @param strata Optional. The name (character string) of the variable in \code{data} indicating strata.
#' @param data A data frame containing the variables in the model.
#' @param inip Optional. Initial values for the baseline hazard and regression parameters.
#' @param iniFpar Optional. Initial value(s) for the frailty distribution parameter(s).
#' @param dist A character string specifying the baseline hazard distribution. 
#' One of: \code{"weibull"}, \code{"inweibull"}, \code{"frechet"}, 
#' \code{"exponential"}, \code{"gompertz"}, \code{"loglogistic"}, 
#' \code{"lognormal"}, \code{"logskewnormal"}.
#' @param frailty A character string specifying the frailty distribution. 
#' One of: \code{"none"}, \code{"gamma"}, \code{"ingau"}, 
#' \code{"possta"}, \code{"lognormal"}, \code{"loglogistic"}.
#' @param method Optimization method passed to \code{optimx}.
#' @param maxitparfm Maximum number of iterations for internal optimization.
#' @param maxit Maximum number of iterations for optimization.
#' @param tolS Convergence tolerance for S.
#' @param tolll Convergence tolerance for log-likelihood.
#' @param Fparscale Scaling value for the frailty parameter in optimisation. 
#' The algorithm optimizes \code{Fpar / Fparscale}.
#' @param showtime Logical. Should the execution time be printed?
#' @param correct Used only with \code{frailty = "possta"}. Correction factor for likelihood 
#' calculation in cases of many events per cluster.
#' @param c Integer. Number of clusters.
#' @param k Integer. Number of nearest neighbors used in the similarity matrix.
#' @param gammapar Parameter for l1 penalization on the similarity matrix.
#' @param lambda0 Optional. Initial value for the lambda penalty term.
#' @param S0 Optional. Initial similarity matrix. Defaults to uniform.
#' @param trace Logical. Should progress and convergence information be printed?
#' @param transform Logical. Whether to apply transformation during distance computation.
#'
#' @return A list containing:
#' \item{estimate}{Named vector of parameter estimates (frailty, baseline hazard, regression).}
#' \item{F}{Low-rank embedding matrix computed from Laplacian eigenvectors.}
#' \item{S}{Estimated similarity matrix.}
#' \item{loglik}{Log-likelihood value at the optimum.}
#' \item{iterations}{Number of optimization iterations.}
#' \item{convergence}{Logical. Whether the algorithm converged.}
#' \item{lambda}{Final value of the penalty parameter lambda.}
#' \item{time}{Execution time, if \code{showtime = TRUE}.}
#'
#' @details The function uses maximum likelihood estimation to fit a parametric frailty model.
#' If \code{frailty = "none"}, then the model is fitted without random effects.
#' Stratification allows different baseline hazards across strata.
#'
#' @examples
#' \dontrun{
#' data(kidney, package = "survival")
#' parfm(Surv(time, status) ~ age + sex, cluster = "id", data = kidney,
#'       dist = "weibull", frailty = "gamma")
#' }
#'
#' @importFrom survival Surv
#' @importFrom stats model.matrix optimHess
#' @importFrom utils data
#' @importFrom optimx optimx
#' @importFrom igraph graph_from_adjacency_matrix components
#' @export
HazClust <- function(formula,
                      cluster   = NULL,
                      strata    = NULL,
                      data,
                      inip      = NULL,
                      iniFpar   = NULL,
                      dist      =  c("weibull", "inweibull", "frechet", "exponential", 
                                   "gompertz", "loglogistic", "lognormal",
                                   "logskewnormal"),
                      frailty   = c("none", "gamma", "ingau", "possta",
                                  "lognormal", "loglogistic"),
                      method    = "nlminb",
                      maxitparfm = 500,
                      maxit     = 500,
                      tolS      = 1e-4,
                      tolll     = 1e-3,
                      Fparscale = 1,
                      showtime  = FALSE,
                      correct   = 0,
                      transform = TRUE,
                      c, # clusters
                      k, # knn
                      gammapar, # parameter for l1 penalization on similarity
                      lambda0 = NULL, # initial value for lambda
                      S0 = NULL, # initial value for S
                      trace = TRUE
                      ){
  
  
  if (missing(data)) {
    data <- eval(parse(text = paste("data.frame(", 
                                    paste(all.vars(formula), collapse = ", "),
                                    ")")))
  }
  
  #----- Check the baseline hazard and the frailty distribution ---------------#
  dist <- tolower(dist)
  frailty <- tolower(frailty)
  if (frailty == "none" &&  !is.null(cluster)) {
    warning(paste0("With frailty='none' the cluster variable '",
                   cluster, "' is not used!"))
  }
  if (frailty == "none" &&  !is.null(iniFpar)) {
    warning("With frailty='none' the argument 'iniFpar' is not used!")
  }
  
  correct   = 0
  #----- 'Correct' is useless except for frailty="possta" -------------------#
  if (frailty == "possta") {  #Do not exaggerate when setting 'correct' !
    if (10 ^ correct == Inf || 10 ^ -correct == 0) {
      stop("'correct' is too large!")
    }
    if (10 ^ correct == 0 || 10 ^ -correct == Inf) {
      stop("'correct' is too small!")
    }
  } else if (correct != 0) {
    warning(paste0("'correct' has no effect when 'frailty = ", frailty, "'"))
  }
  
  
  # Initialization
  iter <- 0
  max_iter <- maxit
  converged <- FALSE
  prev_loglik <- -Inf
  start_time <- Sys.time()
  
  if (!is.null(S0)) {
     S = S0
  } else {
     # S = matrix(1/nrow(data), nrow(data), nrow(data))
     S = matrix(1/(nrow(data)-1), nrow(data), nrow(data))
     diag(S) = 0
  }
  
  if (!is.null(lambda0)) {
      lambda = lambda0
    } else {
      lambda = 1000
  }
  
  if (trace) cat("Starting HazClust optimization...\n")
  
  if(gammapar == 0) {
    x = MYparfm(formula = formula, cluster = cluster, strata = strata, data = data, 
                inip = inip, iniFpar = iniFpar,
                dist = dist, frailty = frailty,
                method = method, maxit = maxitparfm,
                Fparscale = Fparscale,
                showtime  = showtime,
                correct = correct,
                gammapar = gammapar, S = S)
    p <- x[, "ESTIMATE"]
    # se <- x[, "SE"]
    if (trace){
      cat("\nEstimated parameters:\n")
      cat("  p     : ", paste(sprintf("%.4f", p), collapse = ", "), "\n", sep = "")
      # cat("  se    : ", paste(sprintf("%.4f", se), collapse = ", "), "\n", sep = "")
    }
    # add pvalue if available
    if ("p-val" %in% colnames(x)) {
      pval <- x[, "p-val"]
      if (trace) cat("  p-val : ", paste(ifelse(is.na(pval), "NA", 
                                                sprintf("%.4f", pval)), collapse = ", "), "\n", sep = "")
    }
    loglik <- attributes(x)$loglik
    if (trace) cat(sprintf("Optimized MYparfm. Log-likelihood: %.4f\n", loglik))
    
    Fmatr = NULL
    S = NULL
    iterations = 1
    convergence = TRUE
  }
  else {
  
    repeat {
      iter <- iter + 1
      if (trace) cat(sprintf("\n--- Iteration %d ---\n", iter))
      
      # Step 1: Fix S -> Fix S (and w, ψ, θ implicitly), Optimize F
      # here i need S and c
      Fmatr = computeF(S, c) 
        if (trace) cat("Updated F (from Laplacian).\n")
      
      # Step 2a: Fix F , S, Optimize w, ψ, θ
      # here i need eta, gammapar and S
      if (exists("S") && any(!is.finite(S))) stop("S contains non-finite values before MYparfm.")
      if (exists("Fmatr") && any(!is.finite(Fmatr))) stop("Fmatr contains non-finite values before MYparfm.")
      if (!is.finite(lambda)) stop("lambda is non-finite before MYparfm.")
      x = MYparfm(formula = formula, cluster = cluster, strata = strata, data = data, 
                  inip = inip, iniFpar = iniFpar,
                  dist = dist, frailty = frailty,
                  method = method, maxit = maxitparfm,
                  Fparscale = Fparscale,
                  showtime  = showtime,
                  correct = correct,
                  gammapar = gammapar, S = S)
      p <- x[, "ESTIMATE"]
      # se <- x[, "SE"]
      if (trace){
        cat("\nEstimated parameters:\n")
        cat("  p     : ", paste(sprintf("%.4f", p), collapse = ", "), "\n", sep = "")
        ##. cat("  se    : ", paste(sprintf("%.4f", se), collapse = ", "), "\n", sep = "")
      }
      # add pvalue if available
      if ("p-val" %in% colnames(x)) {
        pval <- x[, "p-val"]
        if (trace) cat("  p-val : ", paste(ifelse(is.na(pval), "NA", 
                                                  sprintf("%.4f", pval)), collapse = ", "), "\n", sep = "")
      }
      loglik <- attributes(x)$loglik
      if (trace) cat(sprintf("Optimized MYparfm. Log-likelihood: %.4f\n", loglik))
      
      # Step 2b: Fix F , w, ψ, θ, Optimize S
      # here i need F, w, ψ, θ, k
      D = computeD(formula = formula, cluster = cluster, strata = strata, data = data, 
                   dist = dist, frailty = frailty, 
                   p = p,
                   Fmatr = Fmatr, lambda = lambda,
                   transform = transform, correct = correct)
      D = as.matrix(D)
      
      # CHANGED HERE
      # res = safe_computeS(D = D, k = k, S_prev = S, lambda = lambda, trace = trace)
      # res = computeS(D, k)
      res <- safe_computeS(D, k)
      if (is.null(res) || !is.list(res) || !("S" %in% names(res))) {
        stop("safe_computeS returned invalid result.")
      }
      S_new <- res$S
      if (any(!is.finite(S_new))) {
        stop("safe_computeS produced non-finite entries in S_new.")
      }
      diff_S <- sum(abs(S_new - S)) / length(S)
      S <- S_new
      if (trace) cat(sprintf("Updated S. Mean abs change: %.6f\n", diff_S))
      
      # Convergence check
      delta_ll <- abs(loglik - prev_loglik)
      if (trace) cat(sprintf("Loglik improvement: %.6f\n", delta_ll))
      
      G <- graph_from_adjacency_matrix(S, mode = "undirected", weighted = TRUE)
      membership <- components(G)$membership
      n_components <- components(G)$no           # number of connected components
      if (trace) cat(sprintf("Updated number of connected components: %.6f\n", n_components))
      
      
      if ((diff_S < tolS && delta_ll < tolll && n_components == c) || iter >= max_iter) {
        converged <- TRUE
        if (trace) {
          if (iter >= max_iter) {
            cat("Reached max iterations.\n")
          } else {
            cat("Converged.\n")
          }
        }
        break
      }
      
      
      if(iter == 3) {
        lambda = res$mu_global
      } else {
        
        
        S_sym <- (S + t(S)) / 2
        d <- rowSums(S_sym) # degree matrix
        D <- diag(d)
        L <- D - S_sym # Laplacian
        # Compute the nullity (dimension of the kernel)
        n_components <- ncol(L) - qr(L)$rank
        
        # if seen through a graph
        # G <- graph_from_adjacency_matrix(S, mode = "undirected", weighted = TRUE)
        # membership <- components(G)$membership
        # n_components <- components(G)$no           # number of connected components
        
        if(n_components < c) {
          lambda = lambda * 1.2 # see https://ojs.aaai.org/index.php/AAAI/article/view/10909
        } else {
          lambda = lambda / 1.2
        }
      } 
      
      prev_loglik <- loglik
    }
  } # end else gammapar != 0  
  
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "secs")

  
  return(list(
    estimate    = p,
    # se          = se,
    pval       = if (exists("pval")) pval else NULL,
    MYparfmx    = x,
    Fmatr       = Fmatr,
    S           = S,
    loglik      = loglik,
    iterations  = iter,
    convergence = converged,
    lambda      = lambda,
    time        = if (showtime) total_time else NULL
  ))
  

} 






  



  
  