#' Compute Similarity Matrix S from Distance Matrix D
#'
#' Constructs a symmetric similarity matrix \code{S} from a distance matrix \code{D}
#' using a local scaling method based on the k-nearest neighbors.
#'
#' @param D A square numeric distance matrix.
#' @param k An integer. The number of nearest neighbors to consider for each point, excluding itself.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{S}{A symmetric \eqn{N \times N} similarity matrix with local sparsity structure,
#'            computed based on the top-k neighbors (excluding self-distances) for each point.}
#'   \item{mu_global}{The global average of the local \eqn{\mu_r} scaling parameters across rows.}
#' }
#'
#' @details
#' For each row \eqn{r}, the function excludes the diagonal entry \eqn{D_{rr}}, then identifies
#' the k nearest neighbors from the remaining distances. It computes a local scaling factor
#' \eqn{\mu_r} and offset \eqn{\alpha_r}, and defines a local similarity score using:
#' \deqn{S_{rj} = \max\left(0, -D_{rj} / (2 \mu_r) + \alpha_r\right)}
#' Only the k-nearest entries in each row are initially retained. The similarity matrix \code{S}
#' is finally symmetrized by averaging: \code{S <- (S + t(S)) / 2}.
#' The diagonal of \code{S} remains zero, as self-similarity is not computed.
#'
#' @examples
#' D <- matrix(runif(100), nrow = 10)
#' D <- (D + t(D)) / 2  # make symmetric
#' diag(D) <- 0
#' res <- HazClust:::computeS(D, k = 3)
#' str(res$S)
#'
#' @keywords internal
#' 
computeS <- function(D, k) {
  if (!is.matrix(D) && !is.data.frame(D)) stop("Input 'D' must be a matrix or data frame.")
  if (!is.numeric(k) || length(k) != 1 || k <= 0 || k != round(k)) stop("Input 'k' must be a single positive integer.")
  if (k + 1 > ncol(D)) stop(sprintf("k+1 (=%d) is greater than number of columns in D (=%d).", k+1, ncol(D)))
  
  D <- as.matrix(D)
  N <- nrow(D)
  S <- matrix(0, nrow = N, ncol = N)
  
  mu <- numeric(N)
  alpha <- numeric(N)
  
  eps <- 1e-12  # safeguard for division by zero or negative mu
  
  for (r in 1:N) {
    # Exclude self-distance
    indices <- setdiff(1:N, r)
    dists <- D[r, indices]
    
    # Sort neighbors by distance
    ord <- order(dists)
    sorted_indices <- indices[ord]
    sorted_dists   <- dists[ord]
    
    # Top-k neighbors and (k+1)-th distance
    top_k <- sorted_dists[1:k]
    d_kplus1 <- sorted_dists[k + 1]
    
    # Compute mu and alpha (with numerical safeguard)
    mu[r] <- max((k/2) * d_kplus1 - (1/2) * sum(top_k), eps)
    alpha[r] <- (1/k) + (1/(2 * k * mu[r])) * sum(top_k)
    
    # Compute similarity weights for top-k neighbors
    s_values <- pmax(0, alpha[r] - top_k / (2 * mu[r]))
    
    # Assign to similarity matrix
    S[r, sorted_indices[1:k]] <- s_values
  }
  
  # Global mu
  mu_global <- mean(mu)
  
  return(list(S = S, mu_global = mu_global, mu = mu, alpha = alpha))
}





# -------------- robust computeS wrapper (drop-in) -----------------
# Use instead of direct: res = computeS(D, k); S_new = res$S

# -------------- robust computeS wrapper (drop-in) -----------------
# Use instead of direct: res = computeS(D, k); S_new = res$S

#.  safe_computeS <- function(D, k, S_prev = NULL, lambda = NULL,
#.                            max_attempts = 5, densify_eps = 1e-4,
#.                            cap_S = 1.0, min_row_sum = 1e-8, trace = TRUE) {
#.    # Pre-check D
#.    if (!is.matrix(D)) D <- as.matrix(D)
#.    if (any(is.na(D) | is.nan(D))) {
#.      if (trace) cat("Warning: D contains NA/NaN -> replacing with large finite values\n")
#.      # Replace NA with large finite value (so similarity becomes small)
#.      finite_max <- max(D[is.finite(D)], na.rm = TRUE)
#.      if (!is.finite(finite_max)) finite_max <- 1e6
#.      D[is.na(D) | is.nan(D)] <- finite_max * 10
#.    }
#.    if (any(is.infinite(D))) {
#.      if (trace) cat("Warning: D contains Inf -> replacing with large finite values\n")
#.      finite_max <- max(D[is.finite(D)], na.rm = TRUE)
#.      if (!is.finite(finite_max)) finite_max <- 1e6
#.      D[is.infinite(D)] <- finite_max * 10
#.    }
#.    # clip extreme values in D to a reasonable range to avoid under/overflow
#.    D[is.na(D)] <- max(1e6, max(D[is.finite(D)], na.rm = TRUE))
#.    D <- pmax(D, 0)  # ensure non-negative distances
#.    
#.    attempt <- 1
#.    res <- NULL
#.    while (attempt <= max_attempts) {
#.      # call your computeS safely
#.      res_try <- tryCatch({
#.        computeS(D, k)
#.      }, error = function(e) {
#.        if (trace) cat(sprintf("computeS error (attempt %d): %s\n", attempt, e$message))
#.        return(structure(NULL, error = TRUE, message = e$message))
#.      })
#.      
#.      if (!is.null(attr(res_try, "error")) && attr(res_try, "error") == TRUE) {
#.        # fallback: densify D slightly and retry (reduces identical rows issues)
#.        if (trace) cat(sprintf("computeS failed; densifying D and retrying (attempt %d)\n", attempt))
#.        D <- D + matrix(runif(length(D), 0, densify_eps), nrow(D), ncol(D))
#.        attempt <- attempt + 1
#.        next
#.      }
#.      if (is.null(res_try) || is.null(res_try$S)) {
#.        if (trace) cat(sprintf("computeS returned NULL/invalid (attempt %d). Densifying and retry.\n", attempt))
#.        D <- D + matrix(runif(length(D), 0, densify_eps), nrow(D), ncol(D))
#.        attempt <- attempt + 1
#.        next
#.      }
#.      # Got something: validate S
#.      S_try <- res_try$S
#.      if (!is.matrix(S_try)) S_try <- as.matrix(S_try)
#.      # sanitize numeric issues
#.      S_try[!is.finite(S_try)] <- 0
#.      S_try[S_try < 0] <- 0
#.      # cap extremely large similarities
#.      S_try <- pmin(S_try, cap_S)
#.      # zero diagonal and symmetry
#.      diag(S_try) <- 0
#.      S_try <- (S_try + t(S_try)) / 2
#.      
#.      # avoid empty rows: if a row sums to near-zero, mix with previous S or uniform
#.      row_sums <- rowSums(S_try)
#.      idx_zero <- which(row_sums < min_row_sum)
#.      if (length(idx_zero) > 0) {
#.        if (trace) cat(sprintf("computeS produced %d near-zero rows; repairing.\n", length(idx_zero)))
#.        if (!is.null(S_prev)) {
#.          # mix those rows with S_prev rows to keep some structure
#.          S_try[idx_zero, ] <- 0.5 * S_try[idx_zero, ] + 0.5 * S_prev[idx_zero, ]
#.          S_try[, idx_zero] <- (S_try[, idx_zero] + t(S_try[idx_zero, ])) / 2
#.        } else {
#.          # fallback: connect those rows weakly to all others
#.          S_try[idx_zero, ] <- densify_eps
#.          S_try[, idx_zero] <- densify_eps
#.        }
#.        diag(S_try) <- 0
#.        S_try <- (S_try + t(S_try)) / 2
#.      }
#.      
#.      # final normalization (optional): scale each row so that average off-diagonal weight
#.      # matches initialization scale 1/(n-1)
#.      n <- nrow(S_try)
#.      target_sum <- 1 / (n - 1)
#.      row_sums <- rowSums(S_try)
#.      # avoid divide-by-zero
#.      row_sums[row_sums == 0] <- 1
#.      S_norm <- S_try / row_sums * target_sum
#.      
#.      # Re-symmetrize and ensure zero diag
#.      S_norm <- (S_norm + t(S_norm)) / 2
#.      diag(S_norm) <- 0
#.      
#.      # Good enough -> return
#.      res_try$S <- S_norm
#.      return(res_try)
#.    } # end while attempts
#.    
#.    # If we exit loop, return safest fallback: previous S or uniform S
#.    if (!is.null(S_prev)) {
#.      if (trace) cat("computeS failed after retries: returning previous S (fallback)\n")
#.      return(list(S = S_prev, mu_global = if (!is.null(S_prev)) mean(S_prev) else 0))
#.    } else {
#.      if (trace) cat("computeS failed after retries: returning uniform S (fallback)\n")
#.      n <- nrow(D)
#.      S_uniform <- matrix(1/(n-1), n, n)
#.      diag(S_uniform) <- 0
#.      return(list(S = S_uniform, mu_global = 0))
#.    }
#.  }
# -------------- end safe_computeS -----------------


safe_computeS <- function(D, k) {
  if (!is.matrix(D)) stop("D must be a matrix")
  D <- as.matrix(D)
  N <- nrow(D)
  S <- matrix(0, N, N)
  mu <- numeric(N)
  alpha <- numeric(N)
  eps <- 1e-6   # use bigger safeguard
  
  # Replace bad distances
  D[is.na(D) | is.nan(D) | is.infinite(D)] <- max(D[is.finite(D)], na.rm = TRUE)
  D[D < 0] <- 0
  
  for (r in seq_len(N)) {
    indices <- setdiff(seq_len(N), r)
    dists <- D[r, indices]
    ord <- order(dists)
    sorted_indices <- indices[ord]
    sorted_dists <- dists[ord]
    if (length(sorted_dists) < k + 1) next
    
    top_k <- sorted_dists[1:k]
    d_kplus1 <- sorted_dists[k + 1]
    
    mu_r <- (k / 2) * d_kplus1 - 0.5 * sum(top_k)
    mu[r] <- max(mu_r, eps)
    alpha[r] <- (1 / k) + (sum(top_k) / (2 * k * mu[r]))
    
    s_values <- pmax(0, alpha[r] - top_k / (2 * mu[r]))
    s_values[!is.finite(s_values)] <- 0
    S[r, sorted_indices[1:k]] <- s_values
  }
  
  # Symmetrize and normalize
  S <- (S + t(S)) / 2
  diag(S) <- 0
  S <- S / (max(S) + eps)
  
  mu_global <- mean(mu[is.finite(mu)])
  list(S = S, mu_global = mu_global, mu = mu, alpha = alpha)
}