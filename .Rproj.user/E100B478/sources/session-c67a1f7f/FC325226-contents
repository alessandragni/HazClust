#' Compute Spectral Embedding from a Similarity Matrix
#'
#' Computes an optimal low-dimensional representation \code{F} from a similarity matrix \code{S}
#' by applying spectral decomposition to its unnormalized graph Laplacian.
#' This is commonly used in spectral clustering.
#'
#' @param S A square, symmetric numeric similarity matrix of size \eqn{n \times n}.
#' @param c A positive integer specifying the number of smallest eigenvectors to return,
#' i.e., the number of dimensions (columns) in the resulting matrix \code{F}.
#'
#' @return A numeric matrix \code{F} of size \eqn{n \times c}, whose columns correspond
#' to the eigenvectors associated with the \code{c} smallest eigenvalues of the Laplacian matrix of \code{S}.
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Ensures symmetry of the similarity matrix \code{S}.
#'   \item Constructs the (unnormalized) graph Laplacian \eqn{L = D - S}, where \eqn{D} is the degree matrix.
#'   \item Computes the spectral decomposition of \eqn{L}.
#'   \item Returns the eigenvectors associated with the \code{c} smallest eigenvalues
#'         (excluding the zero eigenvalue if \code{S} is fully connected).
#' }
#' These eigenvectors can be used as features for clustering or dimensionality reduction.
#'
#' @examples
#' set.seed(1)
#' n <- 10
#' S <- matrix(runif(n^2), n, n)
#' S <- (S + t(S)) / 2  # symmetrize
#' diag(S) <- 0         # zero self-similarity
#' 
#' c <- 3
#' F <- HazClust:::computeF(S, c)
#' dim(F)  # Should be 10 x 3
#'
#' @keywords internal

computeF <- function(S, c) {
  # Check if S is symmetric
  S_sym <- (S + t(S)) / 2
  
  # Degree matrix
  d <- rowSums(S_sym)
  D <- diag(d)
  
  # Laplacian
  L <- D - S_sym
  
  # Eigen decomposition
  eig <- eigen(L, symmetric = TRUE)
  
  # Select the c eigenvectors corresponding to the c smallest eigenvalues
  F_opt <- eig$vectors[, (ncol(eig$vectors) - c + 1):ncol(eig$vectors)]
  
  return(F_opt)
}

