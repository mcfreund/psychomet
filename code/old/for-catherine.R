
## function to get cross-distance matrix

pdist2 <- function(A,B) {
  ## https://www.r-bloggers.com/2013/05/pairwise-distances-in-r/
  ## https://blog.smola.org/post/969195661/in-praise-of-the-second-binomial-formula
  
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  
  m = nrow(A)
  n = nrow(B)
  
  tmp = matrix(rep(an, n), nrow=m) 
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  
  tmp - 2 * tcrossprod(A,B)  ## squared euclidean distance
  
}


## params:

n_subj <- 10
n_vert <- 100

## simulate beta contrasts:

B1 <- matrix(
  rnorm(n_subj*n_vert), 
  nrow = n_subj,
  dimnames = list(subj = paste0("subj_", letters[1:10]), vertex = NULL)
  )  ## Hi-Lo contrast coefficients from fold 1 (e.g., )
B2 <- matrix(
  rnorm(n_subj*n_vert), 
  nrow = n_subj,
  dimnames = list(subj = paste0("subj_", letters[1:10]), vertex = NULL)
  )
dimnames(B1)

## estimate cross-distance matrix, D
## diagonal gives within-subject distances, off-diagonals give between-subject distances.
## rows correspond to patterns from fold 1 (e.g., B1, or test), cols correspond to patterns from B2 (e.g., retest)

D <- pdist2(B1, B2)
D


## function to contrast between subject distances (off-diagonals) with within-subject distances (diagonal)

idi <- function(x) {

  lower <- x[lower.tri(x)]
  upper <- x[upper.tri(x)]
  diagonal <- mean(diag(x))
  
  mean(c(lower, upper)) - mean(diagonal)
  
}

idi(D)
