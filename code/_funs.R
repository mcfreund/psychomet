enlist <- function(x) setNames(vector("list", length(x)), x)

pdist2 <- function(A,B) {
  
  ## this function works on matrices A and B.
  ## A and B should be matrices with subjects as rows and vertices as columns.
  ## A and B should be from separate folds (i.e., "test" and "retest" or "run1" and "run2").
  
  ## this function computes the squared euclidean distances between each row of A and each row of B.
  ## the output is therefore a matrix of size nrow(A)*nrow(B), i.e., Sects * Sects.
  
  ## this is an efficient implementation of the pdist::pdist function.
  ## see links for more information:
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




idi <- function(x) {
  
  ## this function works on square distance matrices, and computes the "intersubject discrimination index".
  ## it subtracts the mean diagonal from the mean of the off diagonals.
  
  offdiag <- x[row(x) != col(x)]
  diagonal <- diag(x)
  
  mean(offdiag) - mean(diagonal)
  
}


