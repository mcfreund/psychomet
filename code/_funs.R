enlist <- function(x) setNames(vector("list", length(x)), x)


list2array <- function(l, name_new_dim = NULL) {
	## input: a list of n-dimensional arrays, of identical dimension and dimnames.
	## binds arrays along new dimension.
	## output: a new, (n+1)-dimensional array.
	## 
	## works well for wrangling results of foreach::foreach() call, when n-D array is returned each iteration.
	
	
	## TODO: input validation --- check list, and all dims and dimnames equal
	
	## get info existing dims
	
	dim_existing <- dim(l[[1]]) ## vector with values equal to length of each dim
	dimnames_existing <- dimnames(l[[1]])  ## list of dim names
	dim_new <- length(l)  ## length of new dim
	dimnames_new <- setNames(list(names(l)), name_new_dim)
	
	## make new array:
	
	a <- array(
		NA,
		dim = c(dim_existing, dim_new),
		dimnames = c(dimnames_existing, dimnames_new)
	)
	
	## iteratively add elements of new dim to array:
	
	index_existing <- as.matrix(expand.grid(lapply(dim_existing, seq_len)))
	for (i in seq_len(dim_new)) a[cbind(index_existing, i)] <- l[[i]]
	
	a
	
}



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


