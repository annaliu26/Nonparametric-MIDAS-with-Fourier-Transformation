list.of.packages <- c("Matrix", "matrixcalc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(Matrix, matrixcalc)

## trace of matrix
tr <- function(X) {
  if (!is.numeric(X) || !is.matrix(X) || !nrow(X) == ncol(X)) stop("X must be a square numeric matrix")
  sum(diag(X))
}

## block diagonal matrix function
"adiag" <- function (..., pad = as.integer(0), do.dimnames = TRUE) 
{
  args <- list(...)
  if (length(args) == 1) {
    return(args[[1]])
  }
  if (length(args) > 2) {
    jj <- do.call("Recall", c(args[-1], list(pad = pad)))
    return(do.call("Recall", c(list(args[[1]]), list(jj), 
                               list(pad = pad))))
  }
  a <- args[[1]]
  b <- args[[2]]
  if (is.null(b)) {
    return(a)
  }
  if (is.null(dim(a)) & is.null(dim(b))) {
    dim(a) <- rep(1, 2)
    dim(b) <- rep(1, 2)
  }
  if (is.null(dim(a)) & length(a) == 1) {
    dim(a) <- rep(1, length(dim(b)))
  }
  if (is.null(dim(b)) & length(b) == 1) {
    dim(b) <- rep(1, length(dim(a)))
  }
  if (length(dim.a <- dim(a)) != length(dim.b <- dim(b))) {
    stop("a and b must have identical number of dimensions")
  }
  s <- array(pad, dim.a + dim.b)
  s <- do.call("[<-", c(list(s), lapply(dim.a, seq_len), list(a)))
  ind <- lapply(seq(dim.b), function(i) seq_len(dim.b[[i]]) + 
                  dim.a[[i]])
  out <- do.call("[<-", c(list(s), ind, list(b)))
  n.a <- dimnames(a)
  n.b <- dimnames(b)
  if (do.dimnames & !is.null(n.a) & !is.null(n.b)) {
    dimnames(out) <- mapply(c, n.a, n.b, SIMPLIFY = FALSE)
    names(dimnames(out)) <- names(n.a)
  }
  return(out)
}

## Fourier Transformation function
Fourier = function(m, K, L)
{
  jj = c(1:m)
  S0 = rep(1,m)
  S.L = matrix(rep(NA, L*m), nrow = m)
  for (i in 1:L)
  {
    S.L[,i] = (jj/m)^i
  }
  S.Ksin = S.Kcos = matrix(rep(NA, K*m), nrow = m)
  for (i in 1:K)
  {
    S.Ksin[,i] = sin(2*i*pi*jj/m)
    S.Kcos[,i] = cos(2*i*pi*jj/m)
  }
  SS = cbind(S0,S.L,S.Ksin,S.Kcos)
  return(SS)
}

## BIC function
BIC1.lambda = function(d1, d4, BB, yy, NN, beta, K.hat)
{
  df = (sum(d1)+sum(d4))/NN*K.hat
  bic1 = log(t(yy-BB%*%beta)%*%(yy-BB%*%beta)/NN)+log(NN)*(K.hat+v.len)/NN
  return(bic1)
}

## clustering functions
knots_eq3 <- function(x, k, m){
  #external knots are on boundary
  #return boundary with internal knots only
  #used in bs or bsplineS
  c(min(x), seq(from=min(x), to=max(x), length.out=m+2)[-c(1,m+2)], max(x))
}

create_adjacency <- function(V,n) {
  differences <- apply(V,2,FUN=function(x) {norm(as.matrix(x),'f')})
  connected_ix <- which(differences == 0);
  index = t(combn(n,2));
  i <- index[connected_ix,1]
  j <- index[connected_ix,2]
  A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  A[(j-1)*n + i] <- 1
  return(A)
}


