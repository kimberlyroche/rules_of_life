library(Rcpp)
library(LaplacesDemon)

sourceCpp("Riemann_dist.cpp")

RiemannR <- function(A, B) {
  invA <- solve(A)
  cholInvA <- chol(invA) # upper triangular: U^T U
  ABA <- cholInvA%*%B%*%t(cholInvA)
  eABA <- eigen(ABA)
  D <- diag(length(eABA$values))
  diag(D) <- eABA$values
  diag(D) <- log(diag(D))
  diag(D) <- diag(D)^2
  return(sqrt(sum(diag(D))))
}

RiemannR2 <- function(A, B) {
  U <- chol(A) # upper triangular: U^T U
  invU <- solve(U)
  ABA <- t(invU)%*%B%*%invU
  eABA <- eigen(ABA)
  D <- diag(length(eABA$values))
  diag(D) <- eABA$values
  diag(D) <- log(diag(D))
  diag(D) <- diag(D)^2
  return(sqrt(sum(diag(D))))
}

# these should all be equivalent
for(i in 1:1) {
  P <- 20
  A <- rinvwishart(P+10, diag(P))
  B <- rinvwishart(P+10, diag(P))
  cat("A x B:",RiemannR(A, B),"\n")
  cat("B x A:",RiemannR(B, A),"\n")
  cat("A x B:",RiemannR2(A, B),"\n")
  cat("B x A:",RiemannR2(B, A),"\n")
  cat("A x B:",mat_dist(A, B, TRUE),"\n")
  cat("B x A:",mat_dist(B, A, TRUE),"\n\n")

  invA <- solve(A)
  invB <- solve(B)
  cat("A x B:",RiemannR(invA, invB),"\n")
  cat("B x A:",RiemannR(invB, invA),"\n")
  cat("A x B:",RiemannR2(invA, invB),"\n")
  cat("B x A:",RiemannR2(invB, invA),"\n")
  cat("A x B:",mat_dist(invA, invB, TRUE),"\n")
  cat("B x A:",mat_dist(invB, invA, TRUE),"\n\n")

  X <- matrix(rnorm(P*P), P, P)
  XAXT <- X%*%A%*%t(X)
  XBXT <- X%*%B%*%t(X)
  cat("A x B:",RiemannR(XAXT, XBXT),"\n")
  cat("B x A:",RiemannR(XBXT, XAXT),"\n")
  cat("A x B:",RiemannR2(XAXT, XBXT),"\n")
  cat("B x A:",RiemannR2(XBXT, XAXT),"\n")
  cat("A x B:",mat_dist(XAXT, XBXT, TRUE),"\n")
  cat("B x A:",mat_dist(XBXT, XAXT, TRUE),"\n\n")
}
