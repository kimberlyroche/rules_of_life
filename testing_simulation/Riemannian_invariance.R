# this file generates fake data in order to demonstrate the Riemannian distance is invariant
# to transformation of CLR <-> ALR <-> ILR
#
# just writes results to STDOUT

library(LaplacesDemon)
library(driver)

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

# analogous implementation
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

get.alr <- function(seed_mat) {
  D <- 10
  upsilon <- D+10
  GG <- t(create_alr_base(D=D, d=D))
  Xi <- GG%*%(seed_mat)%*%t(GG)
  Xi <- Xi*(upsilon-D-1)
}

get.ilr <- function(seed_mat) {
  D <- 10
  upsilon <- D+10
  GG <- t(create_default_ilr_base(D))
  Xi <- GG%*%(seed_mat)%*%t(GG)
  Xi <- Xi*(upsilon-D-1)
}

get.clr <- function(seed_mat) {
  D <- 10
  upsilon <- D+10
  GG <- t(create_clr_base(D))
  Xi <- GG%*%(seed_mat)%*%t(GG)
  Xi <- Xi*(upsilon-D-1)
}

seed_mat_A <- rinvwishart(upsilon, diag(D))
seed_mat_B <- rinvwishart(upsilon, diag(D))

A.alr <- get.alr(seed_mat_A)
B.alr <- get.alr(seed_mat_B)

RiemannR(A.alr,B.alr)
RiemannR(B.alr,A.alr)

A.ilr <- get.ilr(seed_mat_A)
B.ilr <- get.ilr(seed_mat_B)

RiemannR(A.ilr,B.ilr)
RiemannR(B.ilr,A.ilr)

A.clr <- get.clr(seed_mat_A)
B.clr <- get.clr(seed_mat_B)

RiemannR(A.clr,B.clr)
RiemannR(B.clr,A.clr)

RiemannR2(A.clr,B.clr)
RiemannR2(B.clr,A.clr)

# note truncation of CLR(Sigma) to give same distance went as:
#   D <- dim(Sigma)[1]
#   N <- dim(Sigma)[3]
#   truncSigma <- array(NA, dim=c(D-1, D-1, N))
#   for(i in 1:N) {
#     truncSigma[,,i] <- Sigma[1:(D-1), 1:(D-1), i]
#   }
#   Sigma <- truncSigma
