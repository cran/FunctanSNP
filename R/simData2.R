##' @title Generate simulated data for \code{\link{SNPinter}} and \code{\link{SNPcvinter}}
##'
##' @description Generate simulated data for users to apply the method SNPinter and SNPcvinter, including response variable y, scalar variable Z and sequence (genotypes) data X.
##'
##' @param n an interger variable specifying the number of samples to be generated.
##' @param m an interger variable specifying the sequence length of each sample.
##' @param seed an integer variable specifyinging the random seed used for random sequence generation.
##'
##' @return An "simData2" object that contains the list of the following items.
##' \itemize{
##' \item{y:}{ a numeric vector representing the response variables.}
##' \item{z:}{ a matrix representing the scalar covariates, with the number of rows equal to the number of samples.}
##' \item{location:}{ a numeric vector defining the sampling sites of the sequence (genotypes) data.}
##' \item{X:}{ a matrix representing the sequence data, with the number of rows equal to the number of samples.}
##' \item{beta:}{ an "fd" object specifying the genetic effect function.}
##' }
##' @seealso See Also as \code{\link{simData1}}, \code{\link{SNPinter}}, \code{\link{SNPcvinter}}.
##'
##' @import fda
##' @import MASS
##' @importFrom stats lm
##' @importFrom stats rnorm
##' @export
##' @examples
##' library(FunctanSNP)
##' n <- 300
##' m <- 30
##' simdata2 <- simData2(n, m, seed = 123)
##'

#library(MASS)
#library(ggplot2)
#library(fda)
simData2 <- function(n, m, seed){
  ### Generate discrete values of Xi(t)
  norder <- 4
  nknots <- 15
  t <- seq(1e-2, 1, length = m)
  breaks <- seq(0, 1, length = nknots)
  basismat <- bsplineS(x = t, breaks = breaks, norder = norder, returnMatrix = FALSE)
  fbasisX <- create.bspline.basis(rangeval = c(0, 1), breaks = breaks, norder = norder)
  nbasisX <- ncol(basismat)
  set.seed(seed = seed + 123); coef <- mvrnorm(n = n, mu = rep(0, nbasisX), Sigma = diag(x = 1, nrow = nbasisX, ncol = nbasisX))
  Rawfvalue <- coef %*% t(basismat)
  fvalue <- Rawfvalue
  funcX <- function(l){
    x <- fvalue[l, ]
    len <- length(x)
    diffmat <- matrix((rep(x, times = 3) - rep(c(0,1,2), each = len))^2, ncol = len, nrow = 3, byrow = TRUE)
    value <- apply(diffmat, 2, which.min) - 1
    return(value)
  }
  dataX <- t(mapply(funcX, 1:n))

  gamma <- c(0.4, 0.8)
  ### Generate Zik
  set.seed(seed = seed + 1234); z <- mvrnorm(n = n, mu = rep(0, 2), Sigma = diag(x = 1, nrow = 2, ncol = 2))
  ### Generate epsiloni
  set.seed(seed = seed + 12345); epsilon <- rnorm(n = n, mean = 0, sd = 0.1)
  ### Generate beta(t)
  region1 <- t[which(t <= 0.3)]
  region2 <- t[which(t > 0.3 & t <= 0.7)]
  region3 <- t[which(t > 0.7)]
  Betapart1 <- 3*(region1 - 1) * sin(2 * pi * (region1 + 0.2))
  Betapart2 <- rep(0, length(region2))
  Betapart3 <- - 3 * region3 * sin(2 * pi * (region3 - 0.2))
  beta0value <- c( Betapart1, Betapart2, Betapart3)
  beta1value <- c(Betapart1, Betapart2, rep(0, length(Betapart3)))
  beta2value <- c(rep(0, length(Betapart1)), Betapart2, Betapart3)

  ### Generate response y
  myfdPar <- fdPar(fbasisX, Lfdobj = 2, lambda = 5e-4)
  beta0fd <- smooth.basis(t, beta0value, myfdPar)$fd
  beta1fd <- smooth.basis(t, beta1value, myfdPar)$fd
  beta2fd <- smooth.basis(t, beta2value, myfdPar)$fd
  basisint0 <- inprod(fdobj1 = fbasisX, fdobj2 = beta0fd, Lfdobj1 = 0, Lfdobj2 = 0)
  basisint1 <- inprod(fdobj1 = fbasisX, fdobj2 = beta1fd, Lfdobj1 = 0, Lfdobj2 = 0)
  basisint2 <- inprod(fdobj1 = fbasisX, fdobj2 = beta2fd, Lfdobj1 = 0, Lfdobj2 = 0)
  basisMatrix <- getbasismatrix(evalarg = t, basisobj = fbasisX, nderiv = 0, returnMatrix = FALSE)
  funcY <- function(i){
    value <- t(z[i,])%*%gamma + t(dataX[i, ])%*%basisMatrix%*%solve(t(basisMatrix) %*% basisMatrix) %*%basisint0 +
             z[i, 1] * (t(dataX[i, ])%*%basisMatrix%*%solve(t(basisMatrix) %*% basisMatrix) %*%basisint1) +
             z[i, 2] * (t(dataX[i, ])%*%basisMatrix%*%solve(t(basisMatrix) %*% basisMatrix) %*%basisint2) + epsilon[i]
    return(value)
    }
  y <- mapply(funcY, 1:n)
  simData <- list(y = y, z = z, location = t, X = dataX ,
                   beta = list(beta0 = beta0fd, beta1 = beta1fd, beta2 = beta2fd))
  class(simData) <- "simData2"
  return(simData)
}

