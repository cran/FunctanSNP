##' @title Cross-validation for \code{\link{SNPinter}}
##'
##' @description Performs K-fold cross validation for the revised partially functional interaction regression analysis over a grid of values for the regularization parameter lambda1 and lambda2.
##'
##' @param y a numeric vector defining the response variables.
##' @param z a matrix defining the scalar covariates, with the number of rows equal to the number of samples.
##' @param location a numeric vector defining the sampling sites of the sequence data.
##' @param X a matrix specifying the sequence (genotypes) data, with the number of rows equal to the number of samples.
##' @param K an integer specifying the number of cross-validation folds, default is 5.
##' @param lambda1 a numeric vector specifying the sparsity penalty parameter to be determined.
##' @param lambda2 a numeric vector specifying the group sparsity penalty parameter to be determined.
##' @param eta a numeric vector specifying the penalty parameter for smoothing analysis.
##' @param type1 a character specifying the type of the basis functions that constitutes the genetic variation function. The options are "Bspline", "Exponential", "Fourier", "Monomial", and "Power".
##' @param nbasis1 an integer specifying the number of basis functions that constitutes the genetic variation function.
##' @param params1 in addition to rangeval1 (a vector of length 2 giving the lower and upper limits of the range of permissible values for the genetic variation function) and nbasis1, all bases have one or two parameters unique to that basis type or shared with one other;
##' \itemize{
##' \item{bspline:}{ Argument norder = the order of the spline, which is one more than the degree of the polynomials used. This defaults to 4, which gives cubic splines.}
##' \item{exponential:}{ Argument ratevec. In fda_2.0.2, this defaulted to 1. In fda_2.0.3, it will default to 0:1.}
##' \item{fourier:}{ Argument period defaults to diff(rangeval).}
##' \item{monomial/power:}{ Argument exponents. Default = 0:(nbasis-1). For monomial bases, exponents must be distinct nonnegative integers. For power bases, they must be distinct real numbers.}
##' }
##' @param Bsplines an integer specifying the number of basis functions that constitutes the genetic effect function.
##' @param norder an integer specifying the order of bsplines that constitutes the genetic effect function, which is one higher than their degree. The default of 4 gives cubic splines.
##' @param intercept should intercept(s) be fitted (TRUE) or set to zero (default = FALSE).
##' @param eps a numeric variable specifying the threshold at which the algorithm terminates, default is 1e-5.
##' @param maxstep a numeric variable specifying the maximum iteration steps, default is 1e5.
##' @param Plot should the estimated genetic effect function beta0(t) and interaction items betak(t) be plotted (TRUE) or not (default = FALSE).
##'
##' @return An "SNPcvinter" object that contains the list of the following items.
##' \itemize{
##' \item{lambda1Select:}{ a numeric value of the sparsity penalty parameter selected by cross validation.}
##' \item{lambda2Select:}{ a numeric value of the group sparsity penalty parameter selected by cross validation.}
##' \item{etaSelect:}{ a numeric value of the smoothing parameter selected by cross validation.}
##' \item{alpha:}{ estimated intercept value..}
##' \item{gamma:}{ estimated coefficients of the scalar covariates.}
##' \item{b:}{ estimated coefficients of the chosen basis functions for the genetic effect function beta0(t) and interaction items betak(t).}
##' \item{betat:}{ an "fd" object, representing the estimated genetic effect function beta(t) and interaction items betak(t).}
##' \item{residuals:}{ the residuals, that is response minus fitted values.}
##' \item{fitted.values: }{ the fitted mean values.}
##' \item{lambda1:}{ a numeric vector specifying the sparsity penalty parameter for cross validation.}
##' \item{lambda2:}{ a numeric vector specifying the group sparsity penalty parameter for cross validation.}
##' \item{eta:}{ a numeric vector specifying the smoothing parameter for cross validation.}
##' \item{CVerror:}{ a numeric vector, containing the mean square errors on testing set during cross validation.}
##' }
##' @seealso See Also as \code{\link{simData2}}, \code{\link{SNPinter}}.
##'
##' @import fda
##' @import glmnet
##' @importFrom caret createFolds
##' @importFrom lava blockdiag
##' @importFrom stats integrate
##' @export
##' @examples
##' library(FunctanSNP)
##' n <- 300
##' m <- 30
##' simdata2 <- simData2(n, m, seed = 123)
##' y <- simdata2$y
##' z <- simdata2$z
##' location <- simdata2$location
##' X <- simdata2$X
##' lambda1 <- c(0.01, 0.05, 0.1)
##' lambda2 <- sqrt(3)*lambda1
##' SNPcvinterres <- SNPcvinter(y, z, location, X, K = 3, lambda1, lambda2, eta = 0,
##'                             type1 = "Bspline", nbasis1 = 5, params1 = 4, Bsplines = 5,
##'                             norder = 4, intercept = TRUE, eps = 1e-2, maxstep = 1e2, Plot = TRUE)
##' SNPcvinterres$lambda1Select
##' SNPcvinterres$lambda2Select
##' SNPcvinterres$alpha
##' SNPcvinterres$gamma
##' SNPcvinterres$b
##'

#library(lava)
#library(glmnet)
#library(caret)
SNPcvinter <- function(y, z, location, X, K, lambda1, lambda2 = NULL, eta = 0,
                       type1, nbasis1, params1, Bsplines = 20, norder = 4,
                       intercept = FALSE, eps = 1e-5, maxstep = 1e5, Plot = FALSE){
  if (!inherits(location, "numeric")) { stop("x should be of vector type.") }
  if (!inherits(y, "numeric")) { stop("y should be of numeric type.") }
  if (!inherits(z, "matrix")) { stop("z should be of matrix type.") }
  if (!inherits(X, "matrix")) { stop("X should be of matrix type.") }
  if (!inherits(lambda1, "numeric")) { stop("lambda1 should be of numeric type, such as a value or vector.")}

  n <- nrow(X)
  m <- ncol(X)
  q <- ncol(z)
  if ( !inherits(lambda2, "numeric") & is.null(lambda2) ) { lambda2 <- sqrt(q+1)*lambda1 }
  if ( !inherits(lambda2, "numeric") & !is.null(lambda2) ) { stop("lambda1 should be of numeric type, such as a value or vector.") }

  foldsList <- createFolds(1:n, K)

  cvinter <- function(k){
    testIndex <- foldsList[[k]]
    trainIndex <- setdiff(1:n, testIndex)
    yTest <- y[testIndex]
    zTest <- z[testIndex, ]
    XTest <- X[testIndex, ]
    yTrain <- y[trainIndex]
    zTrain <- z[trainIndex, ]
    XTrain <- X[trainIndex, ]
    SNPcvtemp <- SNPinter(y = yTrain, z = zTrain, location = location, X = XTrain, lambda1 = lambdaTemp1, lambda2 = lambdaTemp2, eta = etaTemp,
                          type1 = type1, nbasis1 = nbasis1, params1 = params1, Bsplines = Bsplines, norder = norder,
                          intercept = intercept, eps = eps, maxstep = maxstep, Plot = FALSE)
    alpha_temp <- SNPcvtemp$alpha
    gamma_temp <- SNPcvtemp$gamma
    thetalist_temp <- SNPcvtemp$b

    funcXTest <- SNPgvf(location = location, X = XTest, type = type1, nbasis = nbasis1, params = params1, Plot = FALSE)

    if(type1 == "Bspline"){
      fbasis1 <- create.bspline.basis(rangeval = c(min(location), max(location)), nbasis = nbasis1, norder = params1)
    }

    if(type1 == "Exponential"){
      fbasis1 <- create.exponential.basis(rangeval = c(min(location), max(location)), nbasis = nbasis1, ratevec = params1)
    }

    if(type1 == "Fourier"){
      fbasis1 <- create.fourier.basis(rangeval = c(min(location), max(location)), nbasis = nbasis1, period = params1)
    }

    if(type1 == "Monomial"){
      fbasis1 <- create.monomial.basis(rangeval = c(min(location), max(location)), nbasis = nbasis1, exponents = params1)
    }

    if(type1 == "Power"){
      fbasis1 <- create.power.basis(rangeval = c(min(location), max(location)), nbasis = nbasis1, exponents = params1)
    }

    fbasis2 <- create.bspline.basis(rangeval = c(min(location), max(location)), nbasis = Bsplines, norder = norder)

    nTest <- nrow(XTest)
    mTest <- ncol(XTest)
    qTest <- ncol(zTest)

    Mn <- Bsplines - norder + 1
    d <- norder - 1
    funcCoefTest <- t(funcXTest$coefs)
    basisint <- inprod(fdobj1 = fbasis1, fdobj2 = fbasis2, Lfdobj1 = 0, Lfdobj2 = 0)
    funcUTest <- function(i){ funcCoefTest[i, ] %*% basisint }
    UTest <- t(mapply(funcUTest, 1:nTest))
    tildeUTest <- lapply(1:qTest, function(j) sweep(UTest, 1, zTest[, j], "*"))
    RTest <- cbind( UTest, matrix(unlist(tildeUTest), nrow = nTest, byrow = FALSE) )

    if(intercept == FALSE){
      fittedvaluesTest <- zTest%*%gamma_temp + RTest%*%unlist(thetalist_temp)
      residualsTest <- yTest - fittedvaluesTest
      seTest <- sum(residualsTest^2)
    }

    if(intercept == TRUE){
      fittedvaluesTest <- alpha_temp + zTest%*%gamma_temp + RTest%*%unlist(thetalist_temp)
      residualsTest <- yTest - fittedvaluesTest
      seTest <- sum(residualsTest^2)
    }
    return(seTest)
  }

  lambdaGrid <- expand.grid(lambda1, lambda2, eta)
  Num <- nrow(lambdaGrid)
  squareError <- rep(NA, Num)
  for (l in 1:Num) {
    lambdaTemp1 <- lambdaGrid[l, 1]
    lambdaTemp2 <- lambdaGrid[l, 2]
    etaTemp <- lambdaGrid[l, 3]
    squareError[l] <- sum(mapply(cvinter, 1:K)) / n #1:K
  }
  lambdaIndex <- which.min(squareError)[1]
  lambdaSelect1 <- lambdaGrid[lambdaIndex, 1]
  lambdaSelect2 <- lambdaGrid[lambdaIndex, 2]
  etaSelect <- lambdaGrid[lambdaIndex, 3]

  SNPtemp <- SNPinter(y = y, z = z, location = location, X = X, lambda1 = lambdaSelect1, lambda2 = lambdaSelect2, eta = etaSelect,
                        type1 = type1, nbasis1 = nbasis1, params1 = params1, Bsplines = Bsplines, norder = norder,
                        intercept = intercept, eps = eps, maxstep = maxstep, Plot = Plot)

  result <- list(lambda1Select = lambdaSelect1, lambda2Select  = lambdaSelect2, etaSelect = etaSelect,
                 alpha = SNPtemp$alpha, gamma = SNPtemp$gamma, b = SNPtemp$b, betat = SNPtemp$betat,
                 residuals = SNPtemp$residuals, fitted.values = SNPtemp$fittedvalues,
                 lambda1 = lambda1, lambda2 = lambda2, eta = eta, CVerror = squareError)
  class(result) <- "SNPcvinter"
  return(result)
}


