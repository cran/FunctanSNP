##' @title Revised functional linear regression analysis for the sequence (genotypes) data
##'
##' @description This function models the genetic effect of genetic variants by relating the genetic variant function to the phenotype adjusting for covariates.
##'
##' @param y a numeric vector specifying the response variables.
##' @param z a matrix specifying the scalar covariates, with the number of rows equal to the number of samples.
##' @param location a numeric vector defining the sampling sites of the sequence data.
##' @param X a matrix specifying the sequence data, with the number of rows equal to the number of samples.
##' @param type1 a character specifying the type of the basis functions that constitutes the genetic variation function. The options are "Bspline", "Exponential", "Fourier", "Monomial", and "Power".
##' @param type2 a character specifying the type of the basis functions that constitutes the genetic effect function. The options are "Bspline", "Exponential", "Fourier", "Monomial", and "Power".
##' @param nbasis1 an integer specifying the number of basis functions that constitutes the genetic variation function.
##' @param nbasis2 an integer specifying the number of basis functions that constitutes the genetic effect function.
##' @param params1 in addition to rangeval1 (a vector of length 2 giving the lower and upper limits of the range of permissible values for the genetic variation function) and nbasis1, all bases have one or two parameters unique to that basis type or shared with one other;
##' @param params2 in addition to rangeval1 (a vector of length 2 giving the lower and upper limits of the range of permissible values for the genetic effect function) and nbasis1, all bases have one or two parameters unique to that basis type or shared with one other;
##' \itemize{
##' \item{bspline:}{ Argument norder = the order of the spline, which is one more than the degree of the polynomials used. This defaults to 4, which gives cubic splines.}
##' \item{exponential:}{ Argument ratevec. In fda_2.0.2, this defaulted to 1. In fda_2.0.3, it will default to 0:1.}
##' \item{fourier:}{ Argument period defaults to diff(rangeval).}
##' \item{monomial/power:}{ Argument exponents. Default = 0:(nbasis-1). For monomial bases, exponents must be distinct nonnegative integers. For power bases, they must be distinct real numbers.}
##' }
##' @param intercept should intercept(s) be fitted (TRUE) or set to zero (default = FALSE).
##' @param Plot should the estimated genetic effect function beta(t) be plotted (TRUE) or not (default = FALSE).
##'
##' @return An "SNPlm" object that contains the list of the following items.
##' \itemize{
##' \item{alpha:}{ estimated intercept value.}
##' \item{gamma:}{ estimated coefficients of the scalar covariates.}
##' \item{b:}{ estimated coefficients of the chosen basis functions for the genetic effect function beta(t).}
##' \item{betat:}{ an "fd" object, representing the estimated genetic effect function beta(t).}
##' \item{residuals:}{ the residuals, that is response minus fitted values.}
##' \item{fitted.values: }{ the fitted mean values.}
##' \item{terms:}{ the terms object used.}
##' \item{model:}{ if requested (the default), the model frame used.}
##' \item{para:}{ some relevant parameters of "SNPlm" model.}
##' }
##' @seealso See Also as \code{\link{SNPgvf}}.
##'
##' @import fda
##' @export
##' @examples
##' library(FunctanSNP)
##' n <- 300
##' m <- 30
##' simdata1 <- simData1(n, m, seed = 123)
##' y <- simdata1$y
##' z <- simdata1$z
##' location <- simdata1$location
##' X <- simdata1$X
##' SNPlmres <- SNPlm(y, z, location, X, type1 = "Bspline", type2 = "Bspline", nbasis1 = 5,
##'                   nbasis2 = 5, params1 = 4, params2 = 4, intercept = FALSE, Plot = TRUE)
##' SNPlmres$alpha
##' SNPlmres$gamma
##' SNPlmres$b

SNPlm <- function(y, z, location, X, type1, type2, nbasis1, nbasis2, params1, params2, intercept = FALSE, Plot = FALSE){
  if (!inherits(location, "numeric")) { stop("x should be of vector type.") }
  if (!inherits(y, "numeric")) { stop("y should be of numeric type.") }
  if (!inherits(z, "matrix")) { stop("z should be of matrix type.") }
  if (!inherits(X, "matrix")) { stop("X should be of matrix type.") }

  funcX <- SNPgvf(location = location, X = X, type = type1, nbasis = nbasis1, params = params1, Plot = FALSE)

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

  if(type2 == "Bspline"){
    fbasis2 <- create.bspline.basis(rangeval = c(min(location), max(location)), nbasis = nbasis2, norder = params2)
  }

  if(type2 == "Exponential"){
    fbasis2 <- create.exponential.basis(rangeval = c(min(location), max(location)), nbasis = nbasis2, ratevec = params2)
  }

  if(type2 == "Fourier"){
    fbasis2 <- create.fourier.basis(rangeval = c(min(location), max(location)), nbasis = nbasis2, period = params2)
  }

  if(type2 == "Monomial"){
    fbasis2 <- create.monomial.basis(rangeval = c(min(location), max(location)), nbasis = nbasis2, exponents = params2)
  }

  if(type2 == "Power"){
    fbasis2 <- create.power.basis(rangeval = c(min(location), max(location)), nbasis = nbasis2, exponents = params2)
  }

  n <- nrow(X)
  m <- ncol(X)
  Coef <- t(funcX$coefs)
  basisint <- inprod(fdobj1 = fbasis1, fdobj2 = fbasis2, Lfdobj1 = 0, Lfdobj2 = 0)
  funcH <- function(i){ Coef[i, ] %*% basisint  }
  H <- t(mapply(funcH, 1:n))

  lm.x <- cbind(z, H)
  q <- ncol(z)
  colName <- NULL
  for (j in 1:(q + nbasis2)) {
    if(j <= q){colName[j] = paste("Z", j, sep = "")}
    if(j > q){colName[j] = paste("H", j-q, sep = "")}
  }
  colnames(lm.x) <- colName
  wholeData <- data.frame(y, lm.x)
  if(intercept == TRUE){
    lm.model <- lm(y ~ ., data = wholeData)
    alpha <- lm.model$coefficients[1]
    zCoef <- lm.model$coefficients[2:(q+1)]
    basisCoef <- lm.model$coefficients[-c(1:(q+1))]
  }
  if(intercept == FALSE){
    lm.model <- lm(y ~ .-1, data = wholeData)
    alpha <- 0
    zCoef <- lm.model$coefficients[1:q]
    basisCoef <- lm.model$coefficients[-c(1:q)]
  }

  residuals <- lm.model$residuals
  fitted.values <- lm.model$fitted.values
  terms <- lm.model$terms
  model <- lm.model$model
  betat <- fd(coef = basisCoef, basisobj = fbasis2)
  if(Plot == TRUE){ plot(betat, xlab = "Location", ylab = "X") }

  para <- list(location = location, type1 = type1, type2 = type2, nbasis1 = nbasis1,
               nbasis2 = nbasis2, params1 = params1, params2 = params2)
  result <- list(alpha = alpha, gamma = zCoef, b = basisCoef, betat = betat, residuals = residuals,
                 fitted.values = fitted.values, terms = terms, model = model, para = para)
  class(result) <- "SNPlm"
  return(result)
}


