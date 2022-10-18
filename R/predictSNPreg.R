##' @title Make predictions from an SNPlm or SNPinter object
##'
##' @description This functions predicts fitted values for newly entered data from a fitted "SNPlm" or "SNPinter" object.
##'
##' @param x an "SNPlm"/"SNPinter" object obtained by \code{\link{SNPlm}}/\code{\link{SNPinter}}.
##' @param newz matrix of new values for z at which predictions are to be made.
##' @param newX matrix of new values for X at which predictions are to be made.
##'
##' @return
##' \itemize{
##' \item{predicted.values: }{ the predicted mean values.}
##' }
##' @seealso See Also as \code{\link{SNPlm}}, \code{\link{SNPinter}}.
##'
##' @import fda
##' @export
##' @examples
##' library(FunctanSNP)
##' n <- 300
##' m <- 30
##' simdata1 <- simData1(n, m, seed = 123)
##' SNPlmres <- SNPlm(y = simdata1$y[1:200], z = simdata1$z[1:200, ],
##'                   location = simdata1$location, X = simdata1$X[1:200, ],
##'                   type1 = "Bspline", type2 = "Bspline", nbasis1 = 5,
##'                   nbasis2 = 5, params1 = 4, params2 = 4, intercept = FALSE,
##'                   Plot = FALSE)
##' predict.values1 <- predictSNPreg(x = SNPlmres, newz = simdata1$z[-(1:200), ],
##'                                  newX = simdata1$X[-(1:200), ])
##'
##' simdata2 <- simData2(n, m, seed = 123)
##' lambda1 <- 0.05
##' lambda2 <- sqrt(3)*lambda1
##' eta <- 0
##' SNPinterres <- SNPinter(y = simdata2$y[1:200], z = simdata2$z[1:200, ],
##'                         location = simdata2$location, X = simdata2$X[1:200, ],
##'                         lambda1, lambda2, eta, type1 = "Bspline",
##'                         nbasis1 = 5, params1 = 4, Bsplines = 5, norder = 4,
##'                         intercept = TRUE, eps = 1e-2, maxstep = 1e2, Plot = FALSE)
##' predict.values2 <- predictSNPreg(x = SNPinterres, newz = simdata2$z[-(1:200), ],
##'                                  newX = simdata2$X[-(1:200), ])
##'

predictSNPreg <- function(x, newz, newX){
  if (!inherits(x, "SNPlm") & !inherits(x, "SNPinter")) { stop("x should be of 'SNPlm' or 'SNPinter' type.") }
  if (!inherits(newz, "matrix")) { stop("newz should be of matrix type.") }
  if (!inherits(newX, "matrix")) { stop("newX should be of matrix type.") }

  #x <- SNPinterres
  n <- nrow(newX)
  m <- ncol(newX)
  q <- ncol(newz)
  location <- x$para$location
  type1 <- x$para$type1
  nbasis1 <- x$para$nbasis1
  params1 <- x$para$params1
  if(inherits(x, "SNPlm")){
    nbasis2 <- x$para$nbasis2
    params2 <- x$para$params2
  }
  if(inherits(x, "SNPinter")){
    nbasis2 <- x$para$Bsplines
    params2<- x$para$norder
  }
  alpha <- x$alpha
  gamma <- x$gamma
  thetalist <- x$b

  funcX <- SNPgvf(location, X = newX, type = type1, nbasis = nbasis1, params = params1, Plot = FALSE)

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

  fbasis2 <- create.bspline.basis(rangeval = c(min(location), max(location)), nbasis = nbasis2, norder = params2)

  funcCoef <- t(funcX$coefs)
  basisint <- inprod(fdobj1 = fbasis1, fdobj2 = fbasis2, Lfdobj1 = 0, Lfdobj2 = 0)
  funcU <- function(i){ funcCoef[i, ] %*% basisint }
  U <- t(mapply(funcU, 1:n))

  if(inherits(x, "SNPlm")){
    fittedvalues <- alpha + newz%*%gamma + U%*%thetalist
  }

  if(inherits(x, "SNPinter")){
    tildeU <- lapply(1:q, function(j) sweep(U, 1, newz[, j], "*"))
    R <- cbind( U, matrix(unlist(tildeU), nrow = n, byrow = FALSE) )
    fittedvalues <- alpha + newz%*%gamma + R%*%unlist(thetalist)
  }
  return(fittedvalues)
}
