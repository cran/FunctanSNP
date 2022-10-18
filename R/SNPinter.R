##' @title Revised partially functional interaction regression analysis for the sequence (genotypes) data
##'
##' @description This function conducts joint analysis, which includes all scalar covariates Z, genetic variant function X(t), and their interactions in a partially functional interaction regression model. In addition, this function identifies the relevant genetic variation sites with a local sparsity penalty-based method.
##'
##' @param y a numeric vector defining the response variables.
##' @param z a matrix defining the scalar covariates, with the number of rows equal to the number of samples.
##' @param location a numeric vector defining the sampling sites of the sequence data.
##' @param X a matrix specifying the sequence (genotypes) data, with the number of rows equal to the number of samples.
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
##' \item{para:}{ some relevant parameters of "SNPlm" model.}
##' }
##'
##' @param Bsplines an integer specifying the number of basis functions that constitutes the genetic effect function.
##' @param norder an integer specifying the order of bsplines that constitutes the genetic effect function, which is one higher than their degree. The default of 4 gives cubic splines.
##' @param intercept should intercept(s) be fitted (TRUE) or set to zero (default = FALSE).
##' @param eps a numeric variable specifying the threshold at which the algorithm terminates, default is 1e-5.
##' @param maxstep a numeric variable specifying the maximum iteration steps, default is 1e5.
##' @param Plot should the estimated genetic effect function beta0(t) and interaction items betak(t) be plotted (TRUE) or not (default = FALSE).
##'
##' @return An "SNPinter" object that contains the list of the following items.
##' \itemize{
##' \item{alpha:}{ estimated intercept value..}
##' \item{gamma:}{ estimated coefficients of the scalar covariates.}
##' \item{b:}{ estimated coefficients of the chosen basis functions for the genetic effect function beta0(t) and interaction items betak(t).}
##' \item{betat:}{ an "fd" object, representing the estimated genetic effect function beta(t) and interaction items betak(t).}
##' \item{residuals:}{ the residuals, that is response minus fitted values.}
##' \item{fitted.values: }{ the fitted mean values.}
##' }
##' @seealso See Also as \code{\link{simData1}}, \code{\link{SNPcvinter}}.
##'
##' @import fda
##' @import glmnet
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
##' lambda1 <- 0.05
##' lambda2 <- sqrt(3)*lambda1
##' eta <- 0
##' SNPinterres <- SNPinter(y, z, location, X, lambda1, lambda2, eta, type1 = "Bspline",
##'                         nbasis1 = 5, params1 = 4, Bsplines = 5, norder = 4, intercept = TRUE,
##'                         eps = 1e-2, maxstep = 1e2, Plot = TRUE)
##' SNPinterres$alpha
##' SNPinterres$gamma
##' SNPinterres$b
##'

#library(lava)
#library(glmnet)
SNPinter <- function(y, z, location, X, lambda1, lambda2 = NULL, eta,
                     type1, nbasis1, params1, Bsplines = 20, norder = 4,
                     intercept = FALSE, eps = 1e-5, maxstep = 1e5, Plot = FALSE){
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

  fbasis2 <- create.bspline.basis(rangeval = c(min(location), max(location)), nbasis = Bsplines, norder = norder)

  ro <- function(x, mu, alpha) {
    if(mu == 0){ r <- 0 }
    if(mu != 0){
      f <- function(x) mu * (1 > x / (mu * alpha)) * (1 - x / (mu * alpha))
      r <- integrate(f, 0, x)
    }
    return(r)
  }

  ro_d1st <- function(x, mu, alpha) {
    if(mu == 0){ r <- 0 }
    if(mu != 0){ r <- mu * (1 > x / (mu * alpha)) * (1 - x / (mu * alpha)) }
    return(r)
  }

  n <- nrow(X)
  m <- ncol(X)
  q <- ncol(z)
  if (is.null(lambda2)) { lambda2 <- sqrt(q+1)*lambda1 }

  Mn <- Bsplines - norder + 1
  d <- norder - 1
  funcCoef <- t(funcX$coefs)
  basisint <- inprod(fdobj1 = fbasis1, fdobj2 = fbasis2, Lfdobj1 = 0, Lfdobj2 = 0)
  funcU <- function(i){ funcCoef[i, ] %*% basisint }
  U <- t(mapply(funcU, 1:n))
  tildeU <- lapply(1:q, function(j) sweep(U, 1, z[, j], "*"))
  R <- cbind( U, matrix(unlist(tildeU), nrow = n, byrow = FALSE) )
  rawD <- cbind(z, R)

  #rawDmean <- apply(rawD, 2, mean)
  #rawDscale <- apply(rawD, 2, sd)
  #D <- scale(rawD, center = FALSE, scale = TRUE)
  D <- rawD

  if(lambda1 == 0 & lambda2 == 0){
    cv.lambda <- 0
  }else{
    cv.model <- cv.glmnet(x = D, y = y)
    cv.lambda <- cv.model$lambda[which.min(cv.model$cvsd)]
  }
  model <- glmnet(x = D, y = y, lambda = cv.lambda, intercept = intercept)
  Coef <- as.numeric(coef.glmnet(model))

  V <- inprod(fdobj1 = fbasis2, fdobj2 = fbasis2, Lfdobj1 = 2, Lfdobj2 = 2)
  knots <- seq(min(location), max(location), length = Mn + 1)
  Wlist <- lapply(1:Mn, function(l) Mn * inprod(fdobj1 = fbasis2, fdobj2 = fbasis2, Lfdobj1 = 0, Lfdobj2 = 0, rng = c(knots[l], knots[l+1])) / (max(location) - min(location)) )

  if(intercept == FALSE){
    Coef <- Coef[-1]
    alpha0 <- 0
    gamma0 <- Coef[1:q]
    thetalist0 <- lapply(1:(q+1), function(l) Coef[((l-1)*Bsplines+1+q):(l*Bsplines+q)])

    omega_s1 <- omega_s <- Coef
    alpha_s1 <- alpha_s <- alpha0
    gamma_s1 <- gamma_s <- gamma0
    thetalist_s1 <- thetalist_s <-thetalist0

    get_penalty_theta_W0 <- function(l){
      W <- Wlist[[l]]
      norm_temp <- mapply(function(k) t(thetalist_s[[k]]) %*% W %*% thetalist_s[[k]], 1:(q+1))
      norm_sqrt <- sqrt(sum(norm_temp))
      norm_pen <- ro_d1st(norm_sqrt, mu = lambda2, alpha = 6) / (norm_sqrt + 1e-10)
      return( norm_pen )
    }
    get_penalty_thetak_Wk <- function(k, l){
      norm_temp <- t(thetalist_s[[k]]) %*% Wlist[[l]] %*% thetalist_s[[k]]
      norm_sqrt <- sqrt(norm_temp)
      norm_pen <- ro_d1st(norm_sqrt, mu = lambda1, alpha = 6) / (norm_sqrt + 1e-10)
      return( norm_pen )
    }
    get_breve_Wk <- function(k){
      prime_lambda1_temp <- prime_lambda1[(Mn*(k-1)+1):(Mn*k)]
      breve_Wk_temp <- matrix(rep(prime_lambda1_temp + prime_lambda2, each = Bsplines*Bsplines) * unlist(Wlist), ncol = Mn, byrow = FALSE)
      breve_Wk <- 0.5 * matrix(apply(breve_Wk_temp, 1, sum), ncol = Bsplines, byrow = FALSE)
      return(breve_Wk)
    }

    s <- 1
    omegaDiff <- 10
    while (omegaDiff > eps & s <= maxstep){
      gamma_s <- omega_s[1:q]
      thetalist_s <- lapply(1:(q+1), function(l) omega_s[((l-1)*Bsplines+1+q):(l*Bsplines+q)])

      prime_lambda1 <- mapply(get_penalty_thetak_Wk, rep(2:(q+1), each = Mn), rep(1:Mn, times = q))
      prime_lambda2 <- mapply(get_penalty_theta_W0, 1:Mn)
      breve_W0_temp <- matrix(rep(prime_lambda2, each = Bsplines*Bsplines) * unlist(Wlist), ncol = Mn, byrow = FALSE)
      breve_W0 <- 0.5* matrix(apply(breve_W0_temp, 1, sum), ncol = Bsplines, byrow = FALSE)

      tilde_W_temp <- blockdiag( diag(0, q), breve_W0)
      for (k in 1:q) { tilde_W_temp <- blockdiag(tilde_W_temp, get_breve_Wk(k = k)) }
      tilde_W_s <- tilde_W_temp

      tilde_V_temp <- diag(0, q) #Without intercept
      for (k in 1:(q+1)) { tilde_V_temp <- blockdiag(tilde_V_temp, V) }
      tilde_V_s <- tilde_V_temp

      solveMatrix <- t(D) %*% D + 2*n*tilde_W_s + 2*n*eta*tilde_V_s
      solveFunc <- function(lambda){ tryCatch(solve(solveMatrix + lambda*diag(1,ncol(tilde_W_s))), error = function(e) { matrix(0, ncol = ncol(solveMatrix), nrow = nrow(solveMatrix)) }) }
      seqList <- c(0, 1e-20, 1e-19, 1e-18, 1e-17, 1e-16,
                   1e-15, 1e-14, 1e-13, 1e-12, 1e-11,
                   1e-10, 1e-9, 1e-8, 1e-7, 1e-6,
                   1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1)
      invList <- lapply(1:length(seqList), function(l) solveFunc(seqList[l]))
      selectLambda <- seqList[ which(lapply(1:length(seqList), function(l) norm(invList[[l]], type = "F") ) != 0 )[1] ]
      omega_s1 <- solveFunc(lambda = selectLambda) %*% t(D) %*% y

      omegaDiff <- sum((omega_s1 - omega_s)^2) / sum(omega_s^2)
      if(omegaDiff <= eps){
        gamma_s1 <- omega_s1[1:q]
        thetalist_s1 <- lapply(1:(q+1), function(l) omega_s1[((l-1)*Bsplines+1+q):(l*Bsplines+q)])
        break
      }
      omega_s <- omega_s1
      s <- s + 1
    }
    alpha_s1 <- alpha_s
    gamma_s1 <- omega_s1[1:q]
    thetalist_temp <- omega_s1[-(1:q)]

    if(lambda1 != 0 & lambda2 != 0){
      thetalist_temp <- ifelse(abs(thetalist_temp) > 1e-5, thetalist_temp, 0)
    }
    if(lambda1 != 0 & lambda2 == 0){
      thetalist_temp[-(1:Bsplines)] <- ifelse(abs(thetalist_temp[-(1:Bsplines)]) > 1e-5, thetalist_temp[-(1:Bsplines)], 0)
    }

    omega_s1[-(1:q)] <- thetalist_temp
    thetalist_s1 <- lapply(1:(q+1), function(l) omega_s1[((l-1)*Bsplines+1+q):(l*Bsplines+q)])

    fittedvalues <- z%*%gamma_s1 + R%*%unlist(thetalist_s1)
    residuals <- y - fittedvalues
    mse <- mean(residuals^2)
  }

  if(intercept == TRUE){
    intcept <- rep(1, length = nrow(D))
    newD <- cbind(intcept, D)
    alpha0 <- Coef[1]
    gamma0 <- Coef[2:(q+1)]
    thetalist0 <- lapply(1:(q+1), function(l) Coef[((l-1)*Bsplines+2+q):(l*Bsplines+1+q)])

    omega_s1 <- omega_s <- Coef
    alpha_s1 <- alpha_s <- alpha0
    gamma_s1 <- gamma_s <- gamma0
    thetalist_s1 <- thetalist_s <-thetalist0

    get_penalty_theta_W0 <- function(l){
      norm_temp <- mapply(function(k) t(thetalist_s[[k]]) %*% Wlist[[l]] %*% thetalist_s[[k]], 1:(q+1))
      norm_sqrt <- sqrt(sum(norm_temp))
      norm_pen <- ro_d1st(norm_sqrt, mu = lambda2, alpha = 6) / (norm_sqrt+1e-10)
      return( norm_pen )
    }
    get_penalty_thetak_Wk <- function(k, l){
      norm_temp <- t(thetalist_s[[k]]) %*% Wlist[[l]] %*% thetalist_s[[k]]
      norm_sqrt <- sqrt(norm_temp)
      norm_pen <- ro_d1st(norm_sqrt, mu = lambda1, alpha = 6) / (norm_sqrt+1e-10)
      return( norm_pen )
    }
    get_breve_Wk <- function(k){
      prime_lambda1_temp <- prime_lambda1[(Mn*(k-1)+1):(Mn*k)]
      breve_Wk_temp <- matrix(rep(prime_lambda1_temp + prime_lambda2, each = Bsplines*Bsplines) * unlist(Wlist), ncol = Mn, byrow = FALSE)
      breve_Wk <- 0.5 * matrix(apply(breve_Wk_temp, 1, sum), ncol = Bsplines, byrow = FALSE)
      return(breve_Wk)
    }
    s <- 1
    omegaDiff <- 10
    while (omegaDiff > eps & s <= maxstep){
      alpha_s <- omega_s[1]
      gamma_s <- omega_s[2:(q+1)]
      thetalist_s <- lapply(1:(q+1), function(l) omega_s[((l-1)*Bsplines+2+q):(l*Bsplines+1+q)])

      prime_lambda1 <- mapply(get_penalty_thetak_Wk, rep(2:(q+1), each = Mn), rep(1:Mn, times = q))
      prime_lambda2 <- mapply(get_penalty_theta_W0, 1:Mn)
      breve_W0_temp <- matrix(rep(prime_lambda2, each = Bsplines*Bsplines) * unlist(Wlist), ncol = Mn, byrow = FALSE)
      breve_W0 <- 0.5* matrix(apply(breve_W0_temp, 1, sum), ncol = Bsplines, byrow = FALSE)

      tilde_W_temp <- blockdiag( diag(0, (q+1)), breve_W0)
      for (k in 1:q) { tilde_W_temp <- blockdiag(tilde_W_temp, get_breve_Wk(k = k)) }
      tilde_W_s <- tilde_W_temp

      tilde_V_temp <- diag(0, (q+1)) #Having intercept
      for (k in 1:(q+1)) { tilde_V_temp <- blockdiag(tilde_V_temp, V) }
      tilde_V_s <- tilde_V_temp

      solveMatrix <- t(newD) %*% newD + 2*n*tilde_W_s + 2*n*eta*tilde_V_s
      solveFunc <- function(lambda){ tryCatch(solve(solveMatrix + lambda*diag(1,ncol(tilde_W_s))), error = function(e) { matrix(0, ncol = ncol(solveMatrix), nrow = nrow(solveMatrix)) }) }
      seqList <- c(0, 1e-20, 1e-19, 1e-18, 1e-17, 1e-16,
                   1e-15, 1e-14, 1e-13, 1e-12, 1e-11,
                   1e-10, 1e-9, 1e-8, 1e-7, 1e-6,
                   1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1)
      invList <- lapply(1:length(seqList), function(l) solveFunc(seqList[l]))
      selectLambda <- seqList[ which(lapply(1:length(seqList), function(l) norm(invList[[l]], type = "F") ) != 0 )[1] ]
      omega_s1 <- solveFunc(lambda = selectLambda) %*% t(newD) %*% y

      omegaDiff <- sum((omega_s1 - omega_s)^2) / sum(omega_s^2)
      if(omegaDiff <= eps){
        alpha_s1 <- omega_s1[1]
        gamma_s1 <- omega_s1[2:(q+1)]
        thetalist_s1 <- lapply(1:(q+1), function(l) omega_s1[((l-1)*Bsplines+2+q):(l*Bsplines+1+q)])
        break
      }
      omega_s <- omega_s1
      s <- s + 1
    }
    alpha_s1 <- omega_s1[1]
    gamma_s1 <- omega_s1[2:(q+1)]
    thetalist_temp <- omega_s1[-(1:(q+1))]

    if(lambda1 != 0 & lambda2 != 0){
      thetalist_temp <- ifelse(abs(thetalist_temp) > 1e-5, thetalist_temp, 0)
    }
    if(lambda1 != 0 & lambda2 == 0){
      thetalist_temp[-(1:Bsplines)] <- ifelse(abs(thetalist_temp[-(1:Bsplines)]) > 1e-5, thetalist_temp[-(1:Bsplines)], 0)
    }

    omega_s1[-(1:(q+1))] <- thetalist_temp
    thetalist_s1 <- lapply(1:(q+1), function(l) omega_s1[((l-1)*Bsplines+2+q):(l*Bsplines+1+q)])

    fittedvalues <- alpha_s1 + z%*%gamma_s1 + R%*%unlist(thetalist_s1)
    residuals <- y - fittedvalues
    mse <- mean(residuals^2)
  }

  #basisCoef <- matrix(unlist(thetalist_s1)/rawDscale[-(1:q)], nrow = (q+1), byrow = TRUE)
  basisCoef <- matrix(unlist(thetalist_s1), nrow = (q+1), byrow = TRUE)
  betat <- lapply(1:(q+1), function(l) fd(coef = basisCoef[l, ], basisobj = fbasis2))
  names(betat) <- paste("beta", 1:(q+1)-1, "(t)", sep = "")
  if(Plot == TRUE){
    for(l in 1:(q+1)){
      plot(betat[[l]], xlab = "Location", ylab = paste("beta", l-1, "(t)", sep = ""))
    }
  }

  names(thetalist_s1) <- paste("b", 1:(q+1)-1, sep = "")
  para <- list(location = location, type1 = type1, nbasis1 = nbasis1, params1 = params1,
               Bsplines = Bsplines, norder = norder)
  result <- list(alpha = alpha_s1, gamma = gamma_s1, b = thetalist_s1, betat = betat,
                 residuals = residuals, fitted.values = fittedvalues, para = para)
  class(result) <- "SNPinter"
  return(result)
}


