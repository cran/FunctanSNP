##' @title Plot coefficient functions and residuals from an SNPlm or SNPinter object
##'
##' @description Produces a residual diagram or a plot of the genetic effect function beta(t) (and interaction items betak(t), if the entry x is an "SNPinter" object) for a fitted "SNPlm" or "SNPinter" object.
##'
##' @param x an "SNPlm"/"SNPinter" object obtained by \code{\link{SNPlm}}/\code{\link{SNPinter}}.
##' @param type character, "beta" or "residuals" type, if "beta", plot the genetic effect function beta0(t) (and interaction items betak(t), if the entry x is an "SNPinter" object), if "residuals", show the residual plot.
##'
##' @return show the residual plot or coefficient functions.
##'
##' @seealso See Also as \code{\link{SNPlm}}, \code{\link{SNPinter}}.
##'
##' @importFrom fda plot.fd
##' @importFrom graphics abline
##' @export
##' @examples
##' library(FunctanSNP)
##' n <- 300
##' m <- 30
##' simdata1 <- simData1(n, m, seed = 123)
##' SNPlmres <- SNPlm(y = simdata1$y, z = simdata1$z,
##'                   location = simdata1$location, X = simdata1$X,
##'                   type1 = "Bspline", type2 = "Bspline", nbasis1 = 5,
##'                   nbasis2 = 5, params1 = 4, params2 = 4,
##'                   intercept = FALSE, Plot = FALSE)
##' plotSNPreg(x = SNPlmres, type = "beta")
##' plotSNPreg(x = SNPlmres, type = "residuals")
##'
##' simdata2 <- simData2(n, m, seed = 123)
##' lambda1 <- 0.05
##' lambda2 <- sqrt(3)*lambda1
##' eta <- 0
##' SNPinterres <- SNPinter(y = simdata2$y, z = simdata2$z,
##'                         location = simdata2$location, X = simdata2$X,
##'                         lambda1, lambda2, eta, type1 = "Bspline", nbasis1 = 5,
##'                         params1 = 4, Bsplines = 5, norder = 4, intercept = TRUE,
##'                         eps = 1e-2, maxstep = 1e2, Plot = FALSE)
##' plotSNPreg(x = SNPinterres, type = "beta")
##' plotSNPreg(x = SNPinterres, type = "residuals")
##'

plotSNPreg <- function(x, type){
  #x <- SNPinterres
  if(type == "beta"){
    if(inherits(x, "SNPlm")){
      plot.fd(x = x$betat, xlab = "location", ylab = "X")
    }
    if(inherits(x, "SNPinter")){
      Plotnum <- length(x$betat)
      for (l in 1:Plotnum) { plot.fd(x = x$betat[[l]], xlab = "Location", ylab = paste("beta", l-1, "(t)", sep = "")) }
    }
  }

  if(type == "residuals"){
    plot(x$residuals, cex = 0.9, xlab = "subjects", ylab = "residuals") + abline(h = 0, lty = 2, lwd = 2, col="red")
  }
}
