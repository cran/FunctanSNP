% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotSNPreg.R
\name{plotSNPreg}
\alias{plotSNPreg}
\title{Plot coefficient functions and residuals from an SNPlm or SNPinter object}
\usage{
plotSNPreg(x, type)
}
\arguments{
\item{x}{an "SNPlm"/"SNPinter" object obtained by \code{\link{SNPlm}}/\code{\link{SNPinter}}.}

\item{type}{character, "beta" or "residuals" type, if "beta", plot the genetic effect function beta0(t) (and interaction items betak(t), if the entry x is an "SNPinter" object), if "residuals", show the residual plot.}
}
\value{
show the residual plot or coefficient functions.
}
\description{
Produces a residual diagram or a plot of the genetic effect function beta(t) (and interaction items betak(t), if the entry x is an "SNPinter" object) for a fitted "SNPlm" or "SNPinter" object.
}
\examples{
library(FunctanSNP)
n <- 300
m <- 30
simdata1 <- simData1(n, m, seed = 123)
SNPlmres <- SNPlm(y = simdata1$y, z = simdata1$z,
                  location = simdata1$location, X = simdata1$X,
                  type1 = "Bspline", type2 = "Bspline", nbasis1 = 5,
                  nbasis2 = 5, params1 = 4, params2 = 4,
                  intercept = FALSE, Plot = FALSE)
plotSNPreg(x = SNPlmres, type = "beta")
plotSNPreg(x = SNPlmres, type = "residuals")

simdata2 <- simData2(n, m, seed = 123)
lambda1 <- 0.05
lambda2 <- sqrt(3)*lambda1
eta <- 0
SNPinterres <- SNPinter(y = simdata2$y, z = simdata2$z,
                        location = simdata2$location, X = simdata2$X,
                        lambda1, lambda2, eta, type1 = "Bspline", nbasis1 = 5,
                        params1 = 4, Bsplines = 5, norder = 4, intercept = TRUE,
                        eps = 1e-2, maxstep = 1e2, Plot = FALSE)
plotSNPreg(x = SNPinterres, type = "beta")
plotSNPreg(x = SNPinterres, type = "residuals")

}
\seealso{
See Also as \code{\link{SNPlm}}, \code{\link{SNPinter}}.
}
