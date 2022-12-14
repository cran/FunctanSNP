% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SNPcvinter.R
\name{SNPcvinter}
\alias{SNPcvinter}
\title{Cross-validation for \code{\link{SNPinter}}}
\usage{
SNPcvinter(y, z, location, X, K, lambda1, lambda2 = NULL, eta = 0,
  type1, nbasis1, params1, Bsplines = 20, norder = 4,
  intercept = FALSE, eps = 1e-05, maxstep = 1e+05, Plot = FALSE)
}
\arguments{
\item{y}{a numeric vector defining the response variables.}

\item{z}{a matrix defining the scalar covariates, with the number of rows equal to the number of samples.}

\item{location}{a numeric vector defining the sampling sites of the sequence data.}

\item{X}{a matrix specifying the sequence (genotypes) data, with the number of rows equal to the number of samples.}

\item{K}{an integer specifying the number of cross-validation folds, default is 5.}

\item{lambda1}{a numeric vector specifying the sparsity penalty parameter to be determined.}

\item{lambda2}{a numeric vector specifying the group sparsity penalty parameter to be determined.}

\item{eta}{a numeric vector specifying the penalty parameter for smoothing analysis.}

\item{type1}{a character specifying the type of the basis functions that constitutes the genetic variation function. The options are "Bspline", "Exponential", "Fourier", "Monomial", and "Power".}

\item{nbasis1}{an integer specifying the number of basis functions that constitutes the genetic variation function.}

\item{params1}{in addition to rangeval1 (a vector of length 2 giving the lower and upper limits of the range of permissible values for the genetic variation function) and nbasis1, all bases have one or two parameters unique to that basis type or shared with one other;
\itemize{
\item{bspline:}{ Argument norder = the order of the spline, which is one more than the degree of the polynomials used. This defaults to 4, which gives cubic splines.}
\item{exponential:}{ Argument ratevec. In fda_2.0.2, this defaulted to 1. In fda_2.0.3, it will default to 0:1.}
\item{fourier:}{ Argument period defaults to diff(rangeval).}
\item{monomial/power:}{ Argument exponents. Default = 0:(nbasis-1). For monomial bases, exponents must be distinct nonnegative integers. For power bases, they must be distinct real numbers.}
}}

\item{Bsplines}{an integer specifying the number of basis functions that constitutes the genetic effect function.}

\item{norder}{an integer specifying the order of bsplines that constitutes the genetic effect function, which is one higher than their degree. The default of 4 gives cubic splines.}

\item{intercept}{should intercept(s) be fitted (TRUE) or set to zero (default = FALSE).}

\item{eps}{a numeric variable specifying the threshold at which the algorithm terminates, default is 1e-5.}

\item{maxstep}{a numeric variable specifying the maximum iteration steps, default is 1e5.}

\item{Plot}{should the estimated genetic effect function beta0(t) and interaction items betak(t) be plotted (TRUE) or not (default = FALSE).}
}
\value{
An "SNPcvinter" object that contains the list of the following items.
\itemize{
\item{lambda1Select:}{ a numeric value of the sparsity penalty parameter selected by cross validation.}
\item{lambda2Select:}{ a numeric value of the group sparsity penalty parameter selected by cross validation.}
\item{etaSelect:}{ a numeric value of the smoothing parameter selected by cross validation.}
\item{alpha:}{ estimated intercept value..}
\item{gamma:}{ estimated coefficients of the scalar covariates.}
\item{b:}{ estimated coefficients of the chosen basis functions for the genetic effect function beta0(t) and interaction items betak(t).}
\item{betat:}{ an "fd" object, representing the estimated genetic effect function beta(t) and interaction items betak(t).}
\item{residuals:}{ the residuals, that is response minus fitted values.}
\item{fitted.values: }{ the fitted mean values.}
\item{lambda1:}{ a numeric vector specifying the sparsity penalty parameter for cross validation.}
\item{lambda2:}{ a numeric vector specifying the group sparsity penalty parameter for cross validation.}
\item{eta:}{ a numeric vector specifying the smoothing parameter for cross validation.}
\item{CVerror:}{ a numeric vector, containing the mean square errors on testing set during cross validation.}
}
}
\description{
Performs K-fold cross validation for the revised partially functional interaction regression analysis over a grid of values for the regularization parameter lambda1 and lambda2.
}
\examples{
library(FunctanSNP)
n <- 300
m <- 30
simdata2 <- simData2(n, m, seed = 123)
y <- simdata2$y
z <- simdata2$z
location <- simdata2$location
X <- simdata2$X
lambda1 <- c(0.01, 0.05, 0.1)
lambda2 <- sqrt(3)*lambda1
SNPcvinterres <- SNPcvinter(y, z, location, X, K = 3, lambda1, lambda2, eta = 0,
                            type1 = "Bspline", nbasis1 = 5, params1 = 4, Bsplines = 5,
                            norder = 4, intercept = TRUE, eps = 1e-2, maxstep = 1e2, Plot = TRUE)
SNPcvinterres$lambda1Select
SNPcvinterres$lambda2Select
SNPcvinterres$alpha
SNPcvinterres$gamma
SNPcvinterres$b

}
\seealso{
See Also as \code{\link{simData2}}, \code{\link{SNPinter}}.
}
