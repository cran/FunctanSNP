##' @title Plot the genetic variant function
##'
##' @description Plot the estimated genetic variant function obtained by SNPgvf or a given functional data object(s).
##'
##' @param x functional data object(s) to be plotted, such as the estimated genetic variant function obtained by SNPgvf.
##' @param y sequence of points at which to evaluate the functions 'x' and plot on the horizontal axis. Defaults to seq(rangex[1], rangex[2], length = nx).
##' @param Lfdobj either a nonnegative integer or a linear differential operator object. If present, the derivative or the value of applying the operator is plotted rather than the functions themselves.
##' @param href a logical variable: If TRUE, add a horizontal reference line at 0.
##' @param titles a vector of strings for identifying curves
##' @param xlim a vector of length 2 containing axis limits for the horizontal axis.
##' @param ylim a vector of length 2 containing axis limits for the vertical axis.
##' @param xlab a label for the horizontal axis.
##' @param ylab a label for the vertical axis.
##' @param ask a logical value: If TRUE, each curve is shown separately, and the plot advances with a mouse click
##' @param nx the number of points to use to define the plot. The default is usually enough, but for a highly variable function more may be required.
##' @param axes Either a logical or a list or NULL.
##' \itemize{
##' \item{logical}{ whether axes should be drawn on the plot}
##' \item{list}{ a list used to create custom axes used to create axes via x$axes[[1]] and x$axes[-1]. The primary example of this uses list("axesIntervals", ...)}
##' }
##' @return plot the estimated genetic variant function.
##'
##' @seealso See Also as \code{\link{SNPgvf}}.
##'
##' @importFrom fda plot.fd
##' @export
##' @examples
##' library(FunctanSNP)
##' n <- 20
##' m <- 50
##' simdata <- simX(n, m, seed = 1, d.ratio = 0)
##' X <- simdata$X
##' location <- simdata$location
##' SNPgvfres <- SNPgvf(location, X, type = "Bspline", nbasis = 5, params = 4, Plot = FALSE)
##' plotGVF(SNPgvfres)
##'

plotGVF <- function(x, y = NULL, Lfdobj = 0, href = TRUE, titles = NULL,
                    xlim = NULL, ylim = NULL, xlab = NULL,
                    ylab = NULL, ask = FALSE, nx = NULL, axes = NULL){
  gvf <- x
  if (!inherits(gvf, "fd")) { stop("x should be of fd type in fda package.") }

  if(is.null(y) == TRUE){plot.fd(x = gvf, Lfdobj = Lfdobj, href = href, titles = titles,
                                 xlim = xlim, ylim = ylim, xlab = xlab,
                                 ylab = ylab, ask = ask, nx = nx, axes = axes)}
  if(is.null(y) == FALSE){plot.fd(x = gvf, y = y, Lfdobj = Lfdobj, href = href, titles = titles,
                                  xlim = xlim, ylim = ylim, xlab = xlab,
                                  ylab = ylab, ask = ask, nx = nx, axes = axes)}
}
