##' @title Generate the sequence (genotypes) data containing only 0, 1 and 2
##'
##' @description This function provides a method for generating sequence data containing only 0, 1 and 2, which can be used to simulate the generation of sequence genotypes.
##'
##' @param n an interger variable specifying the number of samples to be generated.
##' @param m an interger variable specifying the sequence length of each sample.
##' @param seed an integer variable specifying the random seed used for random sequence generation.
##' @param d.ratio a numeric variable between 0 and 1 indicating the deletion ratio of sample sequences, default value is 0.
##'
##' @return An "simX" object that contains the list of the following items.
##' \itemize{
##' \item{location:}{ a numeric vector defining the sampling sites of the sequence data.}
##' \item{X:}{ a matrix with n rows and m columns representing the sequence data.}
##' }
##' @seealso See Also as \code{\link{plotRawdata}}, \code{\link{SNPgvf}}.
##'
##' @export
##' @examples
##' library(FunctanSNP)
##' n <- 2
##' m <- 50
##' simdata1 <- simX(n, m, seed = 1, d.ratio = 0)
##' simdata2 <- simX(n, m, seed = 1, d.ratio = 0.3)
##' plotRawdata(location = simdata1$location, X = simdata1$X)
##' plotRawdata(location = simdata2$location, X = simdata2$X)
##'

simX <- function(n, m, seed = 1, d.ratio = 0){
  norder <- 4
  nknots <- 20
  t <- seq(1e-2, 1, length = m)
  breaks <- seq(0, 1, length = nknots)
  basismat <- bsplineS(x = t, breaks = breaks, norder = norder, returnMatrix = FALSE)
  nbasisX <- ncol(basismat)
  set.seed(seed = seed + 123); coef <- mvrnorm(n = n, mu = rep(0.5, nbasisX), Sigma = diag(x = 1, nrow = nbasisX, ncol = nbasisX))
  Rawfvalue <- coef %*% t(basismat)
  #fmin <- min(Rawfvalue)
  #fmax <- max(Rawfvalue)
  #fvalue <- (Rawfvalue - fmin) * 3.5 / ( fmax - fmin) - 1awfvalue
  fvalue <- Rawfvalue
  funcX <- function(l){
    x <- fvalue[l, ]
    len <- length(x)
    diffmat <- matrix((rep(x, times = 3) - rep(c(0,1,2), each = len))^2, ncol = len, nrow = 3, byrow = TRUE)
    value <- apply(diffmat, 2, which.min) - 1
    return(value)
  }
  dataX <- t(mapply(funcX, 1:n))
  if(d.ratio == 0){ SNPdata <- dataX }
  if(d.ratio != 0 & d.ratio <= 1){
    func <- function(i){
      set.seed(seed = seed + i); index <- sample(1:m, size = round(m*d.ratio))
      dataX[i, index] <- NA
      return(dataX[i, ])
    }
    SNPdata <- t(mapply(func, 1:n))
  }
  object <- list(location = t, X = SNPdata)
  class(object) <- "simX"
  return(object)
}

