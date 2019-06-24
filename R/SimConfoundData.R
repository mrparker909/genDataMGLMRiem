#' @title SimConfoundData
#' @description SimConfoundData simulates a vector of confounds from a given distribution
#' @param n the sample size to use
#' @param D the distribution function to draw from (eg D=rbinom), must take n (sample size) as first parameter
#' @param ... parameters to pass to the distribution function (eg p=0.25 for D=rbinom)
#' @export
#' @examples
#' C1 <- SimConfoundData(n=10, D=rbinom, size=1, prob=0.25)
SimConfoundData <- function(n,D,...) {
  C <- D(n, ...)
  dim(C) <- c(1,length(C))
  return(C)
}

#' @title SimSNPData
#' @description SimSNPData simulates SNP differences. Wrapper for SimConfoundData
#' @param n the sample size
#' @param d the maximum number of SNP differences (eg if d=3, there will be four SNP categoris: 0, 1, 2, or 3 SNP differences)
#' @param p probability weights for each number of SNP differences
#' @examples
#' C1 <- SimSNPData(n=20, d=3, p=c(7,2,1))
#' @export
SimSNPData <- function(n, d, p=NULL) {
  if(is.null(p)) { p <- rep(1, d+1)}

  C1 <- SimConfoundData(n = n, D = function(n,d,p, ...) {
          sample(x=0:d, size=n, replace=T, prob=p, ...)
        }, d=d, p=p)
  return(C1)
}
