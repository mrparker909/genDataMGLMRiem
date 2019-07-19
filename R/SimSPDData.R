#' @title genSPDdata
#' @description Generate random SPD data with covariates (1 covariate per SPD upper diagonal entry). Default covariate coefficient is beta_j=1. Control the signal to noise ratio via SNR. Noise is gaussian and uncorrelated. A dataframe containing N rows, and 1 column for each matrix upper diagonal entry, and 1 column for each corresponding covariate. Generates Y_dimsxdims, y1...y_(n_y), and x1...x_(n_y). Returns a data.frame with 2*n_y columns, and N rows. n_y is dims*(dims/2+1/2), or dims*(dims/2-1/2) depending on whether the main diagonal is included. The y subsript counts along rows (left to right) first, then moves down one row at a time (eg for dims=3 when the main diagonal is excluded, y1 = Y[1,2], y2=Y[1,3], y3=Y[2,3]).
#' @param N       number of observations (N is the number of SPD matrices generated)
#' @param dims    dimension of SPD matrices to generate
#' @param maxDist maximum distance on the SPD manifold for the generated SPD matrices (measured from the identity matrix)
#' @param SNR     Signal to Noise ratio (noise has standard deviation = 1/SNR)
#' @param includeDiagonal If T, will include the diagonal entries of the SPD matrix. Otherwise will only include the uppder diagonal elements.
#' @param beta    vector of coefficients for the generated covariates. Default is rep(1, times=dims). If beta is too short, will be right padded with zeros. If beta is too long, will be right-truncated.
#' @export
genSPDdata <- function(N=500, dims=5, maxDist = 1, SNR=1, includeDiagonal=F, beta=NULL) {
  if(dims < 2) stop("genDat: dims must be at least 2")

  # generate N random SPDs
  Y <-list()
  for(i in 1:N) {Y[[i]] = MGLMRiem::randspd_FAST(n = dims, maxDist = maxDist)}

  # setup data structure
  n_y = ifelse(includeDiagonal, dims*(dims/2+1/2), dims*(dims/2-1/2))
  y_df = data.frame(matrix(ncol=n_y,nrow=N))
  names(y_df) <- paste0("y",1:n_y)
  x_df = data.frame(matrix(ncol=n_y,nrow=N))
  names(x_df) <- paste0("x",1:n_y)

  # extract y components
  i=0
  if(includeDiagonal) {
    for(row in 1:(dims)) {
      for(col in (row):dims) {
        i=i+1
        y_df[,i] = unlist(lapply(Y, FUN = function(Y) {Y[row,col]}))
      }
    }
  } else {
    for(row in 1:(dims-1)) {
      for(col in (row+1):dims) {
        i=i+1
        y_df[,i] = unlist(lapply(Y, FUN = function(Y) {Y[row,col]}))
      }
    }
  }

  # generate covariates
  e = matrix(rnorm(n = n_y*N, mean = 0, sd = 1/SNR),ncol=n_y)

  if(is.null(beta)) { beta = rep(1, times=n_y) }
  if(length(beta) > n_y) { beta = beta[1:n_y] }
  if(length(beta) < n_y) { beta= c(beta, rep(0, times=n_y))[1:n_y] }

  for(j in 1:n_y) {
    if(beta[j]==0) { x_df[,j] = -e[,j] } else {
      x_df[,j] = (y_df[,j] - e[,j])/beta[j]
    }
  }

  return(list(Y=Y,yx=cbind(y_df,x_df)))
}

