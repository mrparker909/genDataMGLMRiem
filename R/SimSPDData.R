#' @title genSPDdata
#' @description Generate random SPD data with covariates (1 covariate per SPD upper diagonal entry). Default covariate coefficient is beta_j=1. Control the signal to noise ratio via SNR. Noise is gaussian and uncorrelated. A dataframe containing N rows, and 1 column for each matrix upper diagonal entry, and 1 column for each corresponding covariate. Generates Y_dimsxdims, y1...y_(n_y), and x1...x_(n_y). The y subsript counts along rows (left to right) first, then moves down one row at a time (eg for dims=3 when the main diagonal is excluded, y1 = Y[1,2], y2=Y[1,3], y3=Y[2,3]).
#' @return If C=NULL (no confounds), returns a data.frame with 2*n_y columns, and N rows. n_y is dims*(dims/2+1/2), or dims*(dims/2-1/2) depending on whether the main diagonal is included. Otherwise, returns a list with elements Y, X, V, P, and yx, where Y is the array of observed SPD matrices, X is the data.frame of k confounds (kxN), V is the list of k SPD covariate matrices, P is the base point on the manifold, and yx is the data.frame of upper diagonal elements of Y (with each row corresponding to one observation Y[[i]]).
#' @param NUM     number of observations (NUM is the number of SPD matrices generated)
#' @param dims    dimension of SPD matrices to generate
#' @param SNR     Signal to Noise ratio (SNR=0.25 means the noise will be up to 4 times the magnitude of the signal).
#' @param includeDiagonal If T, will include the diagonal entries of the SPD matrix. Otherwise will only include the uppder diagonal elements.
#' @param beta    If C is NULL, this will be used to generate data under a element-wise multivariate normal model (y1...yp ~ X*Beta), NOT a Riemann Manifold model. Vector of coefficients for the generated covariates. Default is rep(1, times=dims). If beta is too short, will be right padded with zeros. If beta is too long, will be right-truncated.
#' @param C       If C is not NULL, this will be used to generate data under a Riemann Manifold model (Y ~ C*Vc + X*Vx), NOT a multivariate normal model. A kxN array of confounds (each row is an N-vector corresponding to a set of observations for one confound, each column corresponds to a k-vector of confound observations for the ith individual)
#' @param X       If X is not NULL, this will X will be the covariates used to generate data under a Riemann Manifold model (Y ~ C*Vc + X*Vx).
#' @param P       If C is not NULL, P is the base point on the manifold from which V displaces (Y ~ exp_m(P,CV)), default is P=Identity matrix (dimension dims.
#' @param maxDist If C is NULL, the maximum distance on the SPD manifold between response values. If C is not NULL, maximum distance on the SPD manifold for the generated SPD matrices for the covariate coefficients, V (distance measured from P).
#' @param minDist If C is NULL, the minimum distance on the SPD manifold between response values. If C is not NULL, maximum distance on the SPD manifold for the generated SPD matrices for the covariate coefficients, V (distance measured from P).
#' @param corr    If corr=T, then the SPD matrices generated will be correlation matrices.
#' @examples
#' set.seed(1234)
#' g <- genSPDdata(NUM = 3, dims = 3, SNR = 2, C=matrix(c(1,2,1),ncol=3))
#' g$Y[[1]]
#' MGLMRiem::expmap_spd(g$P,g$X[1]*g$V)
#'
#' X = matrix(rnorm(10), nrow=2)
#' M = matrix(rnorm(15), nrow=3)
#' C = rbind(X,M)
#'
#' dat = genSPDdata(NUM=5, dims=3, SNR=1, C=C, nX=2)
#'
#' @export
genSPDdata <- function(NUM=500, dims=5, maxDist = 1, minDist=0, SNR=1, includeDiagonal=F, beta=NULL, C=NULL, X=NULL, P=NULL, corr=F) {
  if(dims < 2) stop("genDat: dims must be at least 2")

  # setup data structure
  n_y = ifelse(includeDiagonal, dims*(dims/2+1/2), dims*(dims/2-1/2))
  y_df = data.frame(matrix(ncol=n_y,nrow=NUM))
  names(y_df) <- paste0("y",1:n_y)
  x_df = data.frame(matrix(ncol=n_y,nrow=NUM))
  names(x_df) <- paste0("x",1:n_y)

  Y <-list()
  if(is.null(C)) {
    # generate NUM random SPDs
    for(i in 1:NUM) {Y[[i]] = MGLMRiem::randspd_FAST(n = dims, maxDist = maxDist, minDist=minDist)}
    if(corr) {
      for(i in 1:NUM) {Y[[i]] = cov2cor(Y[[i]])}
    }

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
    e = matrix(rnorm(n = n_y*NUM, mean = 0, sd = 1/SNR),ncol=n_y)

    if(is.null(beta)) { beta = rep(1, times=n_y) }
    if(length(beta) > n_y) { beta = beta[1:n_y] }
    if(length(beta) < n_y) { beta= c(beta, rep(0, times=n_y))[1:n_y] }

    for(j in 1:n_y) {
      if(beta[j]==0) { x_df[,j] = -e[,j] } else {
        x_df[,j] = (y_df[,j] - e[,j])/beta[j]
      }
    }

    return(list(Y=Y,yx=cbind(y_df,x_df)))
  } else {
    if(is.null(P)) { P = diag(dims) }
    if(!is.matrix(X)) { X = matrix(X, nrow=1) }
    if(!is.matrix(C)) { C = matrix(C, nrow=1) }
    nX <- nrow(X) # number of covariates
    nC <- nrow(C) # number of confounds

    if(NUM !=ncol(C)) stop("C must be NULL, or number of columns of C must match sample size N.") # sample size

    spdDat = genSPD_SNR(d=dims, X=X, C=C, scale=1, SNR=SNR, bp=P, maxDist = maxDist, minDist=minDist)
    Y = spdDat$Y

    # extract y components
    i=0
    if(includeDiagonal) {
      for(row in 1:(dims)) {
        for(col in (row):dims) {
          i=i+1
          y_df[,i] = unlist(Y[row,col,])
        }
      }
    } else {
      for(row in 1:(dims-1)) {
        for(col in (row+1):dims) {
          i=i+1
          y_df[,i] = unlist(Y[row,col,])
        }
      }
    }

    # Yret = array(0, dim=c(dims,dims,NUM))
    # for(i in 1:NUM) {
    #   Yret[,,i] = Y[[i]]
    # }

    return(list(Y=Y, X=X, C=C, Vx=spdDat$Vx, Vc=spdDat$Vc, P=P, y_upper=cbind(y_df), Ygroundtruth=spdDat$Ygroundtruth))
  }
}

