
#' @title genSPD_SNR
#' @description Generates SPD matrices with noise from covariates X, confound C, with signal to noise ratio SNR.
#' @param d Dimension of SPD matrices to generate.
#' @param X A x by n matrix, x rows indicate an observed variable, n columns an observed individual.
#' @param C A c by n matrix, c rows indicate an observed variable, n columns an observed individual.
#' @param scale scale factor to scale X and C (reduces computational issues from overly large SPD matrices).
#' @param SNR Signal to Noise ratio to be used in generating the noise. SNR must be larger 0.
#' @param bp Base point on SPD manifold (default is Identity).
#' @param maxDist maximum distance from bp for coefficient vectors.
#' @param minDist minimum distance from bp for coefficient vectors.
#' @param maximumSPDValue If not NULL, limits the maximum value that any entry in the SPD matrix can obtain.
#' @examples
#' @export
genSPD_SNR <- function(d, X, C, scale=1, SNR=1, bp = NULL, maxDist=1.25, minDist=0.75, maximumSPDValue=NULL) {
  n = sizeR(X,2)

  if(is.null(bp)) {
    bp = diag(d)
  }

  # x = covariate obs, c = confound obs
  X = scale*X
  C = scale*C

  nX = nrow(X) # number of covariates
  nC = nrow(C) # number of confounds

  # Vx and Vc are coefficient matrices
  Px = MGLMRiem::randspd_FAST(n = d, NUM = nX+1, maxDist = maxDist, minDist = minDist, maximumSPDValue=maximumSPDValue)
  Pc = MGLMRiem::randspd_FAST(n = d, NUM = nC+1, maxDist = maxDist, minDist = minDist, maximumSPDValue=maximumSPDValue)

  for(i in 1:nX) {
    Vx[,,i] = MGLMRiem::logmap_spd(Px[,,i],Px[,,i+1])
  }
  for(i in 1:nC) {
    Vc[,,i] = MGLMRiem::logmap_spd(Pc[,,i],Pc[,,i+1])
  }

  Vx = MGLMRiem::aug3(Vx)
  Vc = MGLMRiem::aug3(Vc)

  # sum_k x_k * Vx_k
  xVx = array(0, dim = c(d,d,n))
  for(i in 1:n) {
    for(j in 1:nX) {
      xVx[,,i] = xVx[,,i] + X[j,i] * Vx[,,j]
    }
  }
  # sum_k c_k * Vc_k
  cVc = array(0, dim = c(d,d,n))
  for(i in 1:n) {
    for(j in 1:nC) {
      cVc[,,i] = cVc[,,i]+ C[j,i] * Vc[,,j]
    }
  }


  # Signal (covariates x)
  S = array(0, dim=c(d,d,n))
  for(i in 1:n) {
    S[,,i] = MGLMRiem::expmap_spd(bp,xVx[,,i])
  }
  Sbar = MGLMRiem::karcher_mean_spd(S, niter=100)
  sig_S = sqrt(MGLMRiem::gsqerr_spd(Sbar, S))

  # Noise (confounds c)
  N = array(0, dim=c(d,d,n))
  for(i in 1:n) {
    N[,,i] = MGLMRiem::expmap_spd(bp,cVc[,,i])
  }
  Nbar = MGLMRiem::karcher_mean_spd(N, niter=100)
  sig_N = sqrt(MGLMRiem::gsqerr_spd(Nbar, N))

  # scale cVc:
  cUc = array(0, dim=c(d,d,n))
  for(i in 1:n) {
    cUc[,,i] = cVc[,,i] * ( sig_S/ sig_N ) * ( nC / (SNR * nX) )
  }

  N = array(0, dim=c(d,d,n))
  for(i in 1:n) {
    N[,,i] = MGLMRiem::expmap_spd(bp,cUc[,,i])
  }

  if(any(unlist(N) > 10^5)) { warning(paste0("In genSPD_SNR(): N contains very large element(s):",unlist(N)[which(unlist(N) > 10^5)], "\n")) }

  Y = array(0, dim=c(d,d,n))
  for(i in 1:n) {
    Y[,,i] = MGLMRiem::expmap_spd(bp, xVx[,,i]+cUc[,,i])
  }

  return(list(Y = Y, Ygroundtruth=S, Vx=Vx, Vc=Vc, bp=bp, SNR=SNR, scale=scale,X=X,C=C))
}

