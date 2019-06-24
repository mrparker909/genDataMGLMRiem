#' @title SimDTData
#' @description Simulate Diffusion Tensor data with correlated confounds. The overall equation of interest is: Yobs = Exp(Exp(P,VC),eps). Returns a list containing: P (true base point for Y), C (confounds), V (true coefficients for C), Y (ground truth DT's, Y=VC), Yobs (observed, noisy DT's)
#' @param C a kxN array of confounds (each row is an N-vector corresponding to aset of observations for one confound, each column to a k-vector of confound observations for the ith individual)
#' @param P a SPD base point (identity matrix by default)
#' @param V a 3x3xN array of coefficients for C (N-stack of symmetric 3x3 matrices, random by default)
#' @param noise noise to add to generated Diffusion Tensor observations
#' @param Yvar proportional to the variance of base point for Ytrue, the ground truth diffusion tensors
#' @export
SimDTData <- function(C,P=NULL,V=NULL,noise=1,Yvar=1) {
  k <- nrow(C) # number of confounds
  N <- ncol(C) # sample size

  if(is.null(P)) { # identity by default
    P <- diag(rep(1,3))
  }

  if(any(dim(P) != c(3,3))) {
    stop("base point P should be a 3x3 matrix or array")
  }

  # generate N Diffusion Tensors (with dependence on confounds C)
  Yp = array(0, dim=c(3,3,k+1))

  Yp[,,1] = P # base point is Yp[,,1]=P
  for(i in 1:k) {
    Yp[,,i+1] = MGLMRiem::randspd(3,Yvar,3*Yvar)
  }

  if(is.null(V)) { # generate random V
    # Tangent vectors, geodesic bases.
    V = array(0, dim=c(3,3,k))
    for(j in 1:k) {
      V[,,j] = MGLMRiem::logmap_spd(Yp[,,1], Yp[,,j+1])
    }
  }

  if(MGLMRiem::sizeR(V,3) != k) {
    stop("V should be 3x3xN (incorrect sample size N)")
  }

  if(any(MGLMRiem::sizeR(V,c(1,2)) != c(3,3))) {
    stop("V should be 3x3xN (incorrect tensor size, should be 3x3)")
  }

  # Generate Ground Truth Data
  Y2 = array(0, dim=c(3,3, N)) #zeros(3,3,size(X,2));
  for(i in 1:N) {
    Vtmp = array(0, dim=c(3,3,1))  #zeros(3,3,1)
    for(j in 1:k) {
      Vtmp = Vtmp + MGLMRiem::aug3(V[,,j]*C[j,i])
    }
    Y2[,,i] = MGLMRiem::expmap_spd(P=Yp[,,1],X=Vtmp)
  }

  # Generate Y Samples
  Ysample = array(0, dim=c(3,3, N)) #zeros(3,3,size(Y,3)*npairs);
  for(j in 1:N) {
    Ysample[,,j] = MGLMRiem::addnoise_spd(Y2[,,j],noise)
  }

  if(MGLMRiem::isspd_mxstack(Ysample) != 1) { stop("spd stack contains non-spd matrix") }
  Ysample

  return(list(Yobs=Ysample, Xobs=C, Ptrue=P, Vtrue=V, Ytrue=Y2))
}

