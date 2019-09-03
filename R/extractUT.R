#' @title extractUT
#' @description Extract the upper triangular elements of an array of SPDs, return data.frame, each row belongs to one observation. Traverses across columns, starting with row 1, moving right from the diagonal.
#' @param Y an array of dims x dims x N SPD matrices.
#' @param dims the dimension of the SPD matrices.
#' @param includeDiagonal if T, will return the diagonal elements as well as the upper triangular elements.
#' @export
extractUT <- function(Y, dims, includeDiagonal=F) {
  Y = MGLMRiem::aug3(Y)
  N = dim(Y)[3]

  n_y = ifelse(includeDiagonal, dims*(dims/2+1/2), dims*(dims/2-1/2))
  y_df = data.frame(matrix(ncol=n_y,nrow=N))
  names(y_df) <- paste0("y",1:n_y)

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
  return(y_df)
}

