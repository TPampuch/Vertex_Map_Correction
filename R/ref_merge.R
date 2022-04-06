#' Merge reference points
#'
#' This function is not recommended.
#' It removes the doubled reference points (usually trees measured from different positions) and calculates the mean of the points.
#' This function should only be used after using the VMC function and only for purpose of plotting a map without overlapping points.
#'
#' @param data data.frame including coordinates and Tree-IDs
#' @param X_Col name for X-coordinates column
#' @param Y_Col name for Y-coordinates column
#' @param Tree_ID name for tree ID column
#' @return data.frame
#' @export

ref_merge <- function(data, X_col = "X.m.", Y_col = "Y.m.", Tree_ID = "TreeID"){
  colnames(data)[which(colnames(data) == X_col)] <- "X.m."
  colnames(data)[which(colnames(data) == Y_col)] <- "Y.m."
  colnames(data)[which(colnames(data) == Tree_ID)] <- "TreeID"
  
  overlaps <- names(which(table(data$TreeID) > 1))
  #overlaps
  
  for(i in 1:length(overlaps)){
    means = colMeans(data[which(data$TreeID == overlaps[i]),c("X.m.","Y.m.")])
    new = data[which(data$TreeID == overlaps[i]),][1,]
    new$X.m. = means[1]
    new$Y.m. = means[2]
    del <- which(data$TreeID == overlaps[i])
    keep <- del[1]
    del <- del[2:length(del)]
    data[keep,] = new  
    data = data[-del,]
  }
  colnames(data)[which(colnames(data) == "X.m.")] <- X_col
  colnames(data)[which(colnames(data) == "Y.m.")] <- Y_col
  colnames(data)[which(colnames(data) == "TreeID")] <- Tree_ID
  
  return(data)
}