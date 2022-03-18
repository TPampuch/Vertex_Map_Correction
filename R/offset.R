#' Mean offset
#'
#' This function calculates the mean offset between pairs of bases
#' based on reference trees
#' 
#' @param data data.frame including vertex output
#' @param X_Col name for X-coordinates column
#' @param Y_Col name for X-coordinates column
#' @param Base_ID name for base ID column
#' @param Tree_ID name for tree ID column
#' 
#' @return data.frame with corrected X and Y values
#' @export

VM.mean.offset <- function(data, X_col = "X.m.", Y_col = "Y.m.", Base_ID = "Base", Tree_ID = "ID"){
    
    # adjusting colnames based on provided information
    colnames(data)[which(colnames(data) == Base_ID)] <- "Base"
    colnames(data)[which(colnames(data) == Tree_ID)] <- "TreeID"
    colnames(data)[which(colnames(data) == X_col)] <- "X.m."
    colnames(data)[which(colnames(data) == Y_col)] <- "Y.m."
    
  overlaps <- names(which(table(data$TreeID) > 1))
  Bases <- max(na.omit(data$Base))
  
  # calculate initial mean square
  ctrl = NA
  i = 1
  means = data.frame(X = NA, Y = NA, Bases = NA)
  for(a in 1:Bases){ # loop 1 of all Bases
    for(b in 1:Bases){ # loop 2 of all Bases
      if(b == a) next # next iteration if a and b are the same Base
      ctrl = paste(sort(c(a,b))[1], "&", sort(c(a,b))[2]) # control line to avoid doubling base-pairs
      if(ctrl %in% means$Bases) next # next iteration if the combination was already calculated
      
      adata <- data[which(data$TreeID %in% overlaps & data$Base == a),]
      bdata <- data[which(data$TreeID %in% overlaps & data$Base == b),]
      
      adata <- adata[which(adata$TreeID %in% bdata$TreeID),]
      bdata <- bdata[which(bdata$TreeID %in% adata$TreeID),]
      
      if(nrow(adata) < 1 | nrow(bdata) < 1) next # skip to next loop iteration if there is no row in a or b
      
      Base12 <- rbind(adata,bdata)
      
      names = unique(Base12$TreeID)
      if(length(names) < 2) next # skip to next loop iteration if the number of "reference tree" is below the provided amount (should be at least (default) 2)
      ms = data.frame(X_offset = NA, Y_offset = NA)
      for(c in 1:length(names)){
        Basems = Base12[which(Base12$TreeID == names[c]),]
        ms[c,"X_offset"] = (Basems$X.m.[1] - Basems$X.m.[2])
        ms[c,"Y_offset"] = (Basems$Y.m.[1] - Basems$Y.m.[2])
      }
      
      head(ms)
      means[i,] <- c(mean(ms$X), mean(ms$Y), as.character(paste(sort(c(a,b))[1], "&", sort(c(a,b))[2])))
      i = i + 1
    }
  }
  return(means)
}
