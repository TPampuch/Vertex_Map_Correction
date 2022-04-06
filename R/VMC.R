#' Vertex Map Correction
#'
#' This function tries to minimize the offset introduced when mapping an area using a vertex laser geo 360.
#' The function checks for reference trees that where measured from several positions (bases) and simply
#' loops through base-pairs with overlapping information. The offset is reduced by moving all data-points
#' measured from one base towards the other base. Since the data cannot fit perfectly it repeats the step with
#' different base-pairs as long as the sum of the squares of all offsets is further reduced.
#'
#' @param data data.frame including coordinates, Tree-IDs and "Base-number"
#' @param TreeDBH name for DBH column (only needed for plotting)
#' @param X_Col name for X-coordinates column
#' @param Y_Col name for Y-coordinates column
#' @param Alt_Col name for Altitude column (only needed for Z correction)
#' @param Base_ID name for base ID column
#' @param Tree_ID name for tree ID column
#' @param limit number of iterations without change before loop stops (default = 10)
#' @param ref_num minimum number of overlapping trees (default = 2)
#' @param plot logical, if true plots two maps for comparison (requires ggplot2 and ggrepel)
#' @param z_corr logical, if true correction of z-column will be included
#' @return data.frame with corrected X and Y values
#' @export

VMC <- function(data, Tree_DBH = "DBH", X_col = "X.m.", Y_col = "Y.m.", Alt_col = "ALTITUDE", Base_ID = "Base", Tree_ID = "TreeID", limit = 10, ref_num = 2, plot = TRUE, z_corr = FALSE){
  
  ####
  # the data should be loaded as a data.frame with the typical Vertex Laser-Geo information
  # ! important ! #
  # - add ID to the data.frame that indicates individual trees (the IDs are used to find overlapping trees)
  # - add a column with the "base", base should be an integer starting from one (also used to find overlapping trees)
  # - make sure to remove the "settings" rows from the original vertex csv table
  # - Tree IDs should be unique (except for reference points of course)
  # - when there is no improvement in the mean square for a given number of times (variable "limit") in a row the loop stops
  
  # adjusting colnames based on provided information
  if(plot == TRUE){colnames(data)[which(colnames(data) == Tree_DBH)] <- "TreeDBH"}
  if(z_corr == TRUE){colnames(data)[which(colnames(data) == Alt_col)] <- "ALTITUDE"}
  colnames(data)[which(colnames(data) == Base_ID)] <- "Base"
  colnames(data)[which(colnames(data) == Tree_ID)] <- "TreeID"
  colnames(data)[which(colnames(data) == X_col)] <- "X.m."
  colnames(data)[which(colnames(data) == Y_col)] <- "Y.m."

  if(plot == TRUE){
    # create "before" plot of initial data
    p1 = ggplot2::ggplot(data, ggplot2::aes(x = Y.m., y = X.m.))+
      ggplot2::geom_point(ggplot2::aes(color = as.factor(Base), size = TreeDBH))+
      ggplot2::theme_light()+
      ggrepel::geom_text_repel(ggplot2::aes(label=TreeID), cex = 2)+
      ggplot2::labs(x = "X", y = "Y", title = "Before")
  }
  
  overlaps <- names(which(table(data$TreeID) > 1))
  Bases <- max(na.omit(data$Base))
  lms <- data.frame(sum = NA, iteration = NA)
  
  # calculate initial mean square
  ctrl = NA
  i = 1
  if(z_corr == TRUE){
    means = data.frame(X = NA, Y = NA, Z = NA, Base_combo = NA)
  }else{
    means = data.frame(X = NA, Y = NA, Base_combo = NA)
  }
  for(a in 1:Bases){ # loop 1 of all Bases
    for(b in 1:Bases){ # loop 2 of all Bases
      if(b == a) next # next iteration if a and b are the same Base
      ctrl = paste(sort(c(a,b))[1], "&", sort(c(a,b))[2]) # control line to avoid doubling base-pairs
      if(ctrl %in% means$Base_combo) next # next iteration if the combination was already calculated
      
      adata <- data[which(data$TreeID %in% overlaps & data$Base == a),]
      bdata <- data[which(data$TreeID %in% overlaps & data$Base == b),]
      
      adata <- adata[which(adata$TreeID %in% bdata$TreeID),]
      bdata <- bdata[which(bdata$TreeID %in% adata$TreeID),]
      
      if(nrow(adata) < 1 | nrow(bdata) < 1) next # skip to next loop iteration if there is no row in a or b
      
      Base12 <- rbind(adata,bdata)
      
      names = unique(Base12$TreeID)
      if(length(names) < ref_num) next # skip to next loop iteration if the number of "reference tree" is below the provided amount (should be at least (default) 2)
      if(z_corr == TRUE){
        ms = data.frame(X = NA, Y = NA, Z = NA) 
      }else{
        ms = data.frame(X = NA, Y = NA)
      }
      for(c in 1:length(names)){
        if(z_corr == TRUE){
          Basems = Base12[which(Base12$TreeID == names[c]),]
          ms[c,"X"] = (Basems$X.m.[1] - Basems$X.m.[2])
          ms[c,"Y"] = (Basems$Y.m.[1] - Basems$Y.m.[2])
          ms[c,"Z"] = (Basems$ALTITUDE[1] - Basems$ALTITUDE[2])
        }else{
          Basems = Base12[which(Base12$TreeID == names[c]),]
          ms[c,"X"] = (Basems$X.m.[1] - Basems$X.m.[2])
          ms[c,"Y"] = (Basems$Y.m.[1] - Basems$Y.m.[2])
        }
      }
      if(z_corr == TRUE){
        means[i,] <- c(mean(ms$X), mean(ms$Y), mean(ms$Z), as.character(paste(sort(c(a,b))[1], "&", sort(c(a,b))[2])))
      }else{
        means[i,] <- c(mean(ms$X), mean(ms$Y), as.character(paste(sort(c(a,b))[1], "&", sort(c(a,b))[2])))
      }
      i = i + 1
    }
  }
  
  k = 1
  if(z_corr == TRUE){
    lms$sum[k] = sum(as.numeric(means$X)^2 + as.numeric(means$Y)^2 + as.numeric(means$Z)^2) # the sum of all mean offsets (squared to remove negatives)
  }else{
    lms$sum[k] = sum(as.numeric(means$X)^2 + as.numeric(means$Y)^2) 
  }
  lms$iteration[k] = k
  #lms$sum[1]
  k = k + 1
  
  combos = means$Base_combo # now select possible combos and loop through
  death = 1
  repeat{
    for(l in 1:length(combos)){
      b1C = as.numeric(strsplit(combos[l], split = " ")[[1]][1])
      b2C = as.numeric(strsplit(combos[l], split = " ")[[1]][3])
      
      # merging overlapping points of 2 Bases
      adataC <- data[which(data$TreeID %in% overlaps & data$Base == b1C),]
      bdataC <- data[which(data$TreeID %in% overlaps & data$Base == b2C),]
      
      adataC <- adataC[which(adataC$TreeID %in% bdataC$TreeID),]
      bdataC <- bdataC[which(bdataC$TreeID %in% adataC$TreeID),]
      
      Base12C <- rbind(adataC,bdataC)
      Base12C
      
      # calc difference for each tree
      names = unique(Base12C$TreeID)

      if(z_corr == TRUE){
        msC = data.frame(X = NA, Y = NA, Z = NA)
      }else{
        msC = data.frame(X = NA, Y = NA)
      }
      for(a in 1:length(names)){
        if(z_corr == TRUE){
          BasemsC = Base12C[which(Base12C$TreeID == names[a]),]
          msC[a,"X"] = (BasemsC$X.m.[1] - BasemsC$X.m.[2])
          msC[a,"Y"] = (BasemsC$Y.m.[1] - BasemsC$Y.m.[2])
          msC[a,"Z"] = (BasemsC$ALTITUDE[1] - BasemsC$ALTITUDE[2])
        }else{
          BasemsC = Base12C[which(Base12C$TreeID == names[a]),]
          msC[a,"X"] = (BasemsC$X.m.[1] - BasemsC$X.m.[2])
          msC[a,"Y"] = (BasemsC$Y.m.[1] - BasemsC$Y.m.[2])
        } 
      }
      
      if(z_corr == TRUE){
        msdataC <- data.frame(X = mean(msC$X), Y = mean(msC$Y), Z = mean(msC$Z)) # calculate mean offset
        
        data$X.m.[which(data$Base == b2C)] = data$X.m.[which(data$Base == b2C)]+msdataC$X # correct entire dataset
        data$Y.m.[which(data$Base == b2C)] = data$Y.m.[which(data$Base == b2C)]+msdataC$Y
        data$ALTITUDE[which(data$Base == b2C)] = data$ALTITUDE[which(data$Base == b2C)]+msdataC$Z
      }else{
        msdataC <- data.frame(X = mean(msC$X), Y = mean(msC$Y)) # calculate mean offset
        
        data$X.m.[which(data$Base == b2C)] = data$X.m.[which(data$Base == b2C)]+msdataC$X # correct entire dataset
        data$Y.m.[which(data$Base == b2C)] = data$Y.m.[which(data$Base == b2C)]+msdataC$Y
      } 
      
      
      i = 1
      for(o in 1:length(combos)){
        b1 = as.numeric(strsplit(combos[o], split = " ")[[1]][1])
        b2 = as.numeric(strsplit(combos[o], split = " ")[[1]][3])
        
        adata <- data[which(data$TreeID %in% overlaps & data$Base == b1),]
        bdata <- data[which(data$TreeID %in% overlaps & data$Base == b2),]
        
        adata <- adata[which(adata$TreeID %in% bdata$TreeID),]
        bdata <- bdata[which(bdata$TreeID %in% adata$TreeID),]
        
        if(nrow(adata) < 1 | nrow(bdata) < 1) next # skip to next loop iteration if there is no row in a or b
        
        Base12 <- rbind(adata,bdata)
        #Base12
        
        names = unique(Base12$TreeID)
        if(length(names) < ref_num) next
        if(z_corr == TRUE){
          ms = data.frame(X = NA, Y = NA, Z = NA)
        }else{
          ms = data.frame(X = NA, Y = NA) 
        }
        for(c in 1:length(names)){
          if(z_corr == TRUE){
            Basems = Base12[which(Base12$TreeID == names[c]),]
            ms[c,"X"] = (Basems$X.m.[1] - Basems$X.m.[2])
            ms[c,"Y"] = (Basems$Y.m.[1] - Basems$Y.m.[2])
            ms[c,"Z"] = (Basems$ALTITUDE[1] - Basems$ALTITUDE[2])
          }else{
            Basems = Base12[which(Base12$TreeID == names[c]),]
            ms[c,"X"] = (Basems$X.m.[1] - Basems$X.m.[2])
            ms[c,"Y"] = (Basems$Y.m.[1] - Basems$Y.m.[2])
          }
        }
        
        if(z_corr == TRUE){
          means[i,] <- c(mean(ms$X), mean(ms$Y), mean(ms$Z), combos[o])
        }else{
          means[i,] <- c(mean(ms$X), mean(ms$Y), combos[o]) 
        }
        i = i + 1
      }
      #means
      
      if(z_corr == TRUE){
        lms[k,1] = sum(as.numeric(means$X)^2 + as.numeric(means$Y)^2 + as.numeric(means$Z)^2)
        lms[k,2] = k
        lms
      }else{
        lms[k,1] = sum(as.numeric(means$X)^2 + as.numeric(means$Y)^2)
        lms[k,2] = k
        lms
      }
      
      k = k + 1
      if(round(lms[k-2,1],4) <= round(lms[k-1,1],4)){ # compare current ms with the previous iteration, 
        # if the current is similar to the one before or even higher it didn't help
        if(z_corr == TRUE){
          data$X.m.[which(data$Base == b2C)] = data$X.m.[which(data$Base == b2C)]-msdataC$X # undo what have been done, since the change did not improve anything
          data$Y.m.[which(data$Base == b2C)] = data$Y.m.[which(data$Base == b2C)]-msdataC$Y
          data$ALTITUDE[which(data$Base == b2C)] = data$ALTITUDE[which(data$Base == b2C)]-msdataC$Z
        }else{
          data$X.m.[which(data$Base == b2C)] = data$X.m.[which(data$Base == b2C)]-msdataC$X # undo what have been done, since the change did not improve anything
          data$Y.m.[which(data$Base == b2C)] = data$Y.m.[which(data$Base == b2C)]-msdataC$Y 
        }
        k = k - 1
        death = death + 1 # add one to the counter
      }else{death = 1} # every time something was improved the counter is set back to 1
    }  
    
    if(death >= limit) break # if the counter hits the limit the repeat loop stops 
  }
  
  if(plot == TRUE){
    # create "after" plot of adjusted data
    p2 = ggplot2::ggplot(data, ggplot2::aes(x = Y.m., y = X.m.))+
      ggplot2::geom_point(ggplot2::aes(color = as.factor(Base), size = TreeDBH))+
      ggplot2::theme_light()+
      ggrepel::geom_text_repel(ggplot2::aes(label=TreeID), cex = 2)+
      ggplot2::labs(x = "X", y = "Y", title = "After")
  }
  
  if(plot == TRUE){colnames(data)[which(colnames(data) == "TreeDBH")] <- Tree_DBH}
  if(z_corr == TRUE){colnames(data)[which(colnames(data) == "ALTITUDE")] <- Alt_col}
  colnames(data)[which(colnames(data) == "Base")] <- Base_ID
  colnames(data)[which(colnames(data) == "TreeID")] <- Tree_ID 
  colnames(data)[which(colnames(data) == "X.m.")] <- X_col
  colnames(data)[which(colnames(data) == "Y.m.")] <- Y_col
  
  if(plot == TRUE){
    return(list(data,  print(gridExtra::grid.arrange(p1, p2, nrow=2)), message(cat(paste("###################################\nAfter k =",k-1,"steps, no improvement for", limit,"iterations in a row \nLast iteration improved the map to a sum of squares of",lms$sum[k-1],"\n###################################")))))
  }else{
    return(list(data, message(cat(paste("###################################\nAfter k =",k-1,"steps, no improvement for", limit,"iterations in a row \nLast iteration improved the map to a sum of squares of",lms$sum[k-1],"\n###################################")))))
  }
}
