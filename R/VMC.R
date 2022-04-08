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

VMC <- function(data, Tree_DBH = "DBH", X_col = "X.m.", Y_col = "Y.m.", Alt_col = "ALTITUDE", Base_ID = "Base", Tree_ID = "TreeID", limit = 10, ref_num = 2, plot = FALSE, z_corr = FALSE){
  
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
  means = data.frame(X = NA, Y = NA, Base_combo = NA)

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
      ms = data.frame(X = NA, Y = NA)
      for(c in 1:length(names)){
        Basems = Base12[which(Base12$TreeID == names[c]),]
        ms[c,"X"] = (Basems$X.m.[1] - Basems$X.m.[2])
        ms[c,"Y"] = (Basems$Y.m.[1] - Basems$Y.m.[2])
      }
      means[i,] <- c(mean(ms$X), mean(ms$Y), as.character(paste(sort(c(a,b))[1], "&", sort(c(a,b))[2])))
      i = i + 1
    }
  }
  
  k = 1
  lms$sum[k] = sum(as.numeric(means$X)^2 + as.numeric(means$Y)^2) 
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

      msC = data.frame(X = NA, Y = NA)
      for(a in 1:length(names)){
        BasemsC = Base12C[which(Base12C$TreeID == names[a]),]
        msC[a,"X"] = (BasemsC$X.m.[1] - BasemsC$X.m.[2])
        msC[a,"Y"] = (BasemsC$Y.m.[1] - BasemsC$Y.m.[2])
      }
      
      msdataC <- data.frame(X = mean(msC$X), Y = mean(msC$Y)) # calculate mean offset
        
      data$X.m.[which(data$Base == b2C)] = data$X.m.[which(data$Base == b2C)]+msdataC$X # correct entire dataset
      data$Y.m.[which(data$Base == b2C)] = data$Y.m.[which(data$Base == b2C)]+msdataC$Y
      
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
        ms = data.frame(X = NA, Y = NA) 
        for(c in 1:length(names)){
          Basems = Base12[which(Base12$TreeID == names[c]),]
          ms[c,"X"] = (Basems$X.m.[1] - Basems$X.m.[2])
          ms[c,"Y"] = (Basems$Y.m.[1] - Basems$Y.m.[2])
        }
        
        means[i,] <- c(mean(ms$X), mean(ms$Y), combos[o]) 
        i = i + 1
      }
      #means
      
      lms[k,1] = sum(as.numeric(means$X)^2 + as.numeric(means$Y)^2)
      lms[k,2] = k
      #lms

      k = k + 1
      if(round(lms[k-2,1],4) <= round(lms[k-1,1],4)){ # compare current ms with the previous iteration, 
        # if the current is similar to the one before or even higher it didn't help
        data$X.m.[which(data$Base == b2C)] = data$X.m.[which(data$Base == b2C)]-msdataC$X # undo what have been done, since the change did not improve anything
        data$Y.m.[which(data$Base == b2C)] = data$Y.m.[which(data$Base == b2C)]-msdataC$Y 
        
        k = k - 1
        death = death + 1 # add one to the counter
      }else{death = 1} # every time something was improved the counter is set back to 1
    }  
    
    if(death >= limit) break # if the counter hits the limit the repeat loop stops 
  }
  
  ################################################################
  # z_correction to adjust altitude
  ################################################################
  if(z_corr == TRUE){
    overlapsZ <- names(which(table(data$TreeID) > 1))
    BasesZ <- max(na.omit(data$Base))
    lmsZ <- data.frame(sum = NA, iteration = NA)
    
    # calculate initial mean square
    ctrlZ = NA
    iZ = 1
    meansZ = data.frame(Z = NA, Base_combo = NA)
    
    for(a in 1:BasesZ){ # loop 1 of all Bases
      for(b in 1:BasesZ){ # loop 2 of all Bases
        if(b == a) next # next iteration if a and b are the same Base
        ctrlZ = paste(sort(c(a,b))[1], "&", sort(c(a,b))[2]) # control line to avoid doubling base-pairs
        if(ctrlZ %in% meansZ$Base_combo) next # next iteration if the combination was already calculated
        
        adataZ <- data[which(data$TreeID %in% overlaps & data$Base == a),]
        bdataZ <- data[which(data$TreeID %in% overlaps & data$Base == b),]
        
        adataZ <- adataZ[which(adataZ$TreeID %in% bdataZ$TreeID),]
        bdataZ <- bdataZ[which(bdataZ$TreeID %in% adataZ$TreeID),]
        
        if(nrow(adataZ) < 1 | nrow(bdataZ) < 1) next # skip to next loop iteration if there is no row in a or b
        
        Base12Z <- rbind(adataZ,bdataZ)
        
        namesZ = unique(Base12Z$TreeID)
        if(length(namesZ) < ref_num) next # skip to next loop iteration if the number of "reference tree" is below the provided amount (should be at least (default) 2)
        msZ = data.frame(Z = NA)
        for(c in 1:length(namesZ)){
          BasemsZ = Base12Z[which(Base12Z$TreeID == namesZ[c]),]
          msZ[c,"Z"] = (BasemsZ$ALTITUDE[1] - BasemsZ$ALTITUDE[2])
        }
        meansZ[iZ,] <- c(mean(msZ$Z), as.character(paste(sort(c(a,b))[1], "&", sort(c(a,b))[2])))
        iZ = iZ + 1
      }
    }
    
    kZ = 1
    lmsZ$sum[kZ] = sum(as.numeric(meansZ$Z)^2) 
    lmsZ$iteration[kZ] = kZ
    #lms$sum[1]
    kZ = kZ + 1
    
    combosZ = meansZ$Base_combo # now select possible combos and loop through
    deathZ = 1
    repeat{
      for(l in 1:length(combosZ)){
        b1CZ = as.numeric(strsplit(combosZ[l], split = " ")[[1]][1])
        b2CZ = as.numeric(strsplit(combosZ[l], split = " ")[[1]][3])
        
        # merging overlapping points of 2 Bases
        adataCZ <- data[which(data$TreeID %in% overlaps & data$Base == b1CZ),]
        bdataCZ <- data[which(data$TreeID %in% overlaps & data$Base == b2CZ),]
        
        adataCZ <- adataCZ[which(adataCZ$TreeID %in% bdataCZ$TreeID),]
        bdataCZ <- bdataCZ[which(bdataCZ$TreeID %in% adataCZ$TreeID),]
        
        Base12CZ <- rbind(adataCZ,bdataCZ)
        Base12CZ
        
        # calc difference for each tree
        namesZ = unique(Base12CZ$TreeID)
        
        msCZ = data.frame(Z = NA)
        for(a in 1:length(namesZ)){
          BasemsCZ = Base12CZ[which(Base12CZ$TreeID == namesZ[a]),]
          msCZ[a,"Z"] = (BasemsCZ$ALTITUDE[1] - BasemsCZ$ALTITUDE[2])
        }
        
        msdataCZ <- data.frame(Z = mean(msCZ$Z)) # calculate mean offset
        
        data$ALTITUDE[which(data$Base == b2CZ)] = data$ALTITUDE[which(data$Base == b2CZ)]+msdataCZ$Z # correct entire dataset
        
        iZ = 1
        for(o in 1:length(combosZ)){
          b1Z = as.numeric(strsplit(combosZ[o], split = " ")[[1]][1])
          b2Z = as.numeric(strsplit(combosZ[o], split = " ")[[1]][3])
          
          adataZ <- data[which(data$TreeID %in% overlaps & data$Base == b1Z),]
          bdataZ <- data[which(data$TreeID %in% overlaps & data$Base == b2Z),]
          
          adataZ <- adataZ[which(adataZ$TreeID %in% bdataZ$TreeID),]
          bdataZ <- bdataZ[which(bdataZ$TreeID %in% adataZ$TreeID),]
          
          if(nrow(adataZ) < 1 | nrow(bdataZ) < 1) next # skip to next loop iteration if there is no row in a or b
          
          Base12Z <- rbind(adataZ,bdataZ)
          #Base12
          
          namesZ = unique(Base12Z$TreeID)
          if(length(namesZ) < ref_num) next
          msZ = data.frame(Z = NA) 
          for(c in 1:length(namesZ)){
            BasemsZ = Base12Z[which(Base12Z$TreeID == namesZ[c]),]
            msZ[c,"Z"] = (BasemsZ$ALTITUDE[1] - BasemsZ$ALTITUDE[2])
          }
          
          meansZ[iZ,] <- c(mean(msZ$Z), combosZ[o]) 
          iZ = iZ + 1
        }
        #meansZ
        
        lmsZ[kZ,1] = sum(as.numeric(meansZ$Z)^2)
        lmsZ[kZ,2] = kZ
        #lmsZ
        
        kZ = kZ + 1
        if(kZ-1 > length(combosZ)+1){ # additional condition: only compare iterations after "correcting" each base combo once
          if(round(lmsZ[kZ-2,1],4) <= round(lmsZ[kZ-1,1],4)){ # compare current ms with the previous iteration, 
            # if the current is similar to the one before or even higher it didn't help
            data$ALTITUDE[which(data$Base == b2CZ)] = data$ALTITUDE[which(data$Base == b2CZ)]-msdataCZ$Z
            
            kZ = kZ - 1
            deathZ = deathZ + 1 # add one to the counter
          }else{deathZ = 1} # every time something was improved the counter is set back to 1
        }
      }  
      
      if(deathZ >= limit) break # if the counter hits the limit the repeat loop stops 
    }
    
    # something to avoid distribution to be parallel is needed
    
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
