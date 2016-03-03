require(spatstat)
require(jpeg)
require(EBImage)
require(mmand)
require(RNiftyReg)
require(xlsx)
require(PET)

# Directories #

# Working directory (where the pictures are)
wd <- "./Input_Photos"
files <- list.files(path=wd, pattern=".JPG", all.files=T, full.names=T)

# CSV files directory (for saving component labels and output record)
csvd <- "./Component_Labelling/Output_Record"

# Filtered images directory (for saving the images after being processed)
areasd <- "./Component_Labelling/Filtered"

# New reduced size of the pictures (resizing to jxj)
j <- 500

# Photo classification (the output)
output <- c("photo", "status")

###############################################################################
# Image Pre-processing

for (file in files) {
  
  filename <- gsub("/","",gsub(wd,"",file))
  img <- resize(readJPEG(file, native = FALSE), j, j)
  
  kernel <- shapeKernel(c(3,3), type="diamond")
  imgL <- dilate(img[,,3],kernel) - erode(img[,,3],kernel)
 
  ###############################################################################
  ## Component Labelling ##
  
  # Preparations
  labels <- matrix(data = NA, nrow = j, ncol = j)  #Labels matrix
  l <- 1        # Counter for setting label (areas)
  areas <- matrix(data = NA, nrow = j*j, ncol = 1)     #Equivalences table (Areas)
 
  v <- 0.0225     # Threshold value to determine if pixels are related
  
  # Two-pass based algorithm #
  for (iy in 1:j){
    for (ix in 1:j) {
      
      if(ix == 1 & iy == 1){
        labels[ix,iy] <- 1
        areas[1] <- 1
      }
      
      #Check West match -> assign West label
      if(iy > 1){
        if(abs( imgL[ix,(iy-1)] - imgL[ix,iy] ) < v){
          labels[ix,iy] <- labels[ix,(iy-1)]
        }
      }
      
      #Check West & North match -> same area & different labels? -> keep min label & update equivalence table
      if(ix > 1 & iy > 1){
        if(abs( imgL[ix,(iy-1)] - imgL[ix,iy] ) < v  & abs( imgL[(ix-1),iy] - imgL[ix,iy] ) < v  & labels[(ix-1),iy] != labels[ix,(iy-1)]){
          labels[ix,iy] <- min( labels[(ix-1),iy], labels[ix,(iy-1)] )
          
          #Updating equivalence table: labelsnNW <- min(min labelNW, current equivalence assigned)
          areas[labels[(ix-1),iy]] <- min(labels[(ix-1),iy], labels[ix,(iy-1)], areas[labels[(ix-1),iy]], na.rm = TRUE)
          areas[labels[ix,(iy-1)]] <- min(labels[(ix-1),iy], labels[ix,(iy-1)], areas[labels[ix,(iy-1)]], na.rm = TRUE)
        }
      }
      
      #Check North match if West didn't exist -> assign North label
      if(iy == 1 & ix > 1){
        if(abs( imgL[(ix-1),iy] - imgL[ix,iy] ) < v){
          labels[ix,iy] <- labels[(ix-1),iy]
        }
      }
      
      #Check North match if West didn't match -> assign North label
      if(iy > 1 & ix > 1){
        if(abs( imgL[ix,(iy-1)] - imgL[ix,iy] ) >= v  &  abs(imgL[(ix-1),iy] - imgL[ix,iy] ) < v){
          labels[ix,iy] <- labels[(ix-1),iy]
        }
      }
      
      #Neither North or West matched (still NA) -> create new label and assign it to the pixel
      if(is.na(labels[ix,iy])){
        l <- l+1
        labels[ix,iy] <- l
        areas[l] <- (j*j)+1
      }
    }
  }
  
  #Complete the default equivalent areas
  for(ix in 1:(j*j)){
    if(!is.na(areas[ix]) & areas[ix] == (j*j)+1){
      areas[ix] <- ix
    }
  }
  
  #Step 2 Equivalence of every single area
  eq <- matrix(data = NA, nrow = j*j, ncol = 1)
  
  for(ix in 1:(j*j)){
    if(!is.na(areas[ix])){
      ifelse(areas[ix] == ix, eq[ix] <- ix, eq[ix] <- eq[areas[ix]])
    }
  }
  
  #Replacing area labels
  for(ix in (j*j):1){
    if(!is.na(eq[ix])){
      labels[which(labels == ix)] <- eq[ix]
    }
  }
  
  # Saving the area labellings into a csv file
  test <- labels
  #display(labels)  
  fcsv <- gsub(".JPG", ".csv", gsub(wd,csvd,file))
  #write.table(labels, fcsv, sep=",", row.names = FALSE, col.names = FALSE)
  
  
  ###############################################################################
  # Recovering the pixel-area labels from csv (optional) #
  
  #tst <-  read.table(fcsv, header=TRUE, sep=",")
  #display(tst)
  #
  #for (ix in 1:j){
  #  for (iy in 1:j) {
  #    test[ix,iy] <- tst[ix,iy]
  #  }
  #}
  #display(test)
  #rm(tst)
  
  ###############################################################################
  ## Area Description ##
  
  # Area catalogue ( cat[area] = size, minrow, maxrow, mincol, maxcol )
  cat <- matrix(0,max(test),5)
  
  #Number of pixels required to tag an area as potentially useful
  minval <- 0.002*j*j
  maxval <- 0.20*j*j
  
  #Calculating area sizes
  for(i in min(test):max(test)){
    ones <- test == i
    countOnes <- colSums(ones)
    size <- sum(countOnes)
    cat[i,1] <- size
    #print(i)    #Just indicating this is still running and to check the number of regions created
    #print(size) #Just indicating this is still running and to check the number of regions created
  }
  
  #Finding minrow, maxrow, mincol, maxcol for every area
  for (ix in 1:j){
    for(iy in 1:j){
      a <- test[ix,iy]
      ifelse(cat[a,2] == 0, cat[a,2] <- ix, cat[a,2] <- min(cat[a,2], ix, na.rm = T))
      ifelse(cat[a,3] == 0, cat[a,3] <- ix, cat[a,3] <- max(cat[a,3], ix, na.rm = T))
      ifelse(cat[a,4] == 0, cat[a,4] <- iy, cat[a,4] <- min(cat[a,4], iy, na.rm = T))
      ifelse(cat[a,5] == 0, cat[a,5] <- iy, cat[a,5] <- max(cat[a,5], iy, na.rm = T))
    }
  }
  
  ###############################################################################
  ## Area Selection / Filtering ##
  
  #Discarding 'useless' areas by total size (# of pixels)
  for (i in 1:nrow(cat)) {
    if(cat[i,1] < minval | cat[i,1] > maxval ){
      test[which(test == i)] <- 0
    }
    else{
      #Discarding 'useless' areas by box size or shape
      
      #Lacking length
      if((cat[i,5] - cat[i,4] + 1) < j/8 ){
        test[which(test == i)] <- 0
      }
      
      #Exceding height
      if((cat[i,3] - cat[i,2] + 1) > j*.4 & (cat[i,3] - cat[i,2] + 1) * (cat[i,5] - cat[i,4] + 1) * .5 < cat[i,1]){
        test[which(test == i)] <- 0
      }
      
      #Vertical areas
      if( (cat[i,3] - cat[i,2] + 1) > (cat[i,5] - cat[i,4] + 1) * 1.5){
        test[which(test == i)] <- 0
      }
      
      #Tall broad areas
      if( (cat[i,3] - cat[i,2] + 1) > j/7 & (cat[i,3] - cat[i,2] + 1) * (cat[i,5] - cat[i,4] + 1) * .5 < cat[i,1] ){
        test[which(test == i)] <- 0
      }
      
      #Northern frame
      if(cat[i,2] == 1 & (cat[i,5] - cat[i,4] + 1) < j/3){
        test[which(test == i)] <- 0
      }
      
      #Southern frame
      if(cat[i,3] == j & (cat[i,5] - cat[i,4] + 1) < j/3){
        test[which(test == i)] <- 0
      }
      
      #Western frame
      if(cat[i,4] == 1 & (cat[i,5] - cat[i,4] + 1) < j/3){
        test[which(test == i)] <- 0
      }
      
      #Eastern frame
      if(cat[i,5] == j & (cat[i,5] - cat[i,4] + 1) < j/3){
        test[which(test == i)] <- 0
      }
      
      #Sky
      if(cat[i,2] == 1 & (cat[i,5] - cat[i,4] + 1) > j*.9){
        test[which(test == i)] <- 0
      }
      
      #Western frame II
      if(cat[i,4] == 1 & (cat[i,3] - cat[i,2] + 1) > j/6){
        test[which(test == i)] <- 0
      }
      
      #Eastern frame II
      if(cat[i,5] == j & (cat[i,3] - cat[i,2] + 1) > j/6){
        test[which(test == i)] <- 0
      }
      
    }
    #print(i)  #Just indicating this is still running and to check the number of regions created
  }
  
  #display(test)
  imgLabel <- img
  imgLabel[test[,]==0] <- 0
  #display(imgLabel)
  
  #Color Thresholding
  imgT <- imgLabel
  
  for (ix in 1:j){
    for (iy in 1:j) {
      
      # Red
      if( imgT[ix,iy,1]>(imgT[ix,iy,2]+.1) & imgT[ix,iy,1]>(imgT[ix,iy,3]+.1) ){
        imgT[ix,iy,1] <- 0
        imgT[ix,iy,2] <- 0
        imgT[ix,iy,3] <- 0
      }
      
      # Green 
      if( imgT[ix,iy,2]>(imgT[ix,iy,1]+.05) & imgT[ix,iy,2]>(imgT[ix,iy,3]+.05) ){
        imgT[ix,iy,1] <- 0
        imgT[ix,iy,2] <- 0
        imgT[ix,iy,3] <- 0
      }
      
      # Blue
      if( imgT[ix,iy,3]>(imgT[ix,iy,1]+.15) & imgT[ix,iy,3]>(imgT[ix,iy,2]+.15) ){
        imgT[ix,iy,1] <- 0
        imgT[ix,iy,2] <- 0
        imgT[ix,iy,3] <- 0
      }
      
      # Dark
      if( imgT[ix,iy,1]<.2 & imgT[ix,iy,2]<.15 & imgT[ix,iy,3]<.15 ){
        imgT[ix,iy,1] <- 0
        imgT[ix,iy,2] <- 0
        imgT[ix,iy,3] <- 0
      }
      
      # Trees
      if( imgT[ix,iy,2]>imgT[ix,iy,1] & imgT[ix,iy,2]>imgT[ix,iy,3] & imgT[ix,iy,2]<.2 ){
        imgT[ix,iy,1] <- 0
        imgT[ix,iy,2] <- 0
        imgT[ix,iy,3] <- 0
      }
      
      # Trees (more)
      if( imgT[ix,iy,1]>imgT[ix,iy,2] & imgT[ix,iy,1]>imgT[ix,iy,3] & imgT[ix,iy,2]<.2 & imgT[ix,iy,1]<.4){
        imgT[ix,iy,1] <- 0
        imgT[ix,iy,2] <- 0
        imgT[ix,iy,3] <- 0
      }
      
      # Yellow (more green areas)
      if ( (imgT[ix,iy,1]-imgT[ix,iy,2] < .05) & imgT[ix,iy,1]>(imgT[ix,iy,3]+.2) & imgT[ix,iy,2]>(imgT[ix,iy,3]+.2) ) {
        imgT[ix,iy,1] <- 0
        imgT[ix,iy,2] <- 0
        imgT[ix,iy,3] <- 0
      }
    }
  }
  
  
  ###############################################################################
  # Second edge detection
  imgL <- dilate(imgT[,,3],kernel) - erode(imgT[,,3],kernel)
  
  ###############################################################################
  ## Second Component Labelling ##
  
  # Preparations
  labels <- matrix(data = NA, nrow = j, ncol = j)  #Labels matrix
  l <- 1        # Counter for setting label (areas)
  areas <- matrix(data = NA, nrow = j*j, ncol = 1)     #Equivalences table (Areas)
  #v <- sd(imgL, na.rm = TRUE)/8     # Initial threshold value to determine if pixels are related, decided to change due to variant SDs
  v <- 0.02     # Threshold value to determine if pixels are related
  
  labels <- matrix(data = NA, nrow = j, ncol = j)  #Labels matrix
  l <- 1        # Counter for setting label (areas)
  areas <- matrix(data = NA, nrow = j*j, ncol = 1)     #1st step equivalences table (Areas)
  
  
  # Two-pass based algorithm #
  for (iy in 1:j){
    for (ix in 1:j) {
      
      if(ix == 1 & iy == 1){
        labels[ix,iy] <- 1
        areas[1] <- 1
      }
      
      #Check West match -> assign West label
      if(iy > 1){
        if(abs( imgL[ix,(iy-1)] - imgL[ix,iy] ) < v){
          labels[ix,iy] <- labels[ix,(iy-1)]
        }
      }
      
      #Check West & North match -> same area & different labels? -> keep min label & update equivalence table
      if(ix > 1 & iy > 1){
        if(abs( imgL[ix,(iy-1)] - imgL[ix,iy] ) < v  & abs( imgL[(ix-1),iy] - imgL[ix,iy] ) < v  & labels[(ix-1),iy] != labels[ix,(iy-1)]){
          labels[ix,iy] <- min( labels[(ix-1),iy], labels[ix,(iy-1)] )
          
          #Updating equivalence table: labelsnNW <- min(min labelNW, current equivalence assigned)
          areas[labels[(ix-1),iy]] <- min(labels[(ix-1),iy], labels[ix,(iy-1)], areas[labels[(ix-1),iy]], na.rm = TRUE)
          areas[labels[ix,(iy-1)]] <- min(labels[(ix-1),iy], labels[ix,(iy-1)], areas[labels[ix,(iy-1)]], na.rm = TRUE)
        }
      }
      
      #Check North match if West didn't exist -> assign North label
      if(iy == 1 & ix > 1){
        if(abs( imgL[(ix-1),iy] - imgL[ix,iy] ) < v){
          labels[ix,iy] <- labels[(ix-1),iy]
        }
      }
      
      #Check North match if West didn't match -> assign North label
      if(iy > 1 & ix > 1){
        if(abs( imgL[ix,(iy-1)] - imgL[ix,iy] ) >= v  &  abs(imgL[(ix-1),iy] - imgL[ix,iy] ) < v){
          labels[ix,iy] <- labels[(ix-1),iy]
        }
      }
      
      #Neither North or West matched (still NA) -> create new label and assign it to the pixel
      if(is.na(labels[ix,iy])){
        l <- l+1
        labels[ix,iy] <- l
        areas[l] <- (j*j)+1
      }
    }
  }
  
  #Complete the default equivalent areas
  for(ix in 1:(j*j)){
    if(!is.na(areas[ix]) & areas[ix] == (j*j)+1){
      areas[ix] <- ix
    }
  }
  
  #Step 2 Equivalence of every single area
  eq <- matrix(data = NA, nrow = j*j, ncol = 1)
  
  for(ix in 1:(j*j)){
    if(!is.na(areas[ix])){
      ifelse(areas[ix] == ix, eq[ix] <- ix, eq[ix] <- eq[areas[ix]])
    }
  }
  
  #Replacing area labels
  for(ix in (j*j):1){
    if(!is.na(eq[ix])){
      labels[which(labels == ix)] <- eq[ix]
    }
  }
  
  test <- labels
  
  ###############################################################################
  ## Second Area Description ##
  
  #Label catalogue ( cat[area] = size, minrow, maxrow, mincol, maxcol )
  cat <- matrix(0,max(test),5)
  
  #Number of pixels required to tag an area as potentially useful
  minval <- 0.002*j*j
  maxval <- 0.20*j*j
  
  #Calculating area sizes
  for(i in min(test):max(test)){
    ones <- test == i
    countOnes <- colSums(ones)
    size <- sum(countOnes)
    cat[i,1] <- size
    print(i)    #Just indicating this is still running and to check the number of regions created
    print(size) #Just indicating this is still running and to check the number of regions created
  }
  
  #Finding minrow, maxrow, mincol, maxcol for every area
  for (ix in 1:j){
    for(iy in 1:j){
      a <- test[ix,iy]
      ifelse(cat[a,2] == 0, cat[a,2] <- ix, cat[a,2] <- min(cat[a,2], ix, na.rm = T))
      ifelse(cat[a,3] == 0, cat[a,3] <- ix, cat[a,3] <- max(cat[a,3], ix, na.rm = T))
      ifelse(cat[a,4] == 0, cat[a,4] <- iy, cat[a,4] <- min(cat[a,4], iy, na.rm = T))
      ifelse(cat[a,5] == 0, cat[a,5] <- iy, cat[a,5] <- max(cat[a,5], iy, na.rm = T))
    }
  }
  
  ###############################################################################
  ## Second Area Selection / Filtering ##
  
  #Discarding 'useless' areas by total size (# of pixels)
  for (i in 1:nrow(cat)) {
    if(cat[i,1] < minval | cat[i,1] > maxval ){
      test[which(test == i)] <- 0
    }
    else{
      #Discarding 'useless' areas by box size or shape
      
      #Lacking length (reduced)
      if((cat[i,5] - cat[i,4] + 1) < j/16 ){
        test[which(test == i)] <- 0
      }
      
      #Vertical areas
      if( (cat[i,3] - cat[i,2] + 1) > (cat[i,5] - cat[i,4] + 1) * 1.5){
        test[which(test == i)] <- 0
      }
      
      #Northern frame
      if(cat[i,2] == 1 & (cat[i,5] - cat[i,4] + 1) < j/3){
        test[which(test == i)] <- 0
      }
      
      #Southern frame
      if(cat[i,3] == j & (cat[i,5] - cat[i,4] + 1) < j/3){
        test[which(test == i)] <- 0
      }
      
      #Western frame
      if(cat[i,4] == 1 & (cat[i,5] - cat[i,4] + 1) < j/3){
        test[which(test == i)] <- 0
      }
      
      #Eastern frame
      if(cat[i,5] == j & (cat[i,5] - cat[i,4] + 1) < j/3){
        test[which(test == i)] <- 0
      }
    }
    print(i)  #Just indicating this is still running and to check the number of regions created
  }
  
  #display(test)
  imgLabel <- img
  imgLabel[test[,]==0] <- 0
  #display(imgLabel)
  
  
  
  ###############################################################################
  # Saving the filtered areas into a JPG (optional)
  
  #fareas <- gsub(wd,areasd,gsub(".JPG", " filtered areas.JPG", file))
  #display(imgLabel)
  #assign(fareas, imgLabel)
  #dev.copy(jpeg, filename = fareas)
  #dev.off()
  
  ###############################################################################
  # Picture Classification
  
  # Creates a copy of the filtered areas to perform the final region filtering
  # It uses a labels copy so all the values on the matrix are the same (label number) for every element enclosed in the same region
  regions <- test
  
  # Variable that classifies the picture
  tag <- "not useful"
  
  # Regions analysis
  while (max(regions) > 0){
    pick <- max(regions)
    
    # Verifies if the picked area is not the background or an irrelevent region (to discard it if gets picked) and then performs the filtering by shape until finds 1 region matching the criteria or finishes discarding all of them
    if (cat[pick,1] > maxval | cat[pick,5] - cat[pick,4] < j/20 ){
      regions[regions == pick] <- 0
    } else {
      # Creates the matrix that will store the column count of those that matched the criteria
      column <- matrix(0,cat[pick,5] - cat[pick,4] + 1,1)
      
      # Scans through columns counting the number of pixels per column enclosed on the selected region
      for(iy in cat[pick,4]:cat[pick,5]){
        countingPixels <- 0
        for(ix in cat[pick,2]:cat[pick,3]){
          countingPixels <- countingPixels + ifelse(regions[ix,iy] == pick,1,0)
        }
        column[(iy + 1 - cat[pick,4]),1] <- countingPixels
      }
      
      # Count number of columns matching the criteria
      countingColumns <- 0
      for (ix in 1:(cat[pick,5] - cat[pick,4] + 1)) {
        countingColumns <- countingColumns + ifelse(column[ix,1] < j/33, 0, 1)
      }
      # If an area has the minimum number of matching columns then  break the loop and classify the picture as potential, else keep looking
      if(countingColumns >= j/10){
        regions[,] <- 0
        tag <- "crossarm"
      } else {
        regions[regions == pick] <- 0
      }
    }
  } 
  
  # Saving photo classification into output variable
  output <- rbind(output , c(filename,tag))
}

# Saving the output into a csv file
fcsv <- gsub(filename, "output.csv", gsub(wd,csvd,file))
write.table(output, fcsv, sep=",",row.names = FALSE, col.names = FALSE)
#print(output)