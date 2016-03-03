folderfunc <- function(f)
{
  picturenames <- list.files(paste0(overall.direc, "\\",foldernames[f]))
  nopics <- length(picturenames)
  kclusters <- matrix(0, nrow = nopics, ncol = 250)
  
  for (v in 1:nopics)
  {
  
    picture <- readJPEG(paste0(overall.direc, "\\", foldernames[f],"\\",picturenames[v] ))
    
    # initialise stuff
    no_of_octaves <- 5
    sqrt_2 <- sqrt(2)
    kern_size <- 7
    keypoints <- c(0,0,0,0)
    no_of_keypoints <- 0
    
    # convert to greyscale and normalise
    pic_grey <- (picture[,,1] + picture[,,2] + picture[,,3]) / 3
    pic_grey <- (pic_grey - mean(pic_grey)) / sd(pic_grey)
    rm(picture)
    
    # flip it so it's the right way up
    pic_greyv <- rev(pic_grey)
    pic_greyv <- matrix(pic_greyv, nrow(pic_grey), ncol(pic_grey))
    pic_greyv <- t(pic_greyv)
    pic_greyv <- apply(pic_greyv, 2, rev)
    pic_grey <- pic_greyv
    rm(pic_greyv)
    
    # expand by 2, calc new number of rows and cols
    #newcol <- (ncol(pic_grey)) * 2 + 1
    #newrow <-(nrow(pic_grey)) * 2 + 1
    #ss <- raster(ncol = newcol, nrow = newrow)
    #ss[] <- ncell(ss)
    #piccyras <- raster(pic_grey)
    #extent(piccyras) <- extent(ss)
    #tempn <- resample(piccyras, ss, method = "bilinear")
    #temp = as.matrix(tempn)
    
    #big_rows <- nrow(temp)
    #big_cols <- ncol(temp)
    
    temp <- pic_grey
    
    big_rows <- nrow(pic_grey)
    big_cols <- ncol(pic_grey)
    
    # initialise first octave
    blurs <- array(0, c(nrow(temp),ncol(temp),6,no_of_octaves))
    
    # FIRST OCTAVE
    
    # store big version in first plane of first octave
    blurs[,,1,1] <- temp
    
    print("initialising for octave 1")
    # intialise stuff for fourier
    ksize = 7
    gausskern <- matrix(0, nrow = ksize, ncol = ksize)
    sigma <- sqrt(2) 
    m = ksize/2
    
    hb <- nrow(blurs[,,1,1])
    wb <- ncol(blurs[,,1,1])
    
    # make X and Y matricies
    X <- matrix(1:ksize,ksize,ksize, byrow = TRUE)
    Y <- matrix(1:ksize,ksize,ksize, byrow = FALSE)
    
    # add space to blurs[,,1,1]
    impad <- padimage(blurs[,,1,1], ksize)
    
    #initialise kernel for blurring fft
    gausskern <- (1/(2*pi*sigma^2))*exp(-((X-m)^2 + (Y-m)^2)/(2*sigma^2))
    gausskernpad <- matrix(0, nrow(impad), ncol(impad))
    gausskernpad[1:ksize,1:ksize] <- gausskern
    gausskernfft <- fft(gausskernpad)
    
    # blur blurs[,,1,1] with sigma, store as face 2
    print("blur number 1...")
    imfft <- fft(impad)
    blurredfft <- gausskernfft* imfft
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    temp <- blurredim[((ksize+1):(ksize+hb)),((ksize+1):(ksize+wb))]
    blurs[,,2,1] <- (temp - min(temp)) * (1 / (max(temp) - min(temp)))
    
    print("blur number 2...")
    # blur A with sqrt 2, store as face 3
    impad <- padimage(blurs[,,2,1], ksize)
    imfft <- fft(impad)
    blurredfft <- gausskernfft* imfft
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    temp <- blurredim[((ksize+1):(ksize+hb)),((ksize+1):(ksize+wb))]
    blurs[,,3,1] <- (temp - min(temp)) * (1 / (max(temp) - min(temp)))
    
    print("DoG")
    # DoG, subtract B from A, store as face 4
    temp <- abs(blurs[,,2,1] - blurs[,,3,1])
    blurs[,,4,1] <- (temp - min(temp)) * (1 / (max(temp) - min(temp)))
    
    print("gradient calculations")
    #gradient calculations
    # x edges
    xkern <- matrix(c(-1,-2,-1,0,0,0,1,2,1), 3,3)
    hb <- nrow(blurs[,,1,1])
    wb <- ncol(blurs[,,1,1])
    
    impad <- padimage(blurs[,,1,1], 3)
    xkernpad <- matrix(0, nrow(impad), ncol(impad))
    xkernpad[1:3,1:3] <- xkern
    
    xkernfft <- fft(xkernpad)
    imfft <- fft(impad)
    blurredfft <- xkernfft* imfft
    
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    
    xs <- blurredim[((4):(3+hb)),((3+1):(3+wb))]
    
    # y edges
    ykern <- matrix(c(-1,0,1,-2,0,2,-1,0,1), 3,3)
    ykernpad <- matrix(0, nrow(impad), ncol(impad))
    ykernpad[1:3,1:3] <- ykern
    ykernfft <- fft(ykernpad)
    
    blurredfft <- ykernfft* imfft
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    
    ys <- blurredim[((4):(3+hb)),((3+1):(3+wb))]
    
    # assign face five mags, six orients
    blurs[,,5,1] <- xs^2 + ys^2
    blurs[,,6,1] <- atan2(ys, xs)* (180/pi)
    
    
    
    #SECOND OCTAVE
    # downsample by 1.5, calc new number of rows and cols
    newcol <- round((ncol(blurs[,,3,1]))/1.5)
    newrow <-round((nrow(blurs[,,3,1]))/1.5)
    
    ss <- raster(ncol = newcol, nrow = newrow)
    ss[] <- ncell(ss)
    piccyras <- raster(blurs[,,3,1])
    extent(piccyras) <- extent(ss)
    tempn <- resample(piccyras, ss, method = "bilinear")
    temp = as.matrix(tempn)
    
    big_rows <- nrow(temp)
    big_cols <- ncol(temp)
    
    
    # store big version in first plane of octave
    blurs[(1:nrow(temp)),(1:ncol(temp)),1,2] <- temp
    
    print("initialising for octave 1")
    # intialise shit for fourier
    
    # make X and Y matricies
    X <- matrix(1:ksize,ksize,ksize, byrow = TRUE)
    Y <- matrix(1:ksize,ksize,ksize, byrow = FALSE)
    
    # add space to blurs[,,1,1]
    impad <- padimage(blurs[,,1,2], ksize)
    
    # blur blurs[,,1,2] with sigma, store as face 2
    print("blur number 1...")
    imfft <- fft(impad)
    blurredfft <- gausskernfft* imfft
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    temp <- blurredim[((ksize+1):(ksize+hb)),((ksize+1):(ksize+wb))]
    blurs[,,2,2] <- (temp - min(temp)) * (1 / (max(temp) - min(temp)))
    
    print("blur number 2...")
    # blur A with sqrt 2, store as face 3
    impad <- padimage(blurs[,,2,2], ksize)
    imfft <- fft(impad)
    blurredfft <- gausskernfft* imfft
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    temp <- blurredim[((ksize+1):(ksize+hb)),((ksize+1):(ksize+wb))]
    blurs[,,3,2] <- (temp - min(temp)) * (1 / (max(temp) - min(temp)))
    
    print("DoG")
    # DoG, subtract B from A, store as face 4
    temp <- abs(blurs[,,2,2] - blurs[,,3,2])
    blurs[,,4,2] <- (temp - min(temp)) * (1 / (max(temp) - min(temp)))
    
    print("gradient calculations")
    #gradient calculations
    # x edges
    
    impad <- padimage(blurs[,,1,2], 3)
    
    imfft <- fft(impad)
    blurredfft <- xkernfft* imfft
    
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    
    xs <- blurredim[((4):(3+hb)),((3+1):(3+wb))]
    
    # y edges
    
    
    blurredfft <- ykernfft* imfft
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    
    ys <- blurredim[((4):(3+hb)),((3+1):(3+wb))]
    
    # assign face five mags, six orients
    blurs[,,5,2] <- xs^2 + ys^2
    blurs[,,6,2] <- atan2(ys, xs)* (180/pi)
    
    
    #OCTAVE THREE
    
    # downsample by 1.5, calc new number of rows and cols
    newcol <- round((ncol(blurs[,,3,1]))/1.5)
    newrow <-round((nrow(blurs[,,3,1]))/1.5)
    
    ss <- raster(ncol = newcol, nrow = newrow)
    ss[] <- ncell(ss)
    piccyras <- raster(blurs[,,3,2])
    extent(piccyras) <- extent(ss)
    tempn <- resample(piccyras, ss, method = "bilinear")
    temp = as.matrix(tempn)
    
    big_rows <- nrow(temp)
    big_cols <- ncol(temp)
    
    
    # store big version in first plane of octave
    blurs[(1:nrow(temp)),(1:ncol(temp)),1,3] <- temp
    
    print("initialising for octave 3")
    # intialise shit for fourier
    
    # make X and Y matricies
    X <- matrix(1:ksize,ksize,ksize, byrow = TRUE)
    Y <- matrix(1:ksize,ksize,ksize, byrow = FALSE)
    
    # add space to blurs[,,1,3]
    impad <- padimage(blurs[,,1,3], ksize)
    
    # blur blurs[,,1,3] with sigma, store as face 2
    print("blur number 1...")
    imfft <- fft(impad)
    blurredfft <- gausskernfft* imfft
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    temp <- blurredim[((ksize+1):(ksize+hb)),((ksize+1):(ksize+wb))]
    blurs[,,2,3] <- (temp - min(temp)) * (1 / (max(temp) - min(temp)))
    
    print("blur number 2...")
    # blur A with sqrt 2, store as face 3
    impad <- padimage(blurs[,,2,3], ksize)
    imfft <- fft(impad)
    blurredfft <- gausskernfft* imfft
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    temp <- blurredim[((ksize+1):(ksize+hb)),((ksize+1):(ksize+wb))]
    blurs[,,3,3] <- (temp - min(temp)) * (1 / (max(temp) - min(temp)))
    
    print("DoG")
    # DoG, subtract B from A, store as face 4
    temp <- abs(blurs[,,2,3] - blurs[,,3,3])
    blurs[,,4,3] <- (temp - min(temp)) * (1 / (max(temp) - min(temp)))
    
    print("gradient calculations")
    #gradient calculations
    # x edges
    
    impad <- padimage(blurs[,,1,3], 3)
    
    imfft <- fft(impad)
    blurredfft <- xkernfft* imfft
    
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    
    xs <- blurredim[((4):(3+hb)),((3+1):(3+wb))]
    
    # y edges
    
    blurredfft <- ykernfft* imfft
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    
    ys <- blurredim[((4):(3+hb)),((3+1):(3+wb))]
    
    # assign face five mags, six orients
    blurs[,,5,3] <- xs^2 + ys^2
    blurs[,,6,3] <- atan2(ys, xs)* (180/pi)
    
    
    
    #OCTAVE FOUR
    
    # downsample by 1.5, calc new number of rows and cols
    newcol <- round((ncol(blurs[,,3,1]))/1.5)
    newrow <-round((nrow(blurs[,,3,1]))/1.5)
    
    ss <- raster(ncol = newcol, nrow = newrow)
    ss[] <- ncell(ss)
    piccyras <- raster(blurs[,,3,3])
    extent(piccyras) <- extent(ss)
    tempn <- resample(piccyras, ss, method = "bilinear")
    temp = as.matrix(tempn)
    
    big_rows <- nrow(temp)
    big_cols <- ncol(temp)
    
    
    # store big version in first plane of octave
    blurs[(1:nrow(temp)),(1:ncol(temp)),1,4] <- temp
    
    print("initialising for octave 3")
    # intialise shit for fourier
    
    # make X and Y matricies
    X <- matrix(1:ksize,ksize,ksize, byrow = TRUE)
    Y <- matrix(1:ksize,ksize,ksize, byrow = FALSE)
    
    # add space to blurs[,,1,3]
    impad <- padimage(blurs[,,1,4], ksize)
    
    # blur blurs[,,1,4] with sigma, store as face 2
    print("blur number 1...")
    imfft <- fft(impad)
    blurredfft <- gausskernfft* imfft
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    temp <- blurredim[((ksize+1):(ksize+hb)),((ksize+1):(ksize+wb))]
    blurs[,,2,4] <- (temp - min(temp)) * (1 / (max(temp) - min(temp)))
    
    print("blur number 2...")
    # blur A with sqrt 2, store as face 3
    impad <- padimage(blurs[,,2,4], ksize)
    imfft <- fft(impad)
    blurredfft <- gausskernfft* imfft
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    temp <- blurredim[((ksize+1):(ksize+hb)),((ksize+1):(ksize+wb))]
    blurs[,,3,4] <- (temp - min(temp)) * (1 / (max(temp) - min(temp)))
    
    print("DoG")
    # DoG, subtract B from A, store as face 4
    temp <- abs(blurs[,,2,4] - blurs[,,3,4])
    blurs[,,4,4] <- (temp - min(temp)) * (1 / (max(temp) - min(temp)))
    
    print("gradient calculations")
    #gradient calculations
    # x edges
    
    impad <- padimage(blurs[,,1,4], 3)
    
    imfft <- fft(impad)
    blurredfft <- xkernfft* imfft
    
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    
    xs <- blurredim[((4):(3+hb)),((3+1):(3+wb))]
    
    # y edges
    
    blurredfft <- ykernfft* imfft
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    
    ys <- blurredim[((4):(3+hb)),((3+1):(3+wb))]
    
    # assign face five mags, six orients
    blurs[,,5,4] <- xs^2 + ys^2
    blurs[,,6,4] <- atan2(ys, xs)* (180/pi)
    
    
    #OCTAVE FIVE
    
    # downsample by 1.5, calc new number of rows and cols
    newcol <- round((ncol(blurs[,,3,1]))/1.5)
    newrow <-round((nrow(blurs[,,3,1]))/1.5)
    
    ss <- raster(ncol = newcol, nrow = newrow)
    ss[] <- ncell(ss)
    piccyras <- raster(blurs[,,3,4])
    extent(piccyras) <- extent(ss)
    tempn <- resample(piccyras, ss, method = "bilinear")
    temp = as.matrix(tempn)
    
    big_rows <- nrow(temp)
    big_cols <- ncol(temp)
    
    
    # store big version in first plane of octave
    blurs[(1:nrow(temp)),(1:ncol(temp)),1,5] <- temp
    
    print("initialising for octave 3")
    # intialise shit for fourier
    
    # make X and Y matricies
    X <- matrix(1:ksize,ksize,ksize, byrow = TRUE)
    Y <- matrix(1:ksize,ksize,ksize, byrow = FALSE)
    
    # add space to blurs[,,1,3]
    impad <- padimage(blurs[,,1,5], ksize)
    
    # blur blurs[,,1,3] with sigma, store as face 2
    print("blur number 1...")
    imfft <- fft(impad)
    blurredfft <- gausskernfft* imfft
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    temp <- blurredim[((ksize+1):(ksize+hb)),((ksize+1):(ksize+wb))]
    blurs[,,2,5] <- (temp - min(temp)) * (1 / (max(temp) - min(temp)))
    
    print("blur number 2...")
    # blur A with sqrt 2, store as face 3
    impad <- padimage(blurs[,,2,5], ksize)
    imfft <- fft(impad)
    blurredfft <- gausskernfft* imfft
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    temp <- blurredim[((ksize+1):(ksize+hb)),((ksize+1):(ksize+wb))]
    blurs[,,3,5] <- (temp - min(temp)) * (1 / (max(temp) - min(temp)))
    
    print("DoG")
    # DoG, subtract B from A, store as face 4
    temp <- abs(blurs[,,2,5] - blurs[,,3,5])
    blurs[,,4,5] <- (temp - min(temp)) * (1 / (max(temp) - min(temp)))
    
    print("gradient calculations")
    #gradient calculations
    # x edges
    
    impad <- padimage(blurs[,,1,5], 3)
    
    imfft <- fft(impad)
    blurredfft <- xkernfft* imfft
    
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    
    xs <- blurredim[((4):(3+hb)),((3+1):(3+wb))]
    
    # y edges
    
    blurredfft <- ykernfft* imfft
    blurredim <- fft(blurredfft, inverse = TRUE)
    blurredim <- Re(blurredim)
    
    ys <- blurredim[((4):(3+hb)),((3+1):(3+wb))]
    
    # assign face five mags, six orients
    blurs[,,5,5] <- xs^2 + ys^2
    blurs[,,6,5] <- atan2(ys, xs) * (180/pi)
    
    rm(gausskernfft)
    rm(gausskern)
    rm(blurredim)
    rm(blurredfft)
    rm(ykernpad)
    rm(ykernfft)
    rm(ys)
    rm(xs)
    rm(impad)
    rm(imfft)
    rm(gausskernpad)
    rm(xkernfft)
    rm(xkernpad)
    rm(piccyras)
    rm(ss)
    rm(tempn)
    
    
    
    # for each DoG plane, find local max and minima, excluding edges
    
    for (i in 1:no_of_octaves)
    {
      #xmax <- big_rows / (1 + (i - 1) * 1.5)
      #ymax <- big_cols / (1 + (i - 1) * 1.5)
      xmax <- big_rows * 1.5/ (1.5^(i- 1))
      ymax <- big_cols * 1.5/ (1.5^(i- 1))
      print(paste0("Finding max and min for octave ", i))
      
      for (x in 3:(xmax - 2))
      {
        for (y in 3:(ymax - 2))
        {
          if (blurs[x,y,4,i] > 0.3)
          {
            if (blurs[x,y,4,i] > blurs[x + 1,y,4,i])
            {
              if (blurs[x,y,4,i] > blurs[x - 1,y,4,i])
              {
                if (blurs[x,y,4,i] > blurs[x,y + 1,4,i])
                {
                  if (blurs[x,y,4,i] > blurs[x,y - 1,4,i])
                  {
                    if (blurs[x,y,4,i] > blurs[x + 1,y + 1,4,i])
                    {
                      if (blurs[x,y,4,i] > blurs[x + 1,y - 1,4,i])
                      {
                        if (blurs[x,y,4,i] > blurs[x - 1,y + 1,4,i])
                        {
                          if (blurs[x,y,4,i] > blurs[x - 1,y - 1,4,i])
                          {
                            lx <- round(x / 1.5)
                            ly <- round(y / 1.5)
                            if ((i == no_of_octaves) ||
                                (blurs[x,y,4,i] > blurs[lx,ly,4,i + 1]))
                            {
                              if ((i == no_of_octaves) ||
                                  ((blurs[x,y,4,i] > blurs[lx + 1,ly,4,i + 1]) &&
                                   (blurs[x,y,4,i] > blurs[lx - 1,ly,4,i + 1])))
                              {
                                if ((i == no_of_octaves) ||
                                    ((blurs[x,y,4,i] > blurs[lx,ly + 1,4,i + 1]) &&
                                     (blurs[x,y,4,i] > blurs[lx,ly - 1,4,i + 1])))
                                {
                                  if ((i == no_of_octaves) ||
                                      ((blurs[x,y,4,i] > blurs[lx + 1,ly + 1,4,i + 1]) &&
                                       (blurs[x,y,4,i] > blurs[lx + 1,ly - 1,4,i + 1])))
                                  {
                                    if ((i == no_of_octaves) ||
                                        ((blurs[x,y,4,i] > blurs[lx - 1,ly + 1,4,i + 1]) &&
                                         (blurs[x,y,4,i] > blurs[lx - 1,ly - 1,4,i + 1])))
                                    {
                                      ux <- round(x * 1.5)
                                      uy <- round(y * 1.5)
                                      if ((i == 1) ||
                                          ((blurs[x,y,4,i] > blurs[ux,uy,4,i - 1]) &&
                                           (blurs[x,y,4,i] > blurs[ux + 1,uy,4,i - 1])))
                                      {
                                        if ((i == 1) ||
                                            ((blurs[x,y,4,i] > blurs[ux - 1,uy,4,i - 1]) &&
                                             (blurs[x,y,4,i] > blurs[ux,uy + 1,4,i - 1])))
                                        {
                                          if ((i == 1) ||
                                              ((blurs[x,y,4,i] > blurs[ux,uy - 1,4,i - 1]) &&
                                               (blurs[x,y,4,i] > blurs[ux + 1,uy + 1,4,i - 1])))
                                          {
                                            if ((i == 1) ||
                                                ((blurs[x,y,4,i] > blurs[ux - 1,uy - 1,4,i - 1]) &&
                                                 (blurs[x,y,4,i] > blurs[ux + 1,uy - 1,4,i - 1])))
                                            {
                                              if ((i == 1) || (blurs[x,y,4,i] > blurs[ux - 1,uy + 1,4,i - 1]))
                                              {
                                                keypoints <- rbind(keypoints,c(x,y,i,blurs[x,y,4,i]))
                                                no_of_keypoints <-
                                                  no_of_keypoints + 1
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    
    for (i in 2:(no_of_keypoints + 1))
    {
      # convert keypoints back to graphable points
      keypoints[i,1] <- (keypoints[i,1] * 1.5 ^ (keypoints[i,3] - 1)) / 2
      keypoints[i,2] <- (keypoints[i,2] * 1.5 ^ (keypoints[i,3] - 1)) / 2
    }
    
    # get rid of ones round edges
    keypoints <- keypoints[which(keypoints[,1] > 10),]
    no_of_keypoints <- nrow(keypoints) - 1
    
    grads <- gradients(pic_grey)
    
    rm(pic_grey)
    
    mags <- grads[,,1]
    orients <- grads[,,2]
    
    # round orientations
    orients <- orients * (180 / pi) # convert to degrees
    orients <- round(orients / 10) * 10
    orients[which(orients == -180)] <- 180
    
    # threshold magnitudes
    mags[which(mags < (max(mags) / 10))] <- 0
    
    canonicals <- 0 # 
    
    for (j in 1:no_of_keypoints)
    {
      # for this one biggest keypoint
      print(paste0(
        "Calculating orientation of keypoint ", j, " of ", no_of_keypoints
      ))
      
      keypoint <- keypoints[j + 1,]
      
      im_patch <- patch(3 * keypoint[3],mags,keypoint[1],keypoint[2])
      temporients <- orients[which(im_patch > 0)]
      tempimpatch <- im_patch[which(im_patch > 0)]
      rm(im_patch)
      
      # histogram
      directions <-
        c(
          -170,-160,-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180
        )
      hist <- matrix(0, 1, length(directions))
      for (i in 1:length(directions))
      {
        direc <- directions[i]
        hist[i] <- sum(tempimpatch[which(temporients == direc)])
      }
      
      if(max(hist) > 0)
      {
        canonical <- directions[which(hist == max(hist))]
      }
      else
      {
        canonical <- 0
      }
      canonicals <- rbind(canonicals, canonical)
      
    }
    
    # add orientations to existing keypoints matrix
    keypoints <- cbind(keypoints, canonicals)
    
    descriptors <- matrix(0, no_of_keypoints, 128)
    
    # for each keypoint, make descriptor
    for (j in 1:no_of_keypoints)
    {
      print(paste0("Calculating descriptor for keypoint ", j, " of ", no_of_keypoints))
      keypoint <- keypoints[j+1,]
      scale <- keypoint[3]
      x_cor <- keypoint[1]*2/(1.5^(scale - 1))
      y_cor <- keypoint[2]*2/(1.5^(scale - 1))
      canonical <- keypoint[5]
      
      x_start <- x_cor - 25
      x_end <- x_cor + 24
      y_start <- y_cor - 25
      y_end <- y_cor + 24
      
      # if any of these run over the edge of the image, forget it
      if ((x_start < 1) || (x_end > (nrow(blurs[,,5,scale])-1)) || (y_start < 1) || (y_end > (ncol(blurs[,,5,scale])-1)))
      {
        descriptors[j,] <- 0
        
      } else {
        ipatch <- matrix(0,50,50)
        patch_mags <- matrix(0,50,50)
        patch_orients <- matrix(0,50,50)
        
        ipatch <- blurs[x_start:x_end,y_start:y_end,2,scale]
        patch_mags <- blurs[x_start:x_end, y_start:y_end, 5 ,scale]
        patch_orients <- blurs[x_start:x_end, y_start:y_end, 6,scale]
        
        # normalise these orientations to the key's orientation
        patch_orients <- patch_orients - canonical
        patch_orients[which(patch_orients < (-170))] <- patch_orients[which((patch_orients) < (-170))] + 360
        patch_orients[which(patch_orients > (180))] <- patch_orients[which((patch_orients) > (180))] - 360
        
        # we want to rotate negatively here - all three patches
        ipatch <- as.matrix(rotate.im(as.im(ipatch), angle = (canonical* -pi/180)))
        patch_mags <- as.matrix(rotate(as.im(patch_mags), angle = (canonical* -pi/180)))
        patch_orients <- as.matrix(rotate(as.im(patch_orients), angle = (canonical* -pi/180)))
        
        patch_orients_blurred <- matrix(0, 16, 16)
        patch_mags_blurred <- matrix(0,16,16)
        
        # blur a smidge
        for (k in 18:33)
        {
          for (l in 18:33)
          {
            patch_orients_blurred[k-17,l-17] <- weighted.mean(patch_orients[((k-1):(k+1)),((l-1):(l+1))],patch_mags[((k-1):(k+1)),((l-1):(l+1))], na.rm = TRUE)
            patch_mags_blurred[k-17,l-17] <- mean(patch_mags[((k-1):(k+1)),((l-1):(l+1))], na.rm = TRUE)
          }
        }
        
        patch_orients_blurred <- round(patch_orients_blurred/45)
        patch_orients_blurred[which(patch_orients_blurred == -4)] <- 4
        
        patch_directions <-
          c(
            -3,-2,-1,0,1,2,3,4
          )
        
        patch_hist <- matrix(0, 1, length(patch_directions))
        
        for (x_quad in 1:4)
        {
          for (y_quad in 1:4)
          {
            for (i in 1:length(patch_directions))
            {
              direc <- patch_directions[i]
              little_patch_mags <- patch_mags_blurred[(((x_quad-1)*4+1):((x_quad-1)*4+4)),(((y_quad - 1)*4+1):((y_quad-1)*4+4))] 
              little_patch_orients <-patch_orients_blurred[(((x_quad-1)*4+1):((x_quad-1)*4+4)),(((y_quad - 1)*4+1):((y_quad-1)*4+4))] 
              patch_hist[i] <- sum(little_patch_mags[which(little_patch_orients == direc)])
            }
            
            descriptors[j,(8*(x_quad-1) + ((y_quad-1)*4*8)+1):(8*(x_quad-1) + ((y_quad-1)*4*8)+8)] <- patch_hist
            
          }
        }
        
        descriptors[j,] <- descriptors[j,]/(sqrt(sum(descriptors[j,] ^ 2)))
        descriptors[j,which(descriptors[j,] < 0.2)] <- 0
        descriptors[j,] <- descriptors[j,]/(sqrt(sum(descriptors[j,] ^ 2)))
        
      }
    }
    
    # get rid of non existent ones
    keypoints <- keypoints[which(rowSums(descriptors) > 0),]
    descriptors <- descriptors[which(rowSums(descriptors)>0),]
    
    outdescriptors <- cbind(keypoints, descriptors)
    
    imtestdata <- as.matrix(outdescriptors)
    
    imtest_descrips <- imtestdata[,(6:133)]
    
    test_clusters <- matrix(0,1,250)
    
    # distance matrix has im test descrips down side, k center descrips across top
    distances <- (rdist(imtest_descrips, toptrain_centers))^2
    
    assignments <- matrix(0, 1, nrow(imtest_descrips))
    
    for (i in 1:nrow(imtest_descrips)) # find min dist and check if valid match
    {
      dist1 <- min(distances[i,])
      index_d1 <- which.min(distances[i,])
      
      assignments[i] <- index_d1
    }
    
    descrips1 <- imtest_descrips[which(assignments == 1),]
    descrips2 <- imtest_descrips[which(assignments == 2),]
    descrips3 <- imtest_descrips[which(assignments == 3),]
    descrips4 <- imtest_descrips[which(assignments == 4),]
    descrips5 <- imtest_descrips[which(assignments == 5),]
    
    distances1 <- (rdist(descrips1, centers1))^2
    distances2 <- (rdist(descrips2, centers1))^2
    distances3 <- (rdist(descrips3, centers1))^2
    distances4 <- (rdist(descrips4, centers1))^2
    distances5 <- (rdist(descrips5, centers1))^2
    
    littletest_clusters <- matrix(0,1,50)
    
    if(length(descrips1) > 128)
    {
      for (i in 1:nrow(descrips1)) # find min dist and check if valid match
      {
        dist1 <- min(distances1[i,])
        index_d1 <- which.min(distances1[i,])
        
        # second smallest
        distances1[i, index_d1] <- max(distances1[i,])
        dist2 <- min(distances1[i,])
        
        if (dist1/dist2 <0.8)
        {
          littletest_clusters[index_d1] <- littletest_clusters[index_d1] +1
        }
      }
    }
    
    test_clusters[,1:50] <- littletest_clusters
    
    
    
    
    littletest_clusters <- matrix(0,1,50)
    
    if((length(descrips2)> 128))
    {
      for (i in 1:nrow(descrips2)) # find min dist and check if valid match
      {
        dist1 <- min(distances2[i,])
        index_d1 <- which.min(distances2[i,])
        
        # second smallest
        distances2[i, index_d1] <- max(distances2[i,])
        dist2 <- min(distances2[i,])
        
        if (dist1/dist2 <0.8)
        {
          littletest_clusters[index_d1] <- littletest_clusters[index_d1] +1
        }
      }
    }
    test_clusters[,51:100] <- littletest_clusters
    
    
    
    littletest_clusters <- matrix(0,1,50)
    if((length(descrips3)>128))
    {
      for (i in 1:nrow(descrips3)) # find min dist and check if valid match
      {
        dist1 <- min(distances3[i,])
        index_d1 <- which.min(distances3[i,])
        
        # second smallest
        distances3[i, index_d1] <- max(distances3[i,])
        dist2 <- min(distances3[i,])
        
        if (dist1/dist2 <0.8)
        {
          littletest_clusters[index_d1] <- littletest_clusters[index_d1] +1
        }
      }
    }
    test_clusters[,101:150] <- littletest_clusters
    
    
    
    littletest_clusters <- matrix(0,1,50)
    
    if((length(descrips4) > 128))
    {
      for (i in 1:nrow(descrips4)) # find min dist and check if valid match
      {
        dist1 <- min(distances4[i,])
        index_d1 <- which.min(distances4[i,])
        
        # second smallest
        distances4[i, index_d1] <- max(distances4[i,])
        dist2 <- min(distances4[i,])
        
        if (dist1/dist2 <0.8)
        {
          littletest_clusters[index_d1] <- littletest_clusters[index_d1] +1
        }
      }
    }
    test_clusters[,151:200] <- littletest_clusters
    
    
    littletest_clusters <- matrix(0,1,50)
    if((length(descrips5) > 128))
    {
      for (i in 1:nrow(descrips5)) # find min dist and check if valid match
      {
        dist1 <- min(distances5[i,])
        index_d1 <- which.min(distances5[i,])
        
        # second smallest
        distances5[i, index_d1] <- max(distances5[i,])
        dist2 <- min(distances5[i,])
        
        if (dist1/dist2 <0.8)
        {
          littletest_clusters[index_d1] <- littletest_clusters[index_d1] +1
        }
      }
    }
    test_clusters[,201:250] <- littletest_clusters
    
    kclusters[v,] <- test_clusters

  }
  
  kclusters <- kclusters[,-removed]
  
  kclusters <- (t(apply(kclusters, 1, function(x) (x - mins)/(maxs - mins))))
  pr.nn <- compute(nn, kclusters)
  percentages <- pr.nn$net.result
  
  colnames(percentages) <- colnames(img_typesexp)
  foldername <- matrix(foldernames[f], nopics,1)
  outputs <- cbind(foldername, picturenames, percentages)
  colnames(outputs)[1] <- "Pole number"
  colnames(outputs)[2] <- "Photo name"
  write.csv(outputs, paste0(output.direc, "\\", foldernames[f], ".csv"))
}