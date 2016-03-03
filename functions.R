library(raster)
library(spatstat)
library(jpeg)
library(parallel)
library(foreach)
library(doSNOW)
library(fields)
library(neuralnet)
library(nnet)

padimage <- function(I, p)
{
  hp <- nrow(I)
  wp <- ncol(I)
  Ipad <- matrix(0, hp + 2*p, wp + 2*p)
  
  #middle
  Ipad[((p+1):(p+hp)), ((p+1):(p+wp))] = I
  
  #top and bottom
  Ipad[(1:p), ((p+1):(p+wp))] = (matrix(I[1,],p,wp,byrow=TRUE))
  Ipad[((p+hp+1):(nrow(Ipad))), ((p+1):(p+wp))] = matrix(I[hp,],p,wp,byrow=TRUE)
  
  #left and right
  Ipad[((p+1):(p+hp)), (1:p)] = matrix(I[,1],hp,p)
  Ipad[((p+1):(p+hp)), ((p+wp+1):(ncol(Ipad)))] = matrix(I[,wp],hp,p) 
  
  #Corners
  Ipad[(1:p), (1:p)] = I[1,1] #Top-left
  Ipad[(1:p), ((p+wp+1):(ncol(Ipad)))] = I[1,wp]; #Top-right
  Ipad[((p+hp+1):(nrow(Ipad))), 1:p] = I[hp,1]; #Bottom-left
  Ipad[((p+hp+1):(nrow(Ipad))),((p+wp+1):(ncol(Ipad)))] = I[hp,wp]; #Bottom-right
  
  return(Ipad)
}

blur_pic <- function(pic, sigma, ksize)
{
  #ksize = 51
  kern <- matrix(0, nrow = ksize, ncol = ksize)
  #S <- 2 # pick a sigma
  m = ksize/2
  print('hi')
  
  hb <- nrow(pic)
  wb <- ncol(pic)
  
  # make X and Y matricies
  X <- matrix(1:ksize,ksize,ksize, byrow = TRUE)
  Y <- matrix(1:ksize,ksize,ksize, byrow = FALSE)
  
  print('sldkfj')
  kern <- (1/(2*pi*sigma^2))*exp(-((X-m)^2 + (Y-m)^2)/(2*sigma^2))
  
  impad <- padimage(pic, ksize)
  
  kernpad <- matrix(0, nrow(impad), ncol(impad))
  
  kernpad[1:ksize,1:ksize] <- kern
  print("3")
  kernfft <- fft(kernpad)
  imfft <- fft(impad)
  blurredfft <- kernfft* imfft
  print("7")
  blurredim <- fft(blurredfft, inverse = TRUE)
  blurredim <- Re(blurredim)
  print("9")
  
  blurredimn <- blurredim[((ksize+1):(ksize+hb)),((ksize+1):(ksize+wb))]
  
  print("4")
  return(blurredimn)
}


pixsample <- function(x,y)
{
  x1 <- floor(x)
  x2 <- ceiling(x)
  y1 <- floor(y)
  y2 <- ceiling(y)
  
  fq11 <- pic_grey[x1,y1]
  fq21 <- pic_grey[x2,y1]
  fq12 <- pic_grey[x1,y2]
  fq22 <- pic_grey[x2,y2]
  
  if(x1 == x2)
  {
    fxy1 <- fq11
    fxy2 <- fq12
  }
  else
  {
    fxy1 <- ((x2 - x)/(x2 - x1))*fq11 + ((x-x1)/(x2-x1))*fq21
    fxy2 <- ((x2-x)/(x2-x1))*fq12 + ((x-x1)/(x2-x1))*fq22
  }
  
  if(y1 == y2)
  {
    fxy <- fxy1
  }
  
  else
  {
    fxy <- ((y2-y)/(y2-y1))*fxy1+((y-y1)/(y2-y1))*fxy2
  }
  
  return(fxy)
}



patch <- function(sigma, pic, x_cor, y_cor)
{
  ksize <- ceiling(sigma *8.5)
  
  kernel <- matrix(0, 1000, 1000)
  
  Y <- matrix((1:ksize), ksize, ksize)
  X <- matrix((1:ksize), ksize, ksize, byrow = TRUE)
  
  m <- ksize/2
  
  kernel <- (1/(2*pi*sigma^2)) * exp(-((X-m)^2+(Y-m)^2)/(2*sigma^2))
  
  # put in bigger matrix
  
  big_kern <- matrix(0, (nrow(pic) + ksize*2), (ncol(pic) + ksize*2))
  
  
  big_kern[((x_cor + ksize - round(ksize/2)) : (x_cor +ksize - round(ksize/2) + ksize-1)),((y_cor + ksize - round(ksize/2)) : (y_cor + ksize - round(ksize/2) + ksize-1))] <- kernel
  
  big_kern <- big_kern[(ksize:(ksize+nrow(pic)-1)),(ksize:(ksize+ncol(pic)-1))]
  
  pic_patch <- big_kern * pic
  
  return(pic_patch)
  
}



x_edge <- function(pic)
{
  kern <- matrix(c(-1,-2,-1,0,0,0,1,2,1), 3,3)
  
  hb <- nrow(pic)
  wb <- ncol(pic)
  
  impad <- padimage(pic, 3)
  
  kernpad <- matrix(0, nrow(impad), ncol(impad))
  
  kernpad[1:3,1:3] <- kern
  print("3")
  kernfft <- fft(kernpad)
  imfft <- fft(impad)
  blurredfft <- kernfft* imfft
  print("7")
  blurredim <- fft(blurredfft, inverse = TRUE)
  blurredim <- Re(blurredim)
  print("9")
  
  blurredimn <- blurredim[((4):(3+hb)),((3+1):(3+wb))]
  
  print("4")
  return(blurredimn)
}

y_edge <- function(pic)
{
  kern <- matrix(c(-1,0,1,-2,0,2,-1,0,1), 3,3)
  
  hb <- nrow(pic)
  wb <- ncol(pic)
  
  impad <- padimage(pic, 3)
  
  kernpad <- matrix(0, nrow(impad), ncol(impad))
  
  kernpad[1:3,1:3] <- kern
  print("3")
  kernfft <- fft(kernpad)
  imfft <- fft(impad)
  blurredfft <- kernfft* imfft
  print("7")
  blurredim <- fft(blurredfft, inverse = TRUE)
  blurredim <- Re(blurredim)
  print("9")
  
  blurredimn <- blurredim[((4):(3+hb)),((3+1):(3+wb))]
  
  print("4")
  return(blurredimn)
}

gradients <- function(pic)
{
  xs <- x_edge(pic)
  ys <- y_edge(pic)
  
  grads <- array(0,c(nrow(pic), ncol(pic), 2))
  
  grads[,,1] <- xs^2 + ys^2
  grads[,,2] <- atan2(ys, xs)
  
  return(grads)
}
