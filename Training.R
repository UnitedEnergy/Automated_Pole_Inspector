
overall.direc <- "F:\\poletop photos"
no_cores <- 3
no_train <- 500

train.direc <- paste0(overall.direc, "\\calculation files\\training")
descriptor.direc <- paste0(overall.direc, "\\calculation files\\descriptors")
output.direc <- paste0(overall.direc, "\\calculation files\\outputs")


filesnames <- list.files(train.direc)

cl <- makeCluster(no_cores, outfile = "debug.txt")
registerDoSNOW(cl)
clusterEvalQ(cl, library(jpeg))
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(spatstat))
clusterExport(cl, "filesnames")
foreach(x = 1:length(filesnames), .export = 'totalfunc') %dopar%
{
  totalfunc(x)
}
stopCluster(cl)


## K CLUSTERING 

names_of_files <- list.files(descriptor.direc)
#write.csv(names_of_files, (paste0(overall.direc, "\\calculation files\\type labels.csv")))

# import image

# import centers
toptrain_centers <- read.csv(paste0(overall.direc, "\\calculation files\\toptrain centers.csv"))
toptrain_centers <- toptrain_centers[,2:129]

centers1 <- read.csv(paste0(overall.direc, "\\calculation files\\centers1.csv"))
centers1 <- centers1[,2:129]

centers2 <- read.csv(paste0(overall.direc, "\\calculation files\\centers2.csv"))
centers2 <- centers2[,2:129]

centers3 <- read.csv(paste0(overall.direc, "\\calculation files\\centers3.csv"))
centers3 <- centers3[,2:129]

centers4 <- read.csv(paste0(overall.direc, "\\calculation files\\centers4.csv"))
centers4 <- centers4[,2:129]

centers5 <- read.csv(paste0(overall.direc, "\\calculation files\\centers5.csv"))
centers5 <- centers5[,2:129]

test_descrips <- matrix(0, 1, 128)

# do first one manually
f <- 1
imtestdata <- read.csv(paste0(descriptor.direc, "\\",names_of_files[f]))

imtestdata <- as.matrix(imtestdata)

imtest_descrips <- imtestdata[,(7:134)]

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

outs <- test_clusters


for (f in 2:length(names_of_files))
{
  imtestdata <- read.csv(paste0(descriptor.direc, "\\",names_of_files[f]))
  
  imtestdata <- as.matrix(imtestdata)
  
  imtest_descrips <- imtestdata[,(7:134)]
  
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
  
  if(length(descrips2)>128)
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
  if(length(descrips3)>128)
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
  if(length(descrips4)>128)
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
  
  if(length(descrips5)>128)
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
  
  outs <- rbind(outs, test_clusters)
  print(f)
}

row.names(outs) <- names_of_files
write.csv(outs, (paste0(overall.direc, "\\calculation files\\trained clusters.csv")))