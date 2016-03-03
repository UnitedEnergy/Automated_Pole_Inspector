foldernames <- list.files(overall.direc)
foldernames <- foldernames[which(foldernames != "calculation files")]
doneones <- list.files(output.direc)
doneones <- substring(doneones,1,nchar(doneones)-4)
needtobedone <- foldernames[! foldernames %in% doneones]

cl <- makeCluster(no_cores, outfile = "debug.txt")
registerDoSNOW(cl)
clusterEvalQ(cl, library(jpeg))
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(spatstat))
clusterEvalQ(cl, library(fields))
clusterEvalQ(cl, library(neuralnet))
clusterExport(cl, "needtobedone")
foreach(x = 1:400, .export = 'folderfunc') %dopar%
{
  folderfunc(x)
}
stopCluster(cl)
