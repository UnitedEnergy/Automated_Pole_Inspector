
# train the neural net

data <- read.csv(paste0(overall.direc, "\\calculation files\\trained clusters.csv"))

rowsum <- apply(data[2:ncol(data)], 2, sum)
removed <- which(rowsum == 0)
removed <- removed + 1
#if (removed > 0){
data <- data[-removed]
#}

coldata <- ncol(data)

type_lookup <- read.csv(paste0(overall.direc, "\\calculation files\\type labels.csv"))
type_lookup <- type_lookup[, 2:3]

img_labels <- data[,1]
imgdata <- data[,(2:coldata)]

img_types <- matrix(0, 1, length(img_labels))
img_types <- type_lookup[match(img_labels, paste0(type_lookup[,1])),2]


index <- sample(1:nrow(data), no_train)

maxs <- apply(imgdata, 2, max)
mins <- apply(imgdata, 2, min)

scaled <- as.data.frame(t(apply(imgdata, 1, function(x) (x - mins)/(maxs - mins))))

img_typesexp <- class.ind(img_types)

data <- data.frame(img_labels, img_typesexp, scaled)

traindata <- data[index,]
testdata <- data[-index,]

ns <- names(traindata)
n <- ns[(ncol(img_typesexp)+2):length(ns)]
f <- as.formula(paste0(paste(colnames(img_typesexp), collapse = " + ") ,"~", paste(n, collapse = " + ")))
d <- traindata

nn <- neuralnet(f,data = d, hidden = c(102), linear.output = FALSE)

pr.nn <- compute(nn, testdata[,(ncol(img_typesexp) + 2):ncol(testdata)])

actualanswers <- testdata[,(2:(ncol(img_typesexp) + 1))]
guessanswers <- pr.nn$net.result

actualindexes <- apply(actualanswers, 1, which.max)
guessindexes <- apply(guessanswers, 1, which.max)

nocorrect <- length(which(actualindexes == guessindexes))
nototal <- length(actualindexes)

percent <- nocorrect/nototal * 100

rm(actualanswers)
rm(assignments)
rm(d)
rm(data)
rm(descrips1)
rm(descrips2)
rm(descrips3)
rm(descrips4)
rm(descrips5)
rm(distances)
rm(distances1)
rm(distances2)
rm(distances3)
rm(distances4)
rm(distances5)
rm(guessanswers)
rm(imgdata)
rm(imtest_descrips)
rm(imtestdata)
rm(littletest_clusters)
rm(outs)
rm(scaled)
rm(test_clusters)
rm(test_descrips)
rm(testdata)
rm(traindata)
rm(type_lookup)
rm(actualindexes)
rm(coldata)
rm(dist1)
rm(dist2)
rm(f)
rm(guessindexes)
rm(i)
rm(img_labels)
rm(img_types)
rm(index)
rm(index_d1)
rm(n)
rm(names_of_files)
rm(nocorrect)
rm(nototal)
rm(ns)

print(paste0(percent, " percent correct"))

rm(percent)
rm(pr.nn)
rm(rowsum)
