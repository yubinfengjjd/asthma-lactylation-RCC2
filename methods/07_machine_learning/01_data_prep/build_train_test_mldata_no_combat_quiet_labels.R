library(limma)
library(sva)

geneFile <- "intersectGenes.txt"
setwd("D:\\哮喘乳酸化单细胞开篇\\百种机器学习\\14.MLdata")

files <- dir()
files <- grep("normalize.txt$", files, value = TRUE)
geneList <- list()

for (file in files) {
	rt <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE)
	geneNames <- as.vector(rt[, 1])
	uniqGene <- unique(geneNames)
	header <- unlist(strsplit(file, "\\.|\\-"))
	geneList[[header[1]]] <- uniqGene
}

interGenes <- Reduce(intersect, geneList)

allTab <- data.frame()
batchType <- c()

for (i in seq_along(files)) {
	inputFile <- files[i]
	header <- unlist(strsplit(inputFile, "\\.|\\-"))

	rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
	rt <- as.matrix(rt)
	rownames(rt) <- rt[, 1]
	exp <- rt[, -1, drop = FALSE]

	dimnames <- list(rownames(exp), colnames(exp))
	data <- matrix(as.numeric(exp), nrow = nrow(exp), dimnames = dimnames)
	rt <- avereps(data)

	colnames(rt) <- paste0(header[1], "_", colnames(rt))

	if (i == 1) {
		allTab <- rt[interGenes, ]
	} else {
		allTab <- cbind(allTab, rt[interGenes, ])
	}

	batchType <- c(batchType, rep(i, ncol(rt)))
}

svaTab = allTab
geneRT <- read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
geneTab <- svaTab[intersect(rownames(svaTab), geneRT[, 1]), ]
geneTab <- t(geneTab)

merge_samples <- grepl("^merge", rownames(geneTab), ignore.case = TRUE)
merge_data <- geneTab[merge_samples, ]

set.seed(123)
train_indices <- sample(1:nrow(merge_data), size = round(0.7 * nrow(merge_data)))
test_indices <- setdiff(1:nrow(merge_data), train_indices)

trainExp <- merge_data[train_indices, , drop = FALSE]
testExp <- merge_data[test_indices, , drop = FALSE]


train <- grepl("^merge", rownames(geneTab), ignore.case = TRUE)
trainExp <- geneTab[train, , drop = FALSE]
testExp <- geneTab[!train, , drop = FALSE]

rownames(trainExp) <- gsub("merge_", "Train.", rownames(trainExp))
trainType <- gsub("(.*)\\_(.*)", "\\2", rownames(trainExp))
testType <- gsub("(.*)\\_(.*)", "\\2", rownames(testExp))
trainType <- ifelse(trainType == "QUIET", 0, 1)
testType <- ifelse(testType == "QUIET", 0, 1)

trainExp <- cbind(trainExp, Type = trainType)
testExp <- cbind(testExp, Type = testType)

trainOut <- rbind(id = colnames(trainExp), trainExp)
write.table(trainOut, file = "data.train.txt", sep = "\t", quote = FALSE, col.names = FALSE)
testOut <- rbind(id = colnames(testExp), testExp)
write.table(testOut, file = "data.test.txt", sep = "\t", quote = FALSE, col.names = FALSE)
