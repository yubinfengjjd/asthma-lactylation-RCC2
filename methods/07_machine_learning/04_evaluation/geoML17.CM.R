
rsFile="model.classMatrix.txt"
method="RF"
setwd("D:\\哮喘乳酸化单细胞开篇\\8.百种机器学习\\混淆矩阵")

riskRT=read.table(rsFile, header=T, sep="\t", check.names=F, row.names=1)
CohortID=gsub("(.*)\\_(.*)\\_(.*)", "\\1", rownames(riskRT))
CohortID=gsub("(.*)\\.(.*)", "\\1", CohortID)
riskRT$Cohort=CohortID

draw_confusion_matrix <- function(cm=null, titleName=null) {
#	layout(matrix(c(1,1,2)))
	par(mar=c(2,2,3,2))
	plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
	title(paste0('CONFUSION MATRIX (', titleName,')'), cex.main=1.5)

	# create the matrix
	rect(150, 430, 240, 370, col='#3F97D0')
	text(195, 435, 'Control', cex=1.2)
	rect(250, 430, 340, 370, col='#F7AD50')
	text(295, 435, 'Treat', cex=1.2)
	text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
	text(245, 450, 'Actual', cex=1.3, font=2)
	rect(150, 305, 240, 365, col='#F7AD50')
	rect(250, 305, 340, 365, col='#3F97D0')
	text(140, 400, 'Control', cex=1.2, srt=90)
	text(140, 335, 'Treat', cex=1.2, srt=90)

	# add in the cm results
	res <- as.numeric(cm)
	text(195, 400, res[1], cex=1.6, font=2, col='white')
	text(195, 335, res[2], cex=1.6, font=2, col='white')
	text(295, 400, res[3], cex=1.6, font=2, col='white')
	text(295, 335, res[4], cex=1.6, font=2, col='white')
}

for(Cohort in unique(riskRT$Cohort)){
	rt=riskRT[riskRT$Cohort==Cohort,]
	y=gsub(".*_(.*)", "\\1", row.names(rt))
	y=ifelse(y=="control", 0, 1)

	result_matrix=table(rt[,method], y)
	pdf(file=paste0("confusion.", Cohort, ".pdf"), width=6, height=5)
	draw_confusion_matrix(cm=result_matrix, titleName=Cohort)
	dev.off()
}
