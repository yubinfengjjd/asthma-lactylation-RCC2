#install.packages("pROC")


library(pROC)

rsFile="model.riskMatrix.txt"
method="RF"
setwd("D:\\哮喘乳酸化单细胞开篇\\8.百种机器学习\\模型ROC")

riskRT=read.table(rsFile, header=T, sep="\t", check.names=F, row.names=1)
CohortID=gsub("(.*)\\_(.*)\\_(.*)", "\\1", rownames(riskRT))
CohortID=gsub("(.*)\\.(.*)", "\\1", CohortID)
riskRT$Cohort=CohortID

for(Cohort in unique(riskRT$Cohort)){
	rt=riskRT[riskRT$Cohort==Cohort,]
	y=gsub(".*_(.*)", "\\1", row.names(rt))
	y=ifelse(y=="control", 0, 1)
	pdf(paste0(Cohort, "_ROC.pdf"), width=7, height=6)
	roc1=roc(y, as.numeric(rt[,method]))
	plot(roc1,
			 print.auc = FALSE,
			 grid = c(0.2, 0.2),
			 grid.col = "lightgray",
			 legacy.axes = TRUE,
			 xlab = "False Positive Rate (1 - Specificity)",
			 ylab = "True Positive Rate (Sensitivity)",
			 main = "ROC Curve with 95% Confidence Intervals",
			 col = "#d92b50",
			 lwd = 2)
	set.seed(123)
	auc_ci <- ci.auc(roc1, method="bootstrap", boot.n=2000)
	sp_ci <- ci.sp(roc1, sensitivities=seq(0,1,0.01), boot.n=2000, progress="none")

	plot(sp_ci,
			 type = "shape",
			 col = alpha("#f3bbc6", 0.3), # 使用透明度
			 border = NA)

	auc_text <- sprintf("AUC = %.3f (95%% CI: %.3f-%.3f)",
											auc(roc1),
											auc_ci[1],
											auc_ci[3])
	text(x = 0.6, y = 0.2,
			 labels = auc_text,
			 col = "#d92b50",
			 cex = 1.2)


	legend("bottomright",
				 legend = c("ROC Curve", "95% Confidence Band"),
				 col = c("#d92b50", alpha("#f3bbc6", 0.3)),
				 lwd = c(2, 10),
				 bty = "n",
				 cex = 0.9)
	dev.off()
}
