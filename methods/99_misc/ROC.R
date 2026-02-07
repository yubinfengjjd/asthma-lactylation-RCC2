library(pROC)
library(scales)

setwd("D:\\哮喘乳酸化单细胞开篇\\11.ROC")
expFiles <- c("./GSE43696.normalize.txt", "./GSE63142.normalize.txt", "./GSE67940.normalize.txt")
if (!dir.exists("ROC")) dir.create("ROC")

geneRT <- "RCC2"
bioCol <- c("#d92b50")  # 专业配色

for (expFile in expFiles) {
	rt <- read.table(file.path("nor", expFile), header=TRUE, sep="\t",
									 check.names=FALSE, row.names=1)
	y <- ifelse(gsub(".*_", "", colnames(rt)) == "control", 0, 1)
	title <- gsub(".normalize.txt", "", basename(expFile))

	pdf(file.path("ROC", paste0(title, "_ROC.pdf")), width=7, height=6)

	par(mar = c(4.5, 4.5, 2, 2), mgp = c(2.5, 1, 0))

	aucText <- c()
	for (i in seq_along(geneRT)) {
		gene <- geneRT[i]
		roc_obj <- roc(y, as.numeric(rt[gene,]))

		plot(roc_obj,
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
		auc_ci <- ci.auc(roc_obj, method="bootstrap", boot.n=2000)
		sp_ci <- ci.sp(roc_obj, sensitivities=seq(0,1,0.01), boot.n=2000, progress="none")

		plot(sp_ci,
				 type = "shape",
				 col = alpha("#f3bbc6", 0.3), # 使用透明度
				 border = NA)

		auc_text <- sprintf("AUC = %.3f (95%% CI: %.3f-%.3f)",
												auc(roc_obj),
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
	dev.off()}
}
