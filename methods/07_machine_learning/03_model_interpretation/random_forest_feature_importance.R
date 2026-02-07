library(limma)
library(randomForest)
library(caret)
library(stringr)
setwd("D:\\哮喘乳酸化单细胞开篇\\8.百种机器学习\\随机森林")
set.seed(123)

rt=read.table("./merge.normalize.txt", header=T, sep="\t", check.names=F, row.names=1)
genes=read.table("./intersectGenes.txt", header=T, sep="\t", check.names=F)[,1]
samegene=intersect(genes,rownames(rt))
data=t(rt)[,samegene]
colnames(data)=str_replace_all(colnames(data),"-","_")
Type=gsub("(.*)\\_(.*)", "\\2", rownames(data))
y=Type
x=as.matrix(cbind(y=y,data))
rf=randomForest(as.factor(y)~., data=x, ntree=1000,nodesize = 5,importance = T,
								proximity = T,
								forest = T)
pdf(file="forest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=4)
dev.off()
optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(y)~., data=x, ntree=optionTrees,importance = T,
								 proximity = T,
								 forest = T)

importance=importance(x=rf2)

pdf(file="geneImportance.pdf", width=6, height=6)
par(mar = c(5, 8, 4, 2), mgp = c(3, 0.5, 0))

imp_scores <- importance[, "MeanDecreaseAccuracy"]
imp_sorted <- sort(imp_scores, decreasing = F)

bar_col <- "#087687"  # 深蓝色（来自RColorBrewer的Blues调色板）

barplot(imp_sorted,
				horiz = TRUE,
				las = 1,
				col = bar_col,
				border = NA,
				space = 0.5,
				xlab = "Mean Decrease Accuracy",
				main = "Feature Importance Ranking\n(Random Forest)",
				cex.names = 0.7,
				xlim = c(-3, max(imp_sorted)*1.1))

grid(nx = NA, ny = NULL,
		 col = "gray90", lty = 2)

legend("bottomright",
			 legend = "Mean Decrease Accuracy",
			 fill = bar_col,
			 bty = "n",
			 cex = 0.8)

dev.off()
