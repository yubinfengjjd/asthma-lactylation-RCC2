library(limma)
library(sva)

geneFile <- "intersectGenes.txt"
setwd("D:\\哮喘乳酸化单细胞开篇\\9.GMM与LR")

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

geneRT <- read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
geneTab <- allTab[intersect(rownames(allTab), geneRT[, 1]), ]
geneTab <- t(geneTab)

train <- grepl("^merge", rownames(geneTab), ignore.case = TRUE)
trainExp <- geneTab[train, , drop = FALSE]
testExp <- geneTab[!train, , drop = FALSE]

rownames(trainExp) <- gsub("merge_", "Train.", rownames(trainExp))
allType <- gsub("(.*)\\_(.*)", "\\2", rownames(geneTab))
trainType <- gsub("(.*)\\_(.*)", "\\2", rownames(trainExp))
testType <- gsub("(.*)\\_(.*)", "\\2", rownames(testExp))
allType <- ifelse(allType == "control", 0, 1)
trainType <- ifelse(trainType == "control", 0, 1)
testType <- ifelse(testType == "control", 0, 1)

allExp <- cbind(geneTab, Type = allType)
trainExp <- cbind(trainExp, Type = trainType)
testExp <- cbind(testExp, Type = testType)

allOut <- rbind(id = colnames(allExp), allExp)
write.table(allOut, file = "data.all.txt", sep = "\t", quote = FALSE, col.names = FALSE)
trainOut <- rbind(id = colnames(trainExp), trainExp)
write.table(trainOut, file = "data.train.txt", sep = "\t", quote = FALSE, col.names = FALSE)
testOut <- rbind(id = colnames(testExp), testExp)
write.table(testOut, file = "data.test.txt", sep = "\t", quote = FALSE, col.names = FALSE)

#####GMM_LR####
library(caret)
library(pROC)
library(doParallel)
setwd("D:\\哮喘乳酸化单细胞开篇\\9.GMM与LR")
rm(list = ls())
data_file <- "data.train.txt"
if (!file.exists(data_file)) {
	stop("Data file does not exist!")
}
surv.expr <- read.table(data_file, sep = "\t", row.names = 1,
												stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
colnames(surv.expr) <- gsub("-", "_", colnames(surv.expr))

response_variable <- "Type"

surv.expr[[response_variable]] <- as.factor(surv.expr[[response_variable]])
levels(surv.expr[[response_variable]]) <- make.names(levels(surv.expr[[response_variable]]))

cat("修复后的响应变量水平:\n")
print(levels(surv.expr[[response_variable]]))

ctrl <- trainControl(
	method = "cv",
	number = 5,
	savePredictions = "final",
	classProbs = TRUE,
	summaryFunction = twoClassSummary,
	allowParallel = TRUE
)
predictors <- colnames(surv.expr)[colnames(surv.expr) != response_variable]
response <- surv.expr[[response_variable]]

if (length(levels(surv.expr[[response_variable]])) != 2) {
	stop("响应变量必须是二分类因子")
}

generate_combinations <- function(predictors) {
	n <- length(predictors)
	all_combinations <- unlist(
		lapply(1:n, function(k) combn(predictors, k, simplify = FALSE)),
		recursive = FALSE
	)
	return(all_combinations)
}

feature_combinations <- generate_combinations(predictors)
n.model <- length(feature_combinations)

display.progress <- function(index, totalN, breakN = 20) {
	if (index %% ceiling(totalN / breakN) == 0) {
		cat(paste(round(index * 100 / totalN), "% ", sep = ""))
	}
}

run_models_caret <- function(feature_combinations, data, response_var) {
	aucDF <- data.frame(model = integer(), auc = numeric(),
											n.mRNA = integer(), stringsAsFactors = FALSE)
	formula_list <- list()
	model_list <- list()

	cl <- makePSOCKcluster(detectCores() - 1)
	registerDoParallel(cl)
	on.exit(stopCluster(cl))

	for (i in 1:length(feature_combinations)) {
		if(exists("display.progress") && exists("n.model")) {
			display.progress(i, n.model)
		}

		current_features <- feature_combinations[[i]]
		formula <- as.formula(paste(response_var, "~", paste(current_features, collapse = " + ")))
		formula_list[[i]] <- formula

		tryCatch({
			set.seed(123)
			model <- train(
				formula,
				data = data,
				method = "glm",
				family = "binomial",
				trControl = ctrl,
				metric = "ROC"
			)

			auc_value <- max(model$results$ROC)
			model_list[[i]] <- model

			aucDF <- rbind(aucDF, data.frame(
				model = i,
				auc = auc_value,
				n.mRNA = length(current_features)
			))
		}, error = function(e) {
			message("\n模型 ", i, " 出错: ", e$message)
			model_list[[i]] <- NULL

			aucDF <- rbind(aucDF, data.frame(
				model = i,
				auc = NA,
				n.mRNA = length(current_features)
			))
		})
	}

	return(list(
		aucDF = aucDF,
		formula_list = formula_list,
		model_list = model_list
	))
}

results <- run_models_caret(feature_combinations, surv.expr, response_variable)

aucDF <- results$aucDF
formula_list <- results$formula_list
model_list <- results$model_list

library(pROC)
library(mclust)
library(lares)
mod <- Mclust(aucDF[, 1:2])
aucDF$mclust <- as.numeric(mod$classification)
table(aucDF$mclust)

write.table(aucDF, "all_models_with_GMM_clusters.txt", sep = "\t", row.names = FALSE, quote = FALSE)

col.cluster <- c('#1C79B7', '#F38329', '#2DA248', '#DC403E',
								 '#976BA6', '#8D574C',"#DC4055","#F38367","#2DA788")

pdf("GMM_cluster_of_all_auc_model_combination.pdf", width = 7, height = 6)
par(bty = "o", mgp = c(2, 0.5, 0), mar = c(3.1, 4.1, 1.1, 1.1), tcl = -.25, las = 1)

plot(mod, what = "classification",
		 colors = col.cluster[1:mod$G],
		 xlab = paste0("Sorted logistic regression models (1:", n.model, ")"),
		 ylab = "")
mtext("Training models", side = 2, line = 2.5, las = 3)

legend(
	"bottomright",
	legend = paste("Cluster", 1:mod$G),
	col = col.cluster[1:mod$G],
	pch = 16,
	pt.cex = 1.5,
	cex = 0.8,
	bty = "n",
	title = "GMM Clusters"
)

dev.off()

library(RColorBrewer)
library(ggplot2)
library(ggridges)
library(RColorBrewer)
aucDF$auc <- gsub("Area under the curve: ", "", aucDF$auc) %>%
	as.numeric()

str(aucDF$auc)
aucDF$mclust = as.factor(aucDF$mclust)
pdf("custom_ridgeline.pdf",width = 6,height = 8)
ggplot(aucDF, aes(x = auc, y = mclust, fill = ..density..)) +
	geom_density_ridges_gradient(scale = 1, rel_min_height = 0, size = 1) +
	scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,"Spectral")))(32)) +
	theme(
		axis.line = element_blank(),
		axis.ticks = element_blank(),
		panel.border = element_rect(color = "black", size = 1.5, fill = NA),
		panel.grid = element_blank()
	)
dev.off()

library(ggplot2)
library(ggridges)
library(patchwork)

plot_cluster_auc_with_density <- function(aucDF, clusters, outdir = ".") {
	if (!dir.exists(outdir)) dir.create(outdir)

	for (cluster_num in clusters) {
		cluster_data <- subset(aucDF, mclust == cluster_num)
		cluster_data$auc <- as.numeric(gsub("Area under the curve: ", "", cluster_data$auc))
		cluster_data <- cluster_data[order(-cluster_data$auc), ]
		cluster_data$rank <- 1:nrow(cluster_data)

		p_auc = ggplot(cluster_data, aes(x = auc, y = rank)) +
			geom_segment(
				aes(
					x = min(cluster_data$auc),
					xend = max(cluster_data$auc),
					y = 0,
					yend = 0
				),size = 2) +
			geom_ribbon(
				aes(ymin = 0, ymax = rank),
				fill = "#ffe7c5",
				alpha = 0.5
			) +
			geom_smooth(se = FALSE, color = "#eba85b", size = 2.5, span = 0.4) +
			geom_point(
				data = cluster_data[which.max(cluster_data$auc), ],
				aes(x = auc, y = 0),
				color = "#a3570b",
				size = 6
			) +
			annotate(
				"text",
				x = cluster_data[which.max(cluster_data$auc), "auc"],
				y = cluster_data[which.max(cluster_data$auc), "rank"]+0.03*nrow(cluster_data),
				label = round(max(cluster_data$auc), 4),
				color = "#a0570d",
				size = 6
			) +
			scale_y_continuous(
				trans = "reverse",
				expand = expansion(mult = c(0.05, 0.05))
			) +
			scale_x_continuous(
				limits = c(min(cluster_data$auc), max(cluster_data$auc) + 0.005),
				expand = c(0, 0)
			) +
			labs(
				x = "AUC Value",
				y = "Model Rank (Descending)",
				title = paste0("AUC Ranking of Cluster ",cluster_num," Models"),
				subtitle = paste0(nrow(cluster_data),"Models Ordered by AUC (Highest to Lowest)")
			) +
			theme_minimal() +
			theme(
				plot.background = element_blank(),
				panel.grid = element_blank(),
				panel.border = element_rect(color = "black", size = 2, fill = NA),
				axis.text = element_text(color = "black"),
				plot.title = element_text(hjust = 0.5, face = "bold"),
				plot.subtitle = element_text(hjust = 0.5, color = "grey40")
			)

		dens <- density(cluster_data$auc)
		ridge_df <- data.frame(x = dens$x, height = dens$y, y = 1)

		p_ridge <- ggplot(ridge_df, aes(x = x, y = y, height = height)) +
			geom_ridgeline(fill = "grey80", color = "grey40", alpha = 0.7, scale = 0.9) +
			annotate("text",
							 x = min(cluster_data$auc), y = 1.2,
							 label = paste0(nrow(cluster_data), " models"),
							 hjust = 0, size = 4.5, fontface = "bold"
			) +
			theme_void() +
			theme(
				plot.margin = margin(0, 0, 0, 0),
				panel.border = element_blank()
			)

		final_plot <- p_ridge + p_auc + plot_layout(ncol = 2, widths = c(2, 4))

		ggsave(
			filename = file.path(outdir, paste0("cluster_auc_C", cluster_num, ".pdf")),
			plot = final_plot,
			width = 6, height = 8
		)
	}
}
plot_cluster_auc_with_density(aucDF, clusters = c(6), outdir = "cluster_auc_plots")

cluster4_models <- subset(aucDF, mclust == 6)
best_model_row <- cluster4_models[which.max(cluster4_models$auc), ]
best_model_index <- best_model_row$model
best_model_formula <- formula_list[[best_model_index]]
all_vars <- all.vars(best_model_formula)
model_mRNA <- setdiff(all_vars, response_variable)
model_mRNA = gsub("_","-",model_mRNA)
write.table(model_mRNA, "output_mRNA_selected.txt", row.names = FALSE, col.names = T,sep = "\t",quote = F)

final_model = model_list[[best_model_index]]

#####ROC####
train_prob <- predict(final_model, newdata = surv.expr, type = "prob")
positive_class_prob <- train_prob[["X1"]]
library(pROC)

roc_obj <- roc(
	response = surv.expr[[response_variable]],
	predictor = positive_class_prob
)

set.seed(2025)
auc_ci <- ci.auc(roc_obj, method = "bootstrap", boot.n = 2000)

sp_ci <- ci.sp(roc_obj, sensitivities = seq(0, 1, 0.01), boot.n = 2000, progress = "none")

pdf("Best_Model_Training_ROC.pdf", width = 7, height = 6)
par(mar = c(4.5, 4.5, 2, 2), mgp = c(2.5, 1, 0))

plot(roc_obj,
		 print.auc = FALSE,
		 grid = c(0.2, 0.2),
		 grid.col = "lightgray",
		 legacy.axes = TRUE,
		 xlab = "False Positive Rate (1 - Specificity)",
		 ylab = "True Positive Rate (Sensitivity)",
		 main = "ROC Curve For Train Set",
		 col = "#2c7fb8",
		 lwd = 2)

plot(sp_ci,
		 type = "shape",
		 col = alpha("#5bd1d7", 0.3), # 使用透明度
		 border = NA)

auc_text <- sprintf("AUC = %.3f (95%% CI: %.3f-%.3f)",
										auc(roc_obj),
										auc_ci[1],
										auc_ci[3])
text(x = 0.6, y = 0.2,
		 labels = auc_text,
		 col = "#2c7fb8",
		 cex = 1.2)


legend("bottomright",
			 legend = c("ROC Curve", "95% Confidence Band"),
			 col = c("#2c7fb8", alpha("#5bd1d7", 0.3)),
			 lwd = c(2, 10),
			 bty = "n",
			 cex = 0.9)

dev.off()

test.expr <- read.table("./data.test.txt", sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
colnames(test.expr) = gsub("-","_",colnames(test.expr))
test_prob <- predict(final_model, type = "prob",newdata = test.expr)
auc_test <- auc(response = test.expr$Type, predictor = test_prob[["X1"]])
library(pROC)

roc_obj <- roc(
	response = test.expr[[response_variable]],
	predictor = test_prob[["X1"]]
)

set.seed(123)
auc_ci <- ci.auc(roc_obj, method = "bootstrap", boot.n = 2000)

sp_ci <- ci.sp(roc_obj, sensitivities = seq(0, 1, 0.01), boot.n = 2000, progress = "none")

pdf("Best_Model_Testing_ROC.pdf", width = 7, height = 6)
par(mar = c(4.5, 4.5, 2, 2), mgp = c(2.5, 1, 0))

plot(roc_obj,
		 print.auc = FALSE,
		 grid = c(0.2, 0.2),
		 grid.col = "lightgray",
		 legacy.axes = TRUE,
		 xlab = "False Positive Rate (1 - Specificity)",
		 ylab = "True Positive Rate (Sensitivity)",
		 main = "ROC Curve For Test Set",
		 col = "#2c7fb8",
		 lwd = 2)

plot(sp_ci,
		 type = "shape",
		 col = alpha("#5bd1d7", 0.3), # 使用透明度
		 border = NA)

auc_text <- sprintf("AUC = %.3f (95%% CI: %.3f-%.3f)",
										auc(roc_obj),
										auc_ci[1],
										auc_ci[3])
text(x = 0.6, y = 0.2,
		 labels = auc_text,
		 col = "#2c7fb8",
		 cex = 1.2)


legend("bottomright",
			 legend = c("ROC Curve", "95% Confidence Band"),
			 col = c("#2c7fb8", alpha("#5bd1d7", 0.3)),
			 lwd = c(2, 10),
			 bty = "n",
			 cex = 0.9)

dev.off()

#####SHAP####
library(caret)
library(pROC)
library(ggplot2)
library(shapviz)
library(kernelshap)

pred_fun <- function(model, newdata) {
	prob_matrix <- predict(model, newdata = newdata, type = "prob")
	return(as.numeric(prob_matrix[, 2]))
}

X <- as.data.frame(test.expr[, model_mRNA])
X <- X[, colnames(X) != "(Intercept)"]

shap_values <- kernelshap(
	object = final_model,
	X = X,
	pred_fun = pred_fun
)

X <- X[, colnames(shap_values$S)]
shap_vis <- shapviz(shap_values, X = X)

feature_importance <- colMeans(abs(shap_values$S))
sorted_features <- names(sort(feature_importance, decreasing = TRUE))

visualization_theme <- theme_minimal() +
	theme(
		plot.title = element_text(face = "bold", size = 14),
		axis.title = element_text(size = 12)
	)

pdf("2_SHAP_Feature_Importance_Barplot.pdf", width = 8, height = 6)
sv_importance(shap_vis, kind = "bar", show_numbers = TRUE) +
	visualization_theme +
	labs(
		title = "Feature Importance (Mean Absolute SHAP)",
		subtitle = "Average impact magnitude of each feature on model predictions",
		x = "Mean |SHAP value|",
		y = "Feature",
		caption = "Bar height indicates feature importance, with value showing mean absolute SHAP"
	)
dev.off()

pdf("3_SHAP_BeeSwarm_Plot.pdf", width = 9, height = 7)
sv_importance(shap_vis, kind = "bee", show_numbers = TRUE) +
	visualization_theme +
	labs(
		title = "SHAP Value Distribution (Bee Swarm)",
		subtitle = "Each point represents one sample, color indicates feature value",
		x = "SHAP value (impact on model output)",
		y = "Feature",
		caption = "Red: high feature values | Blue: low feature values | Horizontal spread: direction of effect"
	)
dev.off()

for (feat in sorted_features[1:3]) {
	pdf(paste0("4_SHAP_Feature_Dependence_fixe_",feat,".PDF"), width = 6, height = 6)
	p <- sv_dependence(
		object = shap_vis,
		feature = feat,
		v = feat,
		jitter_width = 0.2
	) +
		labs(
			title = paste("SHAP Dependence Plot:", feat),
			color = feat
		)
	print(p)
	dev.off()
}

pdf("5_SHAP_Sample_Waterfall_treat.pdf", width=9, height=6)
sv_waterfall(shap_vis, row_id=106) +  # Explaining sample #2
	labs(
		title = "Prediction Breakdown (Waterfall Plot)",
		subtitle = "Cumulative contribution of features to final prediction",
		caption = "E[f(x)] = base value | f(x) = model output\nBar length: feature contribution | Direction: sign of effect"
	)
dev.off()

pdf("5_SHAP_Sample_Waterfall_control.pdf", width=9, height=6)
sv_waterfall(shap_vis, row_id=1) +  # Explaining sample #2
	labs(
		title = "Prediction Breakdown (Waterfall Plot)",
		subtitle = "Cumulative contribution of features to final prediction",
		caption = "E[f(x)] = base value | f(x) = model output\nBar length: feature contribution | Direction: sign of effect"
	)
dev.off()
pdf("6_SHAP_Force_Plot_treat.pdf", width=10, height=6)
sv_force(shap_vis, row_id=106) +  # Analyzing sample #2
	labs(
		title = "Feature Contributions (Force Plot)",
		subtitle = "Visual forces pushing prediction from base value to output"
	)
dev.off()

pdf("6_SHAP_Force_Plot_control.pdf", width=10, height=6)
sv_force(shap_vis, row_id=1) +  # Analyzing sample #2
	labs(
		title = "Feature Contributions (Force Plot)",
		subtitle = "Visual forces pushing prediction from base value to output"
	)
dev.off()

#####LIME####
library(lime)
library(caret)
exp <- lime(test.expr[, model_mRNA], final_model)
explaintion <- lime::explain(test.expr[, model_mRNA][1, , drop = FALSE],
														 explainer = exp,
														 n_labels = 2,
														 n_features = 6)
pdf("6_LIME_Plot_control.pdf", width=12, height=6)
plot_features(explaintion)
dev.off()

explaintion <- lime::explain(test.expr[, model_mRNA][106, , drop = FALSE],
														 explainer = exp,
														 n_labels = 2,
														 n_features = 6)
pdf("6_LIME_Plot_treat.pdf", width=12, height=6)
plot_features(explaintion)
dev.off()
