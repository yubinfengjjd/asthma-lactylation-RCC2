library(dplyr)
library(pROC)
library(ggplot2)
library(KEGGREST)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(limma)
library(readxl)
library(caret)
library(sva)
rm(list = ls())
setwd("D:\\哮喘乳酸化单细胞开篇\\9.GMM与LR\\模型比对")
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

#allTab <- ComBat(allTab,batch = batchType)
allTab = allTab
samples <- colnames(allTab)
dataset_prefixes <- sapply(strsplit(samples, "_"), `[`, 1)
merge_samples <- samples[dataset_prefixes == "merge"]

train_expr_all <- t(allTab[, merge_samples])
train_response_all <- ifelse(gsub(".*_", "", merge_samples) == "control", "control", "treat")

gmt_files <- dir(pattern = "\\.gmt$")
gmt_list <- list()

for (file in gmt_files) {
	gmt_content <- readLines(file, warn = FALSE)

	pathways <- list()
	for (line in gmt_content) {
		parts <- strsplit(line, "\t")[[1]]
		pathway_name <- parts[1]
		pathway_desc <- parts[2]
		genes <- parts[3:length(parts)]

		pathways <- list(
			genes = genes
		)
	}

	file_name <- sub("^GOMF_|^KEGG_|^GOBP_|^KEGG_MEDICUS_REFERENCE_", "", file)
	file_name <- sub("\\.v2024\\.1\\.Hs\\.gmt$", "", file_name)

	gmt_list[[file_name]] <- pathways
}

samples <- colnames(allTab)
dataset_prefixes <- sapply(strsplit(samples, "_"), `[`, 1)
unique_datasets <- unique(dataset_prefixes)
test_datasets <- setdiff(unique_datasets, "merge")

ctrl <- trainControl(
	method = "cv",
	number = 5,
	savePredictions = "final",
	classProbs = TRUE,
	summaryFunction = twoClassSummary,
	allowParallel = F
)
gene_sets = gmt_list
pathway_models <- list()
for (pathway_name in names(gene_sets)) {
	genes <- gene_sets[[pathway_name]]$genes
	pathway_genes <- intersect(genes, colnames(train_expr_all))
	cat(pathway_name, "genes count:", length(pathway_genes), "\n")
	if (length(pathway_genes) < 1) next

	df_train <- data.frame(train_expr_all[, pathway_genes, drop = FALSE])
	df_train$response <- as.factor(train_response_all)

	print(dim(df_train))
	print(table(df_train$response))
	set.seed(123)
	model <- tryCatch(
		train(
			response ~ .,
			data = df_train,
			method = "glm",
			family = binomial(),
			trControl = ctrl,
			metric = "ROC"
		),
		error = function(e) {
			cat("Error in training", pathway_name, ":\n", e$message, "\n")
			return(NULL)
		}
	)

	if (!is.null(model)) {
		pathway_models[[pathway_name]] <- model
	}
}


gene_lmbc <- read.table("./output_mRNA_selected.txt", header = TRUE, sep = "\t", check.names = FALSE)[,1]
target_genes <- intersect(gene_lmbc, colnames(train_expr_all))
df_train_lmbd <- data.frame(train_expr_all[, target_genes, drop = FALSE], response = train_response_all)
formula_lmbd <- as.formula(paste("response ~", paste(target_genes, collapse = "+")))
ctrl_lmbd <- trainControl(
	method = "cv",
	number = 5,
	savePredictions = "final",
	classProbs = TRUE,
	summaryFunction = twoClassSummary,
	allowParallel = F
)
set.seed(123)
model_lmbd <- train(
	formula_lmbd,
	data = df_train_lmbd,
	method = "glm",
	family = "binomial",
	trControl = ctrl_lmbd,
	metric = "ROC")

train_expr <- train_expr_all
train_response <- ifelse(train_response_all == "control", 0, 1)

pathway_auc_train <- data.frame()

for (pathway_name in names(pathway_models)) {
	model <- pathway_models[[pathway_name]]
	model_genes <- setdiff(names(coef(model$finalModel)), "(Intercept)")
	valid_genes <- intersect(model_genes, colnames(train_expr))

	if (length(valid_genes) == 0) next

	df_test <- data.frame(train_expr[, valid_genes, drop = FALSE], response = train_response)
	pred <- predict(model, newdata = df_test, type = "prob")
	auc_value <- auc(roc(train_response, pred[["treat"]]))

	pathway_auc_train <- rbind(pathway_auc_train, data.frame(
		Pathway = pathway_name,
		AUC = auc_value,
		Type = "Pathway"
	))
}

df_test_lmbd <- data.frame(train_expr[, target_genes, drop = FALSE], response = train_response)
pred_lmbd <- predict(model_lmbd, newdata = df_test_lmbd, type = "prob")
auc_lmbd <- auc(roc(train_response, pred_lmbd[["treat"]]))

pathway_auc_train <- rbind(pathway_auc_train, data.frame(
	Pathway = "LMBD",
	AUC = auc_lmbd,
	Type = "LMBD"
))

plot_data <- pathway_auc_train
plot_data$Group <- ifelse(plot_data$Pathway == "LMBD", "LMBD", "Pathway")

pdf("AUC_Comparison_TrainSet.pdf", width = 8, height = 12)
P1 = ggplot(plot_data, aes(x = reorder(Pathway, AUC), y = AUC, fill = Group)) +
	geom_bar(stat = "identity", width = 0.6, colour = "black", size = 0.5) +
	coord_flip(clip = "off") +
	scale_fill_manual(values = c("Pathway" = "#dcdcdc", "LMBD" = "#ff6400")) +
	scale_y_continuous(
		expand = expansion(mult = c(0, 0.05)),
		limits = c(0, max(plot_data$AUC) * 1.05)
	) +
	geom_segment(
		aes(x = reorder(Pathway, AUC), xend = reorder(Pathway, AUC), y = AUC - max(plot_data$AUC) * 0.03, yend = AUC + max(plot_data$AUC) * 0.03),
		size = 0.8,
		color = "black"
	) +
	geom_text(
		aes(label = sprintf("%.3f", AUC)),
		position = position_nudge(y = max(plot_data$AUC) * 0.04),
		hjust = 0,
		size = 3.5,
		color = "black"
	) +
	labs(
		title = "Train Set",
		x = NULL,
		y = "AUC Value"
	) +
	theme_minimal(base_size = 12) +
	theme(
		legend.position = "bottom",
		legend.title = element_blank(),
		axis.text.y = element_text(size = 12, color = "black"),
		axis.line.x = element_line(color = "black", size = 0.8),
		axis.line.y = element_line(color = "black", size = 0.8),
		panel.grid = element_blank(),
		plot.title = element_text(
			size = 16,
			face = "bold",
			hjust = 0.5,
			margin = margin(b = 15)
		)
	)
print(P1)
dev.off()

roc_obj_train <- roc(train_response, pred_lmbd[["treat"]])
auc_ci_train <- ci.auc(roc_obj_train, method = "bootstrap", boot.n = 2000)
sp_ci_train <- ci.sp(roc_obj_train, sensitivities = seq(0, 1, 0.01), boot.n = 2000, progress = "none")

pdf("ROC_TrainSet_LMBD.pdf", width = 7, height = 6)
par(mar = c(4.5, 4.5, 2, 2), mgp = c(2.5, 1, 0))
plot(roc_obj_train,
		 print.auc = FALSE,
		 grid = c(0.2, 0.2),
		 grid.col = "lightgray",
		 legacy.axes = TRUE,
		 xlab = "False Positive Rate (1 - Specificity)",
		 ylab = "True Positive Rate (Sensitivity)",
		 main = "ROC Curve For Train Set (LMBD)",
		 col = "#2c7fb8",
		 lwd = 2)
plot(sp_ci_train,
		 type = "shape",
		 col = alpha("#5bd1d7", 0.3),
		 border = NA)
auc_text <- sprintf("AUC = %.3f (95%% CI: %.3f-%.3f)",
										auc(roc_obj_train),
										auc_ci_train[1],
										auc_ci_train[3])
text(x = 0.6, y = 0.2,
		 labels = auc_text,
		 col = "#2c7fb8",
		 cex = 1.2)
legend("bottomright",
			 legend = c("ROC Curve ", "95% Confidence Band"),
			 col = c("#2c7fb8", alpha("#5bd1d7", 0.3)),
			 lwd = c(2, 10),
			 bty = "n",
			 cex = 0.9)
dev.off()

test_datasets <- c(setdiff(unique(dataset_prefixes), "merge"))

for (test_dataset in test_datasets) {
	if (test_dataset == "merge_test") {
		test_samples <- test_samples_merge
	} else {
		test_samples <- samples[dataset_prefixes == test_dataset]
	}

	test_expr <- t(allTab[, test_samples])
	test_response <- ifelse(gsub(".*_", "", test_samples) == "control", 0, 1)

	pathway_auc <- data.frame()
	for (pathway_name in names(pathway_models)) {
		model <- pathway_models[[pathway_name]]
		model_genes <- setdiff(names(coef(model$finalModel)), "(Intercept)")
		valid_genes <- intersect(model_genes, colnames(test_expr))

		if (length(valid_genes) == 0) next

		df_test <- data.frame(test_expr[, valid_genes, drop = FALSE], response = test_response)
		pred <- predict(model, newdata = df_test, type = "prob")
		auc_value <- auc(roc(test_response, pred[["treat"]]))

		pathway_auc <- rbind(pathway_auc, data.frame(
			Pathway = pathway_name,
			AUC = auc_value,
			Type = "Pathway"
		))
	}

	df_test_lmbd <- data.frame(test_expr[, target_genes, drop = FALSE], response = test_response)
	pred_lmbd <- predict(model_lmbd, newdata = df_test_lmbd, type = "prob")
	auc_lmbd <- auc(roc(test_response, pred_lmbd[["treat"]]))

	pathway_auc <- rbind(pathway_auc, data.frame(
		Pathway = "LMBD",
		AUC = auc_lmbd,
		Type = "LMBD"
	))

	roc_obj <- roc(
		response = df_test_lmbd$response,
		predictor = pred_lmbd[["treat"]]
	)

	set.seed(123)
	auc_ci <- ci.auc(roc_obj, method = "bootstrap", boot.n = 2000)

	sp_ci <- ci.sp(roc_obj, sensitivities = seq(0, 1, 0.01), boot.n = 2000, progress = "none")

	pdf(paste0("ROC_", test_dataset, ".pdf"), width = 7, height = 6)
	par(mar = c(4.5, 4.5, 2, 2), mgp = c(2.5, 1, 0))

	plot(roc_obj,
			 print.auc = FALSE,
			 grid = c(0.2, 0.2),
			 grid.col = "lightgray",
			 legacy.axes = TRUE,
			 xlab = "False Positive Rate (1 - Specificity)",
			 ylab = "True Positive Rate (Sensitivity)",
			 main = paste0("ROC Curve For ",test_dataset),
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
				 legend = c("ROC Curve ", "95% Confidence Band"),
				 col = c("#2c7fb8", alpha("#5bd1d7", 0.3)),
				 lwd = c(2, 10),
				 bty = "n",
				 cex = 0.9)

	dev.off()
	plot_data <- pathway_auc

	plot_data$Group <- ifelse(plot_data$Pathway == "LMBD", "LMBD", "Pathway")
	library(ggplot2)
	pdf(paste0("AUC_Comparison_", test_dataset, ".pdf"), width=8, height=12)
	P1=ggplot(plot_data, aes(x = reorder(Pathway, AUC), y = AUC, fill = Group)) +
		geom_bar(stat = "identity", width = 0.6, colour = "black", size = 0.5) +
		coord_flip(clip = "off") +
		scale_fill_manual(values = c("Pathway" = "#dcdcdc", "LMBD" = "#ff6400")) +
		scale_y_continuous(
			expand = expansion(mult = c(0, 0.05)),
			limits = c(0, max(plot_data$AUC) * 1.05)
		) +
		geom_segment(
			aes(x = reorder(Pathway, AUC), xend = reorder(Pathway, AUC), y = AUC - max(plot_data$AUC) * 0.03, yend = AUC + max(plot_data$AUC) * 0.03),
			size = 0.8,
			color = "black"
		) +
		geom_text(
			aes(label = sprintf("%.3f", AUC)),
			position = position_nudge(y = max(plot_data$AUC) * 0.04),
			hjust = 0,
			size = 3.5,
			color = "black"
		) +
		labs(
			title = test_dataset,
			x = NULL,
			y = "AUC Value"
		) +
		theme_minimal(base_size = 12) +
		theme(
			legend.position = "bottom",
			legend.title = element_blank(),
			axis.text.y = element_text(size = 12, color = "black"),
			axis.line.x = element_line(color = "black", size = 0.8),
			axis.line.y = element_line(color = "black", size = 0.8),
			panel.grid = element_blank(),
			plot.title = element_text(
				size = 16,
				face = "bold",
				hjust = 0.5,
				margin = margin(b = 15)
			)
		)
	print(P1)
	dev.off()
}
