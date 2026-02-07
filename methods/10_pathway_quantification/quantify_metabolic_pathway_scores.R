library(sva)
library(IOBR)
library(GSVA)
library(ConsensusClusterPlus)
library(Seurat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(GSEABase)
library(AUCell)
library(GSVA)
library(UCell)
library(singscore)
library(tidyverse)
setwd("D:\\哮喘乳酸化单细胞开篇\\13.114种代谢通路的量化")
rm(list = ls())
data(signature_metabolism)
table(oc$celltype)
metabolic_genesets <- lapply(signature_metabolism, function(x) unlist(strsplit(x, split = ",")))
Addsingscore<-function(oc,signatures,seed=123,assay=null,slot="data"){
	set.seed(seed)
	if(is.null(assay)){
		assay=Seurat::DefaultAssay(oc)
	}
	matrix=Seurat::GetAssayData(object=oc,slot=slot,assay=assay)

	markers=signatures
	h.gsets.list=markers %>% purrr::compact()

	singscore.rank=singscore::rankGenes(as.data.frame(matrix))
	#calculate separately
	singscore.scores=list()
	for (i in seq_along(h.gsets.list)){
		if (any(stringr::str_detect(h.gsets.list[[i]],pattern = "\\+$|-$"))){
			h.gsets.list.positive=stringr::str_match(h.gsets.list[[i]],pattern="(.+)\\+")[,2] %>% purrr::discard(is.na)
			h.gsets.list.negative=stringr::str_match(h.gsets.list[[i]],pattern="(.+)-")[,2] %>% purrr::discard(is.na)
			if(length(h.gsets.list.positive)==0){
				singscore.scores[[i]]=singscore::simpleScore(singscore.rank,upSet=h.gsets.list.negative,centerScore=F)
			}
			if(length(h.gsets.list.negative)==0){
				singscore.scores[[i]]=singscore::simpleScore(singscore.rank,upSet=h.gsets.list.positive,centerScore=F)
			}
			if((length(h.gsets.list.positive)!=0)&(length(h.gsets.list.negative)!=0)){
				singscore.scores[[i]]=singscore::simpleScore(singscore.rank,upSet=h.gsets.list.positive,downSet=h.gsets.list.negative,centerScore=F)
			}
		}else{
			singscore.scores[[i]]=singscore::simpleScore(singscore.rank,upSet=h.gsets.list[[i]],centerScore=F)
		}
		TotalScore=NULL
		singscore.scores[[i]]=singscore.scores[[i]] %>%
			dplyr::select(TotalScore) %>%
			magrittr::set_colnames(names(h.gsets.list)[i])}
	names(singscore.scores)=names(h.gsets.list)
	singscore.scores=do.call(cbind,singscore.scores)
	oc=Seurat::AddMetaData(oc,as.data.frame(singscore.scores))
	return(oc)
}
DefaultAssay(oc) <- "RNA"


for (pathway_name in names(metabolic_genesets)) {
	genes <- metabolic_genesets[[pathway_name]]
	genes_vector <- as.character(genes)

	overlap_genes <- intersect(genes_vector, rownames(oc@assays$RNA$data))
	if (length(overlap_genes) == 0) {
		cat("Skipping pathway (no overlap):", pathway_name, "\n")
		next
	}

	cat("Processing pathway:", pathway_name, "\n")
	genes_list <- list(overlap_genes)

	##### AddModuleScore #####
	oc <- AddModuleScore(oc,
											 features = genes_list,
											 name = paste0(pathway_name, "_Add_Score"),
											 slot = "data")
	new_col <- paste0(pathway_name, "_AddModule")
	colnames(oc@meta.data)[ncol(oc@meta.data)] <- new_col

	##### AUCell #####
	cells_rankings <- AUCell_buildRankings(oc@assays$RNA$data, nCores=1)
	cells_AUC <- AUCell_calcAUC(list(genes = overlap_genes), cells_rankings,
															aucMaxRank = nrow(cells_rankings)*0.05)
	oc@meta.data[[paste0(pathway_name, "_AUCell")]] <- getAUC(cells_AUC)["genes", ]

	##### ssGSEA #####
	genes_ssGSEA <- as.data.frame(overlap_genes)
	gene.expr <- as.matrix(oc@assays$RNA$data)
	ssGSEA.result <- ssgseaParam(gene.expr, genes_ssGSEA)
	ssGSEA.result <- gsva(ssGSEA.result)
	oc@meta.data[[paste0(pathway_name, "_ssGSEA")]] <- ssGSEA.result["overlap_genes", ]

	##### Ucell #####
	oc <- AddModuleScore_UCell(oc, features = list(temp = overlap_genes),
														 name = paste0(pathway_name, "_Ucell_"),
														 assay = "RNA", slot = "data")
	colnames(oc@meta.data)[ncol(oc@meta.data)] <- paste0(pathway_name, "_Ucell")

	##### singscore #####
	oc <- Addsingscore(oc = oc, signatures = list(temp = overlap_genes),
										 assay = "RNA", slot = "data")
	colnames(oc@meta.data)[ncol(oc@meta.data)] <- paste0(pathway_name, "_singscore")

	score_columns <- c(
		paste0(pathway_name, "_AddModule"),
		paste0(pathway_name, "_AUCell"),
		paste0(pathway_name, "_ssGSEA"),
		paste0(pathway_name, "_Ucell"),
		paste0(pathway_name, "_singscore"))

	oc@meta.data[[pathway_name]] <- rowSums(
		oc@meta.data[, score_columns], na.rm = TRUE
	)
}


oc@meta.data <- oc@meta.data[, !grepl("(_singscore|_Ucell|_ssGSEA|_AUCell|_AddModule)$", colnames(oc@meta.data))]

save(oc,file = "MBDs_AA.Rdata")

library(data.table)
library(Seurat)
library(AUCell)
library(GSVA)
gene_df <- data.table::fread("./output_mRNA_selected.txt", data.table = FALSE)
genes_raw <- as.character(gene_df[[1]])
genes <- unique(na.omit(genes_raw))

expr <- oc@assays$RNA$data
overlap_genes <- intersect(genes, rownames(expr))
if (length(overlap_genes) == 0) {
	stop("LMBDs: 与表达矩阵无基因交集，无法打分。")
} else {
	cat("LMBDs 基因数（交集）:", length(overlap_genes), "\n")
}

oc <- AddModuleScore(
	object   = oc,
	features = list(overlap_genes),
	name     = "LMBDs_Add_Score",
	slot     = "data"
)
colnames(oc@meta.data)[ncol(oc@meta.data)] <- "LMBDs_AddModule"

## 4) AUCell
cells_rankings <- AUCell_buildRankings(expr, nCores = 1)
cells_AUC <- AUCell_calcAUC(
	geneSets   = list(LMBDs = overlap_genes),
	rankings   = cells_rankings,
	aucMaxRank = nrow(cells_rankings) * 0.05
)
auc_mat <- getAUC(cells_AUC)
oc@meta.data[["LMBDs_AUCell"]] <- as.numeric(auc_mat["LMBDs", rownames(oc@meta.data)])

genes_ssGSEA = as.data.frame(genes)
gene.expr <- as.matrix(oc@assays$RNA$data)
ssGSEA.result <- ssgseaParam(gene.expr, genes_ssGSEA)
ssGSEA.result <- gsva(ssGSEA.result)
oc@meta.data[["LMBDs_ssGSEA"]] <- ssGSEA.result["genes",]

oc <- AddModuleScore_UCell(
	oc,
	features = list(LMBDs = overlap_genes),
	name     = "LMBDs_Ucell_",
	assay    = "RNA",
	slot     = "data"
)
colnames(oc@meta.data)[ncol(oc@meta.data)] <- "LMBDs_Ucell"

oc <- Addsingscore(
	oc = oc,
	signatures = list(LMBDs = overlap_genes),
	assay = "RNA",
	slot  = "data"
)
colnames(oc@meta.data)[ncol(oc@meta.data)] <- "LMBDs_singscore"

score_cols <- c("LMBDs_AddModule", "LMBDs_AUCell", "LMBDs_ssGSEA", "LMBDs_Ucell", "LMBDs_singscore")
oc@meta.data[["LMBDs"]] <- rowSums(oc@meta.data[, score_cols, drop = FALSE], na.rm = TRUE)
oc@meta.data <- oc@meta.data[, !grepl("(_singscore|_Ucell|_ssGSEA|_AUCell|_AddModule)$", colnames(oc@meta.data))]

cat("打分完成：已在 meta.data 中新增列：\n",
		paste(c(score_cols, "LMBDs"), collapse = ", "), "\n")

library(randomForest)
library(tidyverse)
library(dplyr)
library(stringr)
df <- data.frame(oc@meta.data)
selected_features <- gsub("-", ".", names(signature_metabolism))
df <- df %>% dplyr::select(all_of(c("LMBDs", selected_features)))

rf_model <- randomForest(
	LMBDs ~ .,
	data = df,
	ntree = 1000,
	importance = TRUE,
	proximity = TRUE
)

print(rf_model)

importance <- importance(rf_model)
colnames(importance)
library(ggplot2)
library(dplyr)
library(scales)

df_oob <- tibble(
	trees  = seq_along(rf_model$mse),
	OOB_MSE = rf_model$mse,
	OOB_R2  = rf_model$rsq
)

win <- 50
rel_improve <- c(
	rep(NA, win),
	(df_oob$OOB_MSE[(win+1):nrow(df_oob)] - zoo::rollmean(df_oob$OOB_MSE, k = win, align = "right")[1:(nrow(df_oob)-win)]) /
		zoo::rollmean(df_oob$OOB_MSE, k = win, align = "right")[1:(nrow(df_oob)-win)]
)
df_oob$rel_improve <- rel_improve
best_ntree <- df_oob %>%
	filter(!is.na(rel_improve)) %>%
	mutate(flag = abs(rel_improve) < 0.005) %>%
	summarize(best = if (any(flag)) trees[which(flag)[1]] else which.min(OOB_MSE)) %>%
	pull(best)

r2_min <- min(df_oob$OOB_R2, na.rm = TRUE)
r2_max <- max(df_oob$OOB_R2, na.rm = TRUE)
mse_min <- min(df_oob$OOB_MSE, na.rm = TRUE)
mse_max <- max(df_oob$OOB_MSE, na.rm = TRUE)

scale_r2_to_mse <- function(x) {
	(x - r2_min) / (r2_max - r2_min) * (mse_max - mse_min) + mse_min
}
scale_mse_to_r2 <- function(x) {
	(x - mse_min) / (mse_max - mse_min) * (r2_max - r2_min) + r2_min
}

df_oob <- df_oob %>%
	mutate(R2_on_MSE_axis = scale_r2_to_mse(OOB_R2))

p <- ggplot(df_oob, aes(trees)) +
	geom_line(aes(y = OOB_MSE), linewidth = 1) +
	geom_smooth(aes(y = OOB_MSE), method = "loess", se = FALSE, span = 0.15, linewidth = 0.6, linetype = 2) +
	geom_line(aes(y = R2_on_MSE_axis), linewidth = 0.9, alpha = 0.7) +
	geom_vline(xintercept = best_ntree, linewidth = 0.6, linetype = 3) +
	annotate("label",
					 x = best_ntree, y = df_oob$OOB_MSE[best_ntree],
					 label = paste0("拐点 ≈ ", best_ntree, " trees"),
					 hjust = -0.05, vjust = -0.8, size = 3) +
	scale_y_continuous(
		name = "OOB MSE",
		labels = label_number(accuracy = 0.001),
		sec.axis = sec_axis(~ scale_mse_to_r2(.), name = "OOB R²")
	) +
	scale_x_continuous(name = "Number of trees") +
	ggtitle("Random Forest (Regression): OOB Error vs Trees") +
	theme_minimal(base_size = 12) +
	theme(
		plot.title = element_text(face = "bold", hjust = 0.02)
	)

ggsave("OOB_curves_beautified.pdf", p, width = 7, height = 5)
ggsave("OOB_curves_beautified.png", p, width = 7, height = 5, dpi = 300)

p

pdf(file="geneImportance.pdf", width=12, height=8)
par(mar = c(6, 6, 4, 4) + 0.1, mgp = c(4, 1, 0),
		bg = "white", fg = "black", col.axis = "black", col.lab = "black")
imp_scores <- importance[, "%IncMSE"]
imp_sorted <- sort(imp_scores, decreasing = FALSE)

top_index <- length(imp_sorted)
bar_colors <- rep("#1F78B4", length(imp_sorted))  # 全部先设为蓝色
bar_colors[top_index] <- "#FF7F00"  # 最右侧的条形设为橙色

bp <- barplot(imp_sorted,
							horiz = FALSE,
							las = 1,
							col = bar_colors,
							border = "black",
							space = 0.3,
							ylab = "IncNodePurity",
							xlab = "Features (Ordered by Increasing Importance)",
							main = "Feature Importance Ranking\n(Random Forest)",
							ylim = c(0, max(imp_sorted) * 1.4),
							names.arg = rep("", length(imp_sorted)),
							cex.axis = 1.0,
							cex.lab = 1.2,
							cex.main = 1.3)

box(lwd = 2, col = "black")

grid(nx = NA, ny = NULL, col = "gray90", lty = 2)

legend("topleft",
			 inset = c(0.02, 0.02),
			 legend = c("Most Important Feature", "Other Features"),
			 fill = c("#FF7F00", "#1F78B4"),
			 border = "black",
			 box.lwd = 1,
			 bg = "white",
			 cex = 0.9)

top_feature <- names(imp_sorted)[top_index]
top_value <- imp_sorted[top_index]

bar_x <- bp[top_index]
bar_y <- top_value * 1.05

label_x <- mean(par("usr")[1:2])
label_y <- max(imp_sorted) * 1.25

segments(
	x0 = label_x, y0 = label_y,
	x1 = bar_x, y1 = bar_y,
	col = "gray30", lwd = 1.5, lty = 2
)
arrows(
	x0 = bar_x, y0 = bar_y,
	x1 = bar_x, y1 = top_value * 0.98,
	length = 0.1, col = "gray30", lwd = 1.5
)

text(
	x = label_x, y = label_y,
	labels = paste0("Most Important Feature:\n", top_feature),
	pos = 3, cex = 1.0, col = "black", font = 2
)

dev.off()

library(ggplot2)
library(ggpubr)
library(gghalves)
library(dplyr)
library(ggsci)
output_dir <- "MetaAnalysis_Results"
if (!dir.exists(output_dir)) dir.create(output_dir)

plot_data <- oc@meta.data %>%
	rownames_to_column("CellID") %>%
	filter(phenotype %in% c("AA", "ANA")) %>%
	dplyr::select(CellID, celltype, Retinoic_Acid_Metabolism) %>%
	mutate(celltype = factor(celltype))

cell_types <- levels(plot_data$celltype)
color_scheme <- scales::hue_pal()(length(cell_types))

color_scheme <- pal_npg("nrc")(length(cell_types))

p <- ggplot(plot_data, aes(x = celltype, y = Retinoic_Acid_Metabolism, fill = celltype)) +

	gghalves::geom_half_point(
		side = "l",
		size = 1.8,
		alpha = 0.6,
		color = "black",
		shape = 21,
		transformation = position_jitter(width = 0.05, height = 0),
		position = position_nudge(x = -0.05)
	) +

	gghalves::geom_half_boxplot(
		side = "l",
		notch = TRUE,
		outlier.shape = NA,
		width = 0.25,
		alpha = 0.8,
		color = "black",
		lwd = 0.6,
		position = position_nudge(x = 0)
	) +

	gghalves::geom_half_violin(
		side = "r",
		color = "black",
		lwd = 0.6,
		alpha = 0.7,
		trim = TRUE,
		scale = "width",
		width = 0.25,
		position = position_nudge(x = 0.05)
	) +

	scale_fill_manual(values = color_scheme) +

	labs(
		x = "Cell Type",
		y = "Retinoic_Acid_Metabolism Pathway Score",
		title = "Retinoic_Acid_Metabolism Pathway Activity Across Cell Types"
	) +

	theme_minimal(base_size = 14) +
	theme(
		plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
		axis.title.x = element_text(size = 14, face = "bold"),
		axis.title.y = element_text(size = 14, face = "bold"),
		axis.text.x = element_text(size = 12, face = "bold", color = "black", angle = 30, hjust = 1),
		axis.text.y = element_text(size = 12, face = "bold", color = "black"),
		axis.line = element_line(color = "black", size = 0.8),
		axis.ticks = element_line(color = "black", size = 0.8),
		panel.grid.major = element_line(color = "grey90", size = 0.3),
		panel.grid.minor = element_blank(),
		legend.position = "none",
		plot.background = element_rect(fill = "white", color = NA),
		panel.background = element_rect(fill = "white", color = NA)
	) +

	coord_cartesian(ylim = c(min(plot_data$Retinoic_Acid_Metabolism) * 0.9, max(plot_data$Retinoic_Acid_Metabolism) * 1.1))

ggsave(
	file.path(output_dir, "Retinoic_Acid_Metabolism_CellType_Comparison_SingleGroup.pdf"),
	plot = p,
	width = 14,
	height = 8,
	dpi = 300
)

print(p)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(viridis)
df <- data.frame(oc@meta.data, oc@reductions[["umap"]]@cell.embeddings)

class_avg <- df %>%
	group_by(celltype) %>%
	summarise(umap_1 = median(umap_1),
						umap_2 = median(umap_2))
colnames(df)
score_methods <- c("Retinoic_Acid_Metabolism")

for (method in score_methods) {
	pdf(file = paste0(method, "_umap.pdf"), width = 7, height = 6)

	p <- ggplot(df, aes(x = umap_1, y = umap_2)) +
		geom_point(aes_string(color = method), size = 0.5, alpha = 0.8) +
		ggrepel::geom_label_repel(
			aes(label = celltype),
			data = class_avg,
			size = 3,
			box.padding = 0.35,
			segment.color = NA,
			label.size = 0.25,
			color = "black",
			fill = "white"
		) +
		scale_color_viridis(
			option = "H",
			name = "Score",
			guide = guide_colorbar(
				barwidth = 0.8,
				barheight = 5,
				title.position = "top",
				title.hjust = 0.5
			)
		) +
		theme_classic() +
		labs(x = "UMAP1", y = "UMAP2", title = paste(method, "Score")) +
		theme(
			axis.title = element_text(size = 12, face = "bold"),
			axis.text = element_text(color = "black"),
			plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
			panel.border = element_rect(color = "#256a53", size = 1.5, fill = NA),
			panel.background = element_rect(fill = "#f5f5f5"),
			legend.position = "right"
		)

	print(p)
	dev.off()
}
