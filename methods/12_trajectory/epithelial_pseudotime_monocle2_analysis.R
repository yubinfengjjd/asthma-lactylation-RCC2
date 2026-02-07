library(cowplot)
library(dplyr)
library(data.table)
library(tidyverse)
library(paletteer)
library(vcd)
library(limma)
library(gplots)
library(Seurat)
library(clustree)
library(ggplot2)
library(reshape2)
library(dplyr)
library(monocle)
library(RColorBrewer)
library(ggpubr)
library(ggsignif)
library(patchwork)
library(tidydr)
library(ggforce)
library(ggrastr)
library(viridis)
library(gridExtra)
library(ggnewscale)
library(org.Mm.eg.db)
library(ComplexHeatmap)
library(ClusterGVis)
library(monocle)
library(Seurat)
rm(list = ls())
setwd("D:\\哮喘乳酸化单细胞开篇\\上皮细胞拟时序")
load(file.path("D:\\哮喘乳酸化单细胞开篇\\4.乳酸化评分","Lactylation.Rdata"))
oc_subset = subset(oc,celltype == "AEC") %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 1500)
seurat_to_monocle <- function(otherCDS, assay, slot, lowerDetectionLimit = 0, import_all = FALSE) {
	if(class(otherCDS)[1] == 'Seurat') {
		requireNamespace("Seurat")
		data <- GetAssayData(otherCDS, assay = assay, slot = slot)
		data <- data[rowSums(as.matrix(data)) != 0,]
		pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
		fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
		fd <- new("AnnotatedDataFrame", data = fData)
		if(all(data == floor(data))) {
			expressionFamily <- negbinomial.size()
		} else if(any(data < 0)){
			expressionFamily <- uninormal()
		} else {
			expressionFamily <- tobit()
		}
		valid_data <- data[, row.names(pd)]
		monocle_cds <- newCellDataSet(data,
																	phenoData = pd,
																	featureData = fd,
																	lowerDetectionLimit=lowerDetectionLimit,
																	expressionFamily=expressionFamily)
		if(import_all) {
			if("Monocle" %in% names(otherCDS@misc)) {
				otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
				otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
				monocle_cds <- otherCDS@misc$Monocle
				mist_list <- otherCDS
			} else {
				mist_list <- otherCDS
			}
		} else {
			mist_list <- list()
		}
	}
	return(monocle_cds)
}
cds_seurat <- seurat_to_monocle(oc_subset, assay = "RNA", slot = "counts")
cds_seurat <- estimateSizeFactors(cds_seurat)
cds_seurat <- estimateDispersions(cds_seurat)

cds_seurat <- detectGenes(cds_seurat, min_expr = 0.1)
print(head(fData(cds_seurat)))
expressed_genes <- row.names(subset(fData(cds_seurat), num_cells_expressed >= 10))
pData(cds_seurat)$Total_mRNAs <- Matrix::colSums(exprs(cds_seurat))
cds_seurat <- cds_seurat[,pData(cds_seurat)$Total_mRNAs < 1e6]

cds_seurat <- cds_seurat
#seurat_variable_genes <- gene[gene%in%rownames(oc_subset@assays$RNA$data)]
seurat_variable_genes <- VariableFeatures(oc_subset)
cds_seurat <- setOrderingFilter(cds_seurat, seurat_variable_genes)
plot_ordering_genes(cds_seurat)

cds_seurat <- reduceDimension(cds_seurat, max_components = 2,reduction_method = 'DDRTree')
cds_seurat <- orderCells(cds_seurat)

save(cds_seurat,file="cds_seurat.Rdata")

pdf(file = "plot_cell_trajectory.PDF",width = 6,height = 6)
plot_cell_trajectory(cds_seurat, markers = "MNDA",use_color_gradient=T,cell_size = 1,cell_link_size = 1.5)
dev.off()

pdf(file = "plot_genes_in_pseudotime.PDF",width = 6,height = 6)
cds_subset=cds_seurat[c("RCC2"),]
plot_genes_in_pseudotime(cds_subset,color_by = "State")
dev.off()
pdf(file = "celltype.PDF",width = 6,height = 6)
plot_cell_trajectory(cds_seurat,color_by="State")
dev.off()

#####SCP####
library(SCP)
library(DDRTree)
library(monocle)
oc_subset[["RNA"]] <- as(object = oc_subset[["RNA"]], Class = "Assay")
DefaultAssay(oc_subset) <- "RNA"
oc_subset=RunMonocle2(oc_subset,features = VariableFeatures(oc_subset))
oc_subset@tools$Monocle2$cds=cds_seurat
save(oc_subset,file = "oc_subset_monocle.Rdata")
trajectory <- oc_subset@tools$Monocle2$trajectory
pdf("State.PDF", width = 6.5, height = 5.5)
CellDimPlot(oc_subset, group.by = "Monocle2_State", reduction = "DDRTree", label = TRUE, theme_use = "theme_blank") + trajectory
dev.off()
pdf("celltype.PDF", width = 6.5, height = 5.5)
CellDimPlot(oc_subset, group.by = "celltype", reduction = "DDRTree", label = TRUE, theme_use = "theme_blank")
dev.off()
pdf("Pseudotime.PDF", width = 6.5, height = 5.5)
FeatureDimPlot(oc_subset, features = "Monocle2_Pseudotime", reduction = "DDRTree", theme_use = "theme_blank")
dev.off()
pdf("RCC2.PDF", width = 6.5, height = 5.5)
FeatureDimPlot(oc_subset, features = "RCC2", reduction = "DDRTree", theme_use = "theme_blank")
dev.off()
pdf("Scoreing.PDF", width = 6.5, height = 5.5)
FeatureDimPlot(oc_subset, features = "Scoreing", reduction = "DDRTree", theme_use = "theme_blank")
dev.off()

pdf("plot_genes_in_pseudotime.PDF", width = 6.5, height = 5.5)

cds_subset <- cds_seurat[c("RCC2"),]

p <- plot_genes_in_pseudotime(
	cds_subset,
	color_by = "State",
	cell_size = 1.8
)

p +
	scale_color_viridis_d(
		name = "State",
		option = "D",
		end = 0.9,
		guide = guide_legend(
			override.aes = list(size = 3),
			nrow = 2
		)
	) +
	theme_minimal(base_size = 12) +
	theme(
		legend.position = "bottom",
		legend.box = "horizontal",
		legend.spacing.x = unit(0.2, "cm"),
		panel.grid.major = element_line(linewidth = 0.2),
		axis.line = element_line(linewidth = 0.3),
		plot.title = element_text(hjust = 0.5, face = "bold")
	) +
	labs(
		x = "Pseudotime",
		y = "Expression Level (log)",
		title = "RCC2 Expression Dynamics"
	)
dev.off()

library(ClusterGVis)
library(dplyr)
diff_test_res <- differentialGeneTest(cds_seurat,
																			fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res1 <- diff_test_res %>%
	filter(num_cells_expressed > 200)
df=plot_pseudotime_heatmap2(cds_seurat[row.names(subset(diff_test_res1,qval<1e-4)),],
														num_clusters=3)
gene=sample(df$wide.res$gene,20,replace=F)

pdf("plot_genes_heatmap.PDF", width=6, height=8, onefile=FALSE)
visCluster(object=df, plot.type="heatmap", markGenes=gene)
dev.off()

library(RcppML)
library(irGSEA)
oc_subset$State <- pData(cds_seurat)$State

DefaultAssay(oc_subset) <- "RNA"
oc_subset <- SeuratObject::UpdateSeuratObject(object = oc_subset)

cds_scored <- irGSEA.score(
	object = oc_subset,
	assay = "RNA",
	slot = "data",
	seeds = 123,
	ncores = 1,
	min.cells = 3,
	min.feature = 0,
	custom = F,
	geneset = NULL,
	msigdb = T,
	species = "Homo sapiens",
	category = "H",
	subcategory = NULL,
	geneid = "symbol",
	method = c("AUCell", "UCell", "singscore", "ssgsea", "viper"),
	aucell.MaxRank = NULL,
	ucell.MaxRank = NULL,
	kcdf = 'Gaussian'
)

result.dge <- irGSEA.integrate(
	object = cds_scored,
	group.by = "State",
	method = c("AUCell", "UCell", "singscore", "ssgsea", "viper")
)

pdf(file = "heatmap_state.pdf", width = 12, height = 6)
irGSEA.heatmap(
	object = result.dge,
	method = "RRA",
	top = 50,
	show.geneset = NULL
)
dev.off()

plot_data <- data.frame(oc_subset$Scoreing, oc_subset$Monocle2_State)
colnames(plot_data) = c("S100A4_expression", "Subtype")
plot_data$Subtype <- as.factor(plot_data$Subtype)

kruskal_p_value <- kruskal.test(S100A4_expression ~ Subtype, data = plot_data)$p.value

sig_label <- ifelse(kruskal_p_value < 0.001, "***",
										ifelse(kruskal_p_value < 0.01, "**",
													 ifelse(kruskal_p_value < 0.05, "*", "")))

library(ggplot2)

ggplot(plot_data, aes(x = Subtype, y = S100A4_expression)) +
	geom_violin(aes(fill = Subtype), alpha = 0.7) +
	geom_boxplot(aes(fill = Subtype), width = 0.15, outlier.shape = NA, color = "black") +
	geom_text(aes(x = 3, y = max(S100A4_expression)*1.1, label = sig_label),
						size = 5, vjust = -0.5) +
	scale_fill_manual(values = c("1"="#acd5ab", "2"="#feadac", "3"="#adeada","4"="#adadad","5"="#adaeef")) +
	labs(x = "State", y = "Lactylation") +
	theme_classic() +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
		legend.position = "top",
		panel.border = element_rect(color = "black", fill = NA),
		panel.grid.major.y = element_line(color = "grey90")
	) +
	coord_cartesian(ylim = c(min(plot_data$S100A4_expression), max(plot_data$S100A4_expression)*1.15))

ggsave("score_state.pdf", width = 6, height = 6)

#####S100A4_State####
plot_data <- data.frame(oc_subset@assays$RNA$data["RCC2",], oc_subset$Monocle2_State)
colnames(plot_data) = c("S100A4_expression", "Subtype")
plot_data$Subtype <- as.factor(plot_data$Subtype)

kruskal_p_value <- kruskal.test(S100A4_expression ~ Subtype, data = plot_data)$p.value

sig_label <- ifelse(kruskal_p_value < 0.001, "***",
										ifelse(kruskal_p_value < 0.01, "**",
													 ifelse(kruskal_p_value < 0.05, "*", "")))

library(ggplot2)

ggplot(plot_data, aes(x = Subtype, y = S100A4_expression)) +
	geom_violin(aes(fill = Subtype), alpha = 0.7) +
	geom_boxplot(aes(fill = Subtype), width = 0.15, outlier.shape = NA, color = "black") +
	geom_text(aes(x = 3, y = max(S100A4_expression)*1.1, label = sig_label),
						size = 5, vjust = -0.5) +
	scale_fill_manual(values = c("1"="#acd5ab", "2"="#feadac", "3"="#adeada","4"="#adadad","5"="#adaeef")) +
	labs(x = "State", y = "RCC2_expression") +
	theme_classic() +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
		legend.position = "top",
		panel.border = element_rect(color = "black", fill = NA),
		panel.grid.major.y = element_line(color = "grey90")
	) +
	coord_cartesian(ylim = c(min(plot_data$S100A4_expression), max(plot_data$S100A4_expression)*1.15))

ggsave("RCC2_expression_state.pdf", width = 6, height = 6)
