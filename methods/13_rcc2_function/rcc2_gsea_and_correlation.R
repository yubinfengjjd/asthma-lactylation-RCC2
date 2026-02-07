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
setwd("D:\\哮喘乳酸化单细胞开篇\\17.RCC2基因功能探究")
load(file.path("D:\\哮喘乳酸化单细胞开篇\\6.一致性聚类","Cluster.Rdata"))
table(oc_subset$celltype)
oc_subset <- subset(oc_subset, subset = phenotype == "AA") %>% NormalizeData() %>%
	FindVariableFeatures(selection.method = "vst" ,nfeatures = 1500 ,verbose = T)
oc_subset$group <- ifelse(
	oc_subset@assays$RNA$counts["RCC2", ] > 0,
	"RCC2+ Mono/Mac",
	"RCC2- Mono/Mac")
table(oc_subset$group)

gmtFile <- "c2.cp.kegg_legacy.v2024.1.Hs.entrez.gmt"
dataL=as.data.frame(oc_subset@assays$RNA$data[,oc_subset$group == "RCC2- Mono/Mac"])
dataH=as.data.frame(oc_subset@assays$RNA$data[,oc_subset$group == "RCC2+ Mono/Mac"])
data1=cbind(dataH,dataL)
Idents(oc_subset)=oc_subset$group
diffSig=FindMarkers(object = oc_subset,ident.1 = "RCC2+ Mono/Mac",ident.2 = "RCC2- Mono/Mac")
diffSig <- diffSig[with(diffSig, (p_val_adj < 0.05)), ]
write.table(cbind(Gene=rownames(diffSig),diffSig),file = "diffSig.txt",sep="\t",quote=F,row.names = F)
library(clusterProfiler)
logFC = diffSig$avg_log2FC
names(logFC)=rownames(diffSig)
logFC=sort(logFC,decreasing = T)
genes_sorted = names(logFC)
gene.id <- bitr(genes_sorted, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
names(logFC)=gene.id$ENTREZID
range(logFC)
KEGG_gseresult <- gseKEGG(logFC,
													organism = "hsa",
													pAdjustMethod = "BH",
													pvalueCutoff=0.05)
KEGG_gseresult=as.data.frame(KEGG_gseresult)
write.table(cbind(Gene=rownames(KEGG_gseresult),KEGG_gseresult),file = "KEGG_gseresult.txt",sep="\t",quote=F,row.names = F)
pdf(file="GSEA_boxplot.PDF", width=6, height=4)
library(ggplot2)
id=KEGG_gseresult$ID
p1 <- ggplot(data = KEGG_gseresult[id,], aes(x = setSize, y = Description, fill = NES)) +
	scale_fill_distiller(palette = "YlOrRd", direction = 1) +
	geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
	labs(x = "Number of Gene",
			 y = "pathway",
			 title = "Enrichment barplot") +
	geom_text(aes(x = 0.03, label = Description), hjust = 0) +
	theme_classic() +
	theme(axis.title = element_text(size = 13),
				axis.text = element_text(size = 8),
				plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
				legend.title = element_text(size = 13),
				legend.text = element_text(size = 11),
				axis.text.y = element_blank())
print(p1)
dev.off()

library(SCP)
table(oc$celltype)
load(file.path("D:\\哮喘乳酸化单细胞开篇\\6.一致性聚类","Cluster.Rdata"))
oc[["RNA"]] <- as(object = oc[["RNA"]], Class = "Assay")
pdf(file="MNDA_DIM.pdf", width=6, height=6)
FeatureDimPlot(oc,
							 features = "RCC2", reduction = "umap", label = TRUE,
							 cells.highlight = colnames(oc)[oc$celltype == "Mono/Mac"]
)
dev.off()

test <- as.data.frame(oc_subset@meta.data)
cell <- row.names(test)[test$celltype == "Mono/Mac"]
test_sub <- test[cell,]
tar <- oc_subset@assays$RNA$counts["RCC2",test$celltype == "Mono/Mac"]
result_sub <- test_sub[,"Scoreing"]
cor_res <- cor(x = tar, y = result_sub, method = 'pearson')

cor_test_result <- cor.test(tar, result_sub, method = 'pearson')
cor_rho <- round(cor_test_result$estimate, 2)
cor_p <- format.pval(cor_test_result$p.value, digits = 2)
n <- length(tar)
pdf(file="cor.plot.pdf", width=6, height=6)

par(bty = "o",
		mgp = c(2, 0.5, 0),
		mar = c(4.1, 4.1, 2.1, 4.1),
		tcl = -.25,
		font.main = 3)

plot(NULL, NULL,
		 ylim = range(tar),
		 xlim = range(result_sub),
		 xlab = "Lactylation",
		 ylab = "RCC2 Expression",
		 col = "white",
		 main = "")

rect(par("usr")[1], par("usr")[3],
		 par("usr")[2], par("usr")[4],
		 col = "#EAE9E9", border = FALSE)

grid(col = "white", lty = 1, lwd = 1.5)

points(x = result_sub, y = tar,
			 pch = 19,
			 col = scales::alpha("#E51718", 0.8),  # 使用单一红色系
			 cex = 1.5)

abline(lm(tar ~ result_sub),
			 lwd = 2,
			 col = "black")

rug(result_sub, side = 3, col = "black", lwd = 1)
rug(tar, side = 4, col = "black", lwd = 1)

text(x = min(result_sub),
		 y = max(tar) * 0.95,
		 adj = 0,
		 labels = bquote("Mon/Mac: N = " ~ .(n) ~
		 									"; " ~ rho ~ " = " ~ .(cor_rho) ~
		 									"; " ~ italic(P) ~ " = " ~ .(cor_p)),
		 col = "black",
		 cex = 0.8)
dev.off()

plot_data <- data.frame(oc_subset@assays$RNA$data["RCC2",],oc_subset$subtypes)
colnames(plot_data)=c("RCC2_expression","Subtype")
plot_data$Subtype <- as.factor(plot_data$Subtype)

p_value <- wilcox.test(RCC2_expression ~ Subtype, data = plot_data)$p.value

sig_label <- ifelse(p_value < 0.001, "***",
										ifelse(p_value < 0.01, "**",
													 ifelse(p_value < 0.05, "*", "")))


ggplot(plot_data, aes(x = Subtype, y = RCC2_expression)) +
	geom_violin(aes(fill = Subtype), position = position_dodge(0.9), alpha = 0.7) +
	geom_boxplot(aes(fill = Subtype),
							 width = 0.15,
							 position = position_dodge(0.9),
							 outlier.shape = NA,
							 color = "black") +
	geom_text(
		aes(x = 1.5, y = max(plot_data$RCC2_expression)*1.1, label = sig_label),
		size = 5, vjust = -0.5
	) +
	scale_fill_manual(values = c("MBC1"="#acd5ab", "MBC2"="#feadac")) +
	labs(x = "Subtype", y = "RCC2_Expression") +
	theme_classic() +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
		plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
		legend.position = "top",
		panel.border = element_rect(color = "black", fill = NA),
		panel.grid.major.y = element_line(color = "grey90")
	) +
	coord_cartesian(ylim = c(min(plot_data$RCC2_expression), max(plot_data$RCC2_expression)*1.15))

ggsave("RCC2_MBDs.pdf", width = 6, height = 6)


plot_data <- data.frame(oc_subset@assays$RNA$data["RCC2",],oc_subset$phenotype)
colnames(plot_data)=c("RCC2_expression","Phenotype")
plot_data$Phenotype <- as.factor(plot_data$Phenotype)

p_value <- wilcox.test(RCC2_expression ~ Phenotype, data = plot_data)$p.value

sig_label <- ifelse(p_value < 0.001, "***",
										ifelse(p_value < 0.01, "**",
													 ifelse(p_value < 0.05, "*", "")))


ggplot(plot_data, aes(x = Phenotype, y = RCC2_expression)) +
	geom_violin(aes(fill = Phenotype), position = position_dodge(0.9), alpha = 0.7) +
	geom_boxplot(aes(fill = Phenotype),
							 width = 0.15,
							 position = position_dodge(0.9),
							 outlier.shape = NA,
							 color = "black") +
	geom_text(
		aes(x = 1.5, y = max(plot_data$RCC2_expression)*1.1, label = sig_label),
		size = 5, vjust = -0.5
	) +
	scale_fill_manual(values = c("ANA"="#acd5ab", "AA"="#feadac")) +
	labs(x = "Phenotype", y = "RCC2_Expression") +
	theme_classic() +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
		plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
		legend.position = "top",
		panel.border = element_rect(color = "black", fill = NA),
		panel.grid.major.y = element_line(color = "grey90")
	) +
	coord_cartesian(ylim = c(min(plot_data$RCC2_expression), max(plot_data$RCC2_expression)*1.15))

ggsave("RCC2_Phenotype.pdf", width = 6, height = 6)


library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggforce)
library(gghalves)
plot_data <- data.frame(oc_subset$group, oc_subset$GLYCOLYSIS_GLUCONEOGENESIS,
												oc_subset$HIF_1_SIGNALING_PATHWAY,oc_subset$LACTATE_DEHYDROGENASE_ACTIVITY,
												oc_subset$LACTATE_METABOLIC_PROCESS)
colnames(plot_data) <- c("subtypes","GLYCOLYSIS_GLUCONEOGENESIS","HIF_1_SIGNALING_PATHWAY",
												 "LACTATE_DEHYDROGENASE_ACTIVITY","LACTATE_METABOLIC_PROCESS")

plot_data <- plot_data %>%
	pivot_longer(
		cols = -subtypes,
		names_to = "Pathway",
		values_to = "Score"
	) %>%
	mutate(
		Subtype = as.factor(subtypes),
		Pathway = factor(Pathway,
										 levels = c("GLYCOLYSIS_GLUCONEOGENESIS", "HIF_1_SIGNALING_PATHWAY",
										 					 "LACTATE_DEHYDROGENASE_ACTIVITY", "LACTATE_METABOLIC_PROCESS"),
										 labels = c("Glycolysis/Gluconeogenesis", "HIF-1 Signaling",
										 					 "LDH Activity (GO)", "Lactate Metabolism (GO)"))
	)

calc_p_value <- function(data) {
	test <- wilcox.test(Score ~ Subtype, data = data)
	p_value <- test$p.value
	sig_label <- case_when(
		p_value < 0.001 ~ "***",
		p_value < 0.01 ~ "**",
		p_value < 0.05 ~ "*",
		TRUE ~ ""
	)
	return(data.frame(p_value = p_value, sig_label = sig_label))
}

p_values <- plot_data %>%
	group_by(Pathway) %>%
	group_modify(~ calc_p_value(.x)) %>%
	ungroup()


ggplot(plot_data, aes(x = Subtype, y = Score, fill = Subtype)) +
	geom_jitter(
		aes(color = Subtype),
		position = position_jitter(width = 0.1, height = 0.1),
		alpha = 0.6,
		size = 1.5,
		show.legend = FALSE
	) +
	geom_half_violin(
		side = "r",
		alpha = 0.7,
		color = NA,
		trim = FALSE
	) +
	geom_half_boxplot(
		side = "r",
		width = 0.15,
		outlier.shape = NA,
		color = "black",
		alpha = 0.7
	) +
	geom_text(
		data = p_values,
		aes(x = 1.5, y = max(plot_data$Score) * 1.1, label = sig_label),
		size = 5, vjust = -0.5, inherit.aes = FALSE
	) +
	facet_wrap(~ Pathway, nrow = 1, scales = "free_y") +
	scale_fill_manual(values = c("RCC2- Mono/Mac" = "#2c5ca0", "RCC2+ Mono/Mac" = "#88b36e")) +
	scale_color_manual(values = c("RCC2- Mono/Mac" = "#2c5cb0", "RCC2+ Mono/Mac" = "#88b34e")) +
	labs(x = "Subtype", y = "GSVA Score") +
	theme_classic() +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
		axis.title = element_text(face = "bold", size = 12),
		strip.text = element_text(face = "bold", size = 12),
		strip.background = element_blank(),
		legend.position = "top",
		panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
		panel.spacing = unit(1, "lines")
	) +
	coord_cartesian(ylim = c(min(plot_data$Score), max(plot_data$Score) * 1.15))

ggsave("Combined_lactytation.pdf", width = 16, height = 8)

library(clusterProfiler)
library(GSVA)
library(GSEABase)
library(org.Hs.eg.db)
library(KEGGREST)
kegg_pathways <- list(
	"Th1_Th2_diff" = "hsa04658",     # Th1 and Th2 cell differentiation
	"NF_KappaB" = "hsa04064",             # NF-kappa B signaling pathway
	"JAK_STAT" = "hsa04630"          # JAK-STAT signaling pathway
)

gene_sets <- list()

for (name in names(kegg_pathways)) {
	pathway_id <- kegg_pathways[[name]]
	pathway <- keggGet(pathway_id)[[1]]
	genes <- pathway$GENE

	symbols <- genes[seq(2, length(genes), by = 2)]
	symbols <- gsub(";.*", "", symbols)
	symbols <- as.character(na.omit(symbols))

	gene_sets[[name]] <- symbols
}


gsc <- GeneSetCollection(lapply(names(gene_sets), function(name) {
	GeneSet(
		geneIds = gene_sets[[name]],
		collectionType = KEGGCollection(),
		setName = name,
		organism = "Homo sapiens"
	)
}))
exp = oc_subset@assays$RNA$data
gsva_results <- ssgseaParam(
	exp,
	gsc,minSize = 5,normalize = F
)

result <- gsva(gsva_results)
result_scaled <- t(scale(t(result)))
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggforce)
library(gghalves)

plot_data <- data.frame(oc_subset$group, t(result_scaled))
colnames(plot_data)[1] <- "subtypes"
plot_data <- plot_data %>%
	dplyr::select(subtypes,NF_KappaB, JAK_STAT,Th1_Th2_diff) %>%
	pivot_longer(
		cols = -subtypes,
		names_to = "Pathway",
		values_to = "Score"
	) %>%
	mutate(
		Subtype = as.factor(subtypes),
		Pathway = factor(Pathway,
										 levels = c("Th1_Th2_diff", "NF_KappaB", "JAK_STAT"),
										 labels = c("Th1_Th2_diff", "NF_KappaB", "JAK-STAT"))
	)

calc_p_value <- function(data) {
	test <- wilcox.test(Score ~ Subtype, data = data)
	p_value <- test$p.value
	sig_label <- case_when(
		p_value < 0.001 ~ "***",
		p_value < 0.01 ~ "**",
		p_value < 0.05 ~ "*",
		TRUE ~ ""
	)
	return(data.frame(p_value = p_value, sig_label = sig_label))
}

p_values <- plot_data %>%
	group_by(Pathway) %>%
	group_modify(~ calc_p_value(.x)) %>%
	ungroup()


ggplot(plot_data, aes(x = Subtype, y = Score, fill = Subtype)) +
	geom_jitter(
		aes(color = Subtype),
		position = position_jitter(width = 0.1, height = 0.1),
		alpha = 0.6,
		size = 1.5,
		show.legend = FALSE
	) +
	geom_half_violin(
		side = "r",
		alpha = 0.7,
		color = NA,
		trim = FALSE
	) +
	geom_half_boxplot(
		side = "r",
		width = 0.15,
		outlier.shape = NA,
		color = "black",
		alpha = 0.7
	) +
	geom_text(
		data = p_values,
		aes(x = 1.5, y = max(plot_data$Score) * 1.1, label = sig_label),
		size = 5, vjust = -0.5, inherit.aes = FALSE
	) +
	facet_wrap(~ Pathway, nrow = 1, scales = "free_y") +
	scale_fill_manual(values = c("RCC2- Mono/Mac" = "#acd5ab", "RCC2+ Mono/Mac" = "#feadac")) +
	scale_color_manual(values = c("RCC2- Mono/Mac" = "#6a9662", "RCC2+ Mono/Mac" = "#e67c73")) +
	labs(x = "Subtype", y = "GSVA Score") +
	theme_classic() +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
		axis.title = element_text(face = "bold", size = 12),
		strip.text = element_text(face = "bold", size = 12),
		strip.background = element_blank(),
		legend.position = "top",
		panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
		panel.spacing = unit(1, "lines")
	) +
	coord_cartesian(ylim = c(min(plot_data$Score), max(plot_data$Score) * 1.15))

ggsave("Subtype_GSVA_Scores_Modified.pdf", width = 10, height = 6)
