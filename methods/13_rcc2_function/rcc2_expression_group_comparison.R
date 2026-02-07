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
setwd("D:\\哮喘乳酸化单细胞开篇\\17.RCC2基因功能探究")
mm <- subset(oc,subset = celltype == "Mono/Mac") %>% NormalizeData() %>%
	FindVariableFeatures(selection.method = "vst" ,nfeatures = 1500 ,verbose = T)
mm$group <- ifelse(
	mm@assays$RNA$counts["RCC2", ] > median(mm@assays$RNA$counts["RCC2", ]),
	"RCC2high Mono/Mac",
	"RCC2low Mono/Mac")
table(mm$group)
library(SCP)
oc[["RNA"]] <- as(object = oc[["RNA"]], Class = "Assay")
pdf(file="RCC2_DIM.pdf", width=6, height=6)
FeatureDimPlot(oc,
							 features = "RCC2", reduction = "umap", label = TRUE,
							 cells.highlight = colnames(oc)[oc$celltype == "Mono/Mac"]
)
dev.off()
pdf(file="RCC2_vio.pdf", width=12, height=6)
FeatureStatPlot(
	srt = oc,
	group.by = "celltype",
	stat.by = c("RCC2"),
	add_box = TRUE)
dev.off()

output_dir <- "MetaAnalysis_Results"
type_groups <- c("AA", "ANA")
color_scheme <- c("#b71b23", "#0c695e")
plot_data <- oc@meta.data %>%
	rownames_to_column("CellID") %>%
	filter(phenotype %in% type_groups) %>%
	mutate(Group = factor(phenotype, levels = type_groups)) %>%
	dplyr::select(CellID, Group, celltype) %>%
	mutate(Expression = oc@assays$RNA$data["RCC2",])

stat_data <- plot_data %>%
	group_by(celltype) %>%
	summarise(
		p_value = wilcox.test(Expression ~ Group)$p.value,
		y_max = max(Expression) * 1.1,
		.groups = "drop"
	) %>%
	mutate(
		p_label = ifelse(p_value < 0.001,
										 "***P < 0.001",
										 paste0("*P = ", signif(p_value, 3)))
	)

mean_sd <- plot_data %>%
	group_by(celltype, Group) %>%
	summarise(
		mean = mean(Expression, na.rm = TRUE),
		sd = sd(Expression, na.rm = TRUE),
		y_pos = quantile(Expression, 0.95),
		.groups = "drop"
	) %>%
	mutate(label = sprintf("Mean: %.2f\nSD: %.2f", mean, sd))

ggplot(plot_data, aes(x = celltype, y = Expression)) +
	geom_violin(aes(fill = Group),
							position = position_dodge(0.8),
							width = 0.7,
							alpha = 0.5,
							trim = FALSE) +
	geom_boxplot(aes(group = interaction(celltype, Group)),
							 position = position_dodge(0.8),
							 width = 0.2,
							 outlier.shape = NA,
							 alpha = 0.8) +
	geom_jitter(aes(color = Group),
							position = position_jitterdodge(jitter.width = 0.2,
																							dodge.width = 0.8),
							size = 0.8,
							alpha = 0.3) +
	geom_text(data = stat_data,
						aes(x = celltype, y = y_max, label = p_label),
						size = 4.5,
						vjust = -0.5,
						inherit.aes = FALSE) +
	geom_text(data = mean_sd,
						aes(x = celltype, y = y_pos,
								label = label, group = Group),
						position = position_dodge(0.8),
						size = 3,
						color = "black",
						vjust = 0) +
	scale_fill_manual(values = color_scheme) +
	scale_color_manual(values = color_scheme) +
	labs(title = "Expression Comparison Across Cell Types",
			 y = "Pathway Combined Score",
			 x = "Cell Types") +
	theme_classic(base_size = 14) +
	theme(
		plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
		axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
		axis.text.y = element_text(size = 11),
		legend.position = "top",
		legend.title = element_blank(),
		panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
		plot.margin = margin(1, 1, 1, 1, "cm")
	) +
	scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))

ggsave(filename = file.path(output_dir, "CellType_Scoreing_Comparison.pdf"),
			 width = 16,
			 height = 7,
			 dpi = 300)

library(IOBR)
data(signature_metabolism)
cor_results <- data.frame(
	PTM = character(),
	Correlation = numeric(),
	P_value = numeric(),
	stringsAsFactors = FALSE
)

for (ptm in names(signature_metabolism)) {
	if (!ptm %in% colnames(mm@meta.data)) {
		warning(paste("Column", ptm, "not found in mm. Skipping."))
		next
	}

	data <- data.frame(
		Expression = mm@assays$RNA$data["RCC2",],
		Score = mm@meta.data[[ptm]],
		Phenotype = mm$phenotype
	)
	colnames(data) <- c("Expression", "Score", "Phenotype")

	data$Phenotype <- as.factor(data$Phenotype)

	cor_test <- cor.test(data$Expression, data$Score, method = "spearman")
	cor_value <- round(cor_test$estimate, 3)
	p_value <- round(cor_test$p.value, 5)

	cor_results <- rbind(cor_results, data.frame(
		PTM = ptm,
		Correlation = cor_value,
		P_value = p_value
	))
}
library(dplyr)
library(ggplot2)

top3_ptms <- cor_results %>%
	filter(!is.na(Correlation)) %>%
	arrange(desc(abs(Correlation))) %>%
	slice_head(n = 3) %>%
	pull(PTM)

gene_symbol <- "RCC2"

for (ptm in top3_ptms) {
	data <- data.frame(
		Expression = as.numeric(mm@assays$RNA$data[gene_symbol, ]),
		Score      = mm@meta.data[[ptm]],
		Phenotype  = mm$phenotype
	)
	colnames(data) <- c("Expression", "Score", "Phenotype")
	data$Phenotype <- as.factor(data$Phenotype)
	data <- data[complete.cases(data), ]
	if (nrow(data) < 3) next

	rrow      <- cor_results[cor_results$PTM == ptm, , drop = FALSE]
	cor_value <- round(rrow$Correlation[1], 3)
	p_value   <- rrow$P_value[1]

	safe_filename <- gsub("[^[:alnum:]_]", "_", ptm)
	pdf_file      <- paste0("cor_Expression_", safe_filename, ".pdf")
	cor_label     <- paste0("cor = ", cor_value, "\n", "p = ", format.pval(p_value, digits = 3))

	ann_x <- min(data$Expression, na.rm = TRUE) + 0.05 * diff(range(data$Expression, na.rm = TRUE))
	ann_y <- max(data$Score,      na.rm = TRUE) - 0.05 * diff(range(data$Score,      na.rm = TRUE))

	p <- ggplot(data, aes(x = Expression, y = Score, color = Phenotype)) +
		geom_point(size = 3.5, alpha = 1) +
		geom_smooth(
			method = "lm",
			color  = "#5b82a0",
			linetype = "solid",
			size   = 3,
			se     = TRUE,
			fill   = "#5a5676",
			alpha  = 0.3
		) +
		annotate("text", x = ann_x, y = ann_y,
						 label = cor_label, hjust = 0, color = "black", size = 5) +
		labs(
			title = paste("Expression vs.", ptm),
			x     = paste("Expression of", gene_symbol),
			y     = ptm
		) +
		theme_bw(base_size = 14) +
		theme(
			legend.position  = "none",
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border     = element_rect(color = "black", size = 1.2),
			plot.title       = element_text(hjust = 0.5, size = 16, face = "bold"),
			axis.title       = element_text(size = 14),
			axis.text        = element_text(size = 12)
		)

	if (length(levels(data$Phenotype)) == 2) {
		p <- p + scale_color_manual(values = c("#0c695e", "#b81b23"))
	}

	ggsave(pdf_file, p, width = 7, height = 7, dpi = 300)
}

write.table(cor_results, file = "cor_results.xls",
						sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggforce)
library(gghalves)
plot_data <- data.frame(mm$group, mm$GLYCOLYSIS_GLUCONEOGENESIS,
												mm$HIF_1_SIGNALING_PATHWAY,mm$LACTATE_DEHYDROGENASE_ACTIVITY,
												mm$LACTATE_METABOLIC_PROCESS)
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
	scale_fill_manual(values = c("RCC2low Mono/Mac" = "#2c5ca0", "RCC2high Mono/Mac" = "#88b36e")) +
	scale_color_manual(values = c("RCC2low Mono/Mac" = "#2c5cb0", "RCC2high Mono/Mac" = "#88b34e")) +
	labs(x = "Subtype", y = "Score") +
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
#####irGSEA####
library(RcppML)
library(irGSEA)
Idents(mm)=mm$group
pbmc3k.final3=mm
DefaultAssay(mm) <- "RNA"
pbmc3k.final3 <- SeuratObject::UpdateSeuratObject(object = mm)
pbmc3k.final3 <- irGSEA.score(object = pbmc3k.final3,assay = "RNA",
															slot = "data", seeds = 123, ncores = 1,
															min.cells = 3, min.feature = 0,
															custom = F, geneset = NULL, msigdb = T,
															species = "Homo sapiens", category = "H",
															subcategory = NULL, geneid = "symbol",
															method = c("AUCell","UCell","singscore","ssgsea","viper"),
															aucell.MaxRank = NULL, ucell.MaxRank = NULL,
															kcdf = 'Gaussian')
result.dge <- irGSEA.integrate(object = pbmc3k.final3,
															 group.by = "group",
															 method = c("AUCell","UCell","singscore","ssgsea","viper"))

pdf(file="heatmap.plot.H.pdf", width=12, height=6)
irGSEA.heatmap(object = result.dge,
							 method = "RRA",
							 show.geneset = NULL)
dev.off()

library(IOBR)
library(GSEABase)
signatures_list <- lapply(signature_metabolism, function(gs) {
	geneIds(gs)
})
names(signatures_list) <- names(signatures_KEGG_metab)

pbmc3k.final3 <- irGSEA.score(object = pbmc3k.final3,assay = "RNA",
															slot = "data", seeds = 123, ncores = 1,
															min.cells = 3, min.feature = 0,
															custom = T, geneset = signature_metabolism,
															subcategory = NULL, geneid = "symbol",
															method = c("AUCell","UCell","singscore","ssgsea","viper"),
															kcdf = 'Gaussian')

result.dge <- irGSEA.integrate(object = pbmc3k.final3,
															 group.by = "group",
															 method = c("AUCell","UCell","singscore","ssgsea", "viper"))
pdf(file="heatmap.plot_metabolism.pdf", width=12, height=6)
irGSEA.heatmap(object = result.dge,
							 method = "RRA",
							 show.geneset = NULL)
dev.off()


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
exp = mm@assays$RNA$data
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

plot_data <- data.frame(mm$group, t(result_scaled))
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
	scale_fill_manual(values = c("RCC2low Mono/Mac" = "#acd5ab", "RCC2high Mono/Mac" = "#feadac")) +
	scale_color_manual(values = c("RCC2low Mono/Mac" = "#6a9662", "RCC2high Mono/Mac" = "#e67c73")) +
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
