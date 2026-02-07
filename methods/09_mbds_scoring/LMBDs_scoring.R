library(dplyr)
library(tidyverse)
library(GSEABase)
library(GSVA)
library(data.table)
library(SummarizedExperiment)
setwd("D:\\哮喘乳酸化单细胞开篇\\10.MBDs评分")
rm(list = ls())
main_gene = data.table::fread("./output_mRNA_selected.txt")
rt = data.table::fread("./MBDS.xls")
rt = rt %>%
	dplyr::select(avg_log2FC,id) %>%
	column_to_rownames("id")
data = rt[main_gene$x,,drop = F]

up_ragulate = rownames(data)[which(data$avg_log2FC > 0)]
down_ragulate = rownames(data)[which(data$avg_log2FC < 0)]
gene_sets = list("PTMD_up" = up_ragulate,
								 "PTMD_down" = down_ragulate)
load(file.path("D:\\哮喘乳酸化单细胞开篇\\6.一致性聚类","Cluster.Rdata"))
exp <- oc_subset@assays$RNA$data

gsc <- GeneSetCollection(lapply(names(gene_sets), function(name) {
	GeneSet(
		geneIds = gene_sets[[name]],
		collectionType = KEGGCollection(),
		setName = name,
		organism = "Homo sapiens"
	)
}))
gsva_results <- ssgseaParam(
	exp,
	gsc,normalize = F
)

result <- gsva(gsva_results)
result <- t(scale(t(result)))
scores <- t(result)
scores_df <- as.data.frame(scores) %>%
	mutate(PTMDs = PTMD_up)
plot_data <- data.frame(scores_df$PTMDs,
												oc_subset$subtypes)
colnames(plot_data)=c("Scoreing","Subtype")
plot_data$Subtype <- as.factor(plot_data$Subtype)

p_value <- wilcox.test(Scoreing ~ Subtype, data = plot_data)$p.value

sig_label <- ifelse(p_value < 0.001, "***",
										ifelse(p_value < 0.01, "**",
													 ifelse(p_value < 0.05, "*", "")))


ggplot(plot_data, aes(x = Subtype, y = Scoreing)) +
	geom_violin(aes(fill = Subtype), position = position_dodge(0.9), alpha = 0.7, color = NA) +
	geom_boxplot(aes(fill = Subtype),
							 width = 0.3,
							 position = position_dodge(0.9),
							 outlier.shape = NA,
							 color = "black",
							 size = 1,       # Bold box outline
							 notch = TRUE) +
	geom_text(
		aes(x = 1.5, y = max(plot_data$Scoreing)*1.1, label = sig_label),
		size = 5, vjust = -0.5
	) +
	scale_fill_manual(values = c("MBC1"="#acd5ab", "MBC2"="#feadac")) +
	labs(x = "Subtype", y = "LMBCS") +
	theme_classic() +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
		plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
		legend.position = "top",
		panel.border = element_rect(color = "black", fill = NA),
		panel.grid.major.y = element_line(color = "grey90")
	) +
	coord_cartesian(ylim = c(min(plot_data$Scoreing), max(plot_data$Scoreing)*1.15))

ggsave("Subtype_MBDS.pdf", width = 6, height = 6)

plot_data <- data.frame(scores_df$PTMDs,
												oc_subset$phenotype)
colnames(plot_data)=c("Scoreing","Subtype")
plot_data$Subtype <- as.factor(plot_data$Subtype)

p_value <- wilcox.test(Scoreing ~ Subtype, data = plot_data)$p.value

sig_label <- ifelse(p_value < 0.001, "***",
										ifelse(p_value < 0.01, "**",
													 ifelse(p_value < 0.05, "*", "")))


ggplot(plot_data, aes(x = Subtype, y = Scoreing)) +
	geom_violin(aes(fill = Subtype), position = position_dodge(0.9), alpha = 0.7, color = NA) +
	geom_boxplot(aes(fill = Subtype),
							 width = 0.3,
							 position = position_dodge(0.9),
							 outlier.shape = NA,
							 color = "black",
							 size = 1,       # Bold box outline
							 notch = TRUE) +
	geom_text(
		aes(x = 1.5, y = max(plot_data$Scoreing)*1.1, label = sig_label),
		size = 5, vjust = -0.5
	) +
	scale_fill_manual(values = c("ANA"="#acd5ab", "AA"="#feadac")) +
	labs(x = "Phenotype", y = "LMBCS") +
	theme_classic() +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
		plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
		legend.position = "top",
		panel.border = element_rect(color = "black", fill = NA),
		panel.grid.major.y = element_line(color = "grey90")
	) +
	coord_cartesian(ylim = c(min(plot_data$Scoreing), max(plot_data$Scoreing)*1.15))

ggsave("Phenotype_MBDS.pdf", width = 6, height = 6)

library(tidyverse)
library(GSVA)
library(purrr)

process_single_file <- function(input_file) {
	exp <- data.table::fread(input_file, data.table = FALSE) %>%
		column_to_rownames("ID")

	gsva_results <- ssgseaParam(
		as.matrix(exp),
		gsc, normalize = FALSE
	)
	result <- gsva(gsva_results)
	result <- t(scale(t(result)))
	scores <- t(result)
	scores_df <- as.data.frame(scores) %>%
		mutate(PTMDs = PTMD_up)
	Type <- gsub("(.*)_(.*)", "\\2", colnames(exp))
	plot_data <- data.frame(scores_df$PTMDs, Type)
	colnames(plot_data) <- c("Scoreing", "Subtype")
	plot_data$Subtype <- as.factor(plot_data$Subtype)

	p_value <- wilcox.test(Scoreing ~ Subtype, data = plot_data)$p.value
	sig_label <- ifelse(p_value < 0.001, "***",
											ifelse(p_value < 0.01, "**",
														 ifelse(p_value < 0.05, "*", "")))

	p <- ggplot(plot_data, aes(x = Subtype, y = Scoreing)) +
		geom_violin(aes(fill = Subtype), position = position_dodge(0.9), alpha = 0.4, color = NA) +
		geom_boxplot(aes(fill = Subtype),
								 width = 0.3,
								 position = position_dodge(0.9),
								 outlier.shape = NA,
								 color = "black",size = 1) +
		geom_text(
			aes(x = 1.5, y = max(plot_data$Scoreing)*1.1, label = sig_label),
			size = 5, vjust = -0.5
		) +
		scale_fill_manual(values = c("control"="#7aa3f2", "treat"="#ff9a24")) +
		labs(x = "Subtype", y = "LMBDs") +
		theme_classic() +
		theme(
			axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
			plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
			legend.position = "top",
			panel.border = element_rect(color = "black", fill = NA),
			panel.grid.major.y = element_line(color = "grey90")
		) +
		coord_cartesian(ylim = c(min(plot_data$Scoreing), max(plot_data$Scoreing)*1.15))

	output_name <- tools::file_path_sans_ext(basename(input_file)) %>%
		paste0(".pdf")

	ggsave(output_name, plot = p, width = 6, height = 6)
}

list.files(path = "./nor",
					 pattern = "\\.txt$",
					 full.names = TRUE) %>%
	walk(process_single_file)

library(clusterProfiler)
library(GSVA)
library(GSEABase)
library(org.Hs.eg.db)
library(KEGGREST)
kegg_pathways <- list(
	"Th1_Th2_diff" = "hsa04658",     # Th1 and Th2 cell differentiation
	"NFKB" = "hsa04064",             # NF-kappa B signaling pathway
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

gsva_results <- ssgseaParam(
	as.matrix(exp),
	gsc,minSize = 5,normalize = F
)

result <- gsva(gsva_results)
result_scaled <- t(scale(t(result)))

plot_with_statistics <- function(data,
																 x_col = "Variable1",
																 y_col = "Variable2",
																 color_palette = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"),
																 show_fit_line = FALSE,
																 text_size = 14,
																 line_size = 0.5) {

	library(ggplot2)
	library(ggpubr)
	library(ggpointdensity)

	p <- ggplot(data = data,
							mapping = aes_string(x = x_col, y = y_col)) +
		geom_pointdensity(size = 4) +
		scale_color_gradientn(colors = color_palette) +
		theme_classic(base_size = text_size) +
		labs(x = x_col, y = y_col) +
		theme(
			axis.title = element_text(size = text_size, face = "bold"),
			axis.text = element_text(size = text_size - 2, face = "bold")
		)

	if (show_fit_line) {
		p <- p +
			geom_smooth(method = "lm", formula = y ~ x, color = "black", size = line_size) +
			stat_cor(
				method = "spearman",
				aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
				label.x = -1,
				size = text_size / 2.8,
				fontface = "bold"
			)
	}

	return(p)
}

if (!dir.exists("cor")) {
	dir.create("cor")
}

for (pathway_name in names(kegg_pathways)) {

	pathway_scores <- result_scaled[pathway_name, ]

	data <- data.frame(
		LMBDs = scores_df$PTMDs,
		Pathway_Score = pathway_scores
	)

	colnames(data) <- c("LMBDs", pathway_name)

	plot <- plot_with_statistics(data,
															 x_col = "LMBDs",
															 y_col = pathway_name,
															 show_fit_line = TRUE)

	ggsave(paste0("cor/cor_", pathway_name, ".pdf"), plot = plot, width = 6, height = 6)
}
