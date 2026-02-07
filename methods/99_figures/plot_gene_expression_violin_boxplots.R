library(ggplot2)
library(tidyr)
library(ggsci)
library(ggpubr)
library(tibble)
library(tidyverse)
setwd("D:\\哮喘乳酸化单细胞开篇\\12.箱式图")
gene_list_file <- "./maingene.txt"
datasets <- c(
	"GSE43696" = "./GSE43696.normalize.txt",
	"GSE63142" = "./GSE63142.normalize.txt"
)
output_dir <- "Results"
control_suffix <- "control"
treat_suffix <- "treat"
color_palette <- c("#2874C5", "#EABF00")  # 颜色设置（对照，处理）

preprocess_data <- function(expFile){
	rt <- read.table(file.path("nor", expFile), sep = "\t", header = T, check.names = F)
	exp <- as.matrix(rt[, -1])
	rownames(exp) <- rt[, 1]
	return(exp)
}

generate_plot <- function(exp_matrix, gene, dataset_name, color_pal){
	dir.create(file.path(output_dir, dataset_name), showWarnings = FALSE, recursive = TRUE)

	my_comparisons <- list(
		c("control","treat")
	)

	# p_val <- wilcox.test(count ~ type, data = exp_matrix)$p.value

	p <- ggplot(exp_matrix, aes(x = type, y = count, fill = type)) +
		scale_fill_manual(values = color_pal) +
		geom_violin(
			alpha = 0.4,
			position = position_dodge(width = .75),
			size = 0.8,
			color = "black"
		) +
		geom_boxplot(
			notch = TRUE,
			outlier.size = -1,
			lwd = 0.8,
			alpha = 0.7,
			color = "black",
			position = position_dodge(width = .75)
		) +
		geom_point(
			aes(color = NULL),
			shape = 21,
			size = 2,
			color = "black",
			position = position_jitterdodge(dodge.width = .75)
		) +
		stat_compare_means(
			method = "wilcox.test",
			hide.ns = FALSE,
			comparisons = my_comparisons,
			label = "p.signif"
		) +
		theme_bw() +
		labs(
			title = paste0(dataset_name, ": ", gene),
			x = "",
			y = "Normalized Expression"
		) +
		theme(
			panel.background     = element_blank(),
			panel.grid           = element_blank(),
			axis.text.x          = element_text(size = 12, color = "black"),
			axis.text.y          = element_text(size = 12, color = "black"),
			axis.title           = element_text(size = 12),
			axis.ticks           = element_line(size = 0.2, color = "black"),
			axis.ticks.length    = unit(0.2, "cm"),
			legend.position      = "none",
			plot.title           = element_text(face = "bold", hjust = 0.5)
		)

	ggsave(
		filename = file.path(output_dir, dataset_name, paste0(gene, "_", dataset_name, ".pdf")),
		plot     = p,
		width    = 8,
		height   = 6
	)
}


suppressPackageStartupMessages({
	library(tidyr)
	library(ggplot2)
	library(dplyr)
})

gene_list <- read.table(gene_list_file)$V1 %>% as.character()

for(dataset_name in names(datasets)){
	exp_data <- preprocess_data(datasets[dataset_name])

	target_genes <- intersect(gene_list, rownames(exp_data))

	plot_data <- exp_data[target_genes, , drop = F] %>%
		as.data.frame() %>%
		tibble::rownames_to_column("gene") %>%
		pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
		mutate(
			type = case_when(
				grepl(paste0("_", control_suffix, "$"), sample) ~ "control",
				grepl(paste0("_", treat_suffix, "$"), sample) ~ "treat",
				TRUE ~ NA_character_
			)
		) %>%
		filter(!is.na(type))

	for(gene in unique(plot_data$gene)){
		gene_data <- filter(plot_data, gene == !!gene)
		generate_plot(gene_data, gene, dataset_name, color_palette)
	}
}

message("\n[运行状态] 分析已完成！结果保存至：",
				normalizePath(output_dir), "\n")
