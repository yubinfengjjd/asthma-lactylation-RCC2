library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
setwd("D:\\哮喘乳酸化单细胞开篇\\5.hdWGCNA")
theme_set(theme_cowplot())
oc <- SetupForWGCNA(oc,  gene_select = "fraction",
														fraction = 0.05,
														wgcna_name = "tutorial"
														)

oc <- MetacellsByGroups(
	seurat_obj = oc,
	group.by = c("celltype","orig.ident"),
	k = 25,
	reduction = "harmony",
	slot = 'counts',
	ident.group = 'celltype'
	)

oc <- NormalizeMetacells(oc)
metacell_obj <- GetMetacellObject(oc)
table(metacell_obj$celltype,metacell_obj$orig.ident)

oc <- SetDatExpr(oc,group_name = "Mono/Mac",
												 group.by = 'celltype',
												 assay = 'RNA',
								         use_metacells = TRUE,
												 slot = 'data'
												 )

oc <- TestSoftPowers(oc,
										 powers = c(seq(1, 10, by = 1),
														 seq(12, 30, by = 2))
														 )
power_table <- GetPowerTable(oc)
pdf("power1.pdf", 6, 6)

plot(power_table[,1], -sign(power_table[,3]) * power_table[,2],
		 xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
		 main = paste("Scale independence"))
points(power_table[,1], -sign(power_table[,3]) * power_table[,2], pch = 19, col = "#7AC5CD", cex = 1.5)
abline(h = 0.8, col = "red")
grid()

dev.off()

pdf("power2.pdf", 6, 6)

plot(power_table[,1], power_table[,5],
		 xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
		 main = paste("Mean connectivity"))
points(power_table[,1], power_table[,5], pch = 19, col = "#7AC5CD", cex = 1.5)
grid()

dev.off()
oc <- ConstructNetwork(oc,soft_power = 5,
															 setDatExpr = FALSE,
															 corType = "pearson",
															 networkType = "signed",
															 TOMType = "signed",
															 minModuleSize = 50,
											         detectCutHeight=0.995,
															 mergeCutHeight = 0.25,
															 tom_outdir = "TOM",
															 tom_name = 'Epi'
															 )
pdf(file = "Dendrogram_tree.PDF",width = 8,height = 6)
PlotDendrogram(oc, main='Epi hdWGCNA Dendrogram')
dev.off()
oc <- ScaleData(oc, features=VariableFeatures(oc))
oc <- ModuleEigengenes(oc,scale.model.use = "linear",
															 assay = NULL,
															 pc_dim = 1,
															 group.by.vars = "orig.ident"
															 )
hMEs <- GetMEs(oc)
MEs <- GetMEs(oc, harmonized = FALSE)
colnames(oc@misc[["tutorial"]][["hMEs"]])
seurat_obj <- ModuleConnectivity(oc,
																 group.by = 'celltype',
																 corFnc = "bicor",
																 corOptions = "use='p'",
																 harmonized = TRUE,
																 assay = NULL,
																 slot = "data",
																 group_name = 'Mono/Mac'
																 )

seurat_obj <- ResetModuleNames(seurat_obj,new_name = "Epi")
pdf(file = "KMEs.PDF",width = 12,height = 6)
PlotKMEs(seurat_obj, ncol=5)#by kME for each module
dev.off()
modules <- GetModules(seurat_obj)
p <- PlotKMEs(seurat_obj,ncol = 3,
							n_hubs = 20,
							text_size = 2,
							plot_widths = c(3, 2)
							)
p
seurat_obj$Scoreing<-as.numeric(seurat_obj$Scoreing)
seurat_obj$phenotype <- ifelse(seurat_obj$phenotype == "AA", 1, 0)
seurat_obj$phenotype = as.factor(seurat_obj$phenotype)
cur_traits<-c("Scoreing","phenotype")

str(seurat_obj@meta.data[,cur_traits])
seurat_obj<-ModuleTraitCorrelation(
	seurat_obj,
	traits=cur_traits,
	features="hMEs",
	cor_method="pearson",
	group.by='celltype'
)
mt_cor<-GetModuleTraitCorrelation(seurat_obj)

library(ggplot2)
library(reshape2)
library(dplyr)
library(ggnewscale)

moduleTraitCor <- t(mt_cor[["cor"]][["Mono/Mac"]])
moduleTraitPvalue <- t(mt_cor[["pval"]][["Mono/Mac"]])

cor_long <- melt(moduleTraitCor, varnames = c("Module", "Group"))
pvalue_long <- melt(moduleTraitPvalue, varnames = c("Module", "Group"))

data_long <- merge(cor_long, pvalue_long, by = c("Module", "Group"))
names(data_long) <- c("Module", "Group", "Correlation", "Pvalue")

data_long$Group_num <- as.numeric(factor(data_long$Group))
data_long$Module_num <- as.numeric(factor(data_long$Module))

triangle <- function(pairs, type = "upper") {
	x = c(0, 0, 1)
	y = c(0, 1, 1)
	if (type == "lower") {
		x = c(0, 1, 1)
		y = c(0, 0, 1)
	}
	mat = do.call(
		rbind,
		apply(pairs, 1, function (row) {
			a = row[1]
			b = row[2]
			data.frame(
				x = x + a - 0.5,
				y = y + b - 0.5,
				group = paste(a, b, sep = "-")
			)
		}))
	return(mat)
}

pairs <- expand.grid(Group_num = unique(data_long$Group_num), Module_num = unique(data_long$Module_num))

upper <- triangle(pairs, type = "upper")
lower <- triangle(pairs, type = "lower")

upper_lower <- cbind(upper, lower[, c("x", "y")])
colnames(upper_lower) <- c("upper.x", "upper.y", "group", "lower.x", "lower.y")

data_long <- data_long %>%
	mutate(group = paste(Group_num, Module_num, sep = "-")) %>%
	left_join(upper_lower, by = "group")

data_long$text_x <- data_long$Group_num
data_long$text_y <- data_long$Module_num
data_long$cor_x <- data_long$text_x - 0.2
data_long$cor_y <- data_long$text_y + 0.2
data_long$pval_x <- data_long$text_x + 0.2
data_long$pval_y <- data_long$text_y - 0.2

p <- ggplot(data_long, aes(x = Group_num, y = Module_num)) +
	geom_polygon(aes(x = upper.x, y = upper.y, fill = Correlation, group = group), color = "black") +
	scale_fill_gradient(low = "#c79cc3", high = "#ffa288", name = "Correlation", limits = c(-1, 1)) +
	new_scale_fill() +
	geom_polygon(aes(x = lower.x, y = lower.y, fill = Pvalue, group = group), color = "black") +
	scale_fill_gradient(low = "#FFFAFA", high = "#0475b7", name = "P-value", limits = c(0, 0.5)) +
	geom_text(aes(x = cor_x, y = cor_y, label = round(Correlation, 2)), vjust = 0.5, hjust = 0.5, size = 3, color = "black") +
	geom_text(aes(x = pval_x, y = pval_y, label = formatC(Pvalue, format = "e", digits = 2)), vjust = 0.5, hjust = 0.5, size = 3, color = "black") +
	scale_x_continuous(breaks = unique(data_long$Group_num), labels = unique(data_long$Group)) +
	scale_y_continuous(breaks = unique(data_long$Module_num), labels = unique(data_long$Module)) +
	labs(x = "Group", y = "Module") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 10),
				axis.text.y = element_text(size = 10),
				legend.position = "right") +
	ggtitle("Module-trait relationships")

ggsave("Heatmap.pdf", plot = p, width = 6, height = 6)

library(UCell)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(colorspace)

modules_score <- modules %>%
	filter(!grepl("grey", module, ignore.case = TRUE)) %>%
	group_by(module, color) %>%
	summarise() %>%
	ungroup()

color_hex <- c(
	turquoise = "#1F77B4",
	blue = "#0000FF",
	yellow = "#FFD700",
	brown = "#8B4513",
	green = "#00FF00",
	grey = "#808080"
)

modules_score <- modules_score %>%
	mutate(
		color_hex = color_hex[tolower(color)],
		module_label = paste0(module, " (", color, ")")
	)

seurat_obj <- ModuleExprScore(seurat_obj, n_genes = 25, method='UCell')

df <- data.frame(
	seurat_obj@misc$tutorial$module_scores,
	seurat_obj@reductions[["umap"]]@cell.embeddings,
	seurat_obj@meta.data
)

valid_modules <- intersect(modules_score$module, colnames(df))
if (length(valid_modules) == 0) {
	stop("模块名称与评分列不匹配，请检查以下可能原因：
       1. ModuleExprScore生成的列名是否包含模块编号（如Epi1）
       2. modules数据框中的module列是否对应这些编号")
}

color_mapping <- setNames(modules_score$color_hex, modules_score$module)

class_avg <- df %>%
	group_by(celltype) %>%
	summarise(umap_1 = median(umap_1),
						umap_2 = median(umap_2))

for (method in valid_modules) {
	current_module <- modules_score %>% filter(module == method)
	current_color <- current_module$color_hex

	tryCatch({
		gradient_colors <- c("#f0f0f0", lighten(current_color, 0.4), current_color)
	}, error = function(e) {
		message("颜色渐变生成失败，使用默认渐变")
		gradient_colors <- c("#f0f0f0", "#C0C0C0", current_color)
	})

	pdf(file = paste0(method, "_umap.pdf"), width = 7, height = 6)

	p <- ggplot(df, aes(x = umap_1, y = umap_2)) +
		geom_point(
			aes_string(color = method),
			size = 0.8,
			alpha = 0.9,
			shape = 16,
			stroke = 0
		) +
		ggrepel::geom_label_repel(
			aes(label = celltype),
			data = class_avg,
			size = 3.5,
			box.padding = 0.4,
			label.padding = 0.1,
			segment.color = "gray60",
			label.size = 0.3,
			color = "black",
			fill = alpha("white", 0.8)
		) +
		scale_color_gradientn(
			colors = gradient_colors,
			breaks = c(0, 0.5, 1),
			labels = scales::percent_format(),
			name = "Module Score",
			guide = guide_colorbar(
				barwidth = 1,
				barheight = 8,
				title.position = "top",
				title.hjust = 0.5,
				frame.colour = "black",
				ticks.colour = "black"
			)
		) +
		theme_classic(base_size = 12) +
		labs(
			x = "UMAP1",
			y = "UMAP2",
			title = current_module$module_label,
			subtitle = paste("Signature Score |", nrow(filter(modules, module == method)), "genes")
		) +
		theme(
			axis.title = element_text(face = "bold", colour = "black"),
			axis.text = element_text(color = "gray40"),
			plot.title = element_text(hjust = 0.5, face = "bold", size = 14, colour = current_color),
			plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
			panel.border = element_rect(color = "gray30", fill = NA, size = 0.8),
			plot.background = element_rect(fill = "white"),
			legend.background = element_rect(fill = "white")
		)

	print(p)
	dev.off()
}

hdWGCNAgene=modules$gene_name[modules$module == "Epi1"]
write.table(hdWGCNAgene, file = "hdWGCNAgene.txt", sep = "\t", quote = F,
						col.names = F,row.names = F)
saveRDS(oc, "hdwgcna.rds")
