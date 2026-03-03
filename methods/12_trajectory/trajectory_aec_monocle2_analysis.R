setwd("D:\\哮喘乳酸化单细胞开篇\\16.细胞轨迹")
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
library(tidyr)
library(ggforce)
library(ggrastr)
library(viridis)
library(gridExtra)
library(ggnewscale)
library(org.Mm.eg.db)
library(ComplexHeatmap)
library(ClusterGVis)
library(monocle)
rm(list = ls())
##### Epithelial cell re-clustering #####
library(Seurat)
library(sctransform)
library(dplyr)
library(harmony)
library(clustree)
# Get metadata
meta_data <- oc@meta.data

# Group by cell type and compute per-type median
meta_data <- meta_data %>%
  group_by(celltype) %>%
  mutate(
    median = median(Retinoic_Acid_Metabolism, na.rm = TRUE),  # per-celltype median
    group = if_else(Retinoic_Acid_Metabolism > median,        # > median => high, else low
                    paste0("high_", celltype),
                    paste0("low_", celltype))
  ) %>%
  ungroup()

# Add grouping back to the Seurat object
oc$celltype <- meta_data$group
# Inspect group labels
table(oc$celltype)
mm <- subset(oc, subset = celltype == "high_AEC") %>% NormalizeData() %>%
	FindVariableFeatures(selection.method = "vst" ,nfeatures = 1500 ,verbose = T)

# Dimensionality reduction and graph construction
mm <- RunPCA(mm, npcs = 30, verbose = FALSE)
# Optionally inspect elbow plot to choose PCs
# ElbowPlot(mm, ndims = 50)
dims_use <- 1:20

if ("orig.ident" %in% colnames(mm@meta.data)) {
   mm <- RunHarmony(mm, group.by.vars = "orig.ident", dims.use = dims_use)
   emb_use <- "harmony"
 } else {
   emb_use <- "pca"
 }

# emb_use <- "pca"  # if Harmony is not used

mm <- RunUMAP(mm, reduction = emb_use, dims = dims_use)
mm <- FindNeighbors(mm, reduction = emb_use, dims = dims_use)

# Resolution grid clustering (choose based on markers)
for (res in c(0.2, 0.4, 0.6, 0.8, 1.0)) {
	mm <- FindClusters(mm, resolution = res, verbose = FALSE)
}

library(SCP)
library(DDRTree)
library(monocle)
mm[["RNA"]] <- as(object = mm[["RNA"]], Class = "Assay")
DefaultAssay(mm) <- "RNA"
Idents(mm) = mm$phenotype
deg = FindMarkers(object = mm,ident.1 = "AA",ident.2 = "ANA")
deg = deg %>% dplyr::filter(p_val_adj < 0.05)
mm=RunMonocle2(mm,features = VariableFeatures(mm))

pdf("State_umap.PDF", width = 6, height = 6)
CellDimPlot(
	srt = mm, group.by = c("Monocle2_State"), reduction = "umap", theme_use = "theme_blank"
)
dev.off()
# Assume the Seurat object contains Monocle2 pseudotime
# Genes to visualize
source("./DynamicPlot.R")
feats <- c("RCC2")   # edit as needed

# 1) Use a unified brown line / ribbon
brown_line <- "#A6793D"   # brown
line_palcolor <- setNames(rep(brown_line, length(feats)), feats)

# 2) Point colors by group.by (extend as needed)
grp <- "Monocle2_State"          # edit as needed
lv  <- levels(factor(mm@meta.data[[grp]]))
tableau_like <- c(
	"#4E79A7",  # 蓝  (blue)
	"#59A14F",  # 绿  (green)
	"#F28E2B",  # 橙  (orange)
	"#E15759",  # 红  (red)
	"#76B7B2",  # 青  (teal)
	"#EDC948",  # yellow
	"#B07AA1",  # purple
	"#FF9DA7",  # pink
	"#9C755F",  # brown
	"#BAB0AC"   # grey
)
point_palcolor <- setNames(tableau_like[seq_along(lv)], lv)

# 3) 出图（和示意图一致：棕色曲线、半透明带状区、散点按分组上色，底部 rug）
p <- DynamicPlotMonocle2(
	srt = mm,
	pseudotime = "Monocle2_Pseudotime",
	features = feats,
	group.by = grp,
	slot = "data",
	exp_method = "raw",
	add_interval = TRUE,
	add_point = TRUE,
	pt.size = 1,
	line.size = 1.2,
	legend.position = "bottom",
	# 关键：传入上面两组颜色
	line_palcolor  = line_palcolor,
	point_palcolor = point_palcolor,
	# 可选：整体风格
	theme_use  = theme_scp,
	theme_args = list(base_size = 12)
)

# 保存
ggplot2::ggsave("pseudotime_style_like_example.pdf", p, width = 8, height = 6)
#####SCP_拟时序####
mm@tools$Monocle2$cds
save(mm,file = "monocle.Rdata")
trajectory <- mm@tools$Monocle2$trajectory
pdf("State.PDF", width = 6.5, height = 5.5)
CellDimPlot(mm, group.by = "Monocle2_State", reduction = "DDRTree", label = TRUE, theme_use = "theme_blank") + trajectory
dev.off()
pdf("Pseudotime.PDF", width = 6.5, height = 5.5)
FeatureDimPlot(mm, features = "Monocle2_Pseudotime", reduction = "DDRTree", theme_use = "theme_blank")
dev.off()
pdf("RCC2.PDF", width = 6.5, height = 5.5)
FeatureDimPlot(mm, features = "RCC2", reduction = "DDRTree", theme_use = "theme_blank")
dev.off()
pdf("LMBDs.PDF", width = 6.5, height = 5.5)
FeatureDimPlot(mm, features = "LMBDs", reduction = "DDRTree", theme_use = "theme_blank")
dev.off()

mm$LMBDs
cds_seurat = mm@tools$Monocle2$cds

pdf("plot_genes_in_pseudotime.PDF", width = 6, height = 6)

cds_subset <- cds_seurat[c("RCC2"), ]

# 按State分组计算拟时序中位数，并排序
state_order <- pData(cds_subset) %>%
	group_by(State) %>%
	summarise(median_pseudotime = median(Pseudotime)) %>%
	arrange(median_pseudotime) %>%
	pull(State)

# 将State转换为有序因子
pData(cds_subset)$State <- factor(pData(cds_subset)$State, 
																	levels = state_order)

# 创建基础绘图对象
p <- plot_genes_in_pseudotime(
	cds_subset,
	color_by = "State",
	cell_size = 1.8
) 

# 添加美化层
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
		y = "Expression",
		title = "RCC2 Expression Dynamics"
	)

dev.off()
#####相关性#####
library(IOBR)
library(ggplot2)
library(dplyr)
library(ggsci)
data(signature_metabolism)

# 初始化结果数据框
cor_results <- data.frame(
	PTM = character(),
	Correlation = numeric(),
	P_value = numeric(),
	stringsAsFactors = FALSE
)

# 第一阶段：计算所有通路的相关性
for (ptm in names(signature_metabolism)) {
	if (!ptm %in% colnames(mm@meta.data)) {
		warning(paste("Column", ptm, "not found in mm. Skipping."))
		next
	}
	
	# 计算Spearman相关性
	cor_test <- cor.test(mm$Monocle2_Pseudotime, mm@meta.data[[ptm]], method = "spearman")
	
	# 添加到结果数据框
	cor_results <- rbind(
		cor_results,
		data.frame(
			PTM = ptm,
			Correlation = cor_test$estimate,
			P_value = cor_test$p.value
		)
	)
}

# 保存所有相关性结果
write.csv(cor_results, "pseudotime_correlation_results.csv", row.names = FALSE)

# 第二阶段：选择最显著的正负相关通路
top_n <- 1  # 选择最强的正负相关各5个

# 选择正相关最强的通路
top_positive <- cor_results %>%
	filter(Correlation > 0) %>%
	arrange(desc(Correlation)) %>%
	slice_head(n = top_n)

# 选择负相关最强的通路
top_negative <- cor_results %>%
	filter(Correlation < 0) %>%
	arrange(Correlation) %>%
	slice_head(n = top_n)

# 合并要绘图的通路
selected_ptms <- rbind(top_positive, top_negative)

# 第三阶段：仅为选定的通路绘图
for (i in 1:nrow(selected_ptms)) {
	ptm <- selected_ptms$PTM[i]
	cor_value <- selected_ptms$Correlation[i]
	p_value <- selected_ptms$P_value[i]
	
	# 创建数据框
	data <- data.frame(
		Pseudotime = mm$Monocle2_Pseudotime,
		Score = mm@meta.data[[ptm]],
		State = mm$Monocle2_State
	)
	colnames(data) <- c("Pseudotime", "Score", "State")
	
	# 确保State是因子
	data$State <- as.factor(data$State)
	
	# 创建安全文件名
	safe_filename <- gsub("[^[:alnum:]_]", "_", ptm)
	pdf_file <- paste0("cor_Pseudotime_", safe_filename, ".pdf")
	
	# 生成标注文本
	cor_label <- paste0("cor = ", round(cor_value, 3), "\n", "p = ", format.pval(p_value, digits = 3))
	
	# 创建美观的绘图
	p <- ggplot(data, aes(x = Pseudotime, y = Score, color = State)) +
		geom_point(
			alpha = 0.7,
			size = 3.5,
			shape = 21,
			aes(fill = State),
			color = "black",
			stroke = 0.3
		) +
		geom_smooth(
			aes(group = 1),
			method = "loess",
			formula = y ~ x,
			color = ifelse(cor_value > 0, "#e41a1c", "#377eb8"),
			se = TRUE,
			fill = ifelse(cor_value > 0, "#fbb4ae", "#b3cde3"),
			linewidth = 1.8,
			alpha = 0.2
		) +
		scale_fill_npg() +
		scale_color_npg() +
		theme_classic(base_size = 14) +
		labs(
			x = "Pseudotime",
			y = paste(ptm, "Activity"),
			title = ptm,
			subtitle = paste(ifelse(cor_value > 0, "Positive", "Negative"), "correlation with pseudotime"),
			color = "State",
			fill = "State"
		) +
		theme(
			plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
			plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
			axis.line = element_line(linewidth = 0.8),
			axis.ticks = element_line(linewidth = 0.8),
			axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
			axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
			axis.text = element_text(size = 12, face = "bold"),
			legend.title = element_text(face = "bold", size = 12),
			legend.text = element_text(size = 11),
			panel.background = element_rect(fill = "white", color = NA),
			plot.background = element_rect(fill = "white", color = NA),
			plot.margin = margin(15, 15, 15, 15)
		) +
		annotate(
			geom = "label",
			x = Inf, y = ifelse(cor_value > 0, Inf, -Inf),
			label = cor_label,
			hjust = 1.1, vjust = ifelse(cor_value > 0, 1.1, -0.1),
			size = 5,
			color = "white",
			fill = ifelse(cor_value > 0, "#e41a1c", "#377eb8"),
			fontface = "bold",
			label.padding = unit(0.5, "lines")
		)
	
	# 保存高质量图形
	ggsave(
		pdf_file,
		plot = p,
		width = 8,
		height = 6.5,
		dpi = 300
	)
}

# 保存相关性结果
write.csv(cor_results, "pseudotime_correlation_results.csv", row.names = FALSE)
#####提取state下游分析####
library(RcppML)
library(irGSEA)
mm$Monocle2_State

# 设置默认 assay 为 RNA
DefaultAssay(mm) <- "RNA"
mm <- SeuratObject::UpdateSeuratObject(object = mm)

# 基因集打分分析
cds_scored <- irGSEA.score(
	object = mm,
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
	category = "H",  # 使用 Hallmark 基因集
	subcategory = NULL,
	geneid = "symbol",
	method = c("AUCell", "UCell", "singscore", "ssgsea", "viper"),
	aucell.MaxRank = NULL,
	ucell.MaxRank = NULL,
	kcdf = 'Gaussian'
)

# 整合分析结果按 State 分组
result.dge <- irGSEA.integrate(
	object = cds_scored,
	group.by = "Monocle2_State",  # 按细胞状态分组
	method = c("AUCell", "UCell", "singscore", "ssgsea", "viper")
)

# 绘制热图
pdf(file = "heatmap_state.pdf", width = 12, height = 6)
irGSEA.heatmap(
	object = result.dge,
	method = "RRA",  # 使用 Robust Rank Aggregation 方法
	top = 20,        # 显示前 50 个显著基因集
	show.geneset = NULL,
)
dev.off()

#####阶段小提琴####
# ---- packages ----
plot_state_violin <- function(
		data, colname, state_col = "Monocle2_State",
		reference_group = NULL,  # 新增参考组参数
		output_prefix = "state_violin",
		base_width = 6, base_height = 6
) {
	# 构建绘图数据
	plot_data <- data.frame(
		Expression = data[[colname]],
		Subtype = factor(data[[state_col]])
	)
	colnames(plot_data) <- c("Expression", "Subtype")
	
	state_types <- levels(plot_data$Subtype)
	n_groups <- length(state_types)
	color_scheme <- pal_npg("nrc")(n_groups)
	
	# 确定参考组（新增逻辑）
	if(is.null(reference_group)) {
		reference_group <- tail(state_types, 1)  # 默认使用最后一组
	} else if(!reference_group %in% state_types) {
		stop("指定的参考组 '", reference_group, "' 不存在于数据中。可用组别: ", 
				 paste(state_types, collapse = ", "))
	}
	
	# 动态统计和图边距设置
	y_range <- range(plot_data$Expression, na.rm = TRUE)
	y_min <- y_range[1]; y_max <- y_range[2]; range_diff <- diff(y_range)
	
	# 生成比较组合：所有非参考组 vs 参考组
	other_groups <- setdiff(state_types, reference_group)
	comps <- lapply(other_groups, function(x) c(x, reference_group))
	n_comparisons <- length(comps)
	
	label_y_step <- range_diff * 0.08
	label_y_start <- y_max + label_y_step
	label_y <- label_y_start + label_y_step * (0:(n_comparisons - 1))
	plot_upper_limit <- y_max + (n_comparisons + 1) * label_y_step
	
	# 绘图
	p <- ggplot(plot_data, aes(x = Subtype, y = Expression, fill = Subtype)) +
		gghalves::geom_half_boxplot(side = "l", notch = TRUE,
																outlier.shape = NA, width = 0.35,
																alpha = 0.8, color = "black", lwd = 0.6, center = TRUE) +
		gghalves::geom_half_point(side = "l", shape = 21,
															color = "black", stroke = 0.25,
															size = 1.4, alpha = 0.6,
															position = position_nudge(x = 0.09),
															transformation = position_jitter(width = 0.05, height = 0),
															range_scale = 0.9) +
		gghalves::geom_half_violin(side = "r", color = "black", lwd = 0.6, alpha = 0.7,
															 trim = TRUE, scale = "width", width = 0.35,
															 position = position_nudge(x = 0.10)) +
		stat_summary(fun = median, geom = "errorbar",
								 aes(ymax = ..y.., ymin = ..y..),
								 width = 0.25, color = "black", linetype = "solid", size = 1.2,
								 position = position_nudge(x = 0.08)) +
		ggpubr::stat_compare_means(comparisons = comps,
															 method = "wilcox.test",
															 label = "p.signif",
															 tip.length = 0.01,
															 size = 4.5,
															 vjust = 0.4,
															 label.y = label_y) +
		scale_fill_manual(values = color_scheme) +
		labs(
			x = "State",
			y = colname,
			title = paste(colname, "across cell states"),
			subtitle = paste("Reference group:", reference_group)  # 显示参考组
		) +
		theme_minimal(base_size = 14) +
		theme(
			plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
			plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
			axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
			axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
			axis.text.x = element_text(size = 12, face = "bold", color = "black"),
			axis.text.y = element_text(size = 12, face = "bold", color = "black"),
			axis.line = element_line(color = "black", size = 0.8),
			axis.ticks = element_line(color = "black", size = 0.8),
			panel.grid.major = element_line(color = "grey90", size = 0.3),
			panel.grid.minor = element_blank(),
			legend.position = "none",
			plot.background = element_rect(fill = "white", color = NA),
			panel.background = element_rect(fill = "white", color = NA),
			plot.margin = margin(t = 1 + 0.2 * n_comparisons,
													 r = 1, b = 1, l = 1, unit = "lines")
		) +
		coord_cartesian(ylim = c(y_min, plot_upper_limit), clip = "off")
	
	# 保存图
	filename <- paste0(output_prefix, "_", colname, ".pdf")
	ggsave(filename, plot = p,
				 width = base_width + 0.5 * n_groups,
				 height = base_height + 0.3 * n_comparisons,
				 dpi = 300)
	
	return(p)
}

# ---- 使用示例 ----
# 指定参考组为 "3"
columns_to_plot = c("LMBDs","Glycolysis","Scoreing")
plots_list <- lapply(columns_to_plot, function(col) {
	plot_state_violin(
		data = mm@meta.data,
		colname = col,
		state_col = "Monocle2_State",
		reference_group = "2",  # 指定参考组
		output_prefix = "Custom_Reference_State"
	)
})
#####RCC2_State####
library(ggplot2)
library(ggpubr)
library(gghalves)
library(ggsci)

# ---- 基因表达量绘图函数 ----
plot_gene_expression_state <- function(
		seurat_obj, 
		gene_name,
		state_col = "Monocle2_State",
		reference_group = NULL,
		output_file = NULL,
		base_width = 6,
		base_height = 6
) {
	# 检查基因是否存在于对象中
	if (!gene_name %in% rownames(seurat_obj)) {
		stop("基因 '", gene_name, "' 不存在于Seurat对象中")
	}
	
	# 提取基因表达数据
	expr_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")[gene_name, ]
	
	# 构建绘图数据框
	plot_data <- data.frame(
		Expression = as.numeric(expr_data),
		Subtype = factor(seurat_obj@meta.data[[state_col]])
	)
	
	# 获取所有存在的state类型
	state_types <- levels(plot_data$Subtype)
	n_groups <- length(state_types)
	
	# 检查参考组有效性
	if (is.null(reference_group)) {
		reference_group <- tail(state_types, 1)  # 默认使用最后一组
	} else if (!reference_group %in% state_types) {
		stop("指定的参考组 '", reference_group, "' 不存在于数据中。可用组别: ", 
				 paste(state_types, collapse = ", "))
	}
	
	# 使用科研友好配色
	color_scheme <- pal_npg("nrc")(n_groups)
	
	# ---- 统计检验 ----
	# Kruskal-Wallis检验
	kruskal_p_value <- kruskal.test(Expression ~ Subtype, data = plot_data)$p.value
	sig_label <- ifelse(kruskal_p_value < 0.001, "***",
											ifelse(kruskal_p_value < 0.01, "**",
														 ifelse(kruskal_p_value < 0.05, "*", "")))
	
	# 生成两两比较组合（所有组与参考组比较）
	other_groups <- setdiff(state_types, reference_group)
	my_comparisons <- lapply(other_groups, function(x) c(x, reference_group))
	n_comparisons <- length(my_comparisons)
	
	# ---- 动态调整图形参数 ----
	y_range <- range(plot_data$Expression, na.rm = TRUE)
	y_min <- y_range[1]
	y_max <- y_range[2]
	range_diff <- diff(y_range)
	
	# 动态调整显著性标记位置
	label_y_step <- range_diff * 0.08
	label_y_start <- y_max + label_y_step
	label_y <- label_y_start + label_y_step * (0:(n_comparisons - 1))
	
	# 调整图形上边界
	plot_upper_limit <- y_max + (n_comparisons + 1) * label_y_step
	
	# ---- 绘图 ----
	set.seed(1)  # 让抖动可复现
	
	p <- ggplot(plot_data, aes(x = Subtype, y = Expression, fill = Subtype)) +
		# 左侧：带缺口的半箱线图
		gghalves::geom_half_boxplot(
			side = "l",
			notch = TRUE,
			outlier.shape = NA,
			width = 0.35,
			alpha = 0.8,
			color = "black",
			lwd = 0.6,
			center = TRUE
		) +
		# 左侧：抖动点
		gghalves::geom_half_point(
			side = "l",
			shape = 21,
			color = "black",
			stroke = 0.25,
			size = 1.4,
			alpha = 0.6,
			position = position_nudge(x = 0.09),
			transformation = position_jitter(width = 0.05, height = 0)
		) +
		# 右侧：半边小提琴
		gghalves::geom_half_violin(
			side = "r",
			color = "black",
			lwd = 0.6,
			alpha = 0.7,
			trim = TRUE,
			scale = "width",
			width = 0.35,
			position = position_nudge(x = 0.1)
		) +
		# 在小提琴上加中位线
		stat_summary(
			fun = median,
			geom = "errorbar",
			aes(ymax = ..y.., ymin = ..y..),
			width = 0.25,
			color = "black",
			linetype = "solid",
			size = 1.2,
			position = position_nudge(x = 0.08)
		) +
		# 两两比较显著性（动态位置）
		ggpubr::stat_compare_means(
			comparisons = my_comparisons,
			method = "wilcox.test",
			label = "p.signif",
			tip.length = 0.01,
			size = 4.5,
			vjust = 0.4,
			label.y = label_y
		) +
		# 配色
		scale_fill_manual(values = color_scheme) +
		# 坐标轴与标题
		labs(
			x = "State",
			y = paste0(gene_name, " Expression"),
			title = paste0(gene_name, " Expression Across Cell States"),
			subtitle = paste("Reference group:", reference_group)
		) +
		# 主题设置
		theme_minimal(base_size = 14) +
		theme(
			plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
			plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
			axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
			axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
			axis.text.x  = element_text(size = 12, face = "bold", color = "black"),
			axis.text.y  = element_text(size = 12, face = "bold", color = "black"),
			axis.line    = element_line(color = "black", size = 0.8),
			axis.ticks   = element_line(color = "black", size = 0.8),
			panel.grid.major = element_line(color = "grey90", size = 0.3),
			panel.grid.minor = element_blank(),
			legend.position  = "none",
			plot.background  = element_rect(fill = "white", color = NA),
			panel.background = element_rect(fill = "white", color = NA),
			# 动态上边距（按"行数/lines"为单位）
			plot.margin = margin(t = 1 + 0.2 * n_comparisons, r = 1, b = 1, l = 1, unit = "lines")
		) +
		coord_cartesian(ylim = c(y_min, plot_upper_limit), clip = "off")
	
	# ---- 保存 ----
	if (!is.null(output_file)) {
		ggsave(
			filename = output_file,
			plot = p,
			width = base_width + 0.5 * n_groups,     # 根据组数调整宽度
			height = base_height + 0.3 * n_comparisons, # 根据比较数量调整高度
			dpi = 300
		)
		message("已保存: ", output_file)
	}
	
	return(p)
}

# ---- 使用示例 ----
# 单个基因绘图
plot_gene_expression_state(
	seurat_obj = mm,
	gene_name = "RCC2",
	reference_group = "2",  # 指定参考组为State 3
	output_file = "RCC2_Expression_State_Ref3.pdf"
)

#####表达情况聚类####
library(Mfuzz)
library(dplyr)
library(tidyr)
library(tidyverse)

# 提取mm中包含在signature_metabolism中的列
ptm_cols <- names(signature_metabolism)  # 获取所有PTM通路名称
available_ptms <- ptm_cols[ptm_cols %in% colnames(mm@meta.data)]  # 确保列存在于mm中

# 创建包含所需列的数据框
ptm_data <- mm@meta.data %>%
	dplyr::select(Monocle2_State, all_of(available_ptms)) %>%
	mutate(Monocle2_State = as.factor(Monocle2_State))  # 确保状态是因子

# 计算每个状态的平均值
state_means <- ptm_data %>%
	group_by(Monocle2_State) %>%
	summarise(across(all_of(available_ptms), mean, na.rm = TRUE))

# 转换数据框格式
result_df <- state_means %>%
	pivot_longer(cols = -Monocle2_State, 
							 names_to = "Pathway", 
							 values_to = "Mean_Value") %>%
	pivot_wider(names_from = Monocle2_State, 
							values_from = Mean_Value)

# 设置行名为通路名称
result_df <- result_df %>% column_to_rownames("Pathway")

#构建对象
mfuzz_class <- new('ExpressionSet',exprs = as.matrix(result_df))

#预处理缺失值或者异常值
mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
mfuzz_class <- filter.std(mfuzz_class, min.std = 0)

#标准化数据
mfuzz_class <- standardise(mfuzz_class)

#Mfuzz 基于 fuzzy c-means 的算法进行聚类，详情 ?mfuzz
#需手动定义目标聚类群的个数，例如这里我们为了重现原作者的结果，设定为 10，即期望获得 10 组聚类群
#需要设定随机数种子，以避免再次运行时获得不同的结果
set.seed(123)
cluster_num <- 6
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))

library(RColorBrewer)
library(ggplot2)

# 创建自定义美化版的 Mfuzz 绘图函数
mfuzz_plot <- function(eset, cl, mfrow = c(1,1), colo, min.mem = 0, time.labels, new.window = TRUE) {
	# 创建更美观的颜色方案
	if (missing(colo)) {
		colo <- colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(100)
	}
	
	# 设置绘图布局
	opar <- par(no.readonly = TRUE)
	on.exit(par(opar))
	par(mfrow = mfrow, mar = c(5, 5, 4, 2) + 0.1)
	
	# 获取聚类中心
	centers <- cl$centers
	ncluster <- nrow(centers)
	
	# 设置时间点标签
	if (missing(time.labels)) {
		time.labels <- colnames(eset)
	}
	
	# 绘制每个聚类
	for (i in 1:ncluster) {
		# 获取当前聚类中心
		center <- centers[i, ]
		
		# 创建空白绘图区域
		plot(1, type = "n", 
				 xlim = c(1, length(center)), 
				 ylim = range(centers), 
				 xlab = "", ylab = "", 
				 axes = FALSE,
				 main = paste0("Cluster ", i, " (", sum(cl$cluster == i), " pathways)"),
				 cex.main = 1.8, font.main = 2, col.main = "#333333")
		
		# 添加坐标轴
		axis(1, at = 1:length(center), labels = time.labels, 
				 cex.axis = 1.5, col.axis = "#333333", col = "gray70", lwd = 1.5)
		axis(2, cex.axis = 1.5, col.axis = "#333333", col = "gray70", lwd = 1.5)
		
		# 添加网格线
		grid(NA, NULL, col = "gray90", lty = 1, lwd = 0.8)
		abline(h = 0, col = "gray70", lwd = 1.5)
		
		# 添加聚类中心线
		lines(1:length(center), center, 
					col = "#2c3e50", lwd = 3.5, type = "o", pch = 19, cex = 1.5)
		
		# 添加聚类成员通路（透明度较低）
		cluster_members <- exprs(eset)[cl$cluster == i, , drop = FALSE]
		for (j in 1:nrow(cluster_members)) {
			lines(1:length(center), cluster_members[j, ], 
						col = adjustcolor("#3498db", alpha.f = 0.5), 
						lwd = 1.2)
		}
		
		# 添加坐标轴标签
		mtext("Cell State", side = 1, line = 3.5, cex = 1.2, font = 2, col = "#333333")
		mtext("Standardized Activity", side = 2, line = 3.5, cex = 1.2, font = 2, col = "#333333")
		
		# 添加背景色
		rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
				 col = adjustcolor("#f8f9fa", alpha.f = 0.5), border = NA)
		
		# 添加边框
		box(lwd = 2, col = "gray70")
	}
}

# 使用自定义函数绘制美化版图形
pdf("plot_PTMD_cluster.pdf", width = 14, height = 10)
mfuzz_plot(
	mfuzz_class, 
	mfuzz_cluster, 
	mfrow = c(2, 3),  # 2行3列布局
	colo = colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(100),  # 蓝-黄-红渐变
	time.labels = colnames(mfuzz_class)  # 使用状态作为标签
)
dev.off()
#####拟时序聚类热图####
library(ClusterGVis)
diff_test_res <- differentialGeneTest(cds_seurat,
																			fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res1 <- diff_test_res %>%
	filter(num_cells_expressed > 10)
df=plot_pseudotime_heatmap2(cds_seurat[row.names(subset(diff_test_res,qval<1e-20)),],
														num_clusters=6)
gene=sample(df$wide.res$gene,20,replace=F)

# 创建 PDF 文件
pdf("plot_genes_heatmap.PDF", width=6, height=8, onefile=FALSE)
visCluster(object=df, plot.type="heatmap", markGenes=gene)
dev.off()

#####差异分析####
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# 提取State 1和State 5的细胞
state_subset <- subset(mm, subset = Monocle2_State %in% c(1, 2))

# 获取所有可用的通路名称
available_ptms <- names(signature_metabolism)[names(signature_metabolism) %in% colnames(state_subset@meta.data)]

# 创建差异分析数据框
diff_results <- data.frame(
	Pathway = character(),
	Avg_State1 = numeric(),
	Avg_State5 = numeric(),
	Log2FC = numeric(),
	p_value = numeric(),
	stringsAsFactors = FALSE
)

# 循环分析每个通路
for (ptm in available_ptms) {
	# 提取数据
	state1_data <- state_subset@meta.data %>%
		filter(Monocle2_State == 1) %>%
		pull(!!ptm)
	
	state4_data <- state_subset@meta.data %>%
		filter(Monocle2_State == 2) %>%
		pull(!!ptm)
	
	# 计算平均值
	avg1 <- mean(state1_data, na.rm = TRUE)
	avg4 <- mean(state4_data, na.rm = TRUE)
	
	# 计算log2 fold change (State4 vs State1)
	log2fc <- log2(avg4 / avg1)
	
	# 执行Wilcoxon秩和检验
	p_val <- wilcox.test(state1_data, state4_data)$p.value
	
	# 添加到结果
	diff_results <- rbind(diff_results, data.frame(
		Pathway = ptm,
		Avg_State1 = avg1,
		Avg_State4 = avg4,
		Log2FC = log2fc,
		p_value = p_val
	))
}

# 多重检验校正
diff_results$p_adj <- p.adjust(diff_results$p_value, method = "BH")

# 添加显著性标记
diff_results$Significance <- cut(
	diff_results$p_adj,
	breaks = c(0, 0.001, 0.01, 0.05, 1),
	labels = c("***", "**", "*", "NS"),
	include.lowest = TRUE
)

# 保存结果
write.csv(diff_results, "pathway_diff_state1_vs_state.csv", row.names = FALSE)

# 创建美化版火山图
volcano_plot <- ggplot(diff_results, aes(x = Log2FC, y = -log10(p_adj))) +
	# 添加渐变填充和大小的点
	geom_point(
		aes(fill = -log10(p_adj), size = -log10(p_adj)), 
		shape = 21, 
		color = "black", 
		stroke = 0.5,
		alpha = 0.9
	) +
	# 添加参考线
	geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 0.8) +
	geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", size = 0.8) +
	
	# 添加通路标签
	geom_text_repel(
		data = subset(diff_results, p_adj < 0.05 & abs(Log2FC) > 2),
		aes(label = Pathway),
		size = 4.5,
		fontface = "bold",
		box.padding = 0.5,
		point.padding = 0.2,
		max.overlaps = 30,
		segment.color = "grey40",
		segment.size = 0.4,
		min.segment.length = 0.1,
		direction = "both",
		nudge_y = 0.2
	) +
	
	# 添加顶部箭头标注
	annotate(
		"text",
		x = -Inf, y = Inf,                      # 面板左上角
		label = "Enriched in State 1",
		size = 5.5, color = "#4575b4", fontface = "bold",
		hjust = -0.05, vjust = -0.6             # 往左/上各“探出”一点
	) +
	annotate(
		"text",
		x =  Inf, y = Inf,                      # 面板右上角
		label = "Enriched in State 2",
		size = 5.5, color = "#d73027", fontface = "bold",
		hjust =  1.05, vjust = -0.6             # 往右/上各“探出”一点
	) +
	# 颜色和大小渐变设置
	scale_fill_gradientn(
		name = bquote(atop(-Log[10]~"adjusted",
											 "p-value")),
		colours = c("#4575b4", "#74add1", "#e0f3f8", "#ffffbf", "#fee090", "#fdae61", "#d73027"),
		limits = c(0, max(-log10(diff_results$p_adj), na.rm = TRUE)),
		breaks = scales::pretty_breaks(n = 5),
		guide = guide_colorbar(
			title.position = "top",
			title.hjust = 0.5,
			barwidth = unit(3, "cm"),
			barheight = unit(0.3, "cm")
		)
	) +
	scale_size_continuous(
		name = bquote(atop(-Log[10]~"adjusted",
											 "p-value")),
		range = c(2, 8),
		breaks = scales::pretty_breaks(n = 5),
		guide = guide_legend(
			title.position = "top",
			title.hjust = 0.5,
			override.aes = list(color = "black")
		)
	) +
	
	# 坐标轴和标签
	labs(
		x = expression(Log[2]~Fold~Change~(State~3~vs~State~1)),
		y = expression(-Log[10]~adjusted~p~value),
		title = "Pathway Activity Differences",
		subtitle = "State 1 vs State 2"
	) +
	
	# 主题设置
	theme_minimal(base_size = 14) +
	theme(
		legend.position = "bottom",
		legend.box = "horizontal",
		legend.direction = "horizontal",
		legend.title = element_text(face = "bold", size = 12, vjust = 0.5),
		legend.text = element_text(size = 10),
		legend.spacing.x = unit(0.5, "cm"),
		panel.border = element_rect(color = "black", fill = NA, size = 1.2),
		axis.text = element_text(color = "black", face = "bold", size = 12),
		axis.title = element_text(face = "bold", size = 14),
		plot.title = element_text(face = "bold", hjust = 0.5, size = 20, margin = margin(b = 10)),
		plot.subtitle = element_text(hjust = 0.5, size = 14, color = "grey40", margin = margin(b = 15)),
		panel.grid.major = element_line(color = "grey92", size = 0.3),
		plot.background = element_rect(fill = "white", color = NA),
		plot.margin = margin(20, 20, 20, 20)
	) +
	
	# 扩展绘图区域以容纳顶部标签
	coord_cartesian(
		ylim = c(0, max(-log10(diff_results$p_adj), na.rm = TRUE) * 1.1),
		clip = "off"
	)

# 保存图形
ggsave("pathway_volcano_state1_vs_state2.pdf", volcano_plot, width = 10, height = 9)

#####小提琴图####
library(ggplot2)
library(ggpubr)
library(gghalves)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

# 路径选择函数：从差异分析结果中选出上调/下调最显著的通路
select_top10_paths <- function(df, top_n = 5, sig_only = TRUE) {
	d <- df %>% filter(is.finite(Log2FC))
	
	pick_side <- function(dat, side = c("state2", "state1")) {
		side <- match.arg(side)
		if (side == "state2") {
			cand <- dat %>% arrange(desc(Log2FC))
			sig  <- cand %>% filter(p_adj < 0.05)
			top  <- if (sig_only) sig else cand
			top  <- top %>% slice_head(n = top_n) %>% pull(Pathway)
			if (length(top) < top_n) {
				top <- unique(c(top, cand$Pathway))[seq_len(min(top_n, nrow(cand)))]
			}
			top
		} else {
			cand <- dat %>% arrange(Log2FC)
			sig  <- cand %>% filter(p_adj < 0.05)
			top  <- if (sig_only) sig else cand
			top  <- top %>% slice_head(n = top_n) %>% pull(Pathway)
			if (length(top) < top_n) {
				top <- unique(c(top, cand$Pathway))[seq_len(min(top_n, nrow(cand)))]
			}
			top
		}
	}
	
	unique(c(pick_side(d, "state2"), pick_side(d, "state1")))
}

# 通路绘图函数：可视化选定通路在两个亚组中的活性
plot_subtype_pathways <- function(plot_df, pathway_list, title,
																	subtype_col = "Monocle2_State",
																	code_prefix = "M",
																	legend_cols = 2,        # 说明面板列数
																	top_expand = 0.20,
																	legend_trunc = 50,      # 说明面板换行宽度（字符）
																	legend_replace_us = TRUE,
																	legend_text_size = 5,   # 说明面板字体
																	star_size = 6,          # 主图星号大小
																	legend_left_pad = 0.05, # ← 左侧"贴边"程度
																	legend_col_spacing = 1.6 # ← 列间水平间距
) {
	# 1) 短码映射
	code_map    <- setNames(paste0(code_prefix, seq_along(pathway_list)), pathway_list)
	levels_code <- unname(code_map[pathway_list])
	
	# 2) 整形为长表
	plot_data <- plot_df %>%
		select(all_of(pathway_list), all_of(subtype_col)) %>%
		tidyr::pivot_longer(cols = -all_of(subtype_col), 
												names_to = "Pathway", 
												values_to = "Score") %>%
		mutate(
			Pathway = factor(Pathway, levels = pathway_list),
			Code    = factor(code_map[Pathway], levels = levels_code),
			Subtype = factor(.data[[subtype_col]])
		)
	
	# 3) 每通路 Wilcoxon + 星号
	p_values <- sapply(pathway_list, function(p) {
		d <- dplyr::filter(plot_data, Pathway == p)
		out <- try(wilcox.test(Score ~ Subtype, data = d)$p.value, silent = TRUE)
		if (inherits(out, "try-error")) NA_real_ else out
	})
	
	sig_labels <- ifelse(p_values < 0.001, "***",
											 ifelse(p_values < 0.01,  "**",
											 			 ifelse(p_values < 0.05,  "*",  "")))
	
	sig_df <- data.frame(
		Code  = factor(unname(code_map[pathway_list]), levels = levels_code),
		label = sig_labels, 
		stringsAsFactors = FALSE
	)
	
	# 4) 主图
	y_min <- min(plot_data$Score, na.rm = TRUE)
	y_max <- max(plot_data$Score, na.rm = TRUE)
	y_lim <- c(y_min, y_max * (1 + top_expand))
	
	# 定义颜色方案
	fill_colors <- c("2" = "#E64B35", "1" = "#4DBBD5")
	color_colors <- c("2" = "#B71C1C", "1" = "#006064")
	
	p_main <- ggplot(plot_data, aes(x = Code, y = Score, fill = Subtype)) +
		geom_vline(xintercept = seq(1.5, length(levels_code) - 0.5), 
							 color = "gray70", linetype = "dashed", linewidth = 0.5) +
		gghalves::geom_half_violin(
			aes(fill = Subtype),
			position = position_dodge(0.8),
			side = "l", scale = "width",
			width = 0.7, alpha = 0.5, color = NA
		) +
		geom_boxplot(
			aes(color = Subtype),
			width = 0.3, position = position_dodge(0.8),
			outlier.shape = NA, notch = TRUE, notchwidth = 0.5,
			alpha = 0.8, linewidth = 0.6
		) +
		geom_point(
			aes(color = Subtype),
			position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
			size = 0.4, alpha = 0.4
		) +
		geom_text(
			data = sig_df,
			aes(x = Code, y = y_lim[2] - 0.02 * diff(y_lim), label = label),
			inherit.aes = FALSE, size = star_size, vjust = 0.5, fontface = "bold"
		) +
		scale_fill_manual(values = fill_colors, name = "Group") +
		scale_color_manual(values = color_colors, name = "Group", guide = "none") +
		labs(title = title, x = "", y = "Pathway Activity Score") +
		theme_minimal(base_size = 13) +
		theme(
			axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", 
																 size = 11, color = "#333333"),
			axis.text.y = element_text(face = "bold", size = 11, color = "#333333"),
			axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
			plot.title   = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 8)),
			legend.position = "top",
			legend.title    = element_text(face = "bold", size = 12),
			legend.text     = element_text(size = 11),
			panel.grid.major.y = element_line(color = "grey92", linewidth = 0.6),
			panel.grid.minor   = element_blank(),
			panel.background  = element_rect(fill = "white", color = NA),
			plot.background   = element_rect(fill = "white", color = NA),
			panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.8),
			plot.margin       = margin(10, 15, 5, 15)
		) +
		coord_cartesian(ylim = y_lim, clip = "off")
	
	# 5) 说明面板
	legend_map <- tibble(
		Pathway = pathway_list,
		Code    = unname(code_map[pathway_list])
	) %>%
		mutate(
			PrettyName = if (legend_replace_us) str_replace_all(Pathway, "_", " ") else Pathway,
			LabelFull  = paste0(Code, ": ", PrettyName),
			Label      = str_wrap(LabelFull, width = legend_trunc),
			idx        = row_number(),
			row        = ceiling(idx / legend_cols),
			col        = (idx - 1) %% legend_cols + 1,
			x          = (col - 1) * legend_col_spacing + legend_left_pad
		)
	
	x_right <- (legend_cols - 1) * legend_col_spacing + legend_left_pad + 1
	
	p_legend <- ggplot(legend_map, aes(x = x, y = row, label = Label)) +
		geom_text(hjust = 0, vjust = 1, size = legend_text_size, lineheight = 1.05) +
		scale_x_continuous(limits = c(0, x_right), expand = c(0, 0)) +
		scale_y_reverse(expand = c(0.18, 0.18)) +
		coord_cartesian(clip = "off") +
		theme_void(base_size = 12) +
		theme(
			panel.border = element_rect(color = "grey50", linetype = "dashed", fill = NA),
			plot.margin  = margin(0, 24, 10, 10)
		)
	
	# 组合图形
	ggpubr::ggarrange(p_main, p_legend, ncol = 1, heights = c(1, 0.30))
}

# ===== 使用示例 =====
# 假设我们有：
# 1. diff_results: 差异分析结果数据框，包含Pathway, Log2FC, p_adj列
# 2. state_subset: Seurat对象子集，包含细胞亚型信息和通路活性分数

# 步骤1: 选择top通路
sel_paths <- select_top10_paths(diff_results, top_n = 5, sig_only = TRUE)

# 步骤2: 准备绘图数据
# 从Seurat对象中提取元数据（包含通路分数和亚型信息）
plot_df <- state_subset@meta.data

# 步骤3: 绘制通路图
sig_plot <- plot_subtype_pathways(
	plot_df = plot_df,  # 包含通路分数和亚型的数据框
	pathway_list = sel_paths,
	title = "Top pathways (State2 & State1",
	legend_cols = 2,
	legend_trunc = 56,
	legend_text_size = 4,
	star_size = 7,
	legend_left_pad = 0.02,
	legend_col_spacing = 1.8
)

# 步骤4: 保存图形
ggsave("signature_pathways_top10.pdf", sig_plot, width = 12, height = 7)


tar <- mm$Scoreing
result_sub <- mm$LMBDs
cor_res <- cor(x = tar, y = result_sub, method = 'pearson')

# 计算相关性和p值
cor_test_result <- cor.test(tar, result_sub, method = 'pearson')
cor_rho <- round(cor_test_result$estimate, 2)
cor_p <- format.pval(cor_test_result$p.value, digits = 2)
n <- length(tar)
pdf(file="cor.plot.pdf", width=6, height=6)

# 设置图形参数
par(bty = "o", 
		mgp = c(2, 0.5, 0), 
		mar = c(4.1, 4.1, 2.1, 4.1), 
		tcl = -.25, 
		font.main = 3)

# 创建基础画布
plot(NULL, NULL, 
		 ylim = range(tar), 
		 xlim = range(result_sub), 
		 xlab = "Lactylation", 
		 ylab = "LMBDs", 
		 col = "white", 
		 main = "")

# 设置背景色
rect(par("usr")[1], par("usr")[3], 
		 par("usr")[2], par("usr")[4], 
		 col = "#EAE9E9", border = FALSE)

# 添加网格线
grid(col = "white", lty = 1, lwd = 1.5)

# 添加散点
points(x = result_sub, y = tar, 
			 pch = 19, 
			 col = scales::alpha("#E51718", 0.8),  # 使用单一红色系
			 cex = 1)  # 固定点大小

# 添加回归线
abline(lm(tar ~ result_sub), 
			 lwd = 4, 
			 col = "black")

# 添加边缘rug图
rug(result_sub, side = 3, col = "black", lwd = 1)
rug(tar, side = 4, col = "black", lwd = 1)

# 添加统计注释
text(x = min(result_sub), 
		 y = max(tar) * 0.95,  # 将注释放在右上角
		 adj = 0,
		 labels = bquote("AEC: N = " ~ .(n) ~ 
		 									"; " ~ rho ~ " = " ~ .(cor_rho) ~ 
		 									"; " ~ italic(P) ~ " = " ~ .(cor_p)),
		 col = "black",
		 cex = 0.8)
dev.off()
