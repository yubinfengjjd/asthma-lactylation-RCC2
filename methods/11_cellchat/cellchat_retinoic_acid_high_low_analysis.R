library(Seurat)
library(GSVA)
library(GSEABase)
library(dplyr)
library(tidyverse)
rm(list = ls())
setwd("D:\\哮喘乳酸化单细胞开篇\\14.细胞通讯")

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

table(oc$celltype)
library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
create_cellchat <- function(subset_obj, group_name) {
	# Input expression + group labels
	data_input <- subset_obj@assays$RNA$data
	meta_data <- subset_obj@meta.data[, c("celltype"),drop=F]
	colnames(meta_data) <- c("labels")
	
	# Create CellChat object
	cellchat <- createCellChat(object = data_input, meta = meta_data, group.by = "labels")
	cellchat <- addMeta(cellchat, meta = meta_data)
	cellchat <- setIdent(cellchat, ident.use = "labels")
	
	# Configure CellChat database
	CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
	cellchat@DB <- CellChatDB.use
	
	# Standard CellChat workflow
	cellchat <- subsetData(cellchat)
	cellchat <- identifyOverExpressedGenes(cellchat)
	cellchat <- identifyOverExpressedInteractions(cellchat)
	cellchat <- projectData(cellchat, PPI.human)
	
	# Compute communication probabilities
	cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
	#cellchat <- filterCommunication(cellchat, min.cells = 10)
	cellchat <- computeCommunProbPathway(cellchat)
	cellchat <- aggregateNet(cellchat)
	
	return(cellchat)
}

cellchat <- create_cellchat(oc,"oc")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
###### Visualization ######
library(ggplotify)
library(ggplot2)
groupSize = as.numeric(table(cellchat@idents))

# Wrap base plot into a ggplot object
p_num <- as.ggplot(function() {
	netVisual_circle(
		cellchat@net$count,
		vertex.weight = groupSize,
		weight.scale = TRUE,
		label.edge = FALSE
	)
})

p_wgt <- as.ggplot(function() {
	netVisual_circle(
		cellchat@net$weight,
		vertex.weight = groupSize,
		weight.scale = TRUE,
		label.edge = FALSE
	)
})

# 定义一个统一的标题样式
title_theme <- theme(
	plot.title = element_text(hjust = 0.5, vjust = -0.5, face = "bold", size = 16),
	plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5)  # 上边距缩小
)

# 保存 PDF
ggsave(
	"netVisual_circle_Number.pdf",
	p_num + ggtitle("Number of interactions") + title_theme,
	width = 8, height = 8
)

ggsave(
	"netVisual_circle_weight.pdf",
	p_wgt + ggtitle("Interaction weights/strength") + title_theme,
	width = 8, height = 8
)

library(CellChat)
library(grid)
library(ComplexHeatmap)
library(circlize)   # colorRamp2
# -----
pathways.show <- c("TGFb","CXCL","CCL","COMPLEMENT")
pw.avail <- intersect(pathways.show, cellchat@netP$pathways)

for (pw in pw.avail) {
	df_comm <- tryCatch(subsetCommunication(cellchat, signaling = pw),
											error = function(e) NULL)
	if (is.null(df_comm) || nrow(df_comm) == 0) next
	
	ht <- netVisual_heatmap(cellchat, signaling = pw,
													measure = "weight", color.heatmap = "Reds")
	
	pdf(paste0("CellChat_", pw, "_heatmap.pdf"), width = 8.5, height = 6.5)
	
	grid.newpage()
	# 三行：标题 / 热图 / 色带
	lay <- grid.layout(nrow = 3, ncol = 1,
										 heights = unit.c(unit(6, "mm"), unit(1, "null"), unit(7, "mm")))
	pushViewport(viewport(layout = lay))
	
	## 行1：标题
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
	grid.text(paste0(pw, " signaling network"),
						gp = gpar(fontsize = 14, fontface = "bold"))
	popViewport()
	
	## 行2：热图（不要新开页）
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
	ComplexHeatmap::draw(ht, heatmap_legend_side = "right", newpage = FALSE)
	popViewport()
	
	## 行3：自定义水平色带 + 左右标签
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
	
	# 画一个从红到蓝的水平渐变条（左=高 / 右=低）
	n <- 400
	cols <- colorRampPalette(c("red", "blue"))(n)
	mat_cols <- matrix(cols, nrow = 1, ncol = n, byrow = TRUE)  # 1行n列
	
	grid.raster(mat_cols,
							width = unit(0.75, "npc"), height = unit(3, "mm"),
							x = 0.5, y = 0.55, interpolate = TRUE)
	
	# 左右文字
	grid.text("High",
						x = unit(0.5, "npc") - unit(0.375, "npc"), y = unit(0.55, "npc"),
						just = c("right", "center"),
						gp = gpar(col = "red", fontsize = 10, fontface = "bold"))
	grid.text("Low",
						x = unit(0.5, "npc") + unit(0.375, "npc"), y = unit(0.55, "npc"),
						just = c("left", "center"),
						gp = gpar(col = "blue", fontsize = 10, fontface = "bold"))
	
	popViewport()
	dev.off()
}

library(CellChat)
library(ggplot2)
library(ggalluvial)

# 取簇
clusters       <- levels(cellchat@idents)
high_clusters  <- grep("^high_", clusters, value = TRUE)
low_clusters   <- grep("^low_",  clusters, value = TRUE)
stopifnot(length(high_clusters) > 0, length(low_clusters) > 0)

# 右侧接收者 = low_*
vertex.receiver <- which(clusters %in% high_clusters)

# 想画的通路（可自定义）
pathways.want <- c("TGFb","MIF","IL16")

# 仅保留对象里存在且“有边”的通路
pathways.ok <- intersect(pathways.want, cellchat@netP$pathways)
pathways.use <- Filter(function(pw) {
	df <- subsetCommunication(cellchat, signaling = pw)
	!is.null(df) && nrow(df) > 0
}, pathways.ok)
stopifnot(length(pathways.use) > 0)

for (pw in pathways.use) {
	pdf(paste0("hierarchy_high2low1_", pw, ".pdf"), width = 15, height = 6.5)
	p <- netVisual_aggregate(
		object = cellchat,
		signaling = pw,
		layout = "hierarchy",
		vertex.receiver = vertex.receiver
	)
	dev.off()
}


library(CellChat)

# 提取原始细胞类型标签
original_idents <- as.character(cellchat@idents)

# 创建新的细胞类型标签：
# 1. 保留两种AEC类型不变
# 2. 其他细胞类型移除"high_"和"low_"前缀
new_labels <- ifelse(original_idents %in% c("high_AEC", "low_AEC"),
										 original_idents,
										 gsub("^high_|^low_", "", original_idents))

# 将新标签转换为因子
unique_types <- unique(new_labels)
new_labels <- factor(new_labels, levels = unique_types)

# 创建新的CellChat对象
cellchat_merged <- createCellChat(
	object = cellchat@data.signaling,
	meta = data.frame(labels = new_labels, row.names = names(new_labels)),
	group.by = "labels"
)

# 设置相同的数据库
cellchat_merged@DB <- cellchat@DB

# 重新计算通信概率
cellchat_merged <- subsetData(cellchat_merged)
cellchat_merged <- identifyOverExpressedGenes(cellchat_merged)
cellchat_merged <- identifyOverExpressedInteractions(cellchat_merged)
cellchat_merged <- projectData(cellchat_merged, PPI.human)  # 如果是人类数据
cellchat_merged <- computeCommunProb(cellchat_merged)

# 计算聚合网络
# 关键步骤：计算信号通路水平的通信概率
cellchat_merged <- computeCommunProbPathway(cellchat_merged)
cellchat_merged <- aggregateNet(cellchat_merged)

# 设置信号源为两种AEC
sources.use <- c("high_AEC", "low_AEC")

# 目标细胞类型（除AEC外的所有合并类型）
targets.use <- setdiff(levels(cellchat_merged@idents), sources.use)

# 可选：分别绘制high_AEC和low_AEC的信号
pdf("high_AEC_to_Merged_bubble.pdf", width = 14, height = 6)
bubble_plot = netVisual_bubble(
	cellchat_merged,
	sources.use = "high_AEC",
	targets.use = targets.use,
	remove.isolate = FALSE,
	angle.x = 90
) +
	coord_flip() +
	ggtitle("Signaling from high_AEC to Merged Cell Types")
# 增强气泡大小和可读性
enhanced_plot <- bubble_plot +
	# 增大气泡尺寸
	scale_size_continuous(
		range = c(3, 10),  # 最小和最大气泡尺寸
		name = "Interaction\nStrength",  # 图例标题
		breaks = seq(0, max(bubble_plot$data$prob, na.rm = TRUE), length.out = 5)  # 图例断点
	) +
	# 增强颜色对比度
	scale_color_gradientn(
		colors = c("blue", "yellow", "red"),  # 使用更鲜明的颜色梯度
		name = "Communication\nProbability",
		limits = c(0, 1)
	)
print(enhanced_plot)
dev.off()

pdf("low_AEC_to_Merged_bubble.pdf", width = 14, height = 6)  # 增加画布尺寸以适应更大气泡

# 创建基础气泡图
bubble_plot <- netVisual_bubble(
	cellchat_merged,
	sources.use = "low_AEC",
	targets.use = targets.use,
	remove.isolate = FALSE,
	angle.x = 90
) +
	coord_flip() +
	ggtitle("Signaling from low_AEC to Merged Cell Types")

# 增强气泡大小和可读性
enhanced_plot <- bubble_plot +
	# 增大气泡尺寸
	scale_size_continuous(
		range = c(3, 10),  # 最小和最大气泡尺寸
		name = "Interaction\nStrength",  # 图例标题
		breaks = seq(0, max(bubble_plot$data$prob, na.rm = TRUE), length.out = 5)  # 图例断点
	) +
	# 增强颜色对比度
	scale_color_gradientn(
		colors = c("blue", "yellow", "red"),  # 使用更鲜明的颜色梯度
		name = "Communication\nProbability",
		limits = c(0, 1)
	)
print(enhanced_plot)
dev.off()

## =========================
## A) 在整个对象内：按“每个基础细胞类型”的中位数定义 RA_group
## =========================
library(Seurat); library(dplyr)
md <- oc@meta.data

# 去掉旧的 high_/low_ 前缀，还原基础类型名
base_type <- gsub("^high_|^low_", "", md$celltype)
oc$celltype0 <- base_type

# 计算每个基础类型的 RA 代谢中位数
ra_medians <- tapply(md$Retinoic_Acid_Metabolism, oc$celltype0, median, na.rm = TRUE)
oc$RA_group <- ifelse(md$Retinoic_Acid_Metabolism > ra_medians[oc$celltype0], "RA_high", "RA_low")

## =========================
## B) 仅保留 RA_high 的细胞，且在其中把 AEC 按 RCC2 中位数再分 high/low
## =========================
cells_ra_high <- rownames(oc@meta.data)[oc$RA_group == "RA_high"]
oc_ra <- subset(oc, cells = cells_ra_high)

md_ra <- oc_ra@meta.data
base_type_ra <- oc_ra$celltype0
is_AEC_ra <- base_type_ra == "AEC"     # 若你的命名是 AT1/AT2，请改成 %in% c("AT1","AT2")

# RCC2 表达（log-normalized）
rcc2_gene <- grep("^RCC2$", rownames(oc_ra[["RNA"]]), ignore.case = TRUE, value = TRUE)
stopifnot(length(rcc2_gene) == 1)
rcc2_vec <- as.numeric(GetAssayData(oc_ra, assay = "RNA", slot = "data")[rcc2_gene, ])

rcc2_med_AEC_ra <- median(rcc2_vec[is_AEC_ra], na.rm = TRUE)
oc_ra$RCC2_group <- NA_character_
oc_ra$RCC2_group[is_AEC_ra] <- ifelse(rcc2_vec[is_AEC_ra] > rcc2_med_AEC_ra, "RCC2_high", "RCC2_low")

# 生成绘图/CellChat 用标签：
# - AEC → “RAhigh_AEC_RCC2_high/low”
# - 其他细胞类型（都已是 RA_high） → 用其基础类型名
labels_ra <- ifelse(is_AEC_ra,
                    paste0(oc_ra$RCC2_group,"_AEC"),
                    base_type_ra)
labels_ra <- factor(labels_ra, levels = unique(labels_ra))
table(labels_ra)   # 自查

## =========================
## C) 用 “RA_high 子集 + 新标签” 创建 CellChat 并计算通路级网络
## =========================
library(CellChat)
data_input <- oc_ra@assays$RNA$data
meta_df <- data.frame(labels = labels_ra, row.names = colnames(data_input))

cellchat_ra <- createCellChat(object = data_input, meta = meta_df, group.by = "labels")
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cellchat_ra@DB <- CellChatDB.use

cellchat_ra <- subsetData(cellchat_ra)
cellchat_ra <- identifyOverExpressedGenes(cellchat_ra)
cellchat_ra <- identifyOverExpressedInteractions(cellchat_ra)
cellchat_ra <- projectData(cellchat_ra, PPI.human)

cellchat_ra <- computeCommunProb(cellchat_ra, raw.use = TRUE)
# 可选：若有极小簇，可设最小细胞数过滤
# cellchat_ra <- filterCommunication(cellchat_ra, min.cells = 10)

cellchat_ra <- computeCommunProbPathway(cellchat_ra)  # 关键：得到 netP$prob（源×靶×通路）
cellchat_ra <- aggregateNet(cellchat_ra)

# 设置信号源为两种AEC
sources.use <- c("RCC2_high_AEC", "RCC2_low_AEC")

# 目标细胞类型（除AEC外的所有合并类型）
targets.use <- setdiff(levels(cellchat_ra@idents), sources.use)

# 可选：分别绘制high_AEC和low_AEC的信号
pdf("RCC2_high_AEC_to_Merged_bubble.pdf", width = 14, height = 6)
bubble_plot = netVisual_bubble(
  cellchat_ra,
  sources.use = "RCC2_high_AEC",
  targets.use = targets.use,
  remove.isolate = FALSE,
  angle.x = 90
) +
  coord_flip() +
  ggtitle("Signaling from RCC2_high_AEC to Merged Cell Types")
# 增强气泡大小和可读性
enhanced_plot <- bubble_plot +
  # 增大气泡尺寸
  scale_size_continuous(
    range = c(3, 10),  # 最小和最大气泡尺寸
    name = "Interaction\nStrength",  # 图例标题
    breaks = seq(0, max(bubble_plot$data$prob, na.rm = TRUE), length.out = 5)  # 图例断点
  ) +
  # 增强颜色对比度
  scale_color_gradientn(
    colors = c("blue", "yellow", "red"),  # 使用更鲜明的颜色梯度
    name = "Communication\nProbability",
    limits = c(0, 1)
  )
print(enhanced_plot)
dev.off()

pdf("RCC2_low_AEC_to_Merged_bubble.pdf", width = 14, height = 6)  # 增加画布尺寸以适应更大气泡

# 创建基础气泡图
bubble_plot <- netVisual_bubble(
  cellchat_ra,
  sources.use = "RCC2_low_AEC",
  targets.use = targets.use,
  remove.isolate = FALSE,
  angle.x = 90
) +
  coord_flip() +
  ggtitle("Signaling from RCC2_low_AEC to Merged Cell Types")

# 增强气泡大小和可读性
enhanced_plot <- bubble_plot +
  # 增大气泡尺寸
  scale_size_continuous(
    range = c(3, 10),  # 最小和最大气泡尺寸
    name = "Interaction\nStrength",  # 图例标题
    breaks = seq(0, max(bubble_plot$data$prob, na.rm = TRUE), length.out = 5)  # 图例断点
  ) +
  # 增强颜色对比度
  scale_color_gradientn(
    colors = c("blue", "yellow", "red"),  # 使用更鲜明的颜色梯度
    name = "Communication\nProbability",
    limits = c(0, 1)
  )
print(enhanced_plot)
dev.off()

##计算网络中心评分
cellchat_ra <- netAnalysis_computeCentrality(cellchat_ra, slot.name = "netP")
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pdf("netAnalysis_signalingRole_scatter.pdf",width = 6,height = 6)
netAnalysis_signalingRole_scatter(cellchat_ra)
dev.off()

pdf("netAnalysis_signalingRole_heatmap.pdf",width = 10,height = 6)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_ra, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_ra, pattern = "incoming")
ht1 + ht2
dev.off()

# 循环生成每个通路的热图
pathways.show = c("OSM","VISFATIN")
for (pathway in pathways.show) {
  # 动态生成文件名
  pdf_filename <- paste0("netAnalysis_signalingRole_heatmap_", pathway, ".pdf")
  
  # 打开 PDF 设备
  pdf(pdf_filename, width = 8, height = 6)
  
  # 调用 netAnalysis_signalingRole_network 函数
  netAnalysis_signalingRole_network(cellchat, pathway, width = 12, height = 4, font.size = 6)
  
  # 关闭 PDF 设备
  dev.off()
}
