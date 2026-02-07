library(Seurat)
library(GSVA)
library(GSEABase)
library(dplyr)
library(tidyverse)
rm(list = ls())
setwd("D:\\哮喘乳酸化单细胞开篇\\14.细胞通讯")

meta_data <- oc@meta.data

meta_data <- meta_data %>%
	group_by(celltype) %>%
	mutate(
		median = median(Retinoic_Acid_Metabolism, na.rm = TRUE),
		group = if_else(Retinoic_Acid_Metabolism > median,
										paste0("high_", celltype),
										paste0("low_", celltype))
	) %>%
	ungroup()

oc$celltype <- meta_data$group

table(oc$celltype)
library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
create_cellchat <- function(subset_obj, group_name) {
	data_input <- subset_obj@assays$RNA$data
	meta_data <- subset_obj@meta.data[, c("celltype"),drop=F]
	colnames(meta_data) <- c("labels")

	cellchat <- createCellChat(object = data_input, meta = meta_data, group.by = "labels")
	cellchat <- addMeta(cellchat, meta = meta_data)
	cellchat <- setIdent(cellchat, ident.use = "labels")

	CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
	cellchat@DB <- CellChatDB.use

	cellchat <- subsetData(cellchat)
	cellchat <- identifyOverExpressedGenes(cellchat)
	cellchat <- identifyOverExpressedInteractions(cellchat)
	cellchat <- projectData(cellchat, PPI.human)

	cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
	#cellchat <- filterCommunication(cellchat, min.cells = 10)
	cellchat <- computeCommunProbPathway(cellchat)
	cellchat <- aggregateNet(cellchat)

	return(cellchat)
}

cellchat <- create_cellchat(oc,"oc")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
library(ggplotify)
library(ggplot2)

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

title_theme <- theme(
	plot.title = element_text(hjust = 0.5, vjust = -0.5, face = "bold", size = 16),
	plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5)
)

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
pathways.show <- c("TGFb","MIF","IL16")
pw.avail <- intersect(pathways.show, cellchat@netP$pathways)

for (pw in pw.avail) {
	df_comm <- tryCatch(subsetCommunication(cellchat, signaling = pw),
											error = function(e) NULL)
	if (is.null(df_comm) || nrow(df_comm) == 0) next

	ht <- netVisual_heatmap(cellchat, signaling = pw,
													measure = "weight", color.heatmap = "Reds")

	pdf(paste0("CellChat_", pw, "_heatmap.pdf"), width = 8.5, height = 6.5)

	grid.newpage()
	lay <- grid.layout(nrow = 3, ncol = 1,
										 heights = unit.c(unit(6, "mm"), unit(1, "null"), unit(7, "mm")))
	pushViewport(viewport(layout = lay))

	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
	grid.text(paste0(pw, " signaling network"),
						gp = gpar(fontsize = 14, fontface = "bold"))
	popViewport()

	pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
	ComplexHeatmap::draw(ht, heatmap_legend_side = "right", newpage = FALSE)
	popViewport()

	pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))

	n <- 400
	cols <- colorRampPalette(c("red", "blue"))(n)
	mat_cols <- matrix(cols, nrow = 1, ncol = n, byrow = TRUE)

	grid.raster(mat_cols,
							width = unit(0.75, "npc"), height = unit(3, "mm"),
							x = 0.5, y = 0.55, interpolate = TRUE)

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

clusters       <- levels(cellchat@idents)
high_clusters  <- grep("^high_", clusters, value = TRUE)
low_clusters   <- grep("^low_",  clusters, value = TRUE)
stopifnot(length(high_clusters) > 0, length(low_clusters) > 0)

vertex.receiver <- which(clusters %in% high_clusters)

pathways.want <- c("TGFb","MIF","IL16")

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

original_idents <- as.character(cellchat@idents)

new_labels <- ifelse(original_idents %in% c("high_AEC", "low_AEC"),
										 original_idents,
										 gsub("^high_|^low_", "", original_idents))

unique_types <- unique(new_labels)
new_labels <- factor(new_labels, levels = unique_types)

cellchat_merged <- createCellChat(
	object = cellchat@data.signaling,
	meta = data.frame(labels = new_labels, row.names = names(new_labels)),
	group.by = "labels"
)

cellchat_merged@DB <- cellchat@DB

cellchat_merged <- subsetData(cellchat_merged)
cellchat_merged <- identifyOverExpressedGenes(cellchat_merged)
cellchat_merged <- identifyOverExpressedInteractions(cellchat_merged)
cellchat_merged <- projectData(cellchat_merged, PPI.human)
cellchat_merged <- computeCommunProb(cellchat_merged)

cellchat_merged <- computeCommunProbPathway(cellchat_merged)
cellchat_merged <- aggregateNet(cellchat_merged)

sources.use <- c("high_AEC", "low_AEC")

targets.use <- setdiff(levels(cellchat_merged@idents), sources.use)

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
enhanced_plot <- bubble_plot +
	scale_size_continuous(
		range = c(3, 10),
		name = "Interaction\nStrength",
		breaks = seq(0, max(bubble_plot$data$prob, na.rm = TRUE), length.out = 5)
	) +
	scale_color_gradientn(
		colors = c("blue", "yellow", "red"),
		name = "Communication\nProbability",
		limits = c(0, 1)
	)
print(enhanced_plot)
dev.off()

pdf("low_AEC_to_Merged_bubble.pdf", width = 14, height = 6)

bubble_plot <- netVisual_bubble(
	cellchat_merged,
	sources.use = "low_AEC",
	targets.use = targets.use,
	remove.isolate = FALSE,
	angle.x = 90
) +
	coord_flip() +
	ggtitle("Signaling from low_AEC to Merged Cell Types")

enhanced_plot <- bubble_plot +
	scale_size_continuous(
		range = c(3, 10),
		name = "Interaction\nStrength",
		breaks = seq(0, max(bubble_plot$data$prob, na.rm = TRUE), length.out = 5)
	) +
	scale_color_gradientn(
		colors = c("blue", "yellow", "red"),
		name = "Communication\nProbability",
		limits = c(0, 1)
	)
print(enhanced_plot)
dev.off()

cellchat_merged <- netAnalysis_computeCentrality(cellchat_merged, slot.name = "netP")
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pdf("netAnalysis_signalingRole_scatter.pdf",width = 6,height = 6)
netAnalysis_signalingRole_scatter(cellchat_merged)
dev.off()

pdf("netAnalysis_signalingRole_heatmap.pdf",width = 12,height = 6)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_merged, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_merged, pattern = "incoming")
ht1 + ht2
dev.off()

pathways.show = c("MIF","MK","ANNEXIN")
for (pathway in pathways.show) {
	pdf_filename <- paste0("netAnalysis_signalingRole_heatmap_", pathway, ".pdf")

	pdf(pdf_filename, width = 10, height = 6)

	netAnalysis_signalingRole_network(cellchat, pathway, width = 16, height = 4, font.size = 10)

	dev.off()
}
