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
rm(list = ls())
library(Seurat)
library(sctransform)
library(dplyr)
library(harmony)
library(clustree)

table(oc$celltype)
mm <- subset(oc, subset = celltype == "AEC") %>% NormalizeData() %>%
	FindVariableFeatures(selection.method = "vst" ,nfeatures = 1500 ,verbose = T)

mm <- RunPCA(mm, npcs = 50, verbose = FALSE)
# ElbowPlot(mm, ndims = 50)
dims_use <- 1:30

if ("orig.ident" %in% colnames(mm@meta.data)) {
   mm <- RunHarmony(mm, group.by.vars = "orig.ident", dims.use = dims_use)
   emb_use <- "harmony"
 } else {
   emb_use <- "pca"
 }


mm <- RunUMAP(mm, reduction = emb_use, dims = dims_use,
							n.neighbors = 30, min.dist = 0.2, verbose = FALSE)
mm <- FindNeighbors(mm, reduction = emb_use, dims = dims_use, k.param = 30)

for (res in c(0.2, 0.4, 0.6, 0.8, 1.0)) {
	mm <- FindClusters(mm, resolution = res, verbose = FALSE)
}

library(SCP)
library(DDRTree)
library(monocle)
mm[["RNA"]] <- as(object = mm[["RNA"]], Class = "Assay")
DefaultAssay(mm) <- "RNA"
mm=RunMonocle2(mm)

pdf("State_umap.PDF", width = 6, height = 6)
CellDimPlot(
	srt = mm, group.by = c("Monocle2_State"), reduction = "umap", theme_use = "theme_blank"
)
dev.off()
source("./DynamicPlot.R")
feats <- c("RCC2")

brown_line <- "#A6793D"   # 棕色（近似示意图）
line_palcolor <- setNames(rep(brown_line, length(feats)), feats)

grp <- "Monocle2_State"
lv  <- levels(factor(mm@meta.data[[grp]]))
tableau_like <- c(
	"#4E79A7",  # 蓝  (blue)
	"#59A14F",  # 绿  (green)
	"#F28E2B",  # 橙  (orange)
	"#E15759",  # 红  (red)
	"#76B7B2",  # 青  (teal)
	"#EDC948",  # 黄
	"#B07AA1",  # 紫
	"#FF9DA7",  # 粉
	"#9C755F",  # 棕
	"#BAB0AC"   # 灰
)
point_palcolor <- setNames(tableau_like[seq_along(lv)], lv)

p <- DynamicPlotMonocle2(
	srt = mm,
	pseudotime = "Monocle2_Pseudotime",
	features = feats,
	group.by = grp,
	slot = "count",
	exp_method = "log1p",
	add_interval = TRUE,
	add_point = TRUE,
	pt.size = 1,
	line.size = 1.2,
	legend.position = "bottom",
	line_palcolor  = line_palcolor,
	point_palcolor = point_palcolor,
	theme_use  = theme_scp,
	theme_args = list(base_size = 12)
)

ggplot2::ggsave("pseudotime_style_like_example.pdf", p, width = 6, height = 6)
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

state_order <- pData(cds_subset) %>%
	group_by(State) %>%
	summarise(median_pseudotime = median(Pseudotime)) %>%
	arrange(median_pseudotime) %>%
	pull(State)

pData(cds_subset)$State <- factor(pData(cds_subset)$State,
																	levels = state_order)

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
		y = "Expression",
		title = "RCC2 Expression Dynamics"
	)

dev.off()
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
		Pseudotime = mm$Monocle2_Pseudotime,
		Score = mm@meta.data[[ptm]],
		State = mm$Monocle2_State
	)
	colnames(data) <- c("Pseudotime", "Score", "State")

	data$State <- as.factor(data$State)

	cor_test <- cor.test(data$Pseudotime, data$Score, method = "spearman")
	cor_value <- round(cor_test$estimate, 3)
	p_value <- round(cor_test$p.value, 5)

	cor_results <- rbind(cor_results, data.frame(
		PTM = ptm,
		Correlation = cor_value,
		P_value = p_value
	))

	safe_filename <- gsub("[^[:alnum:]_]", "_", ptm)
	pdf_file <- paste0("cor_Pseudotime_", safe_filename, ".pdf")

	cor_label <- paste0("cor = ", cor_value, "\n", "p = ", format.pval(p_value, digits = 3))

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
			color = "#2c3e50",          # 深蓝色趋势线
			se = TRUE,
			fill = "#3498db",           # 浅蓝色置信区间
			linewidth = 1.8,
			alpha = 0.2
		) +
		scale_fill_brewer(palette = "Set1") +
		scale_color_brewer(palette = "Set1") +

		theme_classic(base_size = 14) +
		labs(
			x = "Pseudotime",
			y = paste(ptm, "Activity"),
			color = "State",
			fill = "State"
		) +
		theme(
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
			geom = "text",
			x = Inf, y = Inf,
			label = cor_label,
			hjust = 1.1, vjust = 1.1,
			size = 6,
			color = "black",
			fontface = "bold"
		)

	ggsave(
		pdf_file,
		plot = p,
		width = 7,
		height = 6,
		dpi = 300
	)
}

write.csv(cor_results, "pseudotime_correlation_results.csv", row.names = FALSE)
library(RcppML)
library(irGSEA)
mm$Monocle2_State

DefaultAssay(mm) <- "RNA"
mm <- SeuratObject::UpdateSeuratObject(object = mm)

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
	group.by = "Monocle2_State",
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

# ---- packages ----
library(ggplot2)
library(ggpubr)
library(gghalves)

plot_data <- data.frame(
	Glycosylation_expression = mm$Glycerophospholipid_Metabolism,
	Subtype = mm$Monocle2_State
)
colnames(plot_data) <- c("Glycosylation_expression", "Subtype")
plot_data$Subtype <- as.factor(plot_data$Subtype)

# Kruskal-Wallis
kruskal_p_value <- kruskal.test(Glycosylation_expression ~ Subtype, data = plot_data)$p.value
sig_label <- ifelse(kruskal_p_value < 0.001, "***",
										ifelse(kruskal_p_value < 0.01, "**",
													 ifelse(kruskal_p_value < 0.05, "*", "")))

my_comparisons <- list(c("1", "5"), c("2", "5"), c("3", "5"), c("4", "5"))

y_top <- max(plot_data$Glycosylation_expression, na.rm = TRUE)
label_y <- c(y_top * 1.03,  # 1 vs 2
						 y_top * 1.09,  # 1 vs 3
						 y_top * 1.15,
						 y_top * 1.21)  # 2 vs 3

y_min <- min(plot_data$Glycosylation_expression, na.rm = TRUE)
y_max <- y_top * 1.25

set.seed(1)

p <- ggplot(plot_data, aes(x = Subtype, y = Glycosylation_expression, fill = Subtype)) +
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
	gghalves::geom_half_point(
		side = "l",
		shape = 21,
		color = "black", stroke = 0.25,
		size = 1.4,
		alpha = 0.6,
		position = position_nudge(x = 0.09),
		transformation = position_jitter(width = 0.05, height = 0),
		range_scale = 0.9
	) +
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
	ggpubr::stat_compare_means(
		comparisons = my_comparisons,
		method = "wilcox.test",
		label = "p.signif",
		tip.length = 0.01,
		size = 4.5,
		vjust = 0.4,
		step.increase = 0,
		label.y = label_y
	) +
	scale_fill_manual(values = c("1" = "#acd5ab", "2" = "#feadac", "3" = "#adeada","4" = "#adeeee","5" = "#adaaaa")) +

	labs(
		x = "State",
		y = "Glycerophospholipid_Metabolism",
		title = "Glycerophospholipid_Metabolism Across Cell States"
	) +

	theme_minimal(base_size = 14) +
	theme(
		plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
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
		panel.background = element_rect(fill = "white", color = NA)
	) +
	coord_cartesian(ylim = c(y_min, y_max), clip = "off")

ggsave("Glycerophospholipid_Metabolism_state.pdf", plot = p, width = 9, height = 8, dpi = 300)

#####RCC2_State####
# ---- packages ----
library(ggplot2)
library(ggpubr)
library(gghalves)

plot_data <- data.frame(
	Glycosylation_expression = mm@assays$RNA$data["RCC2",],
	Subtype = mm$Monocle2_State
)
colnames(plot_data) <- c("Glycosylation_expression", "Subtype")
plot_data$Subtype <- as.factor(plot_data$Subtype)

# Kruskal-Wallis
kruskal_p_value <- kruskal.test(Glycosylation_expression ~ Subtype, data = plot_data)$p.value
sig_label <- ifelse(kruskal_p_value < 0.001, "***",
										ifelse(kruskal_p_value < 0.01, "**",
													 ifelse(kruskal_p_value < 0.05, "*", "")))

my_comparisons <- list(c("1", "5"), c("2", "5"), c("3", "5"), c("4", "5"))

y_top <- max(plot_data$Glycosylation_expression, na.rm = TRUE)
label_y <- c(y_top * 1.03,  # 1 vs 2
						 y_top * 1.09,  # 1 vs 3
						 y_top * 1.15,
						 y_top * 1.21)  # 2 vs 3

y_min <- min(plot_data$Glycosylation_expression, na.rm = TRUE)
y_max <- y_top * 1.25

set.seed(1)

p <- ggplot(plot_data, aes(x = Subtype, y = Glycosylation_expression, fill = Subtype)) +
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
	gghalves::geom_half_point(
		side = "l",
		shape = 21,
		color = "black", stroke = 0.25,
		size = 1.4,
		alpha = 0.6,
		position = position_nudge(x = 0.09),
		transformation = position_jitter(width = 0.05, height = 0),
		range_scale = 0.9
	) +
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
	ggpubr::stat_compare_means(
		comparisons = my_comparisons,
		method = "wilcox.test",
		label = "p.signif",
		tip.length = 0.01,
		size = 4.5,
		vjust = 0.4,
		step.increase = 0,
		label.y = label_y
	) +
	scale_fill_manual(values = c("1" = "#acd5ab", "2" = "#feadac", "3" = "#adeada","4" = "#adeeee","5" = "#adaaaa")) +

	labs(
		x = "State",
		y = "RCC2_Expression",
		title = "RCC2_Expression Across Cell States"
	) +

	theme_minimal(base_size = 14) +
	theme(
		plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
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
		panel.background = element_rect(fill = "white", color = NA)
	) +
	coord_cartesian(ylim = c(y_min, y_max), clip = "off")

ggsave("RCC2_Expression_state.pdf", plot = p, width = 9, height = 8, dpi = 300)

library(Mfuzz)
library(dplyr)
library(tidyr)
library(tidyverse)

ptm_cols <- names(signature_metabolism)
available_ptms <- ptm_cols[ptm_cols %in% colnames(mm@meta.data)]

ptm_data <- mm@meta.data %>%
	dplyr::select(Monocle2_State, all_of(available_ptms)) %>%
	mutate(Monocle2_State = as.factor(Monocle2_State))

state_means <- ptm_data %>%
	group_by(Monocle2_State) %>%
	summarise(across(all_of(available_ptms), mean, na.rm = TRUE))

result_df <- state_means %>%
	pivot_longer(cols = -Monocle2_State,
							 names_to = "Pathway",
							 values_to = "Mean_Value") %>%
	pivot_wider(names_from = Monocle2_State,
							values_from = Mean_Value)

result_df <- result_df %>% column_to_rownames("Pathway")

mfuzz_class <- new('ExpressionSet',exprs = as.matrix(result_df))

mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
mfuzz_class <- filter.std(mfuzz_class, min.std = 0)

mfuzz_class <- standardise(mfuzz_class)

set.seed(123)
cluster_num <- 6
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))

library(RColorBrewer)
library(ggplot2)

mfuzz_plot <- function(eset, cl, mfrow = c(1,1), colo, min.mem = 0, time.labels, new.window = TRUE) {
	if (missing(colo)) {
		colo <- colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(100)
	}

	opar <- par(no.readonly = TRUE)
	on.exit(par(opar))
	par(mfrow = mfrow, mar = c(5, 5, 4, 2) + 0.1)

	centers <- cl$centers
	ncluster <- nrow(centers)

	if (missing(time.labels)) {
		time.labels <- colnames(eset)
	}

	for (i in 1:ncluster) {
		center <- centers[i, ]

		plot(1, type = "n",
				 xlim = c(1, length(center)),
				 ylim = range(centers),
				 xlab = "", ylab = "",
				 axes = FALSE,
				 main = paste0("Cluster ", i, " (", sum(cl$cluster == i), " pathways)"),
				 cex.main = 1.8, font.main = 2, col.main = "#333333")

		axis(1, at = 1:length(center), labels = time.labels,
				 cex.axis = 1.5, col.axis = "#333333", col = "gray70", lwd = 1.5)
		axis(2, cex.axis = 1.5, col.axis = "#333333", col = "gray70", lwd = 1.5)

		grid(NA, NULL, col = "gray90", lty = 1, lwd = 0.8)
		abline(h = 0, col = "gray70", lwd = 1.5)

		lines(1:length(center), center,
					col = "#2c3e50", lwd = 3.5, type = "o", pch = 19, cex = 1.5)

		cluster_members <- exprs(eset)[cl$cluster == i, , drop = FALSE]
		for (j in 1:nrow(cluster_members)) {
			lines(1:length(center), cluster_members[j, ],
						col = adjustcolor("#3498db", alpha.f = 0.5),
						lwd = 1.2)
		}

		mtext("Cell State", side = 1, line = 3.5, cex = 1.2, font = 2, col = "#333333")
		mtext("Standardized Activity", side = 2, line = 3.5, cex = 1.2, font = 2, col = "#333333")

		rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
				 col = adjustcolor("#f8f9fa", alpha.f = 0.5), border = NA)

		box(lwd = 2, col = "gray70")
	}
}

pdf("plot_PTMD_cluster.pdf", width = 14, height = 10)
mfuzz_plot(
	mfuzz_class,
	mfuzz_cluster,
	mfrow = c(2, 3),
	colo = colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(100),  # 蓝-黄-红渐变
	time.labels = colnames(mfuzz_class)
)
dev.off()
library(ClusterGVis)
diff_test_res <- differentialGeneTest(cds_seurat,
																			fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res1 <- diff_test_res %>%
	filter(num_cells_expressed > 100)
df=plot_pseudotime_heatmap2(cds_seurat[row.names(subset(diff_test_res1,qval<1e-4)),],
														num_clusters=6)
gene=sample(df$wide.res$gene,20,replace=F)

pdf("plot_genes_heatmap.PDF", width=6, height=8, onefile=FALSE)
visCluster(object=df, plot.type="heatmap", markGenes=gene)
dev.off()

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

state_subset <- subset(mm, subset = Monocle2_State %in% c(1, 5))

available_ptms <- names(signature_metabolism)[names(signature_metabolism) %in% colnames(state_subset@meta.data)]

diff_results <- data.frame(
	Pathway = character(),
	Avg_State1 = numeric(),
	Avg_State5 = numeric(),
	Log2FC = numeric(),
	p_value = numeric(),
	stringsAsFactors = FALSE
)

for (ptm in available_ptms) {
	state1_data <- state_subset@meta.data %>%
		filter(Monocle2_State == 1) %>%
		pull(!!ptm)

	state5_data <- state_subset@meta.data %>%
		filter(Monocle2_State == 5) %>%
		pull(!!ptm)

	avg1 <- mean(state1_data, na.rm = TRUE)
	avg5 <- mean(state5_data, na.rm = TRUE)

	log2fc <- log2(avg5 / avg1)

	p_val <- wilcox.test(state1_data, state5_data)$p.value

	diff_results <- rbind(diff_results, data.frame(
		Pathway = ptm,
		Avg_State1 = avg1,
		Avg_State5 = avg5,
		Log2FC = log2fc,
		p_value = p_val
	))
}

diff_results$p_adj <- p.adjust(diff_results$p_value, method = "BH")

diff_results$Significance <- cut(
	diff_results$p_adj,
	breaks = c(0, 0.001, 0.01, 0.05, 1),
	labels = c("***", "**", "*", "NS"),
	include.lowest = TRUE
)

write.csv(diff_results, "pathway_diff_state1_vs_state5.csv", row.names = FALSE)

volcano_plot <- ggplot(diff_results, aes(x = Log2FC, y = -log10(p_adj))) +
	geom_point(
		aes(fill = -log10(p_adj), size = -log10(p_adj)),
		shape = 21,
		color = "black",
		stroke = 0.5,
		alpha = 0.9
	) +
	geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 0.8) +
	geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", size = 0.8) +

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

	annotate(
		"text",
		x = -Inf, y = Inf,
		label = "Enriched in State 1",
		size = 5.5, color = "#4575b4", fontface = "bold",
		hjust = -0.05, vjust = -0.6
	) +
	annotate(
		"text",
		x =  Inf, y = Inf,
		label = "Enriched in State 5",
		size = 5.5, color = "#d73027", fontface = "bold",
		hjust =  1.05, vjust = -0.6
	) +
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

	labs(
		x = expression(Log[2]~Fold~Change~(State~5~vs~State~1)),
		y = expression(-Log[10]~adjusted~p~value),
		title = "Pathway Activity Differences",
		subtitle = "State 1 vs State 5"
	) +

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

	coord_cartesian(
		ylim = c(0, max(-log10(diff_results$p_adj), na.rm = TRUE) * 1.1),
		clip = "off"
	)

ggsave("pathway_volcano_state1_vs_state5.pdf", volcano_plot, width = 10, height = 9)

select_top10_paths <- function(df, top_n = 5, sig_only = TRUE) {
	d <- df %>% filter(is.finite(Log2FC))
	pick_side <- function(dat, side = c("state3","state1")) {
		side <- match.arg(side)
		if (side == "state3") {
			cand <- dat %>% arrange(desc(Log2FC))
			sig  <- cand %>% filter(p_adj < 0.05)
			top  <- (if (sig_only) sig else cand) %>% slice_head(n = top_n) %>% pull(Pathway)
			if (length(top) < top_n)
				top <- unique(c(top, cand$Pathway))[seq_len(min(top_n, nrow(cand)))]
			top
		} else {
			cand <- dat %>% arrange(Log2FC)
			sig  <- cand %>% filter(p_adj < 0.05)
			top  <- (if (sig_only) sig else cand) %>% slice_head(n = top_n) %>% pull(Pathway)
			if (length(top) < top_n)
				top <- unique(c(top, cand$Pathway))[seq_len(min(top_n, nrow(cand)))]
			top
		}
	}
	unique(c(pick_side(d, "state3"), pick_side(d, "state1")))
}
plot_subtype_pathways <- function(pathway_list, title,
																	code_prefix = "M",
																	legend_cols = 2,
																	top_expand = 0.20,
																	legend_trunc = 50,
																	legend_replace_us = TRUE,
																	legend_text_size = 5,
																	star_size = 6,
																	legend_left_pad = 0.05,
																	legend_col_spacing = 1.6
) {
	code_map    <- setNames(paste0(code_prefix, seq_along(pathway_list)), pathway_list)
	levels_code <- unname(code_map[pathway_list])

	plot_data <- as.data.frame(state_subset@meta.data[, pathway_list, drop = FALSE]) %>%
		mutate(Subtype = state_subset$Monocle2_State) %>%
		tidyr::pivot_longer(cols = -Subtype, names_to = "Pathway", values_to = "Score") %>%
		mutate(
			Pathway = factor(Pathway, levels = pathway_list),
			Code    = factor(code_map[Pathway], levels = levels_code),
			Subtype = factor(Subtype)
		)

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
		label = sig_labels, stringsAsFactors = FALSE
	)

	y_min <- min(plot_data$Score, na.rm = TRUE)
	y_max <- max(plot_data$Score, na.rm = TRUE)
	y_lim <- c(y_min, y_max * (1 + top_expand))

	p_main <- ggplot(plot_data, aes(x = Code, y = Score, fill = Subtype)) +
		geom_vline(xintercept = seq(1.5, length(levels_code) - 0.5), color = "gray70",
							 linetype = "dashed", linewidth = 0.5) +
		gghalves::geom_half_violin(aes(fill = Subtype),
															 position = position_dodge(0.8),
															 side = "l", scale = "width",
															 width = 0.7, alpha = 0.5, color = NA) +
		geom_boxplot(aes(color = Subtype),
								 width = 0.3, position = position_dodge(0.8),
								 outlier.shape = NA, notch = TRUE, notchwidth = 0.5,
								 alpha = 0.8, linewidth = 0.6) +
		geom_point(aes(color = Subtype),
							 position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
							 size = 0.4, alpha = 0.4) +
		geom_text(data = sig_df,
							aes(x = Code, y = y_lim[2] - 0.02 * diff(y_lim), label = label),
							inherit.aes = FALSE, size = star_size, vjust = 0.5, fontface = "bold") +
		scale_fill_manual(values = c("5" = "#E64B35", "1" = "#4DBBD5"), name = "Group") +
		scale_color_manual(values = c("5" = "#B71C1C", "1" = "#006064"), name = "Group", guide = "none") +
		labs(title = title, x = "", y = "Pathway Activity Score") +
		theme_minimal(base_size = 13) +
		theme(
			axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 11, color = "#333333"),
			axis.text.y = element_text(face = "bold", size = 11, color = "#333333"),
			axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
			plot.title   = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 8)),
			legend.position = "top",
			legend.title    = element_text(face = "bold", size = 12),
			legend.text     = element_text(size = 11),
			panel.grid.major.y = element_line(color = "grey92", linewidth = 0.6),
			panel.grid.minor   = element_blank(),
			panel.background    = element_rect(fill = "white", color = NA),
			plot.background     = element_rect(fill = "white", color = NA),
			panel.border        = element_rect(color = "grey80", fill = NA, linewidth = 0.8),
			plot.margin         = margin(10, 15, 5, 15)
		) +
		coord_cartesian(ylim = y_lim, clip = "off")

	legend_map <- tibble(
		Pathway = pathway_list,
		Code    = unname(code_map[pathway_list])
	) %>%
		mutate(
			PrettyName = if (legend_replace_us) stringr::str_replace_all(Pathway, "_", " ") else Pathway,
			LabelFull  = paste0(Code, ": ", PrettyName),
			Label      = stringr::str_wrap(LabelFull, width = legend_trunc),
			idx        = dplyr::row_number(),
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

	ggpubr::ggarrange(p_main, p_legend, ncol = 1, heights = c(1, 0.30))
}

sel_paths <- select_top10_paths(diff_results, top_n = 5, sig_only = TRUE)
sig_plot <- plot_subtype_pathways(
	pathway_list = sel_paths,
	title        = "Top pathways (State5↑ & State1↑)",
	legend_cols  = 2,
	legend_trunc = 56,
	legend_text_size = 6,
	star_size = 7,
	legend_left_pad = 0.02,
	legend_col_spacing = 1.8
)
ggsave("signature_pathways_top10.pdf", sig_plot, width = 12, height = 7)

# ---- packages ----
library(ggplot2)
library(ggpubr)
library(gghalves)

plot_data <- data.frame(
	Glycosylation_expression = mm$Scoreing,
	Subtype = mm$Monocle2_State
)
colnames(plot_data) <- c("Glycosylation_expression", "Subtype")
plot_data$Subtype <- as.factor(plot_data$Subtype)

# Kruskal-Wallis
kruskal_p_value <- kruskal.test(Glycosylation_expression ~ Subtype, data = plot_data)$p.value
sig_label <- ifelse(kruskal_p_value < 0.001, "***",
										ifelse(kruskal_p_value < 0.01, "**",
													 ifelse(kruskal_p_value < 0.05, "*", "")))

my_comparisons <- list(c("1", "5"), c("2", "5"), c("3", "5"), c("4", "5"))

y_top <- max(plot_data$Glycosylation_expression, na.rm = TRUE)
label_y <- c(y_top * 1.03,  # 1 vs 2
						 y_top * 1.09,  # 1 vs 3
						 y_top * 1.15,
						 y_top * 1.21)  # 2 vs 3

y_min <- min(plot_data$Glycosylation_expression, na.rm = TRUE)
y_max <- y_top * 1.25

set.seed(1)

p <- ggplot(plot_data, aes(x = Subtype, y = Glycosylation_expression, fill = Subtype)) +
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
	gghalves::geom_half_point(
		side = "l",
		shape = 21,
		color = "black", stroke = 0.25,
		size = 1.4,
		alpha = 0.6,
		position = position_nudge(x = 0.09),
		transformation = position_jitter(width = 0.05, height = 0),
		range_scale = 0.9
	) +
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
	ggpubr::stat_compare_means(
		comparisons = my_comparisons,
		method = "wilcox.test",
		label = "p.signif",
		tip.length = 0.01,
		size = 4.5,
		vjust = 0.4,
		step.increase = 0,
		label.y = label_y
	) +
	scale_fill_manual(values = c("1" = "#acd5ab", "2" = "#feadac", "3" = "#adeada","4" = "#adeeee","5" = "#adaaaa")) +

	labs(
		x = "State",
		y = "Scoreing of Lactylation",
		title = "Scoreing Across Cell States"
	) +

	theme_minimal(base_size = 14) +
	theme(
		plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
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
		panel.background = element_rect(fill = "white", color = NA)
	) +
	coord_cartesian(ylim = c(y_min, y_max), clip = "off")

ggsave("Scoreing_state.pdf", plot = p, width = 9, height = 8, dpi = 300)

# ---- packages ----
library(ggplot2)
library(ggpubr)
library(gghalves)

plot_data <- data.frame(
	Glycosylation_expression = mm$LMBDs,
	Subtype = mm$Monocle2_State
)
colnames(plot_data) <- c("Glycosylation_expression", "Subtype")
plot_data$Subtype <- as.factor(plot_data$Subtype)

# Kruskal-Wallis
kruskal_p_value <- kruskal.test(Glycosylation_expression ~ Subtype, data = plot_data)$p.value
sig_label <- ifelse(kruskal_p_value < 0.001, "***",
										ifelse(kruskal_p_value < 0.01, "**",
													 ifelse(kruskal_p_value < 0.05, "*", "")))

my_comparisons <- list(c("1", "5"), c("2", "5"), c("3", "5"), c("4", "5"))

y_top <- max(plot_data$Glycosylation_expression, na.rm = TRUE)
label_y <- c(y_top * 1.03,  # 1 vs 2
						 y_top * 1.09,  # 1 vs 3
						 y_top * 1.15,
						 y_top * 1.21)  # 2 vs 3

y_min <- min(plot_data$Glycosylation_expression, na.rm = TRUE)
y_max <- y_top * 1.25

set.seed(1)

p <- ggplot(plot_data, aes(x = Subtype, y = Glycosylation_expression, fill = Subtype)) +
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
	gghalves::geom_half_point(
		side = "l",
		shape = 21,
		color = "black", stroke = 0.25,
		size = 1.4,
		alpha = 0.6,
		position = position_nudge(x = 0.09),
		transformation = position_jitter(width = 0.05, height = 0),
		range_scale = 0.9
	) +
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
	ggpubr::stat_compare_means(
		comparisons = my_comparisons,
		method = "wilcox.test",
		label = "p.signif",
		tip.length = 0.01,
		size = 4.5,
		vjust = 0.4,
		step.increase = 0,
		label.y = label_y
	) +
	scale_fill_manual(values = c("1" = "#acd5ab", "2" = "#feadac", "3" = "#adeada","4" = "#adeeee","5" = "#adaaaa")) +

	labs(
		x = "State",
		y = "LMBDs",
		title = "LMBDs Across Cell States"
	) +

	theme_minimal(base_size = 14) +
	theme(
		plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
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
		panel.background = element_rect(fill = "white", color = NA)
	) +
	coord_cartesian(ylim = c(y_min, y_max), clip = "off")

ggsave("LMBDs_state.pdf", plot = p, width = 9, height = 8, dpi = 300)


tar <- mm$Scoreing
result_sub <- mm$LMBDs
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
		 ylab = "LMBDs",
		 col = "white",
		 main = "")

rect(par("usr")[1], par("usr")[3],
		 par("usr")[2], par("usr")[4],
		 col = "#EAE9E9", border = FALSE)

grid(col = "white", lty = 1, lwd = 1.5)

points(x = result_sub, y = tar,
			 pch = 19,
			 col = scales::alpha("#E51718", 0.8),  # 使用单一红色系
			 cex = 1)

abline(lm(tar ~ result_sub),
			 lwd = 4,
			 col = "black")

rug(result_sub, side = 3, col = "black", lwd = 1)
rug(tar, side = 4, col = "black", lwd = 1)

text(x = min(result_sub),
		 y = max(tar) * 0.95,
		 adj = 0,
		 labels = bquote("AEC: N = " ~ .(n) ~
		 									"; " ~ rho ~ " = " ~ .(cor_rho) ~
		 									"; " ~ italic(P) ~ " = " ~ .(cor_p)),
		 col = "black",
		 cex = 0.8)
dev.off()
