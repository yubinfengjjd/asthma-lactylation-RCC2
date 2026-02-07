library(sva)
library(IOBR)
library(GSVA)
library(ConsensusClusterPlus)
library(Seurat)
data(signature_metabolism)
setwd("D:\\哮喘乳酸化单细胞开篇\\6.一致性聚类")
#gene=data.table::fread("./hdWGCNAgene.txt",header = F,data.table = F)[,1]
#exp = exp %>% column_to_rownames("ID")
metabolic_genesets <- lapply(signature_metabolism, function(x) unlist(strsplit(x, split = ",")))
oc_subset <- subset(oc, subset = celltype == "Mono/Mac") %>% NormalizeData() %>%
	FindVariableFeatures(selection.method = "vst" ,nfeatures = 1500 ,verbose = T)
exp = as.matrix(oc_subset@assays$RNA$data)
Dim(exp)
################################
################################
gsva_results <- ssgseaParam(
	exp,
  metabolic_genesets,minSize = 5,normalize = F
)

result <- gsva(gsva_results)
result_scaled <- t(scale(t(result)))
input = as.matrix(result_scaled)
summary(as.vector(input))
####################################
####################################
result_cluster <- ConsensusClusterPlus(
	d = input,
	maxK = 3,
	reps = 1000,
	pItem = 0.8,
	pFeature = 1,
	clusterAlg = "pam",
	distance = "euclidean",
	seed = 123456,plot = 'pdf'
)
save(result_cluster,file = "result_cluster.Rdata")
cluster_assignment <- result_cluster[[2]]$consensusClass
oc_subset$subtypes <- ifelse(cluster_assignment == 1, "MBC1", "MBC2")
save(oc_subset,file = "Cluster.Rdata")
consensus_matrix <- result_cluster[[2]]$consensusMatrix
colnames(consensus_matrix) <- rownames(consensus_matrix) <- colnames(input)

order_samples <- order(subtypes)
sorted_consensus <- consensus_matrix[order_samples, order_samples]
sorted_subtypes <- subtypes[order_samples]

library(pheatmap)

annotation_colors <- list(
	Subtype = c(MBC1 = "#FF6F61", MBC2 = "#92A8D1")  # 自定义颜色
)

annotation_df <- data.frame(
	Subtype = sorted_subtypes,
	row.names = colnames(sorted_consensus)
)

pdf(file = "1.pdf",width = 6,height = 6)
pheatmap(
	mat = sorted_consensus,
	color = colorRampPalette(c("white", "#0072B2"))(100),  # 蓝白渐变
	cluster_rows = FALSE,
	cluster_cols = FALSE,
	annotation_col = annotation_df,
	annotation_colors = annotation_colors,
	show_rownames = FALSE,
	show_colnames = FALSE,
	main = "Consensus Matrix for K=2 Clustering",
	gaps_col = cumsum(table(sorted_subtypes)),
	gaps_row = cumsum(table(sorted_subtypes))
)
dev.off()


###############################
###############################

library(vegan)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(multcompView)

normalize=function(x){
	return((x-min(x))/(max(x)-min(x)))}
input_PCOA=normalize(result)
input_PCOA=input_PCOA[apply(input_PCOA,1,sd)>0.01,]
pcoa_input <- t(input_PCOA)
rownames(pcoa_input) <- colnames(input_PCOA)

subtypes <- factor(subtypes)
meta_data <- data.frame(
	sample = colnames(input),
	Subtype = subtypes
)

dist_matrix <- vegdist(pcoa_input, method = "bray")

pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)
points<-data.frame(scores(pcoa_result))
point<-data.frame(sample=row.names(points),points)
pc1<-round((pcoa_result$eig/sum(pcoa_result$eig))*100,2)[1]
pc2<-round((pcoa_result$eig/sum(pcoa_result$eig))*100,2)[2]

colnames(points)[1:2]<-c('dim1','dim2')
plotdata=data.frame(rownames(points),points$dim1,points$dim2,meta_data$Subtype)
colnames(plotdata)=c("sample","dim1","dim2","group")


adonis_result_dis <- adonis2(dist_matrix ~ Subtype, data = meta_data, permutations = 999)
R2 = adonis_result_dis$R2[1]
pvalue = adonis_result_dis$`Pr(>F)`[1]
adonis <- paste("PERMANOVA:\nR2 = ", round(R2, 4), "\nP-value = ", pvalue)
length_samples <- length(unique(as.character(meta_data$sample)))
length_groups <- length(unique(as.character(meta_data$Subtype)))

times1 <- length_samples %/% 8
res1 <- length_samples %% 8
times2 <- length_samples %/% 5
res2 <- length_samples %% 5

mycol <- c("#B2182B", "#E69F00", "#56B4E9", "#009E73",
					 "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
					 "#CC6666", "#9999CC", "#66CC99", "#999999", "#ADD1E5")

col1 <- rep(mycol, times1)
col <- c(col1, mycol[1:res1])

pich1 <- rep(c(15:18, 20, 7:14, 0:6), times2)
pich <- c(pich1, 15:(15 + res2))

if (length_groups > 30) {
	n <- 2
} else {
	n <- 1
}
plot <- ggplot(plotdata, aes(dim1, dim2)) +
	geom_point(aes(colour = group, shape = group), size = 2) +
	scale_colour_manual(values = col) +
	xlab(paste("PCoA axis1 (", pc1, "%)", sep = "")) +
	ylab(paste("PCoA axis2 (", pc2, "%)", sep = "")) +
	theme(
		plot.title = element_text(hjust = 0.5, size = 18, colour = "black", face = "bold"),
		panel.background = element_rect(fill = 'white', colour = 'black'),
		panel.grid = element_blank(),
		axis.title = element_text(color = 'black', size = 20),
		axis.text = element_text(colour = 'black', size = 16, margin = unit(0.6, "lines")),
		axis.ticks = element_line(color = 'black'),
		axis.ticks.length = unit(0.4, "lines"),
		legend.title = element_blank(),
		legend.key = element_blank(),
		legend.position = c(0.9, 0.1),
		legend.background = element_rect(colour = "black")
	) +
	geom_vline(xintercept = 0, linetype = "dotted") +
	geom_hline(yintercept = 0, linetype = "dotted") +
	guides(
		col = guide_legend(ncol = n),
		shape = guide_legend(ncol = n)
	) +
	stat_ellipse(aes(x = dim1, y = dim2, color = group), data = plotdata)

print(plot)
box1 <- ggplot(plotdata, aes(x = meta_data$Subtype, y = dim1, fill = meta_data$Subtype)) +
	geom_boxplot(show.legend = FALSE) +
	stat_boxplot(geom = "errorbar", width = 0.1, size = 0.5) +
	scale_fill_manual(values = col) +
	coord_flip() +
	theme_bw() +
	theme(
		panel.grid = element_blank(),
		axis.title = element_blank(),
		axis.line = element_line(colour = "black"),
		axis.ticks = element_line(color = 'black'),
		axis.text.x = element_blank(),
		axis.text.y = element_text(colour = 'black', size = 16)
	)
box1
box2 <- ggplot(plotdata, aes(x = meta_data$Subtype, y = dim2, fill = meta_data$Subtype)) +
	geom_boxplot(show.legend = FALSE) +
	stat_boxplot(geom = "errorbar", width = 0.1, size = 0.5) +
	scale_fill_manual(values = col) +
	theme_bw() +
	theme(
		panel.grid = element_blank(),
		axis.title = element_blank(),
		axis.line = element_line(colour = "black"),
		axis.ticks = element_line(color = 'black'),
		axis.text.x = element_text(colour = 'black', size = 16, angle = 45, vjust = 1, hjust = 1),
		axis.text.y = element_blank()
	)
box3 <- ggplot(plotdata, aes(dim1, dim2)) +
	geom_text(
		aes(x = -0.5, y = 0.6, label = adonis),
		size = 2
	) +
	theme_bw() +
	xlab("") +
	ylab("") +
	theme(
		panel.grid = element_blank(),
		axis.title = element_blank(),
		axis.line = element_blank(),
		axis.ticks = element_blank(),
		axis.text = element_blank()
	)

p <- box1 + box3 + plot + box2 +
	plot_layout(heights = c(1, 4), widths = c(4, 1), ncol = 2, nrow = 2)
p
pdf(file = "PCoA.pdf",width = 6,height = 6)
print(p)
dev.off()

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

pathway_diff <- apply(result_scaled, 1, function(x) {
	mean(x[subtypes == "MBC2"]) - mean(x[subtypes == "MBC1"])
})

mbc2_up <- pathway_diff %>%
	sort(decreasing = TRUE) %>%
	head(10) %>%
	names()

mbc1_up <- pathway_diff %>%
	sort() %>%
	head(10) %>%
	names()

plot_subtype_pathways <- function(pathway_list, title) {
	plot_data <- as.data.frame(t(result_scaled[pathway_list, ])) %>%
		mutate(Subtype = subtypes) %>%
		pivot_longer(cols = -Subtype,
								 names_to = "Pathway",
								 values_to = "Score") %>%
		mutate(Pathway = factor(Pathway, levels = pathway_list))

	p_values <- sapply(pathway_list, function(p) {
		wilcox.test(Score ~ Subtype, data = filter(plot_data, Pathway == p))$p.value
	})
	sig_labels <- ifelse(p_values < 0.001, "***",
											 ifelse(p_values < 0.01, "**",
											 			 ifelse(p_values < 0.05, "*", "")))

	ggplot(plot_data, aes(x = Pathway, y = Score)) +
		geom_violin(aes(fill = Subtype), position = position_dodge(0.9), alpha = 0.4, color = NA) +
		geom_boxplot(aes(fill = Subtype),
								 width = 0.15,
								 position = position_dodge(0.9),
								 outlier.shape = NA,
								 color = "black") +
		geom_text(
			data = data.frame(Pathway = pathway_list, label = sig_labels),
			aes(x = Pathway, y = max(plot_data$Score)*1.1, label = label),
			size = 5, vjust = -0.5
		) +
		scale_fill_manual(values = c("MBC1"="#FF6F61", "MBC2"="#92A8D1")) +
		labs(title = title, x = "", y = "GSVA Enrichment Score") +
		theme_classic() +
		theme(
			axis.text.x = element_text(angle = 45, hjust = 1, face = "bold",size = 6),
			plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
			legend.position = "top",
			panel.border = element_rect(color = "black", fill = NA),
			panel.grid.major.y = element_line(color = "grey90")
		) +
		coord_cartesian(ylim = c(min(plot_data$Score), max(plot_data$Score)*1.15))
}

p_mbc2 <- plot_subtype_pathways(
	pathway_list = mbc2_up,
	title = "Top 10 Pathways Upregulated in MBC2"
)

p_mbc1 <- plot_subtype_pathways(
	pathway_list = mbc1_up,
	title = "Top 10 Pathways Upregulated in MBC1"
)

ggsave("MBC2_Top10_Upregulated.pdf", p_mbc2, width = 12, height = 6)
ggsave("MBC1_Top10_Upregulated.pdf", p_mbc1, width = 12, height = 6)

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

plot_data <- oc_subset@meta.data %>%
	select(Scoreing, subtypes)
colnames(plot_data)=c("Scoreing","Subtype")
plot_data$Subtype <- as.factor(plot_data$Subtype)

p_value <- wilcox.test(Scoreing ~ Subtype, data = plot_data)$p.value

sig_label <- ifelse(p_value < 0.001, "***",
										ifelse(p_value < 0.01, "**",
													 ifelse(p_value < 0.05, "*", "")))


ggplot(plot_data, aes(x = Subtype, y = Scoreing)) +
	geom_violin(aes(fill = Subtype), position = position_dodge(0.9), alpha = 0.7) +
	geom_boxplot(aes(fill = Subtype),
							 width = 0.15,
							 position = position_dodge(0.9),
							 outlier.shape = NA,
							 color = "black") +
	geom_text(
		aes(x = 1.5, y = max(plot_data$Scoreing)*1.1, label = sig_label),
		size = 5, vjust = -0.5
	) +
	scale_fill_manual(values = c("MBC1"="#acd5ab", "MBC2"="#feadac")) +
	labs(x = "Subtype", y = "Lactylation") +
	theme_classic() +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
		plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
		legend.position = "top",
		panel.border = element_rect(color = "black", fill = NA),
		panel.grid.major.y = element_line(color = "grey90")
	) +
	coord_cartesian(ylim = c(min(plot_data$Scoreing), max(plot_data$Scoreing)*1.15))

ggsave("Subtype_Scoreing.pdf", width = 4, height = 6)

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

plot_data <- data.frame(oc_subset$subtypes, t(result_scaled))
colnames(plot_data)[1] <- "subtypes"
plot_data <- plot_data %>%
	select(subtypes, Th1_Th2_diff, NFKB, JAK_STAT) %>%
	pivot_longer(
		cols = -subtypes,
		names_to = "Pathway",
		values_to = "Score"
	) %>%
	mutate(
		Subtype = as.factor(subtypes),
		Pathway = factor(Pathway,
										 levels = c("Th1_Th2_diff", "NFKB", "JAK_STAT"),
										 labels = c("Th1/Th2", "NF-kappaB", "JAK-STAT"))
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
	scale_fill_manual(values = c("MBC1" = "#acd5ab", "MBC2" = "#feadac")) +
	scale_color_manual(values = c("MBC1" = "#6a9662", "MBC2" = "#e67c73")) +
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
