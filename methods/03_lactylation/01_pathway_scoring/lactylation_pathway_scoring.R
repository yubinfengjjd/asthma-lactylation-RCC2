library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(GSEABase)
library(AUCell)
library(GSVA)
library(UCell)
library(singscore)
setwd("D:\\")
Addsingscore<-function(oc,signatures,seed=1,assay=null,slot="data"){
	set.seed(seed)
	if(is.null(assay)){
		assay=Seurat::DefaultAssay(oc)
	}
	matrix=Seurat::GetAssayData(object=oc,slot=slot,assay=assay)

	markers=signatures
	h.gsets.list=markers %>% purrr::compact()

	singscore.rank=singscore::rankGenes(as.data.frame(matrix))
	#calculate separately
	singscore.scores=list()
	for (i in seq_along(h.gsets.list)){
		if (any(stringr::str_detect(h.gsets.list[[i]],pattern = "\\+$|-$"))){
			h.gsets.list.positive=stringr::str_match(h.gsets.list[[i]],pattern="(.+)\\+")[,2] %>% purrr::discard(is.na)
			h.gsets.list.negative=stringr::str_match(h.gsets.list[[i]],pattern="(.+)-")[,2] %>% purrr::discard(is.na)
			if(length(h.gsets.list.positive)==0){
				singscore.scores[[i]]=singscore::simpleScore(singscore.rank,upSet=h.gsets.list.negative,centerScore=F)
			}
			if(length(h.gsets.list.negative)==0){
				singscore.scores[[i]]=singscore::simpleScore(singscore.rank,upSet=h.gsets.list.positive,centerScore=F)
			}
			if((length(h.gsets.list.positive)!=0)&(length(h.gsets.list.negative)!=0)){
				singscore.scores[[i]]=singscore::simpleScore(singscore.rank,upSet=h.gsets.list.positive,downSet=h.gsets.list.negative,centerScore=F)
			}
		}else{
			singscore.scores[[i]]=singscore::simpleScore(singscore.rank,upSet=h.gsets.list[[i]],centerScore=F)
		}
		TotalScore=NULL
		singscore.scores[[i]]=singscore.scores[[i]] %>%
			dplyr::select(TotalScore) %>%
			magrittr::set_colnames(names(h.gsets.list)[i])}
	names(singscore.scores)=names(h.gsets.list)
	singscore.scores=do.call(cbind,singscore.scores)
	oc=Seurat::AddMetaData(oc,as.data.frame(singscore.scores))
	return(oc)
}
DefaultAssay(oc) <- "RNA"

folders <- c("./GLYCOLYSIS_GLUCONEOGENESIS/", "./HIF_1_SIGNALING_PATHWAY/",
						 "./LACTATE_DEHYDROGENASE_ACTIVITY/", "./LACTATE_METABOLIC_PROCESS/")

read_folder_gmt <- function(folder_path){
	gmt_files <- list.files(folder_path, pattern = "\\.gmt$", full.names = TRUE)
	gene_sets <- list()
	for(f in gmt_files){
		gmt <- getGmt(f)
		set_names <- names(geneIds(gmt))
		genes <- geneIds(gmt)
		for(i in seq_along(set_names)){
			new_name <- paste0(tools::file_path_sans_ext(basename(f)),
												 "_", set_names[i])
			gene_sets[[new_name]] <- unique(genes[[i]])
		}
	}
	return(unique(unlist(gene_sets)))
}

for(folder in folders){
	folder_name <- basename(folder)
	cat("Processing folder:", folder_name, "\n")

	genes <- read_folder_gmt(folder)
	if(length(genes) == 0) next

	genes_list <- list(genes)
	genes_vector <- as.character(genes)

	##### AddModuleScore #####
	oc <- AddModuleScore(oc,
											 features = genes_list,
											 name = paste0(folder_name, "_Add_Score"),
											 slot = "data")
	new_col <- paste0(folder_name, "_AddModule")
	colnames(oc@meta.data)[ncol(oc@meta.data)] <- new_col

	##### AUCell #####
	cells_rankings <- AUCell_buildRankings(oc@assays$RNA$data, nCores=1)
	cells_AUC <- AUCell_calcAUC(list(genes = genes_vector), cells_rankings,
															aucMaxRank = nrow(cells_rankings)*0.05)
	oc@meta.data[[paste0(folder_name, "_AUCell")]] <- getAUC(cells_AUC)["genes", ]

	##### ssGSEA #####
	genes_ssGSEA = as.data.frame(genes)
	gene.expr <- as.matrix(oc@assays$RNA$data)
	ssGSEA.result <- ssgseaParam(gene.expr, genes_ssGSEA)
	ssGSEA.result <- gsva(ssGSEA.result)
	oc@meta.data[[paste0(folder_name, "_ssGSEA")]] <- ssGSEA.result["genes", ]

	##### Ucell #####
	oc <- AddModuleScore_UCell(oc, features = list(temp = genes_vector),
														 name = paste0(folder_name, "_Ucell_"),
														 assay = "RNA", slot = "data")
	colnames(oc@meta.data)[ncol(oc@meta.data)] <- paste0(folder_name, "_Ucell")

	##### singscore #####
	oc <- Addsingscore(oc = oc, signatures = list(temp = genes_vector),
										 assay = "RNA", slot = "data")
	colnames(oc@meta.data)[ncol(oc@meta.data)] <- paste0(folder_name, "_singscore")

	score_columns <- c(paste0(folder_name, "_AddModule"),
										 paste0(folder_name, "_AUCell"),
										 paste0(folder_name, "_ssGSEA"),
										 paste0(folder_name, "_Ucell"),
										 paste0(folder_name, "_singscore"))

	oc@meta.data[[paste0(folder_name)]] <- rowSums(
		oc@meta.data[, score_columns], na.rm = TRUE)
}

saveRDS(oc, "oc_with_pathway_scores.rds")

macrophage_cells <- oc$celltype == "Macrophage"
mnda_expr <- oc@assays$RNA$counts["MNDA", macrophage_cells]

median_mnda <- median(as.vector(mnda_expr))

oc$celltype <- case_when(
	oc@assays$RNA$counts["MNDA", ] > median_mnda & oc$celltype == "Macrophage" ~ "MNDA+ Macrophage",
	oc@assays$RNA$counts["MNDA", ] <= median_mnda & oc$celltype == "Macrophage" ~ "MNDA- Macrophage",
	TRUE ~ oc$celltype
)
Idents(oc)=oc$celltype
source("./1.R")
library(RColorBrewer)
library(viridis)
library(wesanderson)
library(cowplot)
n=30
qual_col_pals=brewer.pal.info[brewer.pal.info$category == "qual",]
col_vector=unlist(mapply(brewer.pal,qual_col_pals$maxcolors,rownames(qual_col_pals)))
pie(rep(1,n),col=sample(col_vector,n))
color=grDevices::colors()[grep('gr(a|e)y',grDevices::colors(),invert=T)]
pie(rep(6,n),col=sample(color,n))
col_vector

col_vector=c(wes_palette("Darjeeling1"),wes_palette("GrandBudapest1"),wes_palette("Cavalcanti1"),wes_palette("GrandBudapest2"),wes_palette("FantasticFox1"))
pal=wes_palette("Zissou1",10,type="continuous")
pal2=wes_palette("Zissou1",5,type="continuous")
pal[3:10]
pdf(file="scoredotplot.pdf", width=6, height=10)
colnames(oc@meta.data)
DotPlotScores(object = oc, scores = c("GLYCOLYSIS_GLUCONEOGENESIS","HIF_1_SIGNALING_PATHWAY",
																			"LACTATE_DEHYDROGENASE_ACTIVITY","LACTATE_METABOLIC_PROCESS"),dot.scale = 6)+
	RotatedAxis()+
	theme(axis.text.x=element_text(angle=90,face="italic",hjust=1),axis.text.y=element_text(face="bold"))+
	scale_colour_gradientn(colours=pal)+theme(legend.position="right")+labs(title="Lactylation",y="",x="")
dev.off()

library(ggplot2)
library(dplyr)
library(ggrepel)
library(viridis)
df <- data.frame(oc@meta.data, oc@reductions[["umap"]]@cell.embeddings)

class_avg <- df %>%
	group_by(celltype) %>%
	summarise(umap_1 = median(umap_1),
						umap_2 = median(umap_2))

score_methods <- c("GLYCOLYSIS_GLUCONEOGENESIS","HIF_1_SIGNALING_PATHWAY",
									 "LACTATE_DEHYDROGENASE_ACTIVITY","LACTATE_METABOLIC_PROCESS")

for (method in score_methods) {
	pdf(file = paste0(method, "_umap.pdf"), width = 7, height = 6)

	p <- ggplot(df, aes(x = umap_1, y = umap_2)) +
		geom_point(aes_string(color = method), size = 0.5, alpha = 0.8) +
		ggrepel::geom_label_repel(
			aes(label = celltype),
			data = class_avg,
			size = 3,
			box.padding = 0.35,
			segment.color = NA,
			label.size = 0.25,
			color = "black",
			fill = "white"
		) +
		scale_color_viridis(
			option = "I",
			name = "Score",
			guide = guide_colorbar(
				barwidth = 0.8,
				barheight = 5,
				title.position = "top",
				title.hjust = 0.5
			)
		) +
		theme_classic() +
		labs(x = "UMAP1", y = "UMAP2", title = paste(method, "Score")) +
		theme(
			axis.title = element_text(size = 12, face = "bold"),
			axis.text = element_text(color = "black"),
			plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
			panel.border = element_rect(color = "#256a53", size = 1.5, fill = NA),
			panel.background = element_rect(fill = "#f5f5f5"),
			legend.position = "right"
		)

	print(p)
	dev.off()
}

library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(tidyverse)
output_dir <- "MetaAnalysis_Results"
target_columns <- c("GLYCOLYSIS_GLUCONEOGENESIS","HIF_1_SIGNALING_PATHWAY",
										"LACTATE_DEHYDROGENASE_ACTIVITY","LACTATE_METABOLIC_PROCESS")
type_groups <- c("AA", "ANA")
color_scheme <- c("#FF7F0E", "#1F77B4")
comparison_label <- "Macrophage Subtypes"

plot_data <- oc@meta.data %>%
	rownames_to_column("CellID") %>%
	filter(phenotype %in% type_groups) %>%
	mutate(Group = factor(phenotype, levels = type_groups)) %>%
	dplyr::select(CellID, Group, all_of(target_columns)) %>%
	pivot_longer(
		cols = all_of(target_columns),
		names_to = "Pathway",
		values_to = "Value"
	) %>%
	mutate(Pathway = case_when(
		Pathway == "GLYCOLYSIS_GLUCONEOGENESIS" ~ "GLYCOLYSIS_GLUCONEOGENESIS",
		Pathway == "HIF_1_SIGNALING_PATHWAY" ~ "HIF_1_SIGNALING_PATHWAY",
		Pathway == "LACTATE_DEHYDROGENASE_ACTIVITY" ~ "LACTATE_DEHYDROGENASE_ACTIVITY",
		Pathway == "LACTATE_METABOLIC_PROCESS" ~ "LACTATE_METABOLIC_PROCESS",
		TRUE ~ Pathway
	)) %>%
	mutate(Pathway = factor(Pathway,
													levels = c("GLYCOLYSIS_GLUCONEOGENESIS", "HIF_1_SIGNALING_PATHWAY",
																		 "LACTATE_DEHYDROGENASE_ACTIVITY", "LACTATE_METABOLIC_PROCESS")))

stat_data <- plot_data %>%
	group_by(Pathway) %>%
	summarise(
		p_value = wilcox.test(Value ~ Group)$p.value,
		y_max = max(Value) * 1.1,
		.groups = "drop"
	) %>%
	mutate(
		p_label = ifelse(p_value < 0.001,
										 "***P < 0.001",
										 paste0("*P = ", signif(p_value, 3)))
	)

mean_sd <- plot_data %>%
	group_by(Pathway, Group) %>%
	summarise(
		mean = mean(Value, na.rm = TRUE),
		sd = sd(Value, na.rm = TRUE),
		y_pos = quantile(Value, 0.95),
		.groups = "drop"
	) %>%
	mutate(
		label = sprintf("Mean: %.2f\nSD: %.2f", mean, sd))

ggplot(plot_data, aes(x = Pathway, y = Value)) +
	geom_violin(aes(fill = Group),
							position = position_dodge(0.8),
							width = 0.7,
							alpha = 0.5,
							trim = FALSE) +
	geom_boxplot(aes(group = interaction(Pathway, Group)),
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
						aes(x = Pathway, y = y_max, label = p_label),
						size = 4.5,
						vjust = -0.5,
						inherit.aes = FALSE) +
	geom_text(data = mean_sd,
						aes(x = Pathway, y = y_pos,
								label = label, group = Group),
						position = position_dodge(0.8),
						size = 3.5,
						color = "black",
						vjust = 0) +
	scale_fill_manual(values = color_scheme) +
	scale_color_manual(values = color_scheme) +
	labs(title = "Pathway Activity Comparison between Macrophage Subtypes",
			 y = "Pathway Score",
			 x = "Signaling Pathways") +
	theme_classic(base_size = 14) +
	theme(
		plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
		axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
		axis.text.y = element_text(size = 11),
		legend.position = "top",
		legend.title = element_blank(),
		panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
		plot.margin = margin(1, 1, 1, 1, "cm")
	) +
	scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))

dir.create(output_dir, showWarnings = FALSE)
ggsave(filename = file.path(output_dir, "Combined_Pathway_Comparison_with_Stats.pdf"),
			 width = 10,
			 height = 7,
			 dpi = 300)
plot_data <- oc@meta.data %>%
	rownames_to_column("CellID") %>%
	filter(phenotype %in% type_groups) %>%
	mutate(Group = factor(phenotype, levels = type_groups)) %>%
	dplyr::select(CellID, Group, celltype, Scoreing)

stat_data <- plot_data %>%
	group_by(celltype) %>%
	summarise(
		p_value = wilcox.test(Scoreing ~ Group)$p.value,
		y_max = max(Scoreing) * 1.1,
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
		mean = mean(Scoreing, na.rm = TRUE),
		sd = sd(Scoreing, na.rm = TRUE),
		y_pos = quantile(Scoreing, 0.95),
		.groups = "drop"
	) %>%
	mutate(label = sprintf("Mean: %.2f\nSD: %.2f", mean, sd))

ggplot(plot_data, aes(x = celltype, y = Scoreing)) +
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
						size = 3.5,
						color = "black",
						vjust = 0) +
	scale_fill_manual(values = color_scheme) +
	scale_color_manual(values = color_scheme) +
	labs(title = "Scoreing Comparison Across Cell Types",
			 y = "Pathway Combined Score",
			 x = "Cell Types") +
	theme_classic(base_size = 14) +
	theme(
		plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
		axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
		axis.text.y = element_text(size = 11),
		legend.position = "top",
		legend.title = element_blank(),
		panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
		plot.margin = margin(1, 1, 1, 1, "cm")
	) +
	scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))

ggsave(filename = file.path(output_dir, "CellType_Scoreing_Comparison.pdf"),
			 width = 14,
			 height = 7,
			 dpi = 300)
save(oc,file = "lactylation.Rdata")
pathways <- c("JAK_STAT", "MARK", "NF_Kappa_B", "NOD_like_receptor")
cell_type <- "MNDA+ Macrophage"
output_dir <- "Correlation_Analysis"

dir.create(output_dir, showWarnings = FALSE)

results_df <- data.frame(
	Pathway = character(),
	Correlation = numeric(),
	P_value = numeric(),
	N = integer(),
	stringsAsFactors = FALSE
)

for (pathway in pathways) {
	target_cells <- rownames(oc@meta.data)[oc@meta.data$celltype == cell_type]

	mnda_exp <- GetAssayData(oc, assay = "RNA", slot = "counts")["MNDA", target_cells]
	pathway_score <- oc@meta.data[target_cells, pathway]

	cor_res <- cor.test(
		x = as.numeric(mnda_exp),
		y = as.numeric(pathway_score),
		method = "pearson"
	)
	cor_rho <- round(cor_res$estimate, 2)
	results_df <- rbind(results_df, data.frame(
		Pathway = pathway,
		Correlation = round(cor_res$estimate, 3),
		P_value = signif(cor_res$p.value, 3),
		N = length(target_cells)
	))

	plot_title <- paste("MNDA vs", pathway, "in", cell_type)
	output_file <- file.path(output_dir, paste0(pathway, "_correlation.pdf"))

	pdf(output_file, width = 6, height = 6)

	par(bty = "o",
			mgp = c(2, 0.5, 0),
			mar = c(4.1, 4.1, 2.1, 4.1),
			tcl = -.25,
			font.main = 3)

	plot(NULL, NULL,
			 xlim = range(pathway_score),
			 ylim = range(mnda_exp),
			 xlab = paste(pathway, "Score"),
			 ylab = "MNDA Expression",
			 main = plot_title,
			 col = "white")

	rect(par("usr")[1], par("usr")[3],
			 par("usr")[2], par("usr")[4],
			 col = "#EAE9E9", border = FALSE)

	grid(col = "white", lty = 1, lwd = 1.5)

	points(x = pathway_score,
				 y = as.numeric(mnda_exp),
				 pch = 19,
				 col = scales::alpha("#E51718", 0.8),
				 cex = 1.5)

	abline(lm(as.numeric(mnda_exp) ~ pathway_score),
				 lwd = 2,
				 col = "black")

	rug(pathway_score, side = 3, col = "black", lwd = 1)
	rug(mnda_exp, side = 4, col = "black", lwd = 1)

	text(x = max(pathway_score)*0.95,
			 y = max(mnda_exp)*0.95,
			 adj = 1,
			 labels = bquote("N = " ~ .(n) ~
			 									"; " ~ rho ~ " = " ~ .(cor_rho) ~
			 									"; " ~ italic(P) ~ " = " ~ .(cor_p)),
			 col = "black",
			 cex = 1)
	# ============================================

	dev.off()
}

write.csv(results_df,
					file = file.path(output_dir, "correlation_results.csv"),
					row.names = FALSE)

message("\n分析完成！结果保存至：", normalizePath(output_dir))
cat("相关性分析结果汇总：\n")
print(results_df)

pathways <- c("JAK_STAT", "MARK", "NF_Kappa_B", "NOD_like_receptor")
cell_type <- "MNDA+ Macrophage"
output_dir <- "Correlation_Analysis_lactylation"

dir.create(output_dir, showWarnings = FALSE)

results_df <- data.frame(
	Pathway = character(),
	Correlation = numeric(),
	P_value = numeric(),
	N = integer(),
	stringsAsFactors = FALSE
)

for (pathway in pathways) {
	target_cells <- rownames(oc@meta.data)[oc@meta.data$celltype == cell_type]

	mnda_exp <- oc@meta.data[target_cells, "All"]
	pathway_score <- oc@meta.data[target_cells, pathway]

	cor_res <- cor.test(
		x = as.numeric(mnda_exp),
		y = as.numeric(pathway_score),
		method = "pearson"
	)
	cor_rho <- round(cor_res$estimate, 2)
	results_df <- rbind(results_df, data.frame(
		Pathway = pathway,
		Correlation = round(cor_res$estimate, 3),
		P_value = signif(cor_res$p.value, 3),
		N = length(target_cells)
	))

	plot_title <- paste("Lactylation vs", pathway, "in", cell_type)
	output_file <- file.path(output_dir, paste0(pathway, "_correlation.pdf"))

	pdf(output_file, width = 6, height = 6)

	par(bty = "o",
			mgp = c(2, 0.5, 0),
			mar = c(4.1, 4.1, 2.1, 4.1),
			tcl = -.25,
			font.main = 3)

	plot(NULL, NULL,
			 xlim = range(pathway_score),
			 ylim = range(mnda_exp),
			 xlab = paste(pathway, "Score"),
			 ylab = "Lactylation Score",
			 main = plot_title,
			 col = "white")

	rect(par("usr")[1], par("usr")[3],
			 par("usr")[2], par("usr")[4],
			 col = "#EAE9E9", border = FALSE)

	grid(col = "white", lty = 1, lwd = 1.5)

	points(x = pathway_score,
				 y = as.numeric(mnda_exp),
				 pch = 19,
				 col = scales::alpha("#E51718", 0.8),
				 cex = 1.5)

	abline(lm(as.numeric(mnda_exp) ~ pathway_score),
				 lwd = 2,
				 col = "black")

	rug(pathway_score, side = 3, col = "black", lwd = 1)
	rug(mnda_exp, side = 4, col = "black", lwd = 1)

	text(x = max(pathway_score)*0.95,
			 y = max(mnda_exp)*0.95,
			 adj = 1,
			 labels = bquote("N = " ~ .(n) ~
			 									"; " ~ rho ~ " = " ~ .(cor_rho) ~
			 									"; " ~ italic(P) ~ " = " ~ .(cor_p)),
			 col = "black",
			 cex = 1)
	# ============================================

	dev.off()
}

write.csv(results_df,
					file = file.path(output_dir, "correlation_results.csv"),
					row.names = FALSE)

message("\n分析完成！结果保存至：", normalizePath(output_dir))
cat("相关性分析结果汇总：\n")
print(results_df)
