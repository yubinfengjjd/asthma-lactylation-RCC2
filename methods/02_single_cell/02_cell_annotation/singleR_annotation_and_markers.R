options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(data.table)
library(paletteer)
library(vcd)
library(gplots)
library(limma)
library(stringr)
library(SingleR)
library(monocle)
library(cowplot)
library(patchwork)
library(metap)
library(dplyr)
library(RColorBrewer)
library(harmony)
library(paletteer)
logFCfilter=1
adjPvalFilter=0.05
setwd("D:\\哮喘乳酸化单细胞开篇\\3.单细胞注释")
DefaultAssay(oc) <- "RNA"
Idents(oc) <- "RNA_snn_res.0.2"
counts<-oc@assays$RNA$counts
clusters<-oc@meta.data$RNA_snn_res.0.2
ref=get(load("ref_Human_all.RData"))
ref=ref_Human_all
singler=SingleR(test=counts, ref =ref,
								labels=ref$label.main, clusters = clusters)
singler2=SingleR(test=counts, ref =ref,
								 labels=ref$label.main)
cellAnn=as.data.frame(singler2)
cellAnn=cbind(id=row.names(cellAnn), cellAnn)
cellAnn=cellAnn[,c("id", "labels")]
write.table(cellAnn, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)
newLabels=singler$labels
names(newLabels)=levels(oc)
oc=RenameIdents(oc, newLabels)
pdf(file="singerR.UMAP.pdf",width=6.5,height=6)
DimPlot(oc, reduction = "umap", label = T)
dev.off()
colnames(oc@meta.data)
pdf(file="UMAP_auto.pdf",width=6.5,height=6)
plotScoreHeatmap(singler)
dev.off()

DefaultAssay(oc) <- "RNA"
Idents(oc) <- "RNA_snn_res.0.8"
table(Idents(oc))
oc.markers=FindAllMarkers(object = oc,
													only.pos = FALSE,
													min.pct = 0.25,
													shold = logFCfilter,assay = "RNA")
sig.cellMarkers=oc.markers[as.numeric(as.vector(oc.markers$avg_log2FC))>logFCfilter & as.numeric(as.vector(oc.markers$p_val_adj))<adjPvalFilter,]
write.table(sig.cellMarkers,file="cellMarkers.txt",sep="\t",row.names=F,quote=F)
top10 <- sig.cellMarkers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC )
write.table(top10,file="cellMarkers_top10.txt",sep="\t",row.names=F,quote=F)
small_marker_dict <- list(
	Fibroblast = c("COL1A1","PDGFRA"),
	ECs = c("CD31", "PECAM1","VWF"),
	EpithelialCells = c("EPCAM", "CDH1"),
	Mast = c("KIT", "TPSAB1","TPSB2"),
	BCells = c("MS4A1", "CD79A","CD19","CD79"),
	Plasma = c("MZB1", "IGHA1", "JCHAIN"),
	TCells = c("CD3E","CD3D"),
	DC = c("CLEC10A", "CD1C", "CLEC4C", "PTCRA", "CCR7", "LAMP3"),
	NK=c("XCL1", "FCGR3A", "KLRD1", "KLRF1"),
	Macrophage=c("LGMN", "CTSB", "CD14", "FCGR3A"),
	Monocyte=c("CTSS", "FCN1", "S100A8", "S100A9", "LYZ", "VCAN")
)

all_markers <- unlist(small_marker_dict)


filtered_markers <- sig.cellMarkers %>%
	filter(gene %in% all_markers)

find_cell_type <- function(gene) {
	cell_types <- names(small_marker_dict)[sapply(small_marker_dict, function(markers) gene %in% markers)]
	if (length(cell_types) > 0) {
		return(cell_types[1])
	} else {
		return("Unknown")
	}
}

filtered_markers$cell_type <- sapply(filtered_markers$gene, find_cell_type)

pdf(file = ("gene.pdf"), width = 6, height = 6)
all_markers <- unlist(markers_plot)
FeaturePlot(oc, features = all_markers, max.cutoff = 3, cols = c("grey", "red"),by.col = 2,slot = "data",combine = F,label = T)

dev.off()

Idents(oc)=oc$RNA_snn_res.0.8
genes <- list("T cells" = c("CD3D", "CD3E"),
							"B cells" = c("MS4A1","CD79A"),
							"NK cells" = c("NKG7", "GNLY"),
							"Mast cells" = c("CPA3","SAMSN1"),
							"DCs" = c("FLT3", "HLA-DRA"),
							"Mono/Mac" = c("S100A8", "S100A9","CD68","CD163"),
							"AEC" = c("EPCAM", "KRT19"))

pdf(file="annotation.pdf",width=12,height=8)
p <- SCpubr::do_DotPlot(sample = oc,
												features = genes,
												dot.scale = 6)
print(p)
dev.off()
celltype_map <- c(
	"0" = "T cells",
	"1" = "T cells",
	"2" = "Mono/Mac",
	"3" = "Mono/Mac",
	"4" = "DCs",
	"5" = "AEC",
	"6" = "B cells",
	"7" = "DCs",
	"8" = "DCs",
	"9" = "T cells",
	"10" = "T cells",
	"11" = "AEC",
	"12" = "NK cells",
	"13" = "Mast cells",
	"14" = "AEC",
	"15" = "Mast cells",
	"16" = "AEC",
	"17" = "T cells",
	"18" = "Mono/Mac",
	"19" = "DCs",
	"20" = "AEC"
)

oc$celltype <- factor(
	as.character(oc$RNA_snn_res.0.8),
	levels = names(celltype_map),
	labels = celltype_map
)

Idents(oc) <- "celltype"
pdf(file="umap_annotation.pdf",width=6,height=6)
SCpubr::do_DimPlot(sample = oc,
									 plot.axes = F,label = TRUE,label.size = 2)
dev.off()
save(oc,file = "oc_exegesis.Rdata")

library(ClusterGVis)
library(org.Hs.eg.db)
library(Seurat)
library(dplyr)
pbmc.markers.all=FindAllMarkers(oc,
																only.pos=T,
																min.pct = 0.25,
																logfc.threshold=0.25)
pbmc.markers=pbmc.markers.all[as.numeric(as.vector(pbmc.markers.all$avg_log2FC))>logFCfilter & as.numeric(as.vector(pbmc.markers.all$p_val_adj))<adjPvalFilter,]
pbmc.markers=pbmc.markers %>%
	group_by(cluster) %>%
	top_n(n=20,wt=avg_log2FC)
st.data=prepareDataFromscRNA(
	oc,
	diffData=pbmc.markers,
	showAverage=T
)
str(st.data)
enrich=enrichCluster(
	st.data,
	OrgDb=org.Hs.eg.db,
	type="BP",
	organism="hsa",
	pvalueCutoff=0.5,
	topn=5,
	seed=5201314
)

markGenes=unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)),40,replace = T)]
visCluster(object = st.data,
					 plot.type = "line")
pdf('sc1.pdf',height = 16,width = 14,onefile = F)
visCluster(object=st.data,
					 plot.type = "both",
					 column_names_rot=35,
					 show_row_dend=F,
					 markGenes=markGenes,
					 markGenes.side="left",
					 annoTerm.data=enrich,
					 line.side="left",
					 cluster.order=c(1:7),
					 go.col=rep(jjAnno::useMyCol("stallion",n=7),each=5),
					 add.bar = T)
dev.off()

library(Seurat)
library(Nebulosa)
library(ggnetwork)
library(dplyr)
library(ggunchull)
markers_plot <- list(
	AT1=c("AGER"),
	AT2=c("LAMP3","SFTPC"),
	EpithelialCells = c("EPCAM", "CDH1"),
	Mast = c("KIT","TPSAB1","TPSB2"),
	BCells = c("MS4A1", "CD79A","CD19"),
	Plasma = c("IGHG1","MZB1"),
	TCells = c("CD3E","CD3D"),
	DCs = c("CD83","CD1C","CLEC10A", "CLEC4C", "PTCRA", "CCR7"),
	`Alveolar-Mac`=c("MARCO", "FABP4", "MCEMP1"),
	Monocyte=c("S100A8", "S100A9"),
	Granulocyte=c("CLC"),
	NK=c("NKG7"),
	pDCs=c("IL3RA","CLEC4C"),
	Macrophage=c("LGMN", "CTSB", "CD14", "FCGR3A")
)

plist <- list()
for (i in 1:length(markers_plot)) {
	y_pos <- 14

	text_data <- data.frame(
		x = 5,
		y = y_pos,
		label = ifelse(length(markers_plot[[i]]) > 1,
									 paste(markers_plot[[i]], collapse = " + "),
									 markers_plot[[i]])
	)

	if (length(markers_plot[[i]]) == 1) {
		result <- plot_density(oc, features = markers_plot[[i]], pal = 'magma',
													 raster = FALSE, size = 0.8, joint = F)+
			theme_blank()
		plist[[i]] <- result
	} else {
		result <- plot_density(oc, features = markers_plot[[i]], pal = 'magma',
													 raster = FALSE, size = 0.8, joint = T)

		p <- result[[3]] +
			theme_blank() +
			theme(
				legend.frame = element_rect(colour = "black"),
				legend.ticks = element_line(colour = "black", linewidth = 0),
				legend.key.width = unit(0.3, "cm"),
				legend.key.height = unit(0.6, "cm"),
				legend.title = element_text(color = 'black', face = "bold", size = 8),
				plot.title = element_blank()
			) +
			geom_text(data = text_data, aes(x = x, y = y, label = label),
								vjust = 1.5, size = 5, color = 'black', angle = 0)

		plist[[i]] <- p
	}
}


library(cowplot)
combined_plot = plot_grid(plotlist = plist, ncol = 4)
ggsave('check_paper_markers_cluster.pdf', plot = combined_plot, height = 12, width = 20)

Idents(oc) <- "celltype"
pdf(file="dot_annotation.pdf",width=12,height=8)
p <- SCpubr::do_DotPlot(sample = oc,
												features = genes,
												dot.scale = 8)
print(p)
dev.off()
sample_table <- as.data.frame(
	table(
		oc@meta.data$celltype,
		oc@meta.data$orig.ident
	)
)

names(sample_table) <- c("celltype", "group", "CellNumber")

colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C")

pdf(file="celltypebar.pdf",width=6,height=6)
p3 <- ggplot(
	sample_table,
	aes(x = group, weight = CellNumber, fill = celltype)
) +
	geom_bar(
		position = "fill",
		width = 0.7,
		size = 0.5,
		colour = '#222222'
	) +
	scale_fill_manual(values = colors) +
	theme(
		panel.grid = element_blank(),
		panel.background = element_rect(fill = "transparent", colour = NA),
		axis.line.x = element_line(colour = "black"),
		axis.line.y = element_line(colour = "black"),
		plot.title = element_text(lineheight = .8, face = "bold", hjust = 0.5, size = 16)
	) +
	labs(y = "Percentage") +
	coord_flip()

print(p3)
dev.off()

library(ggpubr)
pdf(file="geneVlnplot.pdf", width=6, height=8)
Data = data.frame('cell_type'=oc$celltype,
									'Score'=oc@assays$RNA$data["RBM25",],
									"group"=oc$phenotype)
names(Data)=c("cell_type","Score","group")

Data_summary = Data %>%
	group_by(cell_type) %>%
	summarise(across(.cols = Score,.fns = list(
		'mean' = mean, 'sd' =sd
	))) %>%
	dplyr::rename('mean' = Score_mean,
								'sd' = Score_sd)
Data = left_join(Data,Data_summary,by="cell_type")

pdf(file="RBM25_expression.PDF", width=10, height=6)
p=ggplot(Data, aes(x=cell_type, y=Score,fill=group)) +
	geom_violin(trim=FALSE,color="white") +
	geom_boxplot(width=0.2,position=position_dodge(0.9))+
	geom_point(aes(x=cell_type, y=mean),pch=19,position=position_dodge(0.9),size=1.5)+
	geom_errorbar(aes(ymin = mean-sd, ymax=mean+sd),
								width=0.1,
								position=position_dodge(0.9),
								color="black",
								alpha = 0.7,
								size=0.5) +
	scale_fill_manual(values = ggsci::pal_aaas()(2))+
	theme_bw()+
	theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=12),
				axis.text.y=element_text(family="Times",size=16,face="plain"),
				axis.title.y=element_text(family="Times",size = 5,face="plain"),
				panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
				legend.text=element_text(face="italic", family="Times", colour="black",
																 size=8),
				legend.title=element_text(face="italic", family="Times", colour="black",
																	size=18),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank())+
	ylab("")+xlab("")+
	ggtitle("RBM25")+
	stat_compare_means(aes(group=oc@meta.data$phenotype),
										 method="wilcox.test",
										 symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
										 label = "p.signif")
print(p)
dev.off()
oc[["RNA"]] <- as(object = oc[["RNA"]], Class = "Assay")
library(SCP)
pdf(file="RBM25_DIM.pdf", width=6, height=6)
FeatureDimPlot(oc,
							 features = "RBM25", reduction = "umap", label = TRUE,
							 cells.highlight = colnames(oc)[oc$celltype == "Mono/Mac"],
							 pt.size = 0.1,split.by = "phenotype"
)
dev.off()

DEGs <- oc@tools$DEtest_celltype$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]
VolcanoPlot(srt = RA, group_by = "celltype",ncol = 4)
# Annotate features with transcription factors and surface proteins
oc <- AnnotateFeatures(oc, species = "Homo_sapiens", db = c("TF", "CSPA"))
ht <- FeatureHeatmap(
	srt = oc, group.by = "celltype", features = DEGs$gene, feature_split = DEGs$group1,
	species = "Homo_sapiens", db = c("GO_BP", "KEGG"), anno_terms = TRUE,
	feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
	height = 5, width = 4
)

library(SCP)
oc[["RNA"]] <- as(object = oc[["RNA"]], Class = "Assay")
DefaultAssay(oc) <- "RNA"
p=FeatureStatPlot(
	oc,
	stat.by = unlist(genes),
	legend.position = "top",
	legend.direction = "horizontal",
	bg.by = "celltype",group.by = "celltype",
	stack = TRUE,
)
print(p)
