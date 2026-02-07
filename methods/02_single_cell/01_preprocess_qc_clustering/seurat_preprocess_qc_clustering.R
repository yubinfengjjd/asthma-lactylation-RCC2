library(SeuratDisk)
library(patchwork)
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
Convert('./GSE193816_all_data_raw_counts.h5ad/GSE193816_all_data_raw_counts.h5ad', "h5seurat",overwrite = TRUE,assay = "RNA")
scRNA=LoadH5Seurat("./GSE193816_all_data_raw_counts.h5ad/GSE193816_all_data_raw_counts.h5seurat", misc = F)
oc <- CreateSeuratObject(counts = GetAssayData(scRNA))
oc@meta.data=oc@meta.data[match(rownames(oc@meta.data),rownames(scRNA@meta.data)),]
oc@meta.data=cbind(oc@meta.data,scRNA@meta.data)
oc=subset(oc,subset = sample%in%c("Ag"))
oc[["percent.mt"]] <- PercentageFeatureSet(oc, pattern = "^MT-")
oc[["percent.rp"]] <- PercentageFeatureSet(oc, pattern = "^RP[SL]")
oc[["percent.hb"]] <- PercentageFeatureSet(oc, pattern = "^HB[^(P)]")
oc <- subset(oc,
						 nFeature_RNA >= 300 & nFeature_RNA <= 6000 &
						 	percent.mt >= 0 & percent.mt <= 15 &
						 	percent.hb >= 0 & percent.hb <= 0.1 &
						 	percent.rp >= 1 & percent.rp <= 100
)
save(oc,file = "oc.Rdata")
vln_plots <- VlnPlot(oc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
										 ncol=3, pt.size=0, group.by="orig.ident", combine = FALSE)

vln_plots <- lapply(vln_plots, function(p) {
	p + theme(axis.text.x = element_text(size = 8),
						axis.text.y = element_text(size = 8),
						axis.title.x = element_text(size = 10),
						axis.title.y = element_text(size = 10),
						legend.position = "none",
						strip.text = element_text(size = 8))
})

combined_plot <- plot_grid(plotlist = vln_plots, ncol = 3)

ggsave("VlnPlot.pdf", plot = combined_plot, width = 12, height = 6)

DefaultAssay(oc) <- "RNA"
oc <- oc %>% NormalizeData(verbose=F) %>%
	FindVariableFeatures(selection.method = "vst" ,nfeatures = 2000 ,verbose = F) %>%
	ScaleData(verbose = F,features = VariableFeatures(oc)) %>%
	RunPCA(npcs=50,verbose = F,features = VariableFeatures(oc))
ElbowPlot(oc,ndims = 50)
p1=DimPlot(oc, reduction = "pca",group.by = "orig.ident")
ggsave("pca点图.pdf", plot = p1, width = 6, height = 6)
colnames(oc@meta.data)
oc <- RunHarmony(oc,
								 group.by.vars="phenotype")
oc <- RunUMAP(oc, reduction = "harmony", dims = 1:30)
oc <- FindNeighbors(oc, reduction = "harmony", dims = 1:30)

for (res in c(0.05, 0.1, 0.2, 0.3, 0.5,0.8,1,1.2)) {
	print(res)
	oc <- FindClusters(oc, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
cluster_umap <- wrap_plots(ncol = 4,
													 DimPlot(oc, reduction = "umap", group.by = "RNA_snn_res.0.05", label = T) & NoAxes(),
													 DimPlot(oc, reduction = "umap", group.by = "RNA_snn_res.0.1", label = T) & NoAxes(),
													 DimPlot(oc, reduction = "umap", group.by = "RNA_snn_res.0.2", label = T)& NoAxes(),
													 DimPlot(oc, reduction = "umap", group.by = "RNA_snn_res.0.3", label = T)& NoAxes(),
													 DimPlot(oc, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T) & NoAxes(),
													 DimPlot(oc, reduction = "umap", group.by = "RNA_snn_res.0.8", label = T) & NoAxes(),
													 DimPlot(oc, reduction = "umap", group.by = "RNA_snn_res.1", label = T) & NoAxes(),
													 DimPlot(oc, reduction = "umap", group.by = "RNA_snn_res.1.2", label = T) & NoAxes()
)
#cluster_umap
ggsave(cluster_umap,filename = "Step4.After_inter.cluster_umap.pdf",
			 width = 20, height = 16)
# Make plot
library(clustree)
p1=clustree(oc@meta.data,prefix="RNA_snn_res.")

ggsave(p1,filename = "聚类树.pdf",
			 width = 12, height = 8)

Idents(object = oc) <- "RNA_snn_res.0.8"

library(SCpubr)
pdf(file="umap.pdf",width=6,height=6)
SCpubr::do_DimPlot(sample = oc,
									 plot.axes = F,label = TRUE)
dev.off()
save(oc,file = "oc_umap.Rdata")
