library(Seurat)
library(ggplot2)
library(dplyr)

setwd("D:\\哮喘乳酸化单细胞开篇\\6.2.Degs与cor")

load(file.path("D:\\哮喘乳酸化单细胞开篇\\6.一致性聚类", "Cluster.Rdata"))
hdwgcna_genes <- data.table::fread("./input/hdWGCNAgene.txt", header = FALSE, data.table = FALSE)[,1]

# Set default assay
DefaultAssay(oc_subset) <- "RNA"

# Define analysis function (corrected)
analyze_gene_set <- function(seurat_obj, genes, score_var = "Scoreing", phenotype_var = "phenotype",
														 p_threshold = 0.05, cor_threshold = 0.2) {
	# Extract expression matrix with drop=FALSE to preserve matrix structure
	expr <- GetAssayData(seurat_obj, slot = "data")[genes, , drop = FALSE]

	# Get scores and phenotype
	scores <- FetchData(seurat_obj, vars = score_var)[,1]
	phenotype <- FetchData(seurat_obj, vars = phenotype_var)[,1]
	phenotype_binary <- ifelse(phenotype == "AA", 1, 0)

	# Initialize results dataframe
	results <- data.frame(
		Gene = rownames(expr),
		Correlation_scores = numeric(nrow(expr)),
		Correlation_phenotype = numeric(nrow(expr)),
		P_Value_scores = numeric(nrow(expr)),
		P_Value_phenotype = numeric(nrow(expr)),
		stringsAsFactors = FALSE
	)

	# Calculate correlations per gene
	for (i in 1:nrow(expr)) {
		gene_expr <- as.numeric(expr[i, ])

		# Skip constant genes
		if (sd(gene_expr) == 0) {
			results[i, -1] <- NA
			next
		}

		# Calculate correlation with scores
		cor_scores <- suppressWarnings(cor.test(gene_expr, scores, method = "spearman", exact = FALSE))
		results$Correlation_scores[i] <- cor_scores$estimate
		results$P_Value_scores[i] <- cor_scores$p.value

		# Calculate correlation with phenotype
		cor_pheno <- suppressWarnings(cor.test(gene_expr, phenotype_binary, method = "spearman", exact = FALSE))
		results$Correlation_phenotype[i] <- cor_pheno$estimate
		results$P_Value_phenotype[i] <- cor_pheno$p.value
	}

	# Filter significant genes
	sig_genes <- results %>%
		filter(P_Value_scores < p_threshold & abs(Correlation_scores) > cor_threshold) %>%
		arrange(desc(abs(Correlation_scores)))

	return(list(results = results, sig_genes = sig_genes))
}

# Define visualization function (corrected)
plot_correlation <- function(results_df, title, filename) {
	# Remove NA values
	results_df <- na.omit(results_df)

	# Add grouping based on significance
	results_df <- results_df %>%
		mutate(group = ifelse(
			abs(Correlation_scores) > 0.2 & P_Value_scores < 0.05,
			"High",
			"Low"
		))

	# Calculate overall Spearman correlation
	cor_test <- suppressWarnings(
		cor.test(
			abs(results_df$Correlation_scores),
			abs(results_df$Correlation_phenotype),
			method = "spearman",
			exact = FALSE
		)
	)

	# Format label
	cor_label <- sprintf(
		"Spearman's rho = %.2f\np = %s",
		cor_test$estimate,
		ifelse(cor_test$p.value < 0.001,
					 "< 0.001",
					 format(cor_test$p.value, digits = 2))
	)

	# Create plot
	p <- ggplot(results_df, aes(x = abs(Correlation_scores), y = abs(Correlation_phenotype), color = group)) +
		geom_point(size = 3.5, alpha = 0.8) +
		geom_smooth(
			method = "lm",
			color = "#5b82a0",
			linetype = "solid",
			size = 1.5,
			se = TRUE,
			fill = "#5a5676",
			alpha = 0.3
		) +
		scale_color_manual(values = c("#0c695e", "#b81b23")) +
		geom_vline(xintercept = 0.2, linetype = "dotted", color = "grey60", size = 1) +
		annotate("text",
						 x = min(abs(results_df$Correlation_scores)) + 0.1,
						 y = max(abs(results_df$Correlation_phenotype)) - 0.05,
						 label = cor_label, hjust = 0, color = "black", size = 5) +
		labs(
			title = title,
			x = "Correlation with Scores",
			y = "Correlation with Phenotype"
		) +
		theme_bw(base_size = 14) +
		theme(
			legend.position = "none",
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_rect(color = "black", size = 1.2),
			plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
			axis.title = element_text(size = 14),
			axis.text = element_text(size = 12)
		)

	# Save plot
	ggsave(paste0(filename, ".pdf"), p, width = 8, height = 7, dpi = 300)
	return(p)
}

# Run analysis and plotting
hdwgcna_results <- analyze_gene_set(oc_subset, hdwgcna_genes)
write.table(
	hdwgcna_results$sig_genes$Gene,
	"hdwgcna_sig_genes.txt",
	sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE
)
plot_correlation(hdwgcna_results$results, "hdWGCNA Genes Correlation", "hdWGCNA_correlation")

Idents(oc_subset) <- "subtypes"
mbc_markers <- FindMarkers(oc_subset, ident.1 = "MBC2", ident.2 = "MBC1")
mbc_genes <- rownames(subset(mbc_markers, p_val_adj < 0.05))

mbc_results <- analyze_gene_set(oc_subset, mbc_genes)
write.table(mbc_results$sig_genes$Gene, "MBDs_sig_genes.txt",
						sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
plot_correlation(mbc_results$results, "MBDs Genes Correlation", "MBDs_correlation")

Idents(oc_subset) <- "phenotype"
pheno_markers <- FindMarkers(oc_subset, ident.1 = "AA", ident.2 = "ANA")
pheno_genes <- rownames(subset(pheno_markers, p_val_adj < 0.05))

pheno_results <- analyze_gene_set(oc_subset, pheno_genes)
write.table(pheno_results$sig_genes$Gene, "Degs_sig_genes.txt",
						sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
plot_correlation(pheno_results$results, "DEGs Correlation", "DEGs_correlation")
library(VennDiagram)
outFile="intersectGenes.txt"
outPic="venn.pdf"
files=dir()
files=grep("txt$",files,value=T)
geneList=list()

for(i in 1:length(files)){
	inputFile=files[i]
	if(inputFile==outFile){next}
	rt=read.table(inputFile,header=F)
	geneNames=as.vector(rt[,1])
	geneNames=gsub("^ | $","",geneNames)
	uniqGene=unique(geneNames)
	header=unlist(strsplit(inputFile,"\\.|\\-"))
	geneList[[header[1]]]=uniqGene
	uniqLength=length(uniqGene)
	print(paste(header[1],uniqLength,sep=" "))
}


venn.plot <- venn.diagram(
	x = geneList,
	filename = NULL,
	output = TRUE,
	imagetype = "PDF",
	fill = c("#E64B3599", "#F39B7F99", "#0099B499", "#0099B488"), # 圈里填充的颜色
	alpha = 0.8,
	color = "black",
	lwd = 3,
	lty = 1,
	cat.cex = 1.5,
	cat.dist = c(0.2, 0.2, 0.1, 0.1),
	cex = 1.2,
	main.cex = 2,
	margin = 0.1
)

pdf("VNN.PDF",width = 6,height = 6)
grid.draw(venn.plot)
dev.off()

intersectGenes=Reduce(intersect,geneList)
intersectGenes_table = MBDs_genes[intersectGenes,]
intersectGenes_table$group = ifelse(intersectGenes_table$avg_log2FC > 0,"up","down")
diffSigOut <- rbind(id = colnames(intersectGenes_table), intersectGenes_table)
write.table(diffSigOut, file = paste0("MBDS", ".xls"), sep = "\t", quote = F, col.names = F)
write.table(intersectGenes, file = "intersectGenes.txt", sep = "\t", quote = F, col.names = T,row.names = F)
