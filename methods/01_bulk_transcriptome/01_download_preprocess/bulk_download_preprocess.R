library(GEOquery)
library(tidyverse)
library(limma)
library(org.Hs.eg.db)
library(patchwork)
library(WGCNA)
library(GSEABase)
library(randomcoloR)
library(AnnoProbe)
library(dplyr)
rm(list = ls())
setwd("D:\\哮喘乳酸化单细胞开篇\\7.转录组数据下载")
files <- list.files("./series文件/", pattern = "*.txt", full.names = TRUE)
names <- gsub("(.*)_(.*)_(.*)\\.gz", "\\1", basename(files))
gsets <- list()
pds <- list()
exprSets <- list()

for (i in 1:length(files)) {
	file_path <- files[i]

	if (!file.exists(file_path)) {
		warning(paste("文件不存在:", file_path))
		next
	}

	gset <- tryCatch({
		getGEO(filename = file_path, getGPL = F)
	}, error = function(e) {
		warning(paste("获取GEO数据集时出错:", file_path, "错误信息:", e$message))
		return(NULL)
	})

	if (is.null(gset)) next

	if (is.list(gset)) {
		gset <- exprs(gset)
	}

	gsets[[names[i]]] <- gset

	pdata <- tryCatch({
		pData(gset)
	}, error = function(e) {
		warning(paste("提取表型数据时出错:", names[i], "错误信息:", e$message))
		return(NULL)
	})

	if (is.null(pdata)) next

	pds[[names[i]]] <- pdata

	exprSet <- tryCatch({
		read.table(file_path,
							 comment.char = "!",
							 stringsAsFactors = FALSE,
							 header = TRUE)
	}, error = function(e) {
		warning(paste("读取表达数据时出错:", file_path, "错误信息:", e$message))
		return(NULL)
	})

	if (is.null(exprSet)) next

	rownames(exprSet) <- exprSet[,1]
	exprSet <- exprSet[,-1]

	ex <- exprSet
	qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
	LogC <- (qx[5] > 100) ||
		(qx[6] - qx[1] > 50 && qx[2] > 0) ||
		(qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

	if (LogC) {
		ex[ex <= 0] <- NaN
		exprSet <- log2(ex+1)
		print("log2 transform finished")
	} else {
		print("LogC is FALSE, no transformation performed")
	}
	exprSet=normalizeBetweenArrays(exprSet)
	exprSet <- as.data.frame(exprSet)
	exprSets[[names[i]]] <- exprSet
}

for (dataset_name in names) {
	assign(paste0("gset_", dataset_name), gsets[[dataset_name]])
	assign(paste0("pdata_", dataset_name), pds[[dataset_name]])
	assign(paste0("exprSet_", dataset_name), exprSets[[dataset_name]])
}

names=names[names%in%names(exprSets)]

for (dataset_name in names) {
	cat("数据集:", dataset_name, "\n")
	if (dataset_name %in% names(pds)) {
		cat("表型数据行数:", nrow(get(paste0("pdata_", dataset_name))), "\n")
	} else {
		cat("表型数据: 未找到\n")
	}
	cat("表达数据行数:", nrow(get(paste0("exprSet_", dataset_name))), ", 列数:", ncol(get(paste0("exprSet_", dataset_name))), "\n\n")
}
index_GSE43696 <- gset_GSE43696@annotation
index_GSE63142 = gset_GSE63142@annotation
index_GSE46171  = gset_GSE46171@annotation

anno_df=data.table::fread("./GPL6480-9577 (3).txt",data.table = F)
anno_df2=data.table::fread("./GPL570-55999.txt",data.table = F,skip = 16)
anno_df <- anno_df %>%
	mutate(GENE_SYMBOL = str_replace_all(GENE_SYMBOL, "---", "")) %>%
	filter(GENE_SYMBOL != "") %>%
	mutate(GENE_SYMBOL = str_replace_all(GENE_SYMBOL, " ", "")) %>%
	dplyr::select(ID, GENE_SYMBOL)
names(anno_df) <- c("probe_id", "symbol")
anno_df$probe_id=as.character(anno_df$probe_id)

anno_df2 <- anno_df2 %>%
	setNames(make.names(names(.), unique = TRUE)) %>%
	mutate(GENE_SYMBOL = str_replace_all(GENE_SYMBOL, "---", "")) %>%
	mutate(GENE_SYMBOL = str_replace_all(GENE_SYMBOL, " ", "")) %>%
	filter(GENE_SYMBOL != "" & !is.na(GENE_SYMBOL)) %>%
	dplyr::select(probe_id = ID, symbol = GENE_SYMBOL) %>%
	mutate(probe_id = as.character(probe_id)) %>%
	distinct(probe_id, symbol, .keep_all = TRUE)
anno_df2 <- anno_df2 %>%
	dplyr::select(ID, SPOT_ID)
biotype <- AnnotationDbi::select(
	org.Hs.eg.db,
	keys = anno_df2[,2],
	columns = c("SYMBOL", "GENETYPE"),
	keytype = "ENTREZID"
)
biotype <- biotype[!is.na(biotype$SYMBOL), ]
anno_df2 = anno_df2 %>% filter(anno_df2$SPOT_ID%in%biotype$ENTREZID)
anno_df2$SPOT_ID = biotype$SYMBOL[match(anno_df2$SPOT_ID,biotype$ENTREZID)]
names(anno_df2) <- c("probe_id", "symbol")
anno_df2$probe_id=as.character(anno_df2$probe_id)
#####\\\#####
anno_df2 <- anno_df2 %>%
	mutate(`Gene Symbol` = str_replace_all(`Gene Symbol`, "---", "")) %>%
	filter(`Gene Symbol` != "") %>%
	separate(`Gene Symbol`, into = c("Gene Symbol","drop"), sep = "///") %>%
	mutate(`Gene Symbol` = str_replace_all(`Gene Symbol`, " ", "")) %>%
	dplyr::select(ID, `Gene Symbol`)
names(anno_df2) <- c("probe_id", "symbol")
anno_df2$probe_id=as.character(anno_df2$probe_id)
#####GPL10558-50081####
anno_df2 <- anno_df2 %>%
	mutate(ILMN_Gene = str_replace_all(ILMN_Gene, "---", "")) %>%
	filter(ILMN_Gene != "") %>%
	mutate(ILMN_Gene = str_replace_all(ILMN_Gene, " ", "")) %>%
	dplyr::select(ID, ILMN_Gene)
names(anno_df2) <- c("probe_id", "symbol")
anno_df2$probe_id=as.character(anno_df2$probe_id)
######exprSet_GSE43696####
exprSet_GSE43696 <- exprSet_GSE43696 %>%
	rownames_to_column("probe_id") %>%
	inner_join(anno_df,by="probe_id") %>%
	dplyr::select(-probe_id) %>%
	dplyr::select(symbol,everything()) %>%
	mutate(rowMean =rowMeans(.[,-1])) %>%
	dplyr::arrange(exprSet_GSE43696$rowMean) %>%
	distinct(symbol,.keep_all = T) %>%
	dplyr::select(-rowMean) %>%
	column_to_rownames("symbol")
anyDuplicated(rownames(exprSet_GSE43696))

######exprSet_GSE63142####
exprSet_GSE63142 <- exprSet_GSE63142 %>%
	rownames_to_column("probe_id") %>%
	inner_join(anno_df,by="probe_id") %>%
	dplyr::select(-probe_id) %>%
	dplyr::select(symbol,everything()) %>%
	mutate(rowMean =rowMeans(.[,-1])) %>%
	dplyr::arrange(exprSet_GSE63142$rowMean) %>%
	distinct(symbol,.keep_all = T) %>%
	dplyr::select(-rowMean) %>%
	column_to_rownames("symbol")
anyDuplicated(rownames(exprSet_GSE63142))

######exprSet_GSE67940####
exprSet_GSE46171 <- exprSet_GSE46171 %>%
	rownames_to_column("probe_id") %>%
	inner_join(anno_df, by = "probe_id") %>%
	dplyr::select(-probe_id) %>%
	dplyr::select(symbol, everything()) %>%
	mutate(rowMean = rowMeans(.[, -1])) %>%
	arrange(desc(rowMean)) %>%
	distinct(symbol, .keep_all = TRUE) %>%
	dplyr::select(-rowMean) %>%
	column_to_rownames("symbol")

exprSet_names <- paste0("exprSet_", names)
pdata_names <- paste0("pdata_", names)
if (!dir.exists("expr_data")) dir.create("expr_data", recursive = TRUE)
if (!dir.exists("pdata")) dir.create("pdata", recursive = TRUE)

for (name in exprSet_names) {
	exprSet <- get(name)
	file_name <- paste0(name, ".txt")
	write.table(cbind(ID=rownames(exprSet),exprSet), file = file.path("expr_data",file_name), sep = "\t", quote = FALSE, row.names = F)
	cat("保存了表达数据:", file_name, "\n")
}

for (name in pdata_names) {
	pdata <- get(name)
	file_name <- paste0(name, ".txt")
	write.table(cbind(ID=rownames(pdata),pdata), file = file.path("pdata",file_name), sep = "\t", quote = FALSE, row.names = F)
	cat("保存了表型数据:", file_name, "\n")
}

rm(list=ls())
if (!dir.exists("nor")) dir.create("nor", recursive = TRUE)
exp_files <- list.files("./expr_data/", pattern = "*.txt", full.names = TRUE)
exp_data_list <- lapply(exp_files, function(file) {
	object_name <- gsub(".*_(.*)\\.txt", "\\1", basename(file))
	assign(object_name, read.table(file, sep="\t", header=TRUE, check.names=FALSE, row.names=1), envir = .GlobalEnv)
})

pdata_files <- list.files("./pdata", pattern = "*.txt", full.names = TRUE)
pdata_data_list <- lapply(pdata_files, function(file) {
	object_name <- gsub(".*_(.*)\\.txt", "\\1", basename(file))
	assign(paste0("geo_cli_", object_name), data.table::fread(file, sep="\t", header=TRUE, check.names=FALSE), envir = .GlobalEnv)
})


group_GSE43696 <- ifelse(str_detect(geo_cli_GSE43696$`disease state:ch1`,"Control"), "Con",
												 "Treat")
group_GSE63142 <- ifelse(str_detect(geo_cli_GSE63142$`individual:ch1`,"control"), "Con",
												 "Treat")
geo_cli_GSE46171 = geo_cli_GSE46171[!str_detect(geo_cli_GSE46171$`group:ch1`,"Allergic Rhinitis"),]
GSE46171 = GSE46171[,geo_cli_GSE46171$ID]
group_GSE46171 <- ifelse(str_detect(geo_cli_GSE46171$`group:ch1`,"healthy"), "Con",
												 "Treat")
exprSets <- exp_data_list
groupLists <- list(group_GSE43696,group_GSE46171,group_GSE63142)
exprSetNames <- gsub(".*_(.*)\\.txt", "\\1", basename(exp_files))

for (i in seq_along(exprSets)) {
	exprSet <- exprSets[[i]]
	groupList <- groupLists[[i]]

	condata <- NULL
	treatdata <- NULL

	if ("Con" %in% groupList) {
		condata <- exprSet[, which(groupList == "Con"), drop = FALSE]
		if(ncol(condata) > 0) {
			colnames(condata) <- paste0(colnames(condata), "_control")
		}
	} else {
		warning(paste("数据集", exprSetNames[i], "缺少Con组（Con）"))
	}

	if ("Treat" %in% groupList) {
		treatdata <- exprSet[, which(groupList == "Treat"), drop = FALSE]
		if(ncol(treatdata) > 0) {
			colnames(treatdata) <- paste0(colnames(treatdata), "_treat")
		}
	} else {
		warning(paste("数据集", exprSetNames[i], "缺少处理组（Treat）"))
	}

	exprSet_combined <- if (!is.null(condata) && !is.null(treatdata)) {
		cbind(condata, treatdata)
	} else if (!is.null(condata)) {
		condata
	} else if (!is.null(treatdata)) {
		treatdata
	} else {
		stop(paste("数据集", exprSetNames[i], "无有效分组"))
	}

	output_file <- file.path("nor", paste0(exprSetNames[i], ".normalize.txt"))
	assign(exprSetNames[i], exprSet_combined)
	write.table(
		cbind(ID = rownames(exprSet_combined), exprSet_combined),
		file = output_file,
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)
}

if("RCC2" %in% rownames(GSE46171)) {
	control_mean <- mean(as.numeric(GSE46171["RCC2", group_GSE46171 == "Con"]))
	treat_mean <- mean(as.numeric(GSE46171["RCC2", group_GSE46171 == "Treat"]))
	cat(paste0("GSE67940 RCC2表达验证:\nControl组均值: ", round(control_mean, 2),
						 "\nTreat组均值: ", round(treat_mean, 2),
						 "\nLog2FC: ", round(treat_mean - control_mean, 2)))
}
