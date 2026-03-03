# -----------------------------------------------
# DynamicPlotMonocle2 - Monocle2版动态曲线绘图（与 DynamicPlot 风格对齐）
# 关键特性：
# 1) GAM拟合 + 预测；若预测全NA，自动 fallback 到 geom_smooth(method="gam")
# 2) 线条/置信带颜色由 line_palcolor/line_palette 控制
# 3) group.by 的散点颜色由 point_palcolor/point_palette 控制（通过 ggnewscale 分离刻度）
# 4) 主题默认使用用户提供的 theme_scp（若不存在则退化为 theme_minimal）
# -----------------------------------------------

DynamicPlotMonocle2 <- function(
		srt, pseudotime, features,
		group.by = NULL, cells = NULL,
		slot = "counts", assay = NULL,
		exp_method = c("log1p", "raw", "zscore", "fc", "log2fc"),
		lib_normalize = identical(slot, "counts"), libsize = NULL,
		compare_features = TRUE,               # 多基因时按 feature 分面
		add_line = TRUE, add_interval = TRUE, line.size = 1, line_palette = "Dark2", line_palcolor = NULL,
		add_point = TRUE, pt.size = 1, point_palette = "Paired", point_palcolor = NULL,
		add_rug = TRUE, flip = FALSE, reverse = FALSE, x_order = c("value", "rank"),
		aspect.ratio = NULL,
		legend.position = "bottom", legend.direction = "horizontal",
		theme_use = NULL, theme_args = list(),   # 默认会在函数内自动设为 theme_scp 或 theme_minimal
		combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, seed = 11
) {
	set.seed(seed)
	
	# ---- 依赖检查 ----
	if (!requireNamespace("mgcv", quietly = TRUE)) stop("Please install 'mgcv'.")
	if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'.")
	if (!requireNamespace("patchwork", quietly = TRUE)) stop("Please install 'patchwork'.")
	if (!requireNamespace("Seurat", quietly = TRUE)) stop("Please install 'Seurat'.")
	if (!requireNamespace("ggnewscale", quietly = TRUE)) stop("Please install 'ggnewscale'.")
	if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("Please install 'RColorBrewer'.")
	
	# ---- 小工具 ----
	`%||%` <- function(a, b) if (is.null(a)) b else a
	
	# 尝试默认主题：优先 theme_scp（若环境里存在），否则 theme_minimal
	if (is.null(theme_use)) {
		if (exists("theme_scp", mode = "function")) {
			theme_use <- get("theme_scp", mode = "function")
		} else {
			theme_use <- ggplot2::theme_minimal
		}
	}
	
	# 颜色工具：优先 palette_scp（若环境里存在），否则 RColorBrewer 兜底
	.palette_vals <- function(vals, palette, palcolor) {
		vals_chr <- unique(as.character(vals))
		if (!is.null(palcolor)) {
			# 若传入具名向量，按名称重排；否则直接按顺序回传
			if (!is.null(names(palcolor)) && length(intersect(names(palcolor), vals_chr)) > 0) {
				pal <- palcolor[vals_chr]
				# 若有缺失，尝试用存在的颜色填充
				pal[is.na(pal)] <- palcolor[1]
				return(pal)
			}
			if (length(palcolor) >= length(vals_chr)) return(palcolor[seq_along(vals_chr)])
			return(grDevices::colorRampPalette(palcolor)(length(vals_chr)))
		}
		if (exists("palette_scp", mode = "function")) {
			return(palette_scp(vals_chr, palette = palette))
		} else {
			n <- max(3, length(vals_chr))
			return(suppressWarnings(RColorBrewer::brewer.pal(min(n, 8), palette)))
		}
	}
	
	# ---- 参数与数据准备 ----
	x_order    <- match.arg(x_order)
	exp_method <- match.arg(exp_method)
	
	if (!pseudotime %in% colnames(srt@meta.data)) stop(pseudotime, " is not in meta.data.")
	
	assay <- assay %||% Seurat::DefaultAssay(srt)
	gene  <- features[features %in% rownames(srt@assays[[assay]])]
	meta  <- features[features %in% colnames(srt@meta.data)]
	features <- c(gene, meta)
	if (length(features) == 0) stop("No requested features found.")
	
	cells <- cells %||% colnames(srt)
	valid_cells <- intersect(cells, colnames(srt))
	if (length(valid_cells) == 0) stop("No valid cells selected.")
	
	# --- pseudotime: 强制数值且剔除 NA
	pseudotime_vec <- srt@meta.data[valid_cells, pseudotime]
	if (!is.numeric(pseudotime_vec)) {
		if (is.factor(pseudotime_vec)) pseudotime_vec <- suppressWarnings(as.numeric(as.character(pseudotime_vec))) else
			pseudotime_vec <- suppressWarnings(as.numeric(pseudotime_vec))
	}
	if (anyNA(pseudotime_vec)) {
		keep <- !is.na(pseudotime_vec)
		warning("NA in pseudotime – removing affected cells.")
		valid_cells    <- valid_cells[keep]
		pseudotime_vec <- pseudotime_vec[keep]
	}
	
	# --- 表达矩阵
	if (length(gene) > 0) {
		raw_gene <- Seurat::GetAssayData(srt, slot = slot, assay = assay)[gene, valid_cells, drop = FALSE]
		raw_gene <- as.matrix(t(raw_gene))
	} else {
		raw_gene <- matrix(nrow = length(valid_cells), ncol = 0)
		colnames(raw_gene) <- character(0)
	}
	if (length(meta) > 0) {
		raw_meta <- as.matrix(srt@meta.data[valid_cells, meta, drop = FALSE])
	} else {
		raw_meta <- matrix(nrow = length(valid_cells), ncol = 0)
		colnames(raw_meta) <- character(0)
	}
	raw_matrix <- cbind(raw_gene, raw_meta)
	
	# --- 文库标准化（仅基因）
	if (isTRUE(lib_normalize) && length(gene) > 0) {
		if (is.null(libsize)) {
			libsize_use <- colSums(Seurat::GetAssayData(srt, slot = "counts", assay = assay)[, valid_cells, drop = FALSE])
			if (any(libsize_use == 0)) { warning("Zero library size detected. Replacing with 1."); libsize_use[libsize_use == 0] <- 1 }
		} else {
			libsize_use <- libsize
			if (length(libsize_use) != length(valid_cells)) stop("'libsize' length must match number of cells.")
		}
		if (!any(raw_matrix[, gene, drop = FALSE] < 0, na.rm = TRUE)) {
			raw_matrix[, gene] <- sweep(raw_matrix[, gene, drop = FALSE], 1, libsize_use, "/") * stats::median(libsize_use)
		}
	}
	
	# --- 变换
	if (exp_method == "zscore") {
		cm <- apply(raw_matrix, 2, function(x) mean(x, na.rm = TRUE))
		cs <- apply(raw_matrix, 2, function(x) stats::sd(x, na.rm = TRUE))
		raw_matrix <- scale(raw_matrix, center = cm, scale = cs)
	} else if (exp_method == "fc") {
		colm <- apply(raw_matrix, 2, function(x) mean(x, na.rm = TRUE))
		raw_matrix <- t(t(raw_matrix) / (colm + 1e-12))
	} else if (exp_method == "log2fc") {
		colm <- apply(raw_matrix, 2, function(x) mean(x, na.rm = TRUE))
		raw_matrix <- t(log2(t(raw_matrix) / (colm + 1e-12) + 1e-6))
	} else if (exp_method == "log1p") {
		raw_matrix <- log1p(raw_matrix)
	}
	raw_matrix[!is.finite(raw_matrix)] <- NA
	
	# --- GAM 拟合（命名 data + fallback 友好）
	ncell <- nrow(raw_matrix); nfeat <- ncol(raw_matrix)
	fitted_matrix <- matrix(NA_real_, nrow = ncell, ncol = nfeat, dimnames = list(NULL, colnames(raw_matrix)))
	upr_matrix    <- matrix(NA_real_, nrow = ncell, ncol = nfeat, dimnames = list(NULL, colnames(raw_matrix)))
	lwr_matrix    <- matrix(NA_real_, nrow = ncell, ncol = nfeat, dimnames = list(NULL, colnames(raw_matrix)))
	
	for (feature in colnames(raw_matrix)) {
		y <- raw_matrix[, feature]
		if (all(is.na(y))) next
		ok <- is.finite(y) & is.finite(pseudotime_vec)
		if (sum(ok) < 5) next
		df_fit <- data.frame(y = y[ok], x = pseudotime_vec[ok])
		fit <- tryCatch(mgcv::gam(y ~ mgcv::s(x, bs = "cs"), data = df_fit, method = "REML"), error = function(e) NULL)
		if (is.null(fit)) next
		df_pred <- data.frame(x = as.numeric(pseudotime_vec))
		pred <- tryCatch(stats::predict(fit, newdata = df_pred, se.fit = TRUE, type = "response"), error = function(e) NULL)
		if (is.null(pred)) next
		fitted_matrix[, feature] <- pred$fit
		upr_matrix[, feature]    <- pred$fit + 1.96 * pred$se.fit
		lwr_matrix[, feature]    <- pred$fit - 1.96 * pred$se.fit
	}
	
	# --- 绘图数据
	plot_data <- data.frame(
		Cell       = rep(valid_cells, times = length(features)),
		Pseudotime = rep(as.numeric(pseudotime_vec), times = length(features)),
		Features   = rep(features, each = length(valid_cells)),
		Raw        = as.vector(raw_matrix[, features, drop = FALSE]),
		Fitted     = as.vector(fitted_matrix[, features, drop = FALSE]),
		Upr        = as.vector(upr_matrix[, features, drop = FALSE]),
		Lwr        = as.vector(lwr_matrix[, features, drop = FALSE]),
		stringsAsFactors = FALSE
	)
	
	if (!is.null(group.by) && group.by %in% colnames(srt@meta.data)) {
		plot_data[[group.by]] <- srt@meta.data[plot_data$Cell, group.by, drop = TRUE]
	} else if (!is.null(group.by)) {
		warning(group.by, " not in meta.data; ignore."); group.by <- NULL
	}
	
	# --- x 轴处理（rank / reverse / flip）
	if (x_order == "rank") {
		plot_data$Pseudotime <- ave(plot_data$Pseudotime, plot_data$Features, FUN = rank, ties.method = "average")
	}
	x_trans <- ifelse(flip, "reverse", "identity")
	x_trans <- ifelse(reverse, setdiff(c("reverse","identity"), x_trans), x_trans)
	
	# --- 构图（每个 feature 一张，主题/调色对齐 DynamicPlot；用 ggnewscale 分离颜色刻度）
	plist <- list()
	for (feat in unique(plot_data$Features)) {
		df <- plot_data[plot_data$Features == feat, , drop = FALSE]
		df <- df[order(df$Pseudotime), , drop = FALSE]
		
		p <- ggplot2::ggplot(df, ggplot2::aes(x = Pseudotime)) +
			ggplot2::scale_x_continuous(trans = x_trans, expand = ggplot2::expansion(c(0, 0))) +
			ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0.1, 0.05)))
		
		# 置信带（按 Features 上色；单基因即单色）
		if (isTRUE(add_interval) && any(is.finite(df$Lwr) & is.finite(df$Upr))) {
			p <- p +
				ggplot2::geom_ribbon(
					ggplot2::aes(ymin = Lwr, ymax = Upr, fill = .data[["Features"]]),
					alpha = 0.4, color = "grey90", na.rm = TRUE
				) +
				ggplot2::scale_fill_manual(
					values = .palette_vals(df$Features, line_palette, line_palcolor),
					guide  = "none"
				) +
				ggnewscale::new_scale_fill()  # 之后若再用 fill，不影响带状区
		}
		
		# 曲线（优先 Fitted；否则 fallback 到 geom_smooth(method="gam")）
		if (isTRUE(add_line)) {
			if (any(is.finite(df$Fitted))) {
				p <- p + ggplot2::geom_line(
					ggplot2::aes(y = Fitted, color = .data[["Features"]]),
					linewidth = line.size, alpha = 0.8, na.rm = TRUE
				)
			} else {
				p <- p + ggplot2::geom_smooth(
					ggplot2::aes(y = Raw, color = .data[["Features"]]),
					method = "gam", formula = y ~ s(x, bs = "cs"),
					se = TRUE, linewidth = line.size
				)
			}
			p <- p + ggplot2::scale_color_manual(
				values = .palette_vals(df$Features, line_palette, line_palcolor),
				guide  = "none"
			) +
				ggnewscale::new_scale_color()  # 切换到另一套颜色刻度，给散点使用
		}
		
		# 原始点（按 group.by 着色时走 point_palette / point_palcolor）
		if (isTRUE(add_point)) {
			if (is.null(group.by)) {
				p <- p + ggplot2::geom_point(ggplot2::aes(y = Raw), size = pt.size, alpha = 0.8, na.rm = TRUE)
			} else {
				p <- p +
					ggplot2::geom_point(ggplot2::aes(y = Raw, color = .data[[group.by]]), size = pt.size, alpha = 0.8, na.rm = TRUE) +
					ggplot2::scale_color_manual(
						values = .palette_vals(df[[group.by]], point_palette, point_palcolor),
						guide  = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 3))
					)
			}
		}
		
		# rug
		if (isTRUE(add_rug)) {
			if (is.null(group.by)) {
				p <- p + ggplot2::geom_rug(sides = "b", alpha = 1, length = grid::unit(0.05, "npc"), show.legend = FALSE)
			} else {
				p <- p + ggplot2::geom_rug(ggplot2::aes(color = .data[[group.by]]), sides = "b",
																	 alpha = 1, length = grid::unit(0.05, "npc"), show.legend = TRUE)
			}
		}
		
		# 主题/图例
		ylab <- if (exp_method == "raw") "Expression" else paste0(exp_method, " expression")
		p <- p +
			ggplot2::labs(
				title = if (isTRUE(compare_features)) feat else NULL,
				x = ifelse(x_order == "rank", "Pseudotime (rank)", "Pseudotime"),
				y = ylab
			) +
			do.call(theme_use, c(list(aspect.ratio = aspect.ratio), theme_args)) +
			ggplot2::theme(
				legend.position = legend.position,
				legend.direction = legend.direction
			)
		
		if (isTRUE(flip)) p <- p + ggplot2::coord_flip()
		
		plist[[feat]] <- p
	}
	
	# 组合
	if (isTRUE(combine)) {
		if (length(plist) > 1) {
			plot <- patchwork::wrap_plots(plist, nrow = nrow, ncol = ncol, byrow = byrow)
		} else {
			plot <- plist[[1]]
		}
		return(plot)
	} else {
		return(plist)
	}
}
