DotPlotScores <- function(
		object,
		scores,  # Score column names (e.g., "Add_score")
		assay = NULL,
		cols = c("lightgrey", "blue"),
		dot.min = 0,
		dot.scale = 6,
		idents = NULL,
		group.by = NULL,
		split.by = NULL,
		scale = TRUE
) {
	assay <- assay %||% DefaultAssay(object = object)
	DefaultAssay(object = object) <- assay
	
	cells <- unlist(x = CellsByIdentities(object = object, cells = colnames(object[[assay]]), idents = idents))
	data.scores <- FetchData(object = object, vars = scores, cells = cells)
	
	# Get cell identities
	data.scores$id <- if (is.null(x = group.by)) {
		Idents(object = object)[cells, drop = TRUE]
	} else {
		object[[group.by, drop = TRUE]][cells, drop = TRUE]
	}
	
	if (!is.factor(x = data.scores$id)) {
		data.scores$id <- factor(x = data.scores$id)
	}
	
	# Compute mean scores by identity
	data.plot <- lapply(
		X = unique(x = data.scores$id),
		FUN = function(ident) {
			data.use <- data.scores[data.scores$id == ident, scores, drop = FALSE]
			avg.score <- colMeans(data.use, na.rm = TRUE)
			return(list(avg.score = avg.score))
		}
	)
	
	names(x = data.plot) <- unique(x = data.scores$id)
	data.plot <- lapply(
		X = names(x = data.plot),
		FUN = function(x) {
			data.use <- as.data.frame(x = data.plot[[x]])
			data.use$features.plot <- rownames(x = data.use)
			data.use$id <- x
			return(data.use)
		}
	)
	data.plot <- do.call(what = 'rbind', args = data.plot)
	
	avg.score.scaled <- sapply(
		X = unique(x = data.plot$features.plot),
		FUN = function(x) {
			data.use <- data.plot[data.plot$features.plot == x, 'avg.score']
			if (scale) {
				data.use <- scale(x = log1p(data.use))
			} else {
				data.use <- log1p(x = data.use)
			}
			return(data.use)
		}
	)
	
	avg.score.scaled <- as.vector(x = t(x = avg.score.scaled))
	data.plot$avg.score.scaled <- avg.score.scaled
	data.plot$features.plot <- factor(
		x = data.plot$features.plot,
		levels = scores
	)
	
	# Bubble size: percent of cells above dot.min
	data.plot$pct.exp <- apply(data.scores[, scores, drop = FALSE], 2, function(x) { sum(x > dot.min) }) / nrow(data.scores) * 100
	
	plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
		geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.score.scaled')) +
		scale_size(range = c(0, dot.scale),limits = c(30, 100)) +
		theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
		labs(
			x = 'Features',
			y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
		) +
		theme_cowplot()
	
	if (length(x = cols) == 1) {
		plot <- plot + scale_color_distiller(palette = cols)
	} else {
		plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
	}
	
	return(plot)
}
