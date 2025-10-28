.PlotCluster <- function(object, reduction = NULL, cells = NULL, outfile = NULL,
				p1.group.by = "orig.ident",      p1.color = NULL, p1.label = FALSE,
				p2.group.by = "seurat_clusters", p2.color = NULL, p2.label = TRUE,
				plot.basic.size = 6, ...){
		if ( is.null(p1.color) && ! is.null(p1.group.by) ) {
				p1.color <- switch(p1.group.by,
						"Groups"          = object@misc[["color.group"]],
						"orig.ident"      = object@misc[["color.sample"]],
						"seurat_clusters" = object@misc[["color.cluster"]])
				if ( ! is.null(cells) ) {
						p1.color <- p1.color[levels(droplevels(object@meta.data[cells, p1.group.by]))]
				}
		}
		if ( is.null(p2.color) && ! is.null(p2.group.by) ) {
				p2.color <- switch(p2.group.by,
						"Groups"          = object@misc[["color.group"]],
						"orig.ident"      = object@misc[["color.sample"]],
						"seurat_clusters" = object@misc[["color.cluster"]])
				if ( ! is.null(cells) ) {
						p2.color <- p2.color[levels(droplevels(object@meta.data[cells, p2.group.by]))]
				}
		}
		p1 <- DimPlot(object, reduction = reduction, cells = cells, group.by = p1.group.by, cols = p1.color, label = p1.label, ...)
		p2 <- DimPlot(object, reduction = reduction, cells = cells, group.by = p2.group.by, cols = p2.color, label = p2.label, ...)
		p1 <- p1 + dot_theme_default() + ggtitle(NULL)
		p2 <- p2 + dot_theme_default() + ggtitle(NULL)
		if ( is.null(p2.group.by) ) {
				p <- p1
				width  <- plot.basic.size * 1.2
				height <- plot.basic.size
		} else if ( is.null(p1.group.by) ) {
				p <- p2
				width  <- plot.basic.size * 1.2
				height <- plot.basic.size
		} else {
				p <- p1 + p2
				width  <- plot.basic.size * 1.2 * 2
				height <- plot.basic.size
		}
		if ( is.null(outfile) ) {
				return(p)
		}else{
				ggsave(p, file = outfile, width = width, height = height, limitsize = FALSE )
		}
}


PlotCluster <- function(object, reduction = 'umap', p1.group.by = "orig.ident", split.by = p1.group.by, p2.group.by = "seurat_clusters", outpref = NULL, ...){
		.PlotCluster(object, reduction = reduction, p1.group.by = p1.group.by, p2.group.by = p2.group.by, outfile = paste0(outpref, ".pdf"), ...)
		for ( i in unique(object@meta.data[[split.by]]) ){
				cells.use <- rownames(object@meta.data)[object@meta.data[[split.by]] == i]
				.PlotCluster(object, reduction = reduction, cells = cells.use, p1.group.by = p1.group.by, p2.group.by = p2.group.by, outfile = paste0(outpref, ".", i, ".pdf"), ...)
		}
		data <- object[[reduction]]@cell.embeddings %>% as.data.frame() %>%
				tibble::rownames_to_column(var = "Cells") %>%
				left_join(.GetMetaData(object, cols = c("Samples" = p1.group.by, "Cluster" = p2.group.by, "Groups")))
		WriteTable(data, file = paste0(outpref, ".plot.data.tmp"))
}

.PlotClusterStat <- function(object, stat.what = "seurat_clusters", group.by = "orig.ident", color.st = NULL, color.gb = NULL, outpref = NULL, ...){
		if ( class(object) == "Seurat" ) {
				metadata <- object@meta.data
		} else {
				metadata <- object
		}
		if ( is.null(color.st) ) {
				if ( "misc" %in% slotNames(object) && exists(stat.what, object@misc) ) {
						color.st <- object@misc[[stat.what]]
				} else {
				color.st <- switch(stat.what,
						"seurat_clusters" = object@misc$color.cluster,
						"orig.ident" = object@misc$color.sample,
						"Groups" = object@misc$color.group
				)
				}
		}
		if ( is.null(color.gb) ) {
				if ( "misc" %in% slotNames(object) && exists(group.by, object@misc) ) {
						color.gb <- object@misc[[group.by]]
				} else {
				color.gb <- switch(group.by, 
						"seurat_clusters" = object@misc$color.cluster,
						"orig.ident" = object@misc$color.sample,
						"Groups" = object@misc$color.group
				)
				}
		}
		name.st <- switch(stat.what, "seurat_clusters" = "Cluster", "orig.ident" = "Samples", stat.what)
		name.gb <- switch(group.by,  "seurat_clusters" = "Cluster", "orig.ident" = "Samples", group.by)
		stat.what <- as.name(stat.what)
		group.by  <- as.name(group.by)

		stat_sample <- metadata %>%
				group_by(!! name.gb := !! group.by, !! name.st := !! stat.what) %>%
				summarise("Number of cells" = n())

		p <- list()
		p[["by"]] <- ggplot(stat_sample, aes_(x = as.name(name.gb), y = ~ `Number of cells`, fill = as.name(name.st)))
		p[["in"]] <- ggplot(stat_sample, aes_(x = as.name(name.st), y = ~ `Number of cells`, fill = as.name(name.gb)))
		if ( ! is.null(color.st) ) p[["by"]] <- p[["by"]] + scale_fill_manual(values = color.st)
		if ( ! is.null(color.gb) ) p[["in"]] <- p[["in"]] + scale_fill_manual(values = color.gb)
		geom_stack <- geom_bar(stat = "identity", position = 'stack')
		geom_fill  <- geom_bar(stat = "identity", position = "fill" )

		if ( is.null(outpref) ) {
				outpref <- paste0( name.st, ".stat")
		}
		for ( i in names(p) ) {
				p[[i]] <- p[[i]] + bar_theme_default()
				ggsave( p[[i]] + geom_stack, file = paste0( outpref, ".", i, name.gb, ".pdf"),     height = 6, width = 8 )
				ggsave( p[[i]] + geom_fill + ylab("Fraction of Cells"),  file = paste0( outpref, ".", i, name.gb, ".pct.pdf"), height = 6, width = 8 )
		}
}
