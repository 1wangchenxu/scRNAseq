PlotHeatmapPlot <- function(object, features = NULL, group.by = "seurat_clusters", is.use.name = TRUE,
							outfile = NULL, group.colors = NULL, color = c("#FF00FF", "#000000", "#FFFF00"), slot = 'scale.data', assay = NULL, ...) {
		if ( is.null(group.colors) ) {
				group.colors <- switch(group.by,
								"seurat_clusters" = object@misc$color.cluster,
								"orig.ident" = object@misc$color.sample,
								"Groups" = object@misc$color.group
								)

		} else {
				group.colors <- group.colors[levels(object@meta.data[[group.by]])]
		}
		if ( slot == 'scale.data' ) {
				if ( ! is.null(features) ) {
						if ( is.null(assay) ) assay <- DefaultAssay(object)
						DefaultAssay(object) <- assay
						if ( ! all(features %in% rownames(slot(object@assays[[assay]], slot))) ) {
								message('Below features arenot in scale.data, recal scale.data:')
								message(paste0(features[! features %in% rownames(slot(object@assays[[assay]], slot))], collapse = ", "))
								object <- ScaleData(object, features = features)
						}
				}
		}
		p <- DoHeatmap( object = object, features = features, cells = NULL,
						group.by = group.by, group.colors = group.colors,
						combine = FALSE, raster = FALSE, slot = slot, assay = assay, ...)
		p <- p[[1]]
		p <- p + theme(legend.title = element_blank())
		p$layers[[2]] <- NULL

		if ( is.use.name ) {
				levels(p$data$Feature) <- FindFeaturesName(object, levels(p$data$Feature))
		}

		if (length(color) == 2) {
				p <- p + scale_fill_gradient(low = color[1], high = color[2], na.value = "white")
		} else if (length(color) == 3) {
				p <- p + scale_fill_gradient2(low = color[1], mid = color[2], high = color[3], na.value = "white")
		} else if (length(color) > 3) {
				p <- p + scale_fill_gradientn(colors = color, na.value = "white")
		}

#		p <- ggplotGrob(p)
#		for (i in grep('strip', p$layout$name)){
#				p$grobs[[i]]$layout$clip <- "off"
#		}

		if ( is.null(outfile) ) {
				return(p)
		} else {
				h <- max(7, length(unique(features)) * 0.11 + 2.5 )
				w <- h * 4 / 3
				p <- p + theme(plot.margin=margin(t=1,r=1,unit="lines"))
				ggsave(p, file = outfile, width = w, height = h, limitsize = FALSE )
		}
}

