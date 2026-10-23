library(ggplot2)

# Proportion heatmap: for each cluster, the fraction of its cells that carry each annotation.
#
# object            Seurat object
# annotation_col    metadata column with the labels (e.g. "celltype.mapped", "SingleR")
# cluster_col       metadata column with the clusters
# title             plot title
# order_annotations place each annotation next to the cluster where it is most common
# exclude           labels to hide, e.g. c("Doublet", "Stripped"); proportions are
#                   then computed from the remaining cells
make_annotation_heatmap <- function(
    object,
    annotation_col,
    cluster_col = "RNA_snn_res.0.3",
    title = NULL,
    order_annotations = TRUE,
    exclude = NULL
) {

  metadata <- data.frame(
    cluster    = as.character(object@meta.data[[cluster_col]]),
    annotation = as.character(object@meta.data[[annotation_col]]),
    stringsAsFactors = FALSE
  )
  metadata <- metadata[!is.na(metadata$cluster) & !is.na(metadata$annotation) &
                         !(metadata$annotation %in% exclude), ]

  # Counts, then proportion within each cluster
  plot_data <- as.data.frame(table(cluster = metadata$cluster,
                                   annotation = metadata$annotation),
                             stringsAsFactors = FALSE)
  names(plot_data)[3] <- "n"
  plot_data <- plot_data[plot_data$n > 0, ]

  cluster_totals <- aggregate(n ~ cluster, data = plot_data, FUN = sum)
  names(cluster_totals)[2] <- "cluster_total"
  plot_data <- merge(plot_data, cluster_totals, by = "cluster")
  plot_data$proportion <- plot_data$n / plot_data$cluster_total

  # Cluster order: numeric if the names are numbers, otherwise alphabetical
  cl <- unique(plot_data$cluster)
  cluster_levels <- if (!anyNA(suppressWarnings(as.numeric(cl)))) {
    as.character(sort(as.numeric(cl)))
  } else sort(cl)
  plot_data$cluster <- factor(plot_data$cluster, levels = rev(cluster_levels))  # cluster 0 on top

  # Annotation order: place each annotation next to the cluster where it is most common
  ann_levels <- sort(unique(plot_data$annotation))
  if (order_annotations) {
    top_cluster <- sapply(ann_levels, function(a) {
      d <- plot_data[plot_data$annotation == a, ]
      match(as.character(d$cluster[which.max(d$proportion)]), cluster_levels)
    })
    ann_levels <- ann_levels[order(top_cluster)]
  }
  plot_data$annotation <- factor(plot_data$annotation, levels = ann_levels)

  p <- ggplot(plot_data, aes(x = annotation, y = cluster, fill = proportion)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_viridis_c(name = "Cells in\ncluster", labels = scales::percent) +
    labs(x = NULL, y = "Cluster", title = title) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
          axis.text.y = element_text(size = 8))

  list(plot = p, data = plot_data)
}
