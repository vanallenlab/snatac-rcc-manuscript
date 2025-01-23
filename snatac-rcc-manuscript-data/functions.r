calculate_signature_scores <- function(seurat_object, sig_peak_set, pct.open, final_score_name) {
  # Average normalized counts for the significant peak set
  sum_per_cell_sig = colSums(seurat_object@assays$ATAC@data[sig_peak_set, ])
  avg_per_cell_sig = sum_per_cell_sig / length(sig_peak_set)

  # Set default assay to ATAC
  DefaultAssay(seurat_object) <- 'ATAC'

  # Find peaks that are open in a significant portion of cells
  twenty_percent = ceiling(nrow(seurat_object@meta.data) * pct.open)

  open_peaks <- AccessiblePeaks(seurat_object, min.cells = twenty_percent)

  # Match the overall GC content in the peak set
  meta.feature <- GetAssayData(seurat_object, assay = "ATAC", slot = "meta.features")

  peaks_matched <- MatchRegionStats(
    meta.feature = meta.feature[open_peaks, ],
    query.feature = meta.feature[sig_peak_set, ]
  )
    

  # Average normalized counts for the reference peak set
  sum_per_cell_ref = colSums(seurat_object@assays$ATAC@data[peaks_matched, ])
  avg_per_cell_ref = sum_per_cell_ref / length(peaks_matched)

  # Calculate final signature scores
  sig_scores = avg_per_cell_sig - avg_per_cell_ref
  sig_scores = data.frame(sig_scores)
  colnames(sig_scores) = final_score_name

  return(sig_scores)
}

create_motif_significance_plot <- function(enriched.motifs, max_logp = 320, n_annotate = 10, specific_substr = NULL, pval_threshold = 0.05, low_color_bound = "lightblue", high_color_bound = "blue", color_gradient_limits = c(0, 2), max_overlaps = Inf) {
  # Transform the data frame
  df <- enriched.motifs %>%
    mutate(neg_log_padjust = ifelse(p.adjust == 0, max_logp, -log10(p.adjust))) %>%
    mutate(neg_log_padjust = pmin(neg_log_padjust, max_logp)) %>%  # Apply the maximum limit
    arrange(-neg_log_padjust) %>%
    mutate(rank = row_number()) %>%
    mutate(significant = p.adjust < pval_threshold) %>% # Determine significance
    mutate(fold_enrichment_color = ifelse(significant, fold.enrichment, NA)) # Create a separate column for coloring based on significance
  
  # Create the plot with conditional coloring
  p <- ggplot(df, aes(x = rank, y = neg_log_padjust)) +
    geom_point(aes(color = fold_enrichment_color), show.legend = TRUE) +
    scale_color_gradient(low = low_color_bound, high = high_color_bound, limits = color_gradient_limits, na.value = "grey", name = "Fold Enrichment") +
    labs(x = "Rank", y = "-log10(p.adjust)", title = "Ranked TF Motif Significance") +
    theme_classic() +
    theme(text = element_text(size = 20))
  
  # Annotate top motifs based on n_annotate
  top_motifs <- head(df, n_annotate)
  
  # Filter and annotate specific motifs if specific_substr is provided
  if (!is.null(specific_substr) && length(specific_substr) > 0) {
    specific_motifs <- df %>%
      filter(grepl(paste(specific_substr, collapse = "|"), motif.name)) %>%
      filter(!motif.name %in% top_motifs$motif.name)  # Exclude motifs already in top_motifs
    
    # Merge specific_motifs with top_motifs to ensure all are included in annotation
    annotated_motifs <- bind_rows(top_motifs, specific_motifs)
  } else {
    annotated_motifs <- top_motifs
  }
  
  # Add labels for annotated motifs
  if (nrow(annotated_motifs) > 0) {
    p <- p + geom_text_repel(data = annotated_motifs, aes(label = motif.name), max.overlaps = max_overlaps,
                             size = 5, nudge_y = 0.5, nudge_x = 0.5, color = "black")
  }
  
  return(p)  # Return the plot
}



plot_signature_by_mutation <- function(data, mutation_column, signature, side='less') {
  # Convert mutation status to a factor
  data[[mutation_column]] <- as.factor(data[[mutation_column]])
  
  # Calculate the median signature for each biopsy
  signature_df <- data %>%
    group_by(biopsy, !!sym(mutation_column)) %>%
    summarise(
      cells_per_biopsy = n(),
      median_sig = median(!!sym(signature)),  # Assuming 'stage' is available in the original 'data'
      .groups = 'drop'
    )

  # Rank the points within each mutation group and select the top and bottom 5
  top_bottom_df <- signature_df %>%
    group_by(!!sym(mutation_column)) %>%
    mutate(rank = row_number(median_sig)) %>%
    filter(rank <= 5 | rank > n() - 5) %>%
    ungroup()

  # Create the box plot with overlaid jittered points
  p <- ggplot(signature_df, aes(x = !!sym(mutation_column), y = median_sig)) +
    geom_boxplot(alpha = 1, aes(fill = !!sym(mutation_column))) +
    geom_jitter(width = 0.2, alpha = 1, aes(size = cells_per_biopsy, fill = !!sym(mutation_column)), stroke = 0.5, shape = 21) +
    labs(
      x = paste(mutation_column, "Mutation Status"),
      y = paste(signature)
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
    ) +
    scale_fill_manual(values = c('#3A3B3C','grey'), 
                      guide = guide_legend(title = mutation_column)) +
    geom_signif(comparisons = list(c('0', '1')), map_signif_level = FALSE, test="wilcox.test", test.args=list(alternative = side))
    p
  # Add text labels for the top and bottom 5 points
  #p + geom_text(data = top_bottom_df, aes(label = biopsy), position = position_jitter(width = 0.2), vjust = -0.5, size = 3)
}

# Example usage:
# plot_signature_by_mutation(data = shared_states, mutation_column = "PBRM1", signature = "some_signature_column")

write_hg19_bed <- function(obj, peaks, output_path) {
    
    granges = Signac:::FindRegion(obj,peaks)
    granges$name = peaks
    path = "data/scripts/hg38ToHg19.over.chain"
    ch = import.chain(path)

    seqlevelsStyle(granges) = "UCSC"  # necessary
    cur19 = liftOver(granges, ch)
    cur19 = unlist(cur19)
    
    bed <- data.frame(seqnames=seqnames(cur19),
    starts=start(cur19)-1,
    ends=end(cur19),
    name = cur19$name)
    write.table(bed, file=output_path, quote=F, sep="\t", row.names=F, col.names=F)
}
