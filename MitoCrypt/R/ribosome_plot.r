#' profile supstarts
#'
#' This function plots ribosome profile and displays supported NTE isoforms as vertical lines
#' @param riboseq riboseq_input.txt from scikit-ribo output
#' @param mitocrypt MitoCrypt table
#' @param query_gene id of gene of interest
#' @param alternative_name alternative name to be displayed as plot title
#' @param print_pdf logical, if TRUE pdf of plot will be printed in working directory
#' @param custom_window distance from first supported isoform to be plotted. In codons.
#' @param mitocrypt_supported MitoCrypt table of supported isoforms (output of isoform_expression_wrapper())
#' @return plot of ribosome profile
#' @keywords ribosome profile
#' @seealso isoform_expression_wrapper()
#' @export
#' @examples
#' 

profile_supstarts <- function(riboseq, mitocrypt, query_gene, alternative_name = NULL, print_pdf = F, custom_window = 300, mitocrypt_supported){
  if (is.null(alternative_name)) {alternative_name <- query_gene}
  my_palette = colorRampPalette(c("blue", "yellow", "red"))(n = 100)
  nte_lengths <- mitocrypt$NTE.Length[which(mitocrypt$Gene.Name == query_gene)]
  mitoprot <- mitocrypt$Mitoprot.Probability[which(mitocrypt$Gene.Name == query_gene)]
  supported_ntes <- mitocrypt_supported$NTE.Length[which(mitocrypt_supported$Gene.Name == query_gene)]
  mitoprot_sup <- mitocrypt_supported$Mitoprot.Probability[which(mitocrypt_supported$Gene.Name == query_gene)]
  feature_length <- riboseq %>% filter(gene == query_gene) %>% nrow()
  threshold <- ifelse(feature_length < custom_window, feature_length, custom_window)
  ribo_count <- ifelse(riboseq$gene_strand[which(riboseq$gene == query_gene)] == "+", riboseq$ribosome_count[which(riboseq$gene == query_gene)][1:threshold], rev(riboseq$ribosome_count[which(riboseq$gene == query_gene)])[1:threshold])
  codon_index <- ifelse(riboseq$gene_strand[which(riboseq$gene == query_gene)] == "+", riboseq$codon_idx[which(riboseq$gene == query_gene)][1:threshold], rev(riboseq$codon_idx[which(riboseq$gene == query_gene)])[1:threshold])
  codon_index <- codon_index - max(mitocrypt_supported$NTE.Length[which(mitocrypt_supported$Gene.Name == query_gene)])
  plot_data <- data_frame(ribo_count,codon_index)
  p <- ggplot(plot_data,aes(x = codon_index, y = ribo_count)) +
    geom_line()+
    theme_bw()+
    ggtitle(alternative_name) + theme(plot.title = element_text(size = 25)) +
    xlab("Codon index") +
    ylab("Ribosome RPKM") + theme(axis.title = element_text(size = 15)) +
    scale_x_continuous(breaks = seq(-75,custom_window, 25))
  for (j in 1:length(nte_lengths)) {
    my_col <- my_palette[ceiling(mitoprot[j] * 100)]
    y_segment <- -max(ribo_count)/40
    p <- p + geom_segment(x = -nte_lengths[j] - 0.25 ,y = y_segment,xend = -nte_lengths[j] + 0.25,yend = y_segment, size = 5, col = 'black')
    p <- p + geom_segment(x = -nte_lengths[j],y = y_segment,xend = ifelse(j == 1, max(codon_index), -nte_lengths[j-1]),yend = y_segment, size = 5, col = my_col)
    #p <- p + geom_vline(xintercept = -nte_lengths[j], col = my_col)
    p <- p + theme(axis.text = element_text(size = 12),
                   axis.title = element_text(size = 14))
  }
  for (j in 1:length(supported_ntes)) {
    my_col <- my_palette[ceiling(mitoprot_sup[j] * 100)]
    #y_segment <- -max(ribo_count)/40
    #p <- p + geom_segment(x = -nte_lengths[j],y = y_segment,xend = ifelse(j == 1, max(codon_index), -nte_lengths[j-1]),yend = y_segment, size = 5, col = my_col)
    p <- p + geom_vline(xintercept = -supported_ntes[j], col = my_col)
  }
  p$layers <- rev(p$layers)
  ggsave(paste(alternative_name, ".pdf", sep = ""), plot = p)
}