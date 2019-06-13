#' Gene NTE assigner
#'
#' This function assigns potential extensions (and truncations) to gene of interest.
#' @param genome BSgenome name, compatible with gtf
#' @param gtf_table gtf data-frame annotation
#' @param start_codons vector of triplets to be searched as start codons.
#' @param stop_codons vector of triplets to be searched as stop codons
#' @param max_nte_length maximum NTE length in aminoacids.
#' @param genetic_code GENETIC_CODE in biostrings genetic code format.
#' @param context_span two element vector, corresponding to start and end of the sequence relative to first position of start codon 
#' @param truncations logical, if TRUE truncations will be searched for.
#' @param truncation_start_codons triplets to be searched for as truncations start codons, "ATG" by default.
#' @param max_truncation_length maximum truncation length in aminoacids.
#' @param gene transcript_id of gene of interest.
#' @return MitoCrypt_table of all isoforms called for particular gene.
#' @keywords NTE extension truncation
#' @seealso assign_NTE_wrapper()
#' @family NTE_assign
#' @export
#' @examples
#'





gene_NTE_assigner <- function(genome, gtf_table, start_codons = c("ACG", "ATA", "ATC", "ATG", "ATT", "CTG","GTG", "TTG"), stop_codons = (c("TAA", "TAG", "TGA")), max_nte_length = 300, genetic_code = GENETIC_CODE, context_span = c(-6,4), truncations = TRUE, truncation_start_codons = "ATG",max_truncation_length = 50, gene) {
  NTE_table <- data.frame()
  gene_gtf <- gtf_table[which(gtf_table$transcript_id == gene),] %>% arrange(start)
  if (gene_gtf$strand[1] == "+") {
    stop_search_threshold <- min(c(gene_gtf$start[1] - 1,max_nte_length * 3 + abs(context_span[1])))
    stop_search_gtf <- gene_gtf
    stop_search_gtf$start[1] <- stop_search_gtf$start[1] - stop_search_threshold
    full_search_sequence <- fetch_sequence(stop_search_gtf, genome, gene, output = "both")
  } else {
    stop_search_threshold <- min(c(seqlengths(genome)[[gene_gtf$seqid[1]]] - gene_gtf$end[nrow(gene_gtf)],max_nte_length * 3 + abs(context_span[1])))
    stop_search_gtf <- gene_gtf
    stop_search_gtf$end[nrow(stop_search_gtf)] <- stop_search_gtf$end[nrow(stop_search_gtf)] + stop_search_threshold
    full_search_sequence  <- fetch_sequence(stop_search_gtf, genome, gene, output = "both")
    full_search_sequence[[2]] <- sort(full_search_sequence[[2]], decreasing = TRUE)
  }
  
  stop_searchspace <- full_search_sequence[[1]][(abs(context_span[1]) + 1):stop_search_threshold]
  stop_matches <- match_multiple_codons(stop_codons, stop_searchspace)
  stop_matches <- stop_matches[which(((stop_matches) %% 3) == ((length(stop_searchspace) + 1) %% 3))]
  if (length(stop_matches) > 0) {
    stop_boundary <- max(stop_matches) + 3 
  } else {stop_boundary <- length(stop_searchspace) %% 3 + 1}
  full_search_sequence[[1]] <- full_search_sequence[[1]][stop_boundary:length(full_search_sequence[[1]])]
  full_search_sequence[[2]] <- full_search_sequence[[2]][stop_boundary:length(full_search_sequence[[2]])]
  canonical_start_codon_index <- length(stop_searchspace) - stop_boundary + abs(context_span[1]) + 2
  
  near_cognate_searchspace <- full_search_sequence[[1]][(abs(context_span[1]) + 1):(canonical_start_codon_index + 2)]
  near_cognate_matches <- match_multiple_codons(start_codons,near_cognate_searchspace)
  near_cognate_matches <- near_cognate_matches[which(near_cognate_matches %% 3 == 1)]
  near_cognate_matches <- sort(near_cognate_matches + abs(context_span[1]))
  start_matches <- near_cognate_matches
  
  
  if (truncations == TRUE) {
    truncation_search_threshold <- min(max_truncation_length * 3, length(full_search_sequence[[1]]) - canonical_start_codon_index - 3 - abs(context_span[1])) 
    truncation_searchspace <- full_search_sequence[[1]][(canonical_start_codon_index + 3) : (canonical_start_codon_index + 3 + truncation_search_threshold)]
    truncation_start_matches <- match_multiple_codons(truncation_start_codons, truncation_searchspace)
    truncation_start_matches <- truncation_start_matches[which(truncation_start_matches %% 3 == 1)]
    truncation_start_matches <- truncation_start_matches + canonical_start_codon_index + 2
    start_matches <- c(start_matches, truncation_start_matches)
  }
  if (length(start_matches) > 0) { 
  for (i in 1:length(start_matches)){
    start_cor <- min(full_search_sequence[[2]][start_matches[i]],full_search_sequence[[2]][length(full_search_sequence[[2]])])
    end_cor <- max(full_search_sequence[[2]][start_matches[i]],full_search_sequence[[2]][length(full_search_sequence[[2]])])
    isoform_gtf <- get_isoform_gtf(start = start_cor, end = end_cor, gtf_table = gene_gtf, transcript_id = gene)
    isoform_block <- block_generate(isoform_gtf, gene)
    NTE_length <- (canonical_start_codon_index - start_matches[i]) / 3
    context <- as.character(full_search_sequence[[1]][(start_matches[i] + context_span[1]):(start_matches[i] + context_span[2] - 1)])
    context <- context_tolower(context,context_span = context_span)
    NTE_sequence <- translate(full_search_sequence[[1]][sort(start_matches[i]:(start_matches[i] + NTE_length*3))][-(abs(NTE_length*3)+1)], genetic.code = genetic_code)
    isoform_sequence <- as.character(translate(full_search_sequence[[1]][start_matches[i]:length(full_search_sequence[[1]])], genetic.code = genetic_code))
    substr(isoform_sequence,1,1) <- "M"
    
    
    NTE_table <-rbind.data.frame(NTE_table,c(isoform_gtf$transcript_id[1],
                                             isoform_gtf$gene_name[1],
                                             isoform_gtf$seqid[1],
                                             start_cor,
                                             end_cor,
                                             isoform_gtf$strand[1],
                                             names(start_matches[i]),
                                             context,
                                             as.character(NTE_sequence),
                                             NTE_length,
                                             as.character(isoform_sequence),
                                             isoform_block$blockCount,
                                             paste(isoform_block$blockSizes, collapse = ","),
                                             paste(isoform_block$blockStarts, collapse = ",")),
                                 stringsAsFactors = FALSE)
    
  } 
  colnames(NTE_table) <- c('Gene.Name',
                           'Standard.Name',
                           'Chromosome',
                           'Start.Coordinate',
                           'End.Coordinate',
                           'Strand',
                           'Start.Codon',
                           'Start.Codon.Context',
                           'NTE.Sequence',
                           'NTE.Length',
                           'Isoform.Sequence',
                           'blockCount',
                           'blockSizes',
                           'blockStarts')
  return(NTE_table)
  } else {return()}
}
