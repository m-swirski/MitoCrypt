#' Isoform expression
#'
#' This function detects expressed NTE isoforms and returns absolute ribosome coverage 
#' @param riboseq riboseq_input.txt from scikit-ribo output
#' @param query gene_id 
#' @param mitocrypt MitoCrypt table
#' @param longest_nte additional MitoCrypt table object, with one entry of longest searched for isoform for each gene.
#' @param p_value p value for permuatational t test comparing consecutive NTEs
#' @param ramp length of canonical ORF to be compared with its proximal isoform. In codons.
#' @param nte_cov minimum coverage of NTE isoform. 
#' @return MitoCrypt table of called isoforms for single gene
#' @keywords isoform profiling detection quantification
#' @seealso isoform_expression_wrapper()
#' @export
#' @examples
#' 


isoform_expression <- function(riboseq, query, mitocrypt, longest_nte, p_value = 0.01, ramp = 50, nte_cov = 1){
  profile <- riboseq %>% filter(gene == query)
  profile <- profile[-1,]
  if (query %in% longest_nte$Gene.Name) {
    max_nte <- longest_nte$NTE.Length[which(longest_nte$Gene.Name == query)]
    ntes <-(mitocrypt%>% filter(Gene.Name == query) %>% filter(NTE.Length <= max_nte))$NTE.Length
  } else {
    ntes <- c(0)
    max_nte <- 0
  }
  # rel_ntes <- ntes
  ntes <- extract_ntes(nte_cov,profile, ntes, max_nte)
  # rel_ntes <- rev(rel_ntes[1:length(ntes)])
  rel_ntes <- -(ntes - max_nte)
  if(length(ntes) == 0) {return()}
  print(rel_ntes)
  print(ntes)
  # profile <- profile[-1,]
  #return(rel_ntes)
  iso_starts <- c()
  iso_occupancy <- c()
  p_values <- c()
  iso_count <- 1
  curstart_no <- 1
  threshold <- ifelse(length(profile$codon_idx[ntes[length(ntes)]:nrow(profile)]) < ramp, length(profile$codon_idx[ntes[length(ntes)]:nrow(profile)]) - 1, ramp)
  
  
  if (length(ntes) > 1){
    for (i in 1:(length(ntes)-1)){
      current_iso <- profile$ribosome_count[(ntes[curstart_no] + 1) : ntes[i+1]]
      candidate_iso <- profile$ribosome_count[(ntes[i+1] + 1) : ifelse(i == (length(ntes) - 1), ntes[i + 1] + threshold,ntes[i+2])]
      permtest <- as.numeric(twotPermutation(current_iso, candidate_iso, nsim = 10000,plotit = FALSE))
      if ((permtest < p_value ) & (mean(current_iso, na.rm = TRUE) < mean(candidate_iso, na.rm = TRUE))) {
        iso_starts[iso_count] <- rel_ntes[curstart_no] 
        iso_occupancy[iso_count] <-  mean(current_iso)
        p_values[iso_count] <- permtest
        iso_count <- iso_count + 1
        curstart_no <- i + 1
        if (i >= (length(ntes) - 1)) {
          iso_starts[iso_count] <- rel_ntes[curstart_no] 
          iso_occupancy[iso_count] <-  mean(candidate_iso) 
        }
      } else if (i >= length(ntes) - 1) {
        iso_starts[iso_count] <- rel_ntes[curstart_no]
        iso_occupancy[iso_count] <- mean(c(current_iso, candidate_iso))
        p_values[iso_count] <- permtest
      }
    }
  }else {
    iso_starts[1] <- min(rel_ntes)
    iso_occupancy[1] <- mean(profile$ribosome_count[1:threshold])
    p_values[1] <- NA
  }
  p_values[which(iso_starts == min(rel_ntes))] <- NA
  return(data.frame(gene = rep(query, length(iso_starts)),
                    actual_nte = rev(iso_starts),
                    mean_occupancy = rev(iso_occupancy) %>% round(.,digits = 3),
                    equal_p_value = rev(p_values) %>% round(.,digits = 4)))
}