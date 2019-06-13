ls
#' Isoform expression wrapper
#'
#' This function is multicore wrapper of isoform_expression()
#' @param riboseq riboseq_input.txt from scikit-ribo output
#' @param mitocrypt MitoCrypt table
#' @param longest_nte additional MitoCrypt table object, with one entry of longest searched for isoform for each gene.
#' @param p_value p value for permuatational t test comparing consecutive NTEs
#' @param ramp length of canonical ORF to be compared with its proximal isoform. In codons.
#' @param nte_cov minimum coverage of NTE isoform. 
#' @param cores number of cores to use. Passes argument to mcmapply() mc.cores argument. 
#' @return MitoCrypt table of called isoforms for all genes present in longest_nte MitoCrypt_table
#' @keywords isoform profiling detection quantification
#' @seealso isoform_expression()
#' @export
#' @examples
#' 



isoform_expression_wrapper <- function(riboseq,mitocrypt,longest_nte,p_value = 0.01, ramp = 50, nte_cov = 1, cores = 1) {
  gene_list <- levels(as.factor(as.character(longest_nte$Gene.Name)))
  MitoCrypt_output <- isoform_expression(riboseq, gene_list[1], mitocrypt = mitocrypt, longest_nte  = longest_nte, p_value = p_value, ramp = ramp, nte_cov = nte_cov)
  list_of_dataframes <- mclapply(gene_list[2:length(gene_list)], function(x) isoform_expression(riboseq = riboseq, query = x, mitocrypt = mitocrypt, longest_nte = longest_nte), mc.cores = cores)
  MitoCrypt_output <- do.call("rbind", list_of_dataframes)
  return(MitoCrypt_output)
}