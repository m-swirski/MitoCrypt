#' A read gtf function wrapper
#'
#' This functions is wrapper of read.gtf function from refGenome.
#' @param file gtf file path.
#' @keywords gtf
#' @export
#' @examples
#' read_gtf_as_data_frame()

read_gtf_as_data_frame <- function(file) {
  ensemble <- ensemblGenome()
  read.gtf(ensemble, file)
  gtf_table <- ensemble@ev$gtf
  if (!("gene_name" %in% colnames(gtf_table))) gtf_table$gene_name <- NA
  gtf_table <- gtf_table %>% mutate(gene_name = ifelse(is.na(gene_name),transcript_id, gene_name))
  return(gtf_table)
}