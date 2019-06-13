#' Fetch gene sequence based on gtf annotation and genome sequence
#'
#' This functions is wrapper of read.gtf function from refGenome.
#' @param gtf_table gtf annotation dataframe
#' @param genome BSgenome name, compatible with gtf
#' @param gene_id target gene 
#' @param output output type. either "DNAstring" biostrings object, indexes, or both as a list.
#' @param fetch_feature gene feature to be fetched. "CDS" by default.
#' @keywords sequence gtf
#' @export
#' @examples
#' fetch_sequence()



fetch_sequence <- function(gtf_table, genome, transcript_id, output = c("DNAString", "indexes", "both"), fetch_feature = "CDS") {
  indexes <- c()
  cur_gtf <- gtf_table[which(gtf_table$transcript_id == transcript_id),]
  cur_gtf <- cur_gtf[which(cur_gtf$feature == fetch_feature),]
  cur_gtf <- arrange(cur_gtf, start)
  for (i in 1:nrow(cur_gtf)) {
    indexes <- c(indexes, cur_gtf$start[i]:cur_gtf$end[i])
  }
  if (cur_gtf$strand[1] == '+') {
    output_sequence <- genome[[cur_gtf$seqid[1]]][indexes] 
  } else {
    output_sequence <- reverseComplement(genome[[cur_gtf$seqid[1]]][indexes])
  }
  if (output[1] == "DNAString") {
    return(output_sequence)
  } else if (output[1] == "indexes") {
    return(indexes)
  } else return(list(output_sequence, indexes))
}
