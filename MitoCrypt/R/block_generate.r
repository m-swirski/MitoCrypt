#' Block Generate
#'
#' This function takes gtf and returns .bed compatible blocks
#' @param gtf_table gtf annotation dataframe
#' @param transcript_id target gene
#' @param fetch_feature gene feature to be fetched. "CDS" by default.
#' @return list of blockCount, blockSizes and blockStarts. See UCSC .bed file documentation for details.
#' @keywords bed block
#' @export
#' @examples
#' 


block_generate <- function(gtf_table, transcript_id, fetch_feature = "CDS"){
  cur_gtf <- gtf_table[which(gtf_table$transcript_id == transcript_id),]
  cur_gtf <- cur_gtf[which(cur_gtf$feature == fetch_feature),]
  Count <- nrow(cur_gtf)
  cur_gtf <- arrange(cur_gtf, start)
  Sizes <- integer()
  Starts <- integer()
  for (i in 1:nrow(cur_gtf)){
    Sizes[i] <- cur_gtf$end[i] - cur_gtf$start[i] + 1
    Starts[i] <- cur_gtf$start[i] - cur_gtf$start[1] + 1
  }
  return(list(blockCount = Count, blockSizes = Sizes, blockStarts = Starts))
}