#' Get start/end modified isoform gtf 
#'
#' This function returns gtf entry modified with start or end location, taking into account gene structure (introns)
#' @param start new start of isoform
#' @param end new end of isoform 
#' @param gtf_table gtf annotation dataframe
#' @param transcript_id target gene
#' @param fetch_sequence gene feature to be fetched. "CDS" by default.
#' @return modified gtf entry
#' @keywords
#' @export
#' @examples
#' 





get_isoform_gtf <- function(start,end,gtf_table, transcript_id, fetch_feature = "CDS") {
  cur_gtf <- gtf_table[which(gtf_table$transcript_id == transcript_id),]
  cur_gtf <- cur_gtf[which(cur_gtf$feature == fetch_feature),]
  cur_gtf <- arrange(cur_gtf, start)
  if (cur_gtf$strand[1] == "+") {
    if(start <= cur_gtf$start[1]){
      cur_gtf$start[1] <- start
      iso_gtf <- cur_gtf
    }else {
      stop <- FALSE
      exon_count <- 1
      while (stop == FALSE) {
        if (start %in% cur_gtf$start[exon_count]:cur_gtf$end[exon_count]) {
          stop <- TRUE
          iso_gtf <- cur_gtf[exon_count:nrow(cur_gtf),]
          iso_gtf$start[1] <- start
        } else {exon_count <- exon_count + 1}
      }
    }
  } else {
    if (end >= cur_gtf$end[nrow(cur_gtf)]) {
      cur_gtf$end[nrow(cur_gtf)] <- end
      iso_gtf <- cur_gtf
    } else {
      stop <- FALSE
      exon_count = nrow(cur_gtf)
      while (stop == FALSE) {
        if (end %in% cur_gtf$start[exon_count]:cur_gtf$end[exon_count]) {
          stop <- TRUE
          iso_gtf <- cur_gtf[1:exon_count,]
          iso_gtf$end[nrow(iso_gtf)] <- end
        } else {exon_count <- exon_count - 1}
      }
    }
  }
  return(iso_gtf)
}