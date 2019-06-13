#' Match multiple codons 
#'
#' This functions is wrapper of matchPattern from biostrings package. Allows for multiple 3-nt long sequences matching.
#' @param codons character vector of searched codons.
#' @param Seq biostrings DNAstring object serving as template for search.
#' @return vector of codon start locations, named with codon. 
#' @keywords codon match
#' @export
#' @examples
#' 





match_multiple_codons <- function(codons, Seq) {
  matches <- c()
  for (codon in codons) {
    cur_matches <- matchPattern(codon, Seq)
    cur_matches <- cur_matches@ranges@start
    names(cur_matches) <- rep(codon, length(cur_matches))
    matches <- c(matches, cur_matches)
  }
  return(matches)
} 
