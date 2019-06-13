#' Context tolower
#'
#' This function takes start character string codon context and returns character string where only start codon is displayed by capital letters.
#' @param context_sequence codon context character string object 
#' @param context_span two element vector, corresponding to start and end of the sequence relative to first position of start codon 
#' @return start codon context character string
#' @keywords context start codon
#' @export
#' @examples context_tolower('AAAAAAATGG', c(-6,5)) 
#'

context_tolower <- function(context_sequence, context_span) {
  if (abs(context_span[1]) > 0) substr(context_sequence,1,abs(context_span[1])) <- tolower(substr(context_sequence,1,abs(context_span[1])))
  if (abs(context_span[1]) + abs(context_span[2]) > abs(context_span[1]) + 3) substr(context_sequence,abs(context_span[1]) + 4, abs(context_span[1]) + abs(context_span[2])) <- tolower(substr(context_sequence,abs(context_span[1]) + 4, abs(context_span[1]) + abs(context_span[2])))
  return(context_sequence)
}