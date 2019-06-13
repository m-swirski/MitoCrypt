#'TargetP runner
#' 
#' This function submits aminoacid sequence query to loacally installed TargetP.
#' @param sequence aminoacid sequence character string 
#' @param targetP_path path to targetP executable. By default "targetp", works when "targetp" is in environment $PATH.
#' @param remove_intermediatory_files logical, if TRUE intermediatory files created in working directory are deleted. Strongly suggested to be set as TRUE, otherwise working directory might become cluttered with several thousands intermediatory files. 
#' @param query_id id of submitted query. Important when using multicore wrapper. Suggested query id is "paste(Gene.Name, NTE.Length, sep = '_')".
#' @param organism either "non-plant" (default) or "plant"
#' @return vector with TargetP.Probability, TargetP.Localization and TargetP.Signal.Peptide.Probability
#' @keywords targetp prediction
#' @seealso targetp_wrapper()
#' @family prediction_software
#' @export
#' @examples
#' 



targetp_runner <- function(sequence, targetP_path = "targetp", query_id, remove_intermediatory_files = TRUE, organism = c('non-plant', 'plant')) {
  if (nchar(sequence) >= 20) {
    codename = paste('targetP_query_code', query_id, sep = '_')
    writeChar(paste('>protein; \n', sequence, sep = ''), codename, eos = NULL)
    system(paste(targetP_path, " ",ifelse(organism == 'plant', '-P ', ''), codename, " | grep protein > ", codename, '.targetp', sep = ''))
    TargetP.Probability <- system(paste('cat ', codename, '.targetp | awk \'{print $3}\'', sep = ''), intern = TRUE)
    TargetP.Localization <- system(paste('cat ', codename, '.targetp | awk \'{print $6}\'', sep = ''), intern = TRUE)
    TargetP.Signal.Peptide.Probability <- system(paste('cat ', codename, '.targetp | awk \'{print $4}\'', sep = ''), intern = TRUE)
    if (remove_intermediatory_files) {
      system(paste('rm ', codename, sep = ''))
      system(paste('rm ', codename, '.targetp', sep = ''))}
    return(c(TargetP.Probability, TargetP.Localization, TargetP.Signal.Peptide.Probability))
  } else return(c(NA,NA,NA))
}
