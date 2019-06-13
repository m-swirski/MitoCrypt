#'TargetP wrapper
#' 
#' This function is multicore wrapper for targetp_runner(). 
#' @param MitoCrypt_table MitoCrypt table dataframe
#' @param targetP_path path to targetP executable. By default "targetp", works when "targetp" is in environment $PATH.
#' @param remove_intermediatory_files logical, if TRUE intermediatory files created in working directory are deleted. Strongly suggested to be set as TRUE, otherwise working directory might become cluttered with several thousands intermediatory files. 
#' @param organism either "non-plant" (default) or "plant"
#' @param cores number of cores to use. Passes argument to mcmapply() mc.cores argument. 
#' @return vector with TargetP.Probability, TargetP.Localization and TargetP.Signal.Peptide.Probability
#' @keywords targetp parallel
#' @seealso targetp_runner()
#' @family prediction_software MitoCrypt
#' @export
#' @examples
#' 

targetp_wrapper <- function(MitoCrypt_table, targetP_path = "targetp", remove_intermediatory_files = TRUE, cores = 1,organism = c('non-plant', 'plant')) {
  ids <- mcmapply(function(x,y) paste(x,y, sep = '_'), as.character(MitoCrypt_table$Gene.Name), as.character(MitoCrypt_table$NTE.Length), mc.cores = cores)
  output <- mcmapply(function(x,y) targetp_runner(x, targetP_path, y, organism = organism), MitoCrypt_table$Isoform.Sequence, ids,mc.cores = cores, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  output <- t(as.data.frame(output))
  MitoCrypt_table$TargetP.Probability <- output[,1]
  MitoCrypt_table$TargetP.Localization <- output[,2]
  MitoCrypt_table$TargetP.Signal.Peptide.Probability <- output[,3]
  return(output)
}