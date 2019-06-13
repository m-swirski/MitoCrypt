#' Mitoprot wrapper
#'
#' This function is multicore wrapper for mitoprot_runner(). 
#' @param MitoCrypt_table MitoCrypt table dataframe
#' @param mitoprot_path path to mitoprot executable. By default "mitoprot", works when "mitoprot" is in environment $PATH.
#' @param remove_intermediatory_files logical, if TRUE intermediatory files created in working directory are deleted. Strongly suggested to be set as TRUE, otherwise working directory might become cluttered with several thousands intermediatory files. 
#' @param cores number of cores to use. Passes argument to mcmapply() mc.cores argument. 
#' @return MitoCrypt_table with two added columns: MitoProt.Cleavage.Site and MitoProt.Probability
#' @seealso mitoprot_runner() mcmapply()
#' @family prediction_software MitoCrypt
#' @keywords mitoprot parallel
#' @export
#' @examples
#' 

mitoprot_wrapper <- function(MitoCrypt_table,mitoprot_path = "mitoprot", remove_intermediatory_files = TRUE, cores = 1) {
  ids <- mcmapply(function(x,y) paste(x,y, sep = '_'), as.character(MitoCrypt_table$Gene.Name), as.character(MitoCrypt_table$NTE.Length), mc.cores = cores)
  output <- mcmapply(function(x,y) mitoprot_runner(x, mitoprot_path, y), MitoCrypt_table$Isoform.Sequence, ids,mc.cores = cores, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  output <- t(as.data.frame(output))
  MitoCrypt_table$MitoProt.Cleavage.Site <- output[,1]
  MitoCrypt_table$MitoProt.Probability <- output[,2]
  return(MitoCrypt_table)
}