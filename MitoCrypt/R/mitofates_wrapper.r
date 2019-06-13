#'Mitofates wrapper
#' 
#' This function is multicore wrapper for mitofates_runner().
#' @param MitoCrypt_table MitoCrypt table dataframe
#' @param mitofates_path path to mitofates executable. By default "MitoFates.pl", works when "MitoFates.pl" is in environment $PATH. Perl dependencies must be in environment as well. Rstudio doesn't inherit global environment by default; "printenv | head -n -1 > .Renviron" should solve the issue.
#' @param remove_intermediatory_files logical, if TRUE intermediatory files created in working directory are deleted. Strongly suggested to be set as TRUE, otherwise working directory might become cluttered with several thousands intermediatory files. 
#' @param organism either "fungi" or "metazoa" or "plant"
#' @param cores number of cores to use. Passes argument to mcmapply() mc.cores argument. 
#' @return vector with Mitofates.Cleavage.Site and Mitofates.Probability
#' @keywords mitofates parallel
#' @seealso mitofates_runner()
#' @family prediction_software MitoCrypt
#' @export
#' @examples
#' 

mitofates_wrapper <- function(MitoCrypt_table, mitofates_path = "MitoFates.pl", remove_intermediatory_files = TRUE, cores = 1, organism = c('fungi', 'metazoa', 'plant')) {
  ids <- mcmapply(function(x,y) paste(x,y, sep = '_'), as.character(MitoCrypt_table$Gene.Name), as.character(MitoCrypt_table$NTE.Length), mc.cores = cores)
  output <- mcmapply(function(x,y) mitofates_runner(x, mitofates_path, y, organism = organism), MitoCrypt_table$Isoform.Sequence, ids,mc.cores = cores, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  output <- t(as.data.frame(output))
  MitoCrypt_table$MitoFates.Cleavage.Site <- output[,1]
  MitoCrypt_table$MitoFates.Probability <- output[,2]
  return(MitoCrypt_table)
}