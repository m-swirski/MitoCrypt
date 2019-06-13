#'Mitofates runner
#' 
#' This function submits aminoacid sequence query to loacally installed MitoFates.
#' @param sequence aminoacid sequence character string 
#' @param mitofates_path path to mitofates executable. By default "MitoFates.pl", works when "MitoFates.pl" is in environment $PATH. Perl dependencies must be in environment as well. Rstudio doesn't inherit global environment by default; "printenv | head -n -1 > .Renviron" should solve the issue.
#' @param remove_intermediatory_files logical, if TRUE intermediatory files created in working directory are deleted. Strongly suggested to be set as TRUE, otherwise working directory might become cluttered with several thousands intermediatory files. 
#' @param query_id id of submitted query. Important when using multicore wrapper. Suggested query id is "paste(Gene.Name, NTE.Length, sep = '_')".
#' @param organism either "fungi" or "metazoa" or "plant"
#' @return vector with Mitofates.Cleavage.Site and Mitofates.Probability
#' @keywords mitofates prediction
#' @seealso mitofates_wrapper()
#' @family prediction_software
#' @export
#' @examples
#' 


mitofates_runner <- function(sequence, mitofates_path = 'MitoFates.pl', query_id, remove_intermediatory_files = TRUE, organism = c('fungi', 'metazoa', 'plant')) {
  if (nchar(sequence) >= 20) {
    codename = paste('mitofates_query_code', query_id, sep = '_')
    writeChar(paste('>protein; \n', sequence, sep = ''), codename, eos = NULL)
    system(paste(mitofates_path, " ", codename, " ", organism," | grep protein > ", codename, '.mitofates', sep = ''))
    Mitofates.Cleavage.Site <- system(paste('cat ', codename, '.mitofates | awk \'{print $6}\'', sep = ''), intern = TRUE)
    Mitofates.Cleavage.Site <- sub('\\(.*', '', Mitofates.Cleavage.Site)
    Mitofates.Probability <- system(paste('cat ', codename, '.mitofates | awk \'{print $2}\'', sep = ''), intern = TRUE)
    if (remove_intermediatory_files) {
      system(paste('rm ', codename, sep = ''))
      system(paste('rm ',codename,'.mitofates', sep = ''))
    }
    return(c(Mitofates.Cleavage.Site, Mitofates.Probability))
  } else return(c(NA,NA))
}