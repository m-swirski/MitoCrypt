#' Mitoprot runner
#' 
#' This function submits aminoacid sequence query to loacally installed mitoprot.
#' @param sequence aminoacid sequence character string 
#' @param mitoprot_path path to mitoprot executable. By default "mitoprot", works when "mitoprot" is in environment $PATH.
#' @param query_id id of submitted query. Important when using multicore wrapper. Suggested query id is "paste(Gene.Name, NTE.Length, sep = '_')".
#' @param Cleavage_site logical, if TRUE cleavage site is predicted, otherwise step is skipped.
#' @param Mitoprot_probability logical, if TRUE mitoprot probability is calculated, otherwise step is skipped.
#' @param remove_intermediatory_files logical, if TRUE intermediatory files created in working directory are deleted. Strongly suggested to be set as TRUE, otherwise working directory might become cluttered with several thousands intermediatory files. 
#' @return vector c(CleavageSite, MitoprotProbability)
#' @seealso mitoprot_wrapper()
#' @family prediction_software
#' @keywords mitoprot prediction 
#' @export
#' @examples
#' 

mitoprot_runner <- function(sequence, mitoprot_path = "mitoprot", query_id,Cleavage_site = TRUE, Mitoprot_probability = TRUE,remove_intermediatory_files = TRUE) {
  if (nchar(sequence) >= 20) {
    codename = paste('mitoprot_query_code',query_id,sep = "_")
    writeChar(sequence, codename, eos = NULL)
    system(paste(mitoprot_path ,' -f a ',codename, sep =""), intern = TRUE, ignore.stdout = TRUE)
    output <- c()
    if (Cleavage_site) output <- system(paste('grep CleavSite ',codename,'.mitoprot | sed \'s/CleavSite *: *//\'', sep = ''),intern = TRUE)
    if (Mitoprot_probability) output  <- c(output, system(paste('grep \'DFM \' ', codename, '.mitoprot | sed \'s/DFM  *: * [[:digit:]].[[:digit:]]* *//\'', sep = ''), intern = TRUE))
    if (remove_intermediatory_files) {system(paste('rm ',codename,'.mitoprot', sep = ''))
      system(paste('rm ', codename, sep = ''))
    }
    return(output)
  } else return(c(NA,NA))
}
