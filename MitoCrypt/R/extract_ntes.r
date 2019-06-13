#' Extract NTEs
#'
#' This function eliminates potential NTEs with ribosome profile coverage below given threshold
#' @param min_mean coverage threshold in RPKM
#' @param profile scikit-ribo A-site profile for processed gene
#' @param ntes ascending vector of NTE isoforms lengths.
#' @return filtered NTEs vector
#' @keywords NTE filter
#' @export
#' @examples
#' 


extract_ntes <- function(min_mean, profile, ntes_input, max_nte){
  ntes_proc <- rev(abs(ntes_input - max(ntes_input,max_nte)))
  excluded <- c()
  excluded_counter <- 1
  if (length(ntes_proc) > 1) {
    for (i in 1:(length(ntes_proc)-1)){
      if (abs(ntes_proc[i] - ntes_proc[i + 1]) <= 3){
        if (profile$ribosome_count[ntes_proc[i] + 1] >= profile$ribosome_count[ntes_proc[i+1] + 1]){
          excluded[excluded_counter] <- ntes_proc[i + 1]
        } else {excluded[excluded_counter] <- ntes_proc[i]}
        excluded_counter <- excluded_counter + 1
      }
      
    }
  }
  ntes_proc <- ntes_proc[which(!(ntes_proc %in% excluded))]
  ntes_proc[length(ntes_proc) + 1] <- ntes_proc[length(ntes_proc)] + min(c(51,nrow(profile) - ntes_proc[length(ntes_proc)]))
  clean <- 1
  turn <- 1
  while (clean != 0 & length(ntes_proc) > 1) {
    turn <- turn + 1
    clean <- 0
    excluded <- c()
    excluded_counter <- 1
    for (k in 1:(length(ntes_proc) - 1)) {
      cur_mean <- mean(profile$ribosome_count[(ntes_proc[k] + 1) : (ntes_proc[k+1])])
      cur_median <- median(profile$ribosome_count[(ntes_proc[k] + 1) : (ntes_proc[k+1])])
      if ((cur_mean < min_mean) | (cur_median == 0)) {
        excluded[excluded_counter] <- ntes_proc[k]
        excluded_counter <- excluded_counter + 1
        clean <- clean + 1
      }
    }
    ntes_proc <- ntes_proc[which(!(ntes_proc %in% excluded))]
  }
  ntes_proc <- ntes_proc[-length(ntes_proc)]
  return(ntes_proc)
}