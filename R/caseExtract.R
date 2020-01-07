#' Extract various information about one case throughtout an MCMC output
#' 
#' @param record MCMC output produced by inferTTree
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @param case String matching one of the case names in the data
#' @return Data frame with rows over the record iterations and columns with the following data: 
#' \describe{
#'   \item{\dQuote{is_index}}{whether the case is the index case in this iteration}
#'   \item{\dQuote{intermediates_from_start}}{the number of intermediates in the transmission chain from this case to the index}
#'   \item{\dQuote{infector}}{the name of the infector of this case if sampled otherwise unsampled}
#'   \item{\dQuote{genTime_infector}}{the time from the infection of the infector to the time of infection of the target case}
#'   \item{\dQuote{time_Infected}}{time the target case was infected}
#'   \item{\dQuote{time_sampled}}{time the target case was sampled}
#'    \item{\dQuote{transmissions_to_sampled}}{the number of transmissions to sampled cases}
#'   \item{\dQuote{transmissions_to_unsampled}}{the number of transmissions to unsampled cases}
#'   \item{\dQuote{posterior}}{the posterior probability of the transmission tree at this iteration}
#' } 
#' @author Matthew Gould 
#' @export
recordCaseExtract <- function(record,case,burnin=0.5){
  record=record[max(1,round(length(record)*burnin)):length(record)] 
  
  caseInfo = t(sapply(1:length(record), function(x) {
    tt=extractTTree(record[[x]]$ctree); 
    
    caseRow <- which(tt$nam == case)
    
    indexRow <- which(tt$ttree[,3] == 0)
    
    is_index <- ifelse(caseRow == indexRow,TRUE,FALSE)
    
    itermediates_from_start <- numInter(indexRow,caseRow,tt$ttree) 
    
    rowInfector <- tt$ttree[caseRow,3]
    infector <- ifelse(is.na(tt$nam[rowInfector]) || is_index,"unsampled",tt$nam[rowInfector])
    
    time_Infected <- tt$ttree[caseRow,1]
    
    genTime_infector <- ifelse(is_index,NA,time_Infected - tt$ttree[rowInfector,1]) 
    
    time_sampled <- tt$ttree[caseRow,2]
    
    rowsTrans <- which(tt$ttree[,3] == caseRow)
    
    transmissions_to_sampled <- sum(ifelse(is.na(tt$ttree[rowsTrans,2]),FALSE,TRUE)) 
    
    transmissions_to_unsampled <- sum(ifelse(is.na(tt$ttree[rowsTrans,2]),TRUE,FALSE)) 
    
    
    return(c(is_index, itermediates_from_start, infector, genTime_infector, time_Infected,
             time_sampled, transmissions_to_sampled, transmissions_to_unsampled, record[[x]]$posterior))
    
  }))
  
  colnames(caseInfo) <- c("is_index", "intermediates_from_start", "infector", "genTime_infector", "time_Infected",
                          "time_sampled", "transmissions_to_sampled", "transmissions_to_unsampled","posterior")
  
  return(caseInfo)
  
}

numInter <- function(hostA,hostB,ttree){
  count=0
  while (hostA!=hostB) {
    if (ttree[hostA,1]>ttree[hostB,1]) hostA=ttree[hostA,3] else hostB=ttree[hostB,3]
    count=count+1
  }
  return(count)
}
