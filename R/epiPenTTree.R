#' Penalize a transmission tree based on epidemiological data
#' @param ttree Transmission tree
#' @param epiData A named list of two data.frames: \cr \cr 
#' \strong{case.info} with row.names matching the case names and columns: 
#' \describe{
#'   \item{\dQuote{start}}{the start of the probable exposure for this case}
#'   \item{\dQuote{end}}{the end of the probable exposure for this case}
#'   \item{\dQuote{location}}{a string name of the location for this case}
#' } 
#' \strong{contact.info} with columns: 
#' \describe{
#'   \item{\dQuote{case.A}}{the name of the first case in this probable transmission pair}
#'   \item{\dQuote{case.B}}{the name of the second case in this probable transmission pair}
#'   \item{\dQuote{start}}{the start of the probable contact for the pair of cases}
#'   \item{\dQuote{end}}{the end of the probable contact for the pair of cases}
#' }  
#' @param penaltyInfo Whether or not to return information about which penalties are broken for the given transmission tree 
#' @return Either a list of the numeric penalties or a list of the numeric penalties and data.frame of information about which events 
#' caused the penalties
#' @author Matthew Gould 
#' @export
epiPenTTree = function(ttree, epiData, penaltyInfo = FALSE){
  # check input format of the epiData 
  penalty <- list("exposurePen" = exposurePen(ttree, epiData, penaltyInfo = penaltyInfo),
                  "contactPen" = contactPen(ttree, epiData, penaltyInfo = penaltyInfo), 
                  "locationPen" = locationPen(ttree, epiData, penaltyInfo = penaltyInfo)) 
  if(!penaltyInfo){
    return(penalty) 
  } else {
    return(list("penalties" = sapply(penalty, function(i){
      if(is.null(nrow(i))){
        i
      } else {
        nrow(i)
      }
    }), "penaltyInfo" = penalty)) 
  }
} 

# if provided case location information in the epiData this returns the number of times transmission between locations occurred 
locationPen <- function(ttree, epiData, penaltyInfo = FALSE){
  # check that there is case.info and check that there is location info 
  if(!is.null(epiData$case.info) && !is.null(epiData$case.info$location)){
    nsam <- sum(!is.na(ttree$ttree[,2]))
    
    wiw <- cbind(ttree$ttree[,3],1:length(ttree$ttree[,1]))
    
    # remove the index infection 
    wiw <- wiw[wiw[,2] != 0 & wiw[,1] != 0,]
    
    # filter to be only between sampled cases 
    sampledCaseRows <- wiw[,2] <= nsam & wiw[,1] <= nsam 
    
    if(any(sampledCaseRows)){
      wiw <- wiw[sampledCaseRows ,]
      
      # extract the locations 
      if(length(wiw) == 2){
        locationInfectors <- epiData$case.info[ttree$nam[wiw],"location"] 
        locationInfectee <- epiData$case.info[ttree$nam[wiw],"location"]
        
      } else {
        locationInfectors <- epiData$case.info[ttree$nam[wiw[,1]],"location"] 
        locationInfectee <- epiData$case.info[ttree$nam[wiw[,2]],"location"]
      }
      
      # penalize all transmission between locations 
      transToPen <- ifelse(is.na(locationInfectee != locationInfectors), FALSE, locationInfectee != locationInfectors)
      pen <- sum(transToPen)
      
      if(penaltyInfo && any(transToPen)){
        wiw <- matrix(wiw, ncol = 2)
        rowsInfectees <- wiw[transToPen,2]
        rowsInfectors <- wiw[transToPen,1]
        return(data.frame(infector = ttree$nam[rowsInfectors], infectee = ttree$nam[rowsInfectees], infector.location = locationInfectors[transToPen], infectee.location = locationInfectee[transToPen]))
      } else {
        return(pen)
      }
    } 
  } 
  return(0)
}

# if provided contact information in the epiData this returns the number of missing transmissions between contact pairs 
# and if a window of contact is provided this also returns the number of times the transmission occured outside that window 
contactPen <- function(ttree, epiData, penaltyInfo = FALSE){
  if(!is.null(epiData$contact.info)){
    nsam <- sum(!is.na(ttree$ttree[,2]))
    
    wiw <- cbind(ttree$ttree[,3],1:length(ttree$ttree[,1]))
    
    # remove the index infection 
    wiw <- wiw[wiw[,2] != 0 & wiw[,1] != 0,]
    
    # filter to be only between sampled cases 
    sampledCaseRows <- wiw[,2] <= nsam & wiw[,1] <= nsam 
    
    if(any(sampledCaseRows)){
      wiw <- wiw[sampledCaseRows ,]
      
      epiLength <- length(epiData$contact.info$case.B)
      timeOfTrans <- vector(length = epiLength)
      presentTrans <- vector(length = epiLength)
      validTimeTrans <-  vector(length = epiLength)
      
      windowStart <- epiData$contact.info$start
      if(is.null(windowStart)) windowStart <- rep(NA,length(epiData$contact.info$case.A)) 
      windowEnd <- epiData$contact.info$end
      if(is.null(windowEnd)) windowEnd <- rep(NA,length(epiData$contact.info$case.A)) 
      
      pen <- 0
      # loop over the provided contact pairs and perform checks to see if the current ttree matches the epi events 
      for(i in 1:epiLength){
        A <- which(epiData$contact.info$case.A[i] == ttree$nam, arr.ind = T)
        B <- which(epiData$contact.info$case.B[i] == ttree$nam, arr.ind = T)
        wiw <- matrix(wiw, ncol = 2)
        if(any(wiw[,1] == A & wiw[,2] == B)){
          presentTrans[i] <- TRUE
          timeOfTrans[i] <- ttree$ttree[B,1]
          if(!is.na(windowStart[i]) && !is.na(timeOfTrans[i])){
            if(timeOfTrans[i] < windowStart[i] || timeOfTrans[i] > windowEnd[i]){
              validTimeTrans[i] <- FALSE
              pen <- pen + 1
            } else {
              validTimeTrans[i] <- TRUE
            }
          } else {
            # no given times so is valid 
            validTimeTrans[i] <- TRUE
          }
        } else if(any(wiw[,1] == B & wiw[,2] == A)) {
          presentTrans[i] <- TRUE
          timeOfTrans[i] <- ttree$ttree[A,1]
          if(!is.na(windowStart[i]) && !is.na(timeOfTrans[i])){
            if(timeOfTrans[i] < windowStart[i] || timeOfTrans[i] > windowEnd[i]){
              validTimeTrans[i] <- FALSE
              pen <- pen + 1
            } else {
              validTimeTrans[i] <- TRUE
            }
          } else {
            validTimeTrans[i] <- TRUE
          }
        } else {
          timeOfTrans[i] <- NA
          presentTrans[i] <- FALSE
          validTimeTrans[i] <- ifelse(!is.na(windowStart[i]) && !is.na(windowEnd[i]),FALSE,TRUE) # if missing is valid 
          pen <- pen + 1
        }
      }
      
      if(penaltyInfo && any(presentTrans)){
        penalizedTrans <- !presentTrans | !validTimeTrans
        res <- data.frame(transmission.present = presentTrans[penalizedTrans],
                          case.A = epiData$contact.info$case.A[penalizedTrans],
                          case.B = epiData$contact.info$case.B[penalizedTrans],
                          epi.start = windowStart[penalizedTrans],
                          epi.end = windowEnd[penalizedTrans],
                          ttree.transmission = timeOfTrans[penalizedTrans],
                          stringsAsFactors = F)
        if(nrow(res) != 0){
          return(res)
        } else {
          return(0)
        }
      } else {
        return(pen)
      }
    } 
  } 
  return(0)
}

# if provided case exposure information in the epiData this returns the number of events in the ttree that disagree with the exposure time in 
# two ways: 
# A) transmission occurred after the exposure end time of the infector 
# B) transmission occurred before the exposure start time of the infectee
exposurePen <- function(ttree, epiData, penaltyInfo = FALSE){
  if(!is.null(epiData$case.info) && !is.null(epiData$case.info$end) && !is.null(epiData$case.info$start)){
    pen <- 0
    nSam <- length(ttree$nam)
    
    infectors <- ttree$ttree[,3]
    infectees <- 1:nrow(ttree$ttree)
    
    infectors[infectors == 0] <- nSam + 1 
    samInfectors <- infectors <= nSam 
    infectors <- ifelse(samInfectors,ttree$nam[infectors],"unsampled")
    
    samInfectees <- infectees <= nSam
    infectees <- ifelse(samInfectees,ttree$nam[infectees],"unsampled")
    
    transmissionEvents <- ttree$ttree[,1]
    
    # scenario A
    samInfectorEndTime <- epiData$case.info[infectors[samInfectors],"end"]
    
    transAfterEnd <- transmissionEvents[samInfectors] > samInfectorEndTime
    transAfterEnd <- ifelse(is.na(transAfterEnd),FALSE,transAfterEnd)
    pen <- pen + sum(transAfterEnd)
    
    # scenario B
    samInfecteesStart <- epiData$case.info[infectees[samInfectees],"start"]
    
    transBeforeStart <- transmissionEvents[samInfectees] < samInfecteesStart
    transBeforeStart <- ifelse(is.na(transBeforeStart),FALSE,transBeforeStart)
    pen <- pen + sum(transBeforeStart)
    
    if(penaltyInfo){
      res <- data.frame(scenario= c(rep("A", sum(transAfterEnd)), rep("B", sum(transBeforeStart))) , 
                        infector= c(infectors[samInfectors][transAfterEnd], infectors[samInfectees][transBeforeStart]),
                        infectee= c(infectees[samInfectors][transAfterEnd], infectees[samInfectees][transBeforeStart]),
                        transmission.time=c(transmissionEvents[samInfectors][transAfterEnd], transmissionEvents[samInfectees][transBeforeStart] ),
                        epi.end.infector = c(epiData$case.info[infectors[samInfectors][transAfterEnd],"end"], epiData$case.info[infectors[samInfectees][transBeforeStart],"end"]), 
                        epi.start.infectee = c(epiData$case.info[infectees[samInfectors][transAfterEnd],"start"], epiData$case.info[infectees[samInfectees][transBeforeStart],"start"]), 
                        stringsAsFactors = F)
      if(nrow(res) != 0){
        return(res)
      } else {
        return(0)
      }
    } else {
      return(pen)
    }
  }
  return(0)
}



