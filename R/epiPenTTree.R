#' Penalize a transmission tree based on epidemiological data
#' @param ttree Transmission tree
#' @param epiData A list containing any of the following named elements: \cr \cr 
#' \strong{exposure} a data.frame with row.names matching the case names and columns: 
#' \describe{
#'   \item{\dQuote{start}}{the start time of the probable exposure for this case}
#'   \item{\dQuote{end}}{the end time of the probable exposure for this case}
#' } 
#' \strong{contact} a data.frame with columns: 
#' \describe{
#'   \item{\dQuote{case.A}}{the name of the first case in this probable transmission pair}
#'   \item{\dQuote{case.B}}{the name of the second case in this probable transmission pair}
#'   \item{\dQuote{start}}{the start of the probable contact time for the pair of cases}
#'   \item{\dQuote{end}}{the end of the probable contact time for the pair of cases}
#' }  
#' \strong{location} a named list of character vectors with vector names matching the case names and each vector listing the locations of that case  \cr  
#' @param penaltyInfo Whether or not to return information about which penalties are broken for the given transmission tree 
#' @return Either a vector of the numeric penalties or a list of the numeric penalties and a data.frame of information about which events 
#' caused the penalties
#' @author Matthew Gould
#' @export
epiPenTTree = function(ttree, epiData, penaltyInfo = FALSE){
  exposureP = 1
  contactP = 1
  locationP = 1
  
  # check input format of the epiData
  if(!is(epiData,"list")){
    stop("epiData must be a list, see ?epiPenTTree for format")
  }
  
  if(is.null(epiData$contact) && is.null(epiData$exposure) && is.null(epiData$location)){
    stop("epiData missing or incorrect format, see ?epiPenTTree for format")
  }
  
  if(!is.null(epiData$contact)){
    if(!is(epiData$contact,"data.frame")){
      stop("contact information must be a data.frame, see ?epiPenTTree for format")
    }
  }
  
  if(!is.null(epiData$exposure)){
    if(!is(epiData$exposure,"data.frame")){
      stop("exposure information must be a data.frame, see ?epiPenTTree for format")
    }
  }
    
  if(!is.null(epiData$location)){
    if(!is(epiData$location,"list")){
      stop("location information must be a list, see ?epiPenTTree for format")
    }
  }
  
  if(!all(names(epiData) %in% c("exposure","contact","location"))){
    stop("names in epiData not valid, see ?epiPenTTree for format")
  } 
  
  penWeights <- c(exposureP,contactP,locationP)
  
  penalty <- list("exposurePen" = exposurePen(ttree, epiData, penaltyInfo = penaltyInfo),
                  "contactPen" = contactPen(ttree, epiData, penaltyInfo = penaltyInfo), 
                  "locationPen" = locationPen(ttree, epiData, penaltyInfo = penaltyInfo)) 
  if(!penaltyInfo){
    res <- sapply(1:3, function(i) penalty[[i]]*penWeights[i])
    names(res) <- names(penalty)
    return(res) 
  } else {
    return(list("penalties" = sapply(1:3, function(i){
      penRes <- penalty[[i]]
      if(is.null(nrow(penRes))){
        penRes
      } else {
        nrow(penRes)*penWeights[i]
      }
    }), "penalty.info" = penalty)) 
  }
} 

# if provided case location information in the epiData this returns the number of times transmission between locations occurred 
locationPen <- function(ttree, epiData, penaltyInfo = FALSE){
  # check that there is location info 
  if(!is.null(epiData$location)){
    nsam <- sum(!is.na(ttree$ttree[,2]))
    
    wiw <- cbind(ttree$ttree[,3],1:length(ttree$ttree[,1]))
    
    # remove the index infection 
    wiw <- wiw[wiw[,2] != 0 & wiw[,1] != 0,]
    
    # handle small wiw matrix case 
    if(length(wiw) == 2) wiw <- matrix(wiw, ncol = 2)
    
    # filter to be only between sampled cases 
    sampledCaseRows <- wiw[,2] <= nsam & wiw[,1] <= nsam  
    
    if(any(sampledCaseRows)){
      wiw <- wiw[sampledCaseRows ,]
      
      if(length(wiw) == 2) wiw <- matrix(wiw, ncol = 2)
      
      # extract the locations 
      locationInfectors <- epiData$location[ttree$nam[wiw[,1]]]
      locationInfectee <- epiData$location[ttree$nam[wiw[,2]]]
      
      # penalize all transmission between locations 
      transToPen <- !(sapply(1:length(locationInfectors), function(i){
        if(!is.null(locationInfectee[[i]]) && !is.null(locationInfectors[[i]])){
          return(any(locationInfectee[[i]] %in% locationInfectors[[i]])) 
        } else {
          TRUE
        }
      }))
      
      pen <- sum(transToPen)
      
      if(penaltyInfo && any(transToPen)){
        wiw <- matrix(wiw, ncol = 2)
        rowsInfectees <- wiw[transToPen,2]
        rowsInfectors <- wiw[transToPen,1]
        info <- data.frame(infector = ttree$nam[rowsInfectors], infectee = ttree$nam[rowsInfectees])
        info["infector.locations"] <-  list(locationInfectors[transToPen])
        info["infectee.locations"] <-  list(locationInfectee[transToPen])
        return(info)
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
  if(!is.null(epiData$contact)){
    nsam <- sum(!is.na(ttree$ttree[,2]))
    
    wiw <- cbind(ttree$ttree[,3],1:length(ttree$ttree[,1]))
    
    # remove the index infection 
    wiw <- wiw[wiw[,2] != 0 & wiw[,1] != 0,]
    
    # filter to be only between sampled cases 
    sampledCaseRows <- wiw[,2] <= nsam & wiw[,1] <= nsam 
    
    if(any(sampledCaseRows)){
      wiw <- wiw[sampledCaseRows ,]
      
      epiLength <- length(epiData$contact$case.B)
      timeOfTrans <- vector(length = epiLength)
      presentTrans <- vector(length = epiLength)
      validTimeTrans <-  vector(length = epiLength)
      
      windowStart <- epiData$contact$start
      if(is.null(windowStart)) windowStart <- rep(NA,length(epiData$contact$case.A)) 
      windowEnd <- epiData$contact$end
      if(is.null(windowEnd)) windowEnd <- rep(NA,length(epiData$contact$case.A)) 
      
      pen <- 0
      # loop over the provided contact pairs and perform checks to see if the current ttree matches the epi events 
      for(i in 1:epiLength){
        A <- which(epiData$contact$case.A[i] == ttree$nam, arr.ind = T)
        B <- which(epiData$contact$case.B[i] == ttree$nam, arr.ind = T)
        wiw <- matrix(wiw, ncol = 2)
        if(any(wiw[,1] == A & wiw[,2] == B)){
          presentTrans[i] <- TRUE
          timeOfTrans[i] <- ttree$ttree[B,1]
          if(!is.na(timeOfTrans[i])){
            if(ifelse(is.na(windowStart[i]),FALSE,timeOfTrans[i] < windowStart[i]) || ifelse(is.na(windowEnd),FALSE,timeOfTrans[i] > windowEnd[i])){
              validTimeTrans[i] <- FALSE
              pen <- pen + 1
            } else {
              validTimeTrans[i] <- TRUE
            }
          } 
        } else if(any(wiw[,1] == B & wiw[,2] == A)) {
          presentTrans[i] <- TRUE
          timeOfTrans[i] <- ttree$ttree[A,1]
          if(!is.na(timeOfTrans[i])){
            if(ifelse(is.na(windowStart[i]),FALSE,timeOfTrans[i] < windowStart[i]) || ifelse(is.na(windowEnd),FALSE,timeOfTrans[i] > windowEnd[i])){
              validTimeTrans[i] <- FALSE
              pen <- pen + 1
            } else {
              validTimeTrans[i] <- TRUE
            }
          }
        } else {
          timeOfTrans[i] <- NA
          presentTrans[i] <- FALSE
          validTimeTrans[i] <- ifelse(!is.na(windowStart[i]) && !is.na(windowEnd[i]),FALSE,TRUE) # if missing is valid 
          pen <- pen + 1
        }
      }
      
      if(penaltyInfo){
        penalizedTrans <- !presentTrans | !validTimeTrans
        res <- data.frame(transmission.present = presentTrans[penalizedTrans],
                          case.A = epiData$contact$case.A[penalizedTrans],
                          case.B = epiData$contact$case.B[penalizedTrans],
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
  if(!is.null(epiData$exposure) && !is.null(epiData$exposure$end) && !is.null(epiData$exposure$start)){
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
    samInfectorEndTime <- epiData$exposure[infectors[samInfectors],"end"]
    
    transAfterEnd <- transmissionEvents[samInfectors] > samInfectorEndTime
    transAfterEnd <- ifelse(is.na(transAfterEnd),FALSE,transAfterEnd)
    pen <- pen + sum(transAfterEnd)
    
    # scenario B
    samInfecteesStart <- epiData$exposure[infectees[samInfectees],"start"]
    
    transBeforeStart <- transmissionEvents[samInfectees] < samInfecteesStart
    transBeforeStart <- ifelse(is.na(transBeforeStart),FALSE,transBeforeStart)
    pen <- pen + sum(transBeforeStart)
    
    if(penaltyInfo){
      res <- data.frame(scenario= c(rep("A", sum(transAfterEnd)), rep("B", sum(transBeforeStart))) , 
                        infector= c(infectors[samInfectors][transAfterEnd], infectors[samInfectees][transBeforeStart]),
                        infectee= c(infectees[samInfectors][transAfterEnd], infectees[samInfectees][transBeforeStart]),
                        transmission.time=c(transmissionEvents[samInfectors][transAfterEnd], transmissionEvents[samInfectees][transBeforeStart] ),
                        epi.end.infector = c(epiData$exposure[infectors[samInfectors][transAfterEnd],"end"], epiData$exposure[infectors[samInfectees][transBeforeStart],"end"]), 
                        epi.start.infectee = c(epiData$exposure[infectees[samInfectors][transAfterEnd],"start"], epiData$exposure[infectees[samInfectees][transBeforeStart],"start"]), 
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



