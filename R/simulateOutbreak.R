#' Simulate an outbreak
#' @param off.r First parameter of the negative binomial distribution for offspring number
#' @param off.p Second parameter of the negative binomial distribution for offspring number
#' @param neg the within-host effective population size (Ne) timesgeneration duration (g)
#' @param nSampled number of sampled infected individuals, or NA for any
#' @param pi probability of sampling an infected individual
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation time
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation time 
#' @param ws.shape Shape parameter of the Gamma probability density function representing the sampling time
#' @param ws.scale Scale parameter of the Gamma probability density function representing the sampling time 
#' @param dateStartOutbreak Date when index case becomes infected
#' @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
#' @return Combined phylogenetic and transmission tree
#' @examples
#' plotCTree(simulateOutbreak())
#' @export
simulateOutbreak = function(off.r=1,off.p=0.5,neg=0.25,nSampled=NA,pi=0.5,w.shape=2,w.scale=1,ws.shape=w.shape,ws.scale=w.scale,dateStartOutbreak=2000,dateT=Inf) {
  #Create a transmission tree with nSampled infected sampled individuals
  nsam<-0
  nh<-0
  rejected=-1
  while (is.na(nSampled)||nsam!=nSampled) {
    ttree=NULL
    while (is.null(ttree)) {
      mtt<-makeTTree(off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT-dateStartOutbreak,nSampled)
      rejected=rejected+1
      ttree<-mtt$ttree
      probttree<-mtt$prob
      if (mtt$pruned>0) {
        dateStartOutbreak=dateStartOutbreak+mtt$pruned
        cat(sprintf('Note that simulated outbreak was pruned: in order to have %d sampled by present date %f, the start date was set to %f\n',nSampled,dateT,dateStartOutbreak))
      }
      }
    nsam<-length(which(!is.na(ttree[,2])))
    nh=nrow(ttree)-nsam
    if (is.na(nSampled)) nSampled=nsam
  }
  if (rejected>0) cat(sprintf('Note that rejection sampling was used %d times to simulate outbreak with %d sampled individuals\n',rejected,nSampled))
  n<-nsam+nh
  
  #Create a within-host phylogenetic tree for each infected host
  wtree<-vector('list',n)
  probwithin=0
  for (i in (1:n)) {
    if (is.na(ttree[i,2])) {times<-c(           ttree[which(ttree[,3]==i),1])-ttree[i,1]}
                      else {times<-c(ttree[i,2],ttree[which(ttree[,3]==i),1])-ttree[i,1]}
    a<-withinhost(times,neg)
    wtree[[i]]=a$nodes
    probwithin=probwithin+a$prob
  }
  
  #Glue these trees together
  truth<-.glueTrees(ttree,wtree)
  truth[,1]<-truth[,1]+dateStartOutbreak
  return(list(ctree=truth,nam=mtt$nam,probttree=probttree,probwithin=probwithin))
}  

#' Simulate exposure data for the given ctree 
#' @param ctree Combined phylogenetic/transmission tree
#' @param min The min time difference when determining how close the exposure is to the ctree times 
#' @param max The max time difference when determining how close the exposure is to the ctree times 
#' @param propCases Proportion of cases to generate exporsure data for
#' @return The case.info data.frame containing the simulated exposure times for each case 
#' @author Matthew Gould
#' @export
simulateExposure <- function(ctree, min = 0.1, max = 1, propCases = 0.3){
  
  ttree <- extractTTree(ctree)
  
  nsam <- length(ttree$nam)
  
  numCasesToPick <- ceiling(propCases*nsam)
  if(!(numCasesToPick %in% 1:nsam)){
    stop("incorrect propCases must be a (non-zero) percentage of cases to generate an exporsure window for")
  }
  
  
  # pick a random sample of sampled cases 
  casesToUse <- sample(ttree$nam, numCasesToPick)
  
  
  startTimes <- sapply(casesToUse, function(i){
    row <- which(ttree$nam == i, arr.ind = T)
    ttree$ttree[row,1] - runif(1,min,max)
  })
  
  endTimes <- sapply(casesToUse, function(i){
    row <- which(ttree$nam == i, arr.ind = T)
    
    infectees <- which(ttree$ttree[,3] == row, arr.ind = T)
    infecteesInfT <- ttree$ttree[infectees,1]
    
    max(c(ttree$ttree[row,2],infecteesInfT)) + runif(1,min,max)
  })
  
  return(data.frame(row.names = casesToUse, start = startTimes, end = endTimes)) 
}


#' Simulate contact data for the given ctree
#' @param ctree Combined phylogenetic/transmission tree
#' @param min The min time difference when determining how close the contact times are to the ctree times 
#' @param max The max time difference when determining how close the contact times are to the ctree times 
#' @param propCases Proportion of the valid transmissions to simulate contact data for 
#' @return The contact.info data.frame containing the simulated contact pairs 
#' @author Matthew Gould
#' @export
simulateContact <- function(ctree, min = 0.1, max = 1, propCases = 0.2){
  
  ttree <- extractTTree(ctree)
  
  nsam <- length(ttree$nam)
  
  numCasesToPick <- ceiling(propCases*nsam)
  if(!(numCasesToPick %in% 1:nsam)){
    stop("incorrect propCases must be a (non-zero) percentage of cases to try and simulate contact info for")
  }
  
  infectees <- 1:nrow(ttree$ttree)
  infectors <- ttree$ttree[,3]
  
  zeroInfectors <- infectors == 0
  infectors[zeroInfectors] <- NA
  unsampled <- is.na(ttree$nam[infectors]) | is.na(ttree$nam[infectees])
  
  infectors <- infectors[!unsampled]
  infectees <- infectees[!unsampled]
  
  numCasesToPick <- ifelse(numCasesToPick>length(infectors),length(infectors),numCasesToPick) 
  
  numToUse <- sample(1:length(infectors),numCasesToPick)
  
  infectors <- infectors[numToUse]
  infectees <- infectees[numToUse]
  
  infectTimes <- ttree$ttree[infectees,1]
  
  data.frame(case.A = ttree$nam[infectors], case.B = ttree$nam[infectees], start = infectTimes-runif(length(infectors), min = min, max = max),
                     end = infectTimes+runif(length(infectors), min = min, max = max)) # epi end time should be after last infection
}





