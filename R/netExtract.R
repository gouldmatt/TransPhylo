#' Extract a node and edge list from a transmission tree 
#' 
#' @param ttree 
#' @param colorMethod method to color the cases: 
#' \describe{
#'   \item{\dQuote{all}}{give each case a unique color}
#'   \item{\dQuote{sampled}}{give the sampled cases a unique color and color the unsampled cases grey}
#'   \item{\dQuote{time}}{give each case a color based on the time of infection}
#' }
#' @return a list of two lists one with the node list and one with the edge list 
#' @author Matthew Gould 
#' @note The format of the node and edge list is ready to be used with the \strong{visnetwork} package
#' @export
netExtractTtree <- function(ttree, colorMethod = c("all","sampled","time")){
  if(!missing(colorMethod) & length(colorMethod)>1) stop("only one 'colorMethod' allowed")
  colorMethod <- match.arg(colorMethod)
  
  nam <- ttree$nam
  tree <- ttree$ttree
  
  numHosts <- nrow(tree) 
  
  infectors <- tree[,3]
  infectees <- 1:numHosts
  
  # remove when index was infected 
  indexInf <- infectors == 0
  infectors <- infectors[!indexInf]
  infectees <- infectees[!indexInf]
  
  edges <- data.frame(from = infectors, to = infectees, arrows="to")
  
  hostNames <- ifelse(is.na(nam[1:numHosts]),"unsampled",nam)
  samOrUnsam <- ifelse(hostNames == "unsampled","unsampled","sampled")
  
  if(colorMethod == "time"){
    colfunc <- colorRampPalette(c("#00BFC0", "#FB737A")) 
    cols <- colfunc(numHosts)[numHosts-rank(tree[,1])+1]
    
  } else if(colorMethod == "all"){
    cols <- hcl(h = seq(10, 375, length = numHosts + 1), l = 65, c = 100)[1:numHosts]
  } else if(colorMethod == "sampled"){
    numSam <- sum(hostNames != "unsampled")
    cols <- hcl(h = seq(10, 375, length = numSam + 1), l = 65, c = 100)[1:numSam]
    cols <- ifelse(hostNames != "unsampled",cols,"grey")
  }
  
  nodes <- data.frame(id = 1:numHosts, color = cols, label = hostNames, title = tree[,1], sampledTime = tree[,2], group = samOrUnsam, stringsAsFactors = F) 
  
  list("nodes" = nodes, "edges" = edges)
} 

#' Extract a node and edge list from a MCMC output
#' 
#' @param record MCMC output produced by inferTTree
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @param colorMethod method to color the cases: 
#' \describe{
#'   \item{\dQuote{all}}{give each case a unique color}
#'   \item{\dQuote{sampled}}{give the sampled cases a unique color and color the unsampled cases grey}
#'   \item{\dQuote{time}}{give each case a color based on the time of infection}
#' }
#' @param baseTree The base tree to use: 
#' \describe{
#'   \item{\dQuote{maxPosterior}}{the transmission tree with the highest posterior probability}
#'   \item{\dQuote{consensus}}{the consensus transmission tree}
#' }
#' @param minimum Minimum probability for inclusion of a partition in the consensus, ignored if not using the consensus \strong{baseTree}
#' @param wex A weight expansion factor to multiply the edge weights by 
#' @return a list of two lists one with the node list and one with the edge list 
#' @author Matthew Gould 
#' @note The format of the node and edge list is ready to be used with the \strong{visnetwork} package
#' @export
netExtractRecord <- function(record, burnin = 0, colorMethod = c("all","sampled","time"), baseTree = c("maxPosterior","consensus"), minimum = 0.2, wex = 1) {
  if(!missing(baseTree) & length(baseTree)>1) stop("only one 'baseTree' allowed")
  baseTree <- match.arg(baseTree)
  
  if(!missing(colorMethod) & length(colorMethod)>1) stop("only one 'colorMethod' allowed")
  colorMethod <- match.arg(colorMethod)
  
  # burnin
  record=record[max(1,round(length(record)*burnin)):length(record)]
  
  ps <- sapply(record,function(x) x$posterior)
  
  # extract network data from correct base tree 
  if(baseTree == "maxPosterior"){
    ttree <- extractTTree(record[[which.max(ps)]]$ctree)
    netData = netExtractTtree(ttree, colorMethod = colorMethod)
  } else if(baseTree == "consensus"){
    ttree <- consTTree(record, burnin = 0, minimum = minimum)
    netData = netExtractTtree(ttree, colorMethod = colorMethod)
  }
  
  numCases <- length(ttree$nam)
  
  # extract direct transmission weights 
  wiw=computeMatWIW(record)
  
  # compute weight for each edge in the edgelist, only using wiw for those transmissions from sampled to sampled cases
  wiwInd = which(netData$edges$from <= numCases & netData$edges$to <= numCases & netData$edges$from>0)
  
  ws=vapply(1:length(wiwInd), function(x) wiw[netData$edges$from[wiwInd[x]], netData$edges$to[wiwInd[x]]], FUN.VALUE = 1)
  
  # add the weights for direct trans
  netData$edges$width <- 1
  netData$edges$width[wiwInd] <- 1 + wex*ws
  
  return(netData)
}



