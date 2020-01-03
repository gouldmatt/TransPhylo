#' Plot both phylogenetic and transmission trees using colors on the phylogeny
#' @param tree Combined phylogenetic/transmission tree
#' @param maxTime Maximum time to show on the x axis
#' @param showLabels Whether or not to show the labels 
#' @param alignLabels Whether to align the labels of the tree using dotted lines from the tips 
#' @param colUnsam Give the unsampled cases colors (TRUE), or have all unsampled cases be grey (FALSE)
#' @param trCol The color for the symbol annotating transmission events on the tree 
#' @param trShape Integer code of the plotting symbol for annotating transmission events on the tree 
#' @param cex Numeric character expansion factor for the tip labels 
#' @examples
#' plotCTree(simulateOutbreak())
#' @export
plotCTree = function(tree, maxTime=NA, showLabels=TRUE, alignLabels=FALSE, colUnsam = FALSE, trCol="black", trShape = 20, cex = 1){
  nam=tree$nam
  tree=tree$ctree
  nsam <- sum(tree[ ,2]+tree[ ,3] == 0) 
  nh <- nrow(tree)-3*nsam+1
  ntot <- nsam+nh
  par(yaxt='n',bty='n')
  plot(0,0,type='l',xlim=c(min(tree[,1]),ifelse(is.na(maxTime),max(tree[,1]),maxTime)),ylim=c(0,nsam+1),xlab='',ylab='')
  host <- tree[ ,4] 
  
  #Determine ys for leaves
  root<-which(host==0)
  ys <- matrix(0, nsam, 1) 
  todo <- cbind(root,0,0.5,1);#Matrix of nodes to do,with associated starting x and y coordinates and scale 
  while (nrow(todo) > 0)  { 
    w <- todo[1,1] 
    x <- todo[1,2] 
    y <- todo[1,3] 
    scale <- todo[1,4] 
    if (tree[w,2] == 0 && tree[w,3] == 0)  { 
      #Leaf node 
      ys[w] <- y 
    } else if (tree[w,3] == 0)  { 
      #Transmission node 
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1],y,scale,deparse.level=0)) 
    } else { 
      #Binary node 
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1],y + scale/2,scale/2,deparse.level=0),cbind(tree[w,3],tree[w,1],y-scale/2,scale/2,deparse.level=0))
    } 
    todo <- rbind(todo[-1,])
  } 
  ys<-rank(ys)
  
  #Determine ys for non-leaves
  for (i in ((nsam+1):nrow(tree))) { 
    children <- c()
    todo <- i
    while (length(todo)>0) {
      children=c(children,todo[1])
      todo=c(todo[-1],setdiff(tree[todo[1],2:3],0))
    }
    ys[i] <- mean(ys[children[which(children<=nsam)]])
  } 
  
  # determine the types of the rows 
  wLeaf <- tree[,2] == 0 & tree[,3] == 0
  wTrans <- tree[,3] == 0 & ! wLeaf
  wBin <- !wLeaf & !wTrans 
  
  wLeaf <- which(wLeaf)
  wTrans <- which(wTrans)
  wBin <- which(wBin) 
  
  # transmission infectors and infectees 
  infectors <- tree[-c(wLeaf,root,wBin),4]
  infectees <- tree[tree[-c(wLeaf,root,wBin),2],4]
 
  # interweave evenly spaced hues  
  hues <- vector(length = ntot)
  if(ntot > 1){
    huesSpaced <- seq(15,375,length = ntot+1)[1:ntot]
    huesSpaced <- c(huesSpaced[seq(1,length(huesSpaced),2)],huesSpaced[seq(2,length(huesSpaced),2)])
    
    # space the hues of transmission chains as far as possible 
    uniInf <- unique(infectors)
    spotForUniInf <- seq(1, ntot, length.out = length(uniInf))
    hues[uniInf]  <- huesSpaced[spotForUniInf]
      
    # give the remaining transmissions the other hues 
    for(i in uniInf){
      secondary <- infectors == i
      alsoInfectee <- hues[infectees[secondary]] != 0
      hues[infectees[secondary]][!alsoInfectee] <- huesSpaced[!(huesSpaced %in% hues)][1:sum(!alsoInfectee)]
    }
  } else {
    hues <- seq(15,375,length = ntot+1)
  }

  grDevices::palette(hcl(h = hues, l = 65, c = 100))
  l <- 0
  todo <- cbind(root,tree[root,1]);
  while (nrow(todo) > 0)  { 
    w <- todo[1,1] 
    x <- todo[1,2] 
    y <- ys[w] 
    col=host[w]
    l <- l + 1 
  
    # check if in a sampled host 
    if( is.na(nam[host[w]] ) && host[w] != 0 && !colUnsam) {
      if (tree[w,2] == 0 && tree[w,3] == 0)  { 
        #Leaf node 
        lines(c(x,tree[w,1]),c(y,y),col= "gray",lwd=2) 
        if (showLabels && !alignLabels) text(tree[w,1] + (max(cbind(tree[ ,1]))-min(cbind(tree[ ,1])))/100,y,nam[w],cex = cex)
      } else if (tree[w,3] == 0)  { 
        #Transmission node 
        lines(c(x,tree[w,1]),c(y,y),col= "gray",lwd=2) 
        points(tree[w,1],y,col = trCol,pch=trShape) 
        todo <- rbind(todo,cbind(tree[w,2],tree[w,1])) 
      } else { 
        #Binary node 
        lines(c(x,tree[w,1]),c(y,y),col= "gray",lwd=2)
        lines(c(tree[w,1],tree[w,1]),cbind(ys[tree[w,2]],ys[tree[w,3]]),col= "gray",lwd=2)
        todo <- rbind(todo,cbind(tree[w,2],tree[w,1]),cbind(tree[w,3],tree[w,1])) 
      } 
      todo <- rbind(todo[-1,]);
    } else {
      if (tree[w,2] == 0 && tree[w,3] == 0)  { 
        #Leaf node 
        lines(c(x,tree[w,1]),c(y,y),col=col,lwd=2) 
        if (showLabels && !alignLabels) text(tree[w,1] + (max(cbind(tree[ ,1]))-min(cbind(tree[ ,1])))/100,y,nam[w],cex = cex)
      } else if (tree[w,3] == 0)  { 
        #Transmission node 
        lines(c(x,tree[w,1]),c(y,y),col=col,lwd=2) 
        points(tree[w,1],y,col = trCol,pch=trShape) 
        todo <- rbind(todo,cbind(tree[w,2],tree[w,1])) 
      } else { 
        #Binary node 
        lines(c(x,tree[w,1]),c(y,y),col=col,lwd=2)
        lines(c(tree[w,1],tree[w,1]),cbind(ys[tree[w,2]],ys[tree[w,3]]),col=col,lwd=2)
        todo <- rbind(todo,cbind(tree[w,2],tree[w,1]),cbind(tree[w,3],tree[w,1])) 
      } 
      todo <- rbind(todo[-1,]);
    }
  } 
  
  if(alignLabels){
    leafYs <- ys[wLeaf]
    leafxs <- tree[wLeaf,1]
    maxX <- max(tree[,1])
    
    for(i in 1:length(leafYs)){
      lines(c(leafxs[i],maxX),c(leafYs[i],leafYs[i]),lwd=2, lty = "dotted", col = "gray")
    } 
  }
  if (showLabels && alignLabels) text(maxX+0.03,leafYs,nam[wLeaf],xpd=NA,cex = cex) 
} 
