rm(list = ls())

plotTitle =""; removeSendingNodes = character(0); removeReceivingNodes = character(0); morespacefor7 = 0; palatte =  c("yellow", "pink", "orange", "lavender", "dark green", "blue", "gold"); leg = T
plotxclust = function(plotTitle ="", removeSendingNodes = character(0), removeReceivingNodes = character(0), morespacefor7 = 0, palatte =  c("yellow", "pink", "orange", "lavender", "dark green", "blue", "gold"), leg = T){
  source('~/Dropbox/clusteringT+K/Convert to PNAS/cele/loadCel.R')
  load('~/Dropbox/clusteringT+K/Convert to PNAS/cele/originalOrder.RData')
  Y = kmallu
  remrow = which(rowSums(A)==0)
  remcol = which(colSums(A)==0)
  outneuron[remrow] = NA
  inneuron[remcol] = NA
  
  if(length(removeSendingNodes) > 0){
    outneuron = outneuron[-removeSendingNodes]
    A = A[-removeSendingNodes, ]
    Y = Y[-removeSendingNodes]
  }
  if(length(removeReceivingNodes) > 0){
    inneuron = inneuron[-removeReceivingNodes]
    A = A[, -removeReceivingNodes]
  }
  
  link = rep(NA,length(outneuron))
  for(i in 1:length(outneuron)){
    if(complete.cases(outneuron[i]) ) if(as.character(outneuron[i])!= ""){
      tmp = which(as.character(outneuron[i]) == inneuron)
      if(length(tmp) > 0) link[i] = tmp
    }
  }
  
  
  linkback = rep(NA,length(inneuron))
  for(i in 1:length(inneuron)){
    if(complete.cases(inneuron[i]) )   if(as.character(inneuron[i])!= ""){
      tmp = which(inneuron[i] == as.character(outneuron))
      if(length(tmp) > 0 )  linkback[i] = tmp
    }
  }
  
  A = log(A+1)
  
  rs = rowSums(A); mrs = mean(rs)
  cs = colSums(A); mcs = mean(cs)
  Dr = diag((rs + mrs)^(-1/2))
  Dc = diag((cs + mcs)^(-1/2))
  L = Dr%*%A%*%Dc
  s = svd(L)
  k = 7
  
  u = t(apply(s$u[,1:k]%*%diag(s$d[1:k]), 1, 
              function(x) return(x/sqrt(sum(x^2) + 10^(-10)))))
  v = t(apply(s$v[,1:k]%*%diag(s$d[1:k]), 1, 
              function(x) return(x/sqrt(sum(x^2) + 10^(-10)))))
  dim(u)
  dim(v)
  
  # these two lines are to ensure zero degree nodes are not included in the kmeans step.
  
  kmall = kmeans(rbind(u[-remrow,],v[-remcol,]), centers = 7, nstart = 100)
  clustlab = kmall$cluster
  
  # after removing the zero degree nodes, we must ensure that the rest of the code is unaffected.  So, put in a bunch of NA's for the 0-deg nodes. 
  ucnew = rep(NA, nrow(u)); vcnew = rep(NA, nrow(v))
  ucnew[-remrow] = clustlab[1:nrow(u[-remrow,])]
  vcnew[-remcol] = clustlab[-(1:nrow(u[-remrow,]))]
  clustlab = c(ucnew,vcnew)
  
  kmallu = clustlab[1:nrow(u)]
  newLabels = lm(Y[-remrow] ~ as.factor(kmallu[-remrow])-1)$coef
  newLabels = rank(newLabels)
  allzeros = c(remrow, nrow(u) + remcol)
  clustlab[-allzeros] = newLabels[clustlab[-allzeros]]
#   
#   tmp = clustlab
#   tmp[clustlab==1] = 2
#   tmp[clustlab== 2] = 1
#   tmp[clustlab==3] = 4
#   tmp[clustlab==4] = 6
#   tmp[clustlab==5] = 7
#   tmp[clustlab==6] = 3
#   tmp[clustlab==7] = 5
#   clustlab = tmp
#   
  kmallu = clustlab[1:nrow(u)]
  kmallv = clustlab[-(1:nrow(u))]
  
  adjm = table(kmallu[1:144], kmallv[1:144])
  # adjm = adjm[-1,]; adjm = adjm[,-1]  # remove NA's
  g = graph.adjacency(adjm, mode = "directed")
  ord = order(page.rank(g)$vec)
  
  tmp = clustlab
  tmp[clustlab==ord[1]] = 1
  tmp[clustlab== ord[2]] = 2
  tmp[clustlab==ord[3]] = 3
  tmp[clustlab==ord[4]] = 4
  tmp[clustlab==ord[5]] = 5
  tmp[clustlab==ord[6]] = 6
  tmp[clustlab==ord[7]] = 7
  clustlab = tmp
  
  kmallu = clustlab[1:nrow(u)]
  kmallv = clustlab[-(1:nrow(u))]
  
  
  
  
  
  
  # there are two types of sending (or receiving) nodes that do not receive (send) edges. 
  # those that are not included in inneuron (outneuron)
  # those that have colSums(A) (rowSums(A)) equal to zero.
  
  
  #these are columns that do not have a corresponding row in A
  nooutedges = match(1:length(inneuron), link)
  nooutedges = which(is.na(nooutedges))  
  
  # these are rows that do not have a corresponding columns
  noinedges = which(is.na(inneuron))
  noinedges = noinedges[which(noinedges < nrow(A))]
  # noinedges = which(is.na(link))
  # noinedges= which(colSums(A[,match(outneuron, inneuron)])==0)   # we can use match here because all outneurons have a corresponding inneuron.
  # nooutedges = which(is.na(match(inneuron, outneuron)))
  # nooutedges = c(nooutedges, which(rowSums(A[match(inneuron, outneuron),])==0))
  
  
  
  
  plot(c(-.5,5), c(.5,7.5 + morespacefor7), col = "white", ylab = "", xlab =  "", main = plotTitle, xaxt="n", yaxt = "n", bty="n")
  # points(rep(.8, length(noinedges)), kmallu[noinedges] + runif(length(noinedges), -.2,.2))
  points(rep(2.2, length(nooutedges)), kmallv[nooutedges]+ runif(length(nooutedges), -.2,.2))
  newu = kmallu + runif(length(kmallu), -.2,.2)
  newv = kmallv + runif(length(kmallv), -.2,.2)
  points(rep(1, length(kmallu) ), newu)
  points(rep(2, length(kmallv) ), newv)
  #how far are the points in the latent space?
  # nodesep = rep(0, sum(!is.na(link)))
  
  
  nodesep = rep(0, length(link))
  tablematch = matrix(0, nrow = nrow(A), ncol = 4)
  for(i in 1:nrow(A)){
    matchi = which(as.character(outneuron[i]) == inneuron)
    if(length(matchi) > 0)   nodesep[i] = sqrt(sum((u[i,] - v[matchi,])^2))
    tmp = .1 + 1.1*(nodesep[i] > .75)
    if(nodesep[i] > .75) color = "black"
    if(nodesep[i]<= .75) color = "grey"
    if(length(matchi) > 0)   lines(c(1,2), c(newu[i], newv[matchi]), lwd = tmp, col = color)
    #tablematch[i,] = c(i, matchi, kmallu[])
  }
  
  
  
  library(graphics)
  
  haveinedges = (1:length(kmallu))
  if(length(noinedges)>0) haveinedges = (1:length(kmallu))[-noinedges]  # these are the sending nodes that also have in edges.
  haveoutedges = (1:length(kmallv))[-nooutedges]  # these are the receiving nodes that also have out edges. 
  
  allGroups = c("A", "B", "R", "L", "D", "E", "F")
  for(clust in 1:7){
    sset = which(kmallu== clust); sset = intersect(sset, haveinedges)  # this is the set of sending nodes, in cluster "clust," that also have in edges.
    rset = which(kmallv== clust); rset = intersect(rset, haveoutedges)  #visa versa
    
    # write down the sending clusters (for the nodes that also receive edges)
    xtick = -.5
    ytick = clust + .3 + .1*(clust==7)  + .1*(clust==6)
    for(nodeid in sset){
      lab = outneuron[nodeid]
      size = (nodesep[nodeid] < .75)/1.6 + (nodesep[nodeid] > .75)/1.2
      grouppp = match(outClusterLabels[nodeid], allGroups)
      colll = palatte[grouppp]
      text(xtick,ytick, paste(lab,",", sep=""), cex = size, col = colll)
      xtick = xtick + .2 #+ .05*(nodesep[nodeid] > .75)
      if(xtick>.6){
        ytick = ytick -.2
        xtick = -.5
      }
    }
    # # write down the sending nodes that don't receive any edges.
    # sset = which(kmallu== clust); sset = intersect(sset, noinedges)  # this is the set of sending nodes, in cluster "clust," that have NO in edges.
    # if(length(sset)>0){
    # xtick = -1.5; ytick = clust + .3 + (clust==2)*.3 
    # size = 1/1.6
    # for(nodeide in sset){
    # lab = outneuron[nodeid]
    # grouppp = match(outClusterLabels[nodeid], allGroups)
    # colll = palatte[grouppp]
    # text(xtick,ytick, paste(lab,",", sep=""), cex = size, col = colll)
    # xtick = xtick + .2
    # if(xtick>-.6){
    # ytick = ytick -.2
    # xtick = -1.5
    # }
    # }
    # }
    
    
    
    # write the names of the receiving nodes that send edges
    xtick = 2.5
    ytick = clust + .3 + .7*(clust ==7) + (clust==6)*.3
    for(nodeid in rset){
      lab = inneuron[nodeid]
      travels = nodesep[linkback[nodeid]]
      size = (travels < .75)/1.6 + (travels > .75)/1.2
      grouppp = match(inClusterLabels[nodeid], allGroups)
      colll = palatte[grouppp]
      text(xtick,ytick, paste(as.character(lab),",", sep=""), cex = size, col = colll)
      xtick = xtick + .2 #+ .05*(travels > .75)
      if(xtick>3.4){
        ytick = ytick -.2
        xtick = 2.5
      }
    }
    # write the names of the receiving nodes that DO NOT send edges  (e.g. muscles)
    rset = which(kmallv== clust); rset = intersect(rset, nooutedges)  
    if(length(rset)>0){
      xtick = 3.9; ytick = clust + .3 +  (clust==7)*(-.3)  #+ .3*(clust==3)
      size = 1/1.6
      for(nodeid in rset){
        lab = inneuron[nodeid]
        grouppp = match(inClusterLabels[nodeid], allGroups)
        colll = palatte[grouppp]
        text(xtick,ytick, paste(lab,",", sep=""), cex = size, col = colll)
        xtick = xtick + .3 
        if(xtick>4.9){
          ytick = ytick -.2
          xtick = 3.9
        }
      }
    }
    
    
    
    
  }
  if(leg == T){
    legend("topright", c("C-Response", "C-Locomotion", "D-R(1-5)A", "E-PVV", "F-Insemination"), text.col = palatte[-(1:2)])
    
  }
}




setwd('~/Dropbox/clusteringT+K/Convert to PNAS/')
pdf(file = "xclust2.pdf", height = 8, width =12)
plotxclust(morespacefor7 = 1)
dev.off()
