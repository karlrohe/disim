rm(list = ls())
more = c(-.10,0,.1,.15,1.1) +.05  # to adjust the text.

  
  # ytick = clust + .3 + .15*(clust==2)   + .4*(clust==5)

# plotting parameters:
morespacefor7 = .5; 
plotTitle =""; removeSendingNodes = character(0); removeReceivingNodes = character(0); 
palatte =  c("yellow", "pink", "orange", "lavender", "dark green", "blue", "gold"); leg = F
library(graphics)

pdf(file = "figures/xclust5.pdf", height = 8, width =12)
k=5
source('makecElegans.R')
kmallu[kmallu==0]=NA

# in case the labeling is off between in and out neurons...
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

# there are two types of sending (or receiving) nodes that do not receive (send) edges. 
# those that are not included in inneuron (outneuron)
# those that have colSums(A) (rowSums(A)) equal to zero.


#these are columns that do not have a corresponding row in A
# nooutedges = match(1:length(inneuron), link)
nooutedges = which(rowSums(A)==0)  

# these are rows that do not have a corresponding columns
noinedges = which(is.na(inneuron))
noinedges = noinedges[which(noinedges < nrow(A))]
# noinedges = which(is.na(link))
# noinedges= which(colSums(A[,match(outneuron, inneuron)])==0)   # we can use match here because all outneurons have a corresponding inneuron.
# nooutedges = which(is.na(match(inneuron, outneuron)))
# nooutedges = c(nooutedges, which(rowSums(A[match(inneuron, outneuron),])==0))




plot(c(-.5,5), c(.5,5.5 + morespacefor7), col = "white", ylab = "", xlab =  "", main = plotTitle, xaxt="n", yaxt = "n", bty="n")
# points(rep(.8, length(noinedges)), kmallu[noinedges] + runif(length(noinedges), -.2,.2))
points(rep(2.2, length(nooutedges)), kmallv[nooutedges]*1.1+ runif(length(nooutedges), -.2,.2))
newu = kmallu*1.1 + runif(length(kmallu), -.2,.2)
newv = kmallv*1.1 + runif(length(kmallv), -.2,.2)
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





haveinedges = (1:length(kmallu))
if(length(noinedges)>0) haveinedges = (1:length(kmallu))[-noinedges]  # these are the sending nodes that also have in edges.
haveoutedges = (1:length(kmallv))[-nooutedges]  # these are the receiving nodes that also have out edges. 

allGroups = c("A", "B", "R", "L", "D", "E", "F")
for(clust in 1:7){
  sset = which(kmallu== clust); sset = intersect(sset, haveinedges)  # this is the set of sending nodes, in cluster "clust," that also have in edges.
  rset = which(kmallv== clust); rset = intersect(rset, haveoutedges)  #visa versa
  
  # write down the sending clusters (for the nodes that also receive edges)
  xtick = -.5
  ytick = clust*1.1 + .3 + more[clust]
  # ytick = clust + .3 + -.15*(clust==2)  + -.05*(clust==4) -.2*(clust == 6)
  # ytick = clust + .3 + .1*(clust==7)  + .1*(clust==6)
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
  ytick = clust*1.1 + .3 + more[clust]
  # -.1*(clust==4) +.1*(clust == 5) +.1*(clust == 6) +.3*(clust == 7)
  # ytick = clust + .3 + .7*(clust ==7) + (clust==6)*.3
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
    xtick = 3.9; 
    ytick = clust*1.1 + .3 + more[clust]
    # ytick = clust + .3 + .15*(clust==2)   + .4*(clust==5)
    # -.1*(clust==4) +.1*(clust == 5) +.1*(clust == 6) +.3*(clust == 7)
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
  legend("topright", c("C-Response", "C-Locomotion", "D-R(1-5)A", "E-PVV", "F-Insemination"), text.col = palatte[-(1:2)], cex = .7)
  
}




dev.off()
