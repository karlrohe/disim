rm(list = ls())
library(fields)
library(igraph)
# setwd("/Users/karlrohe/Dropbox/DirectedClustering/jrssb2/code/celegans")
setwd('~/Dropbox/clusteringT+K/Convert to PNAS/cele/')

plotTitle =""; removeSendingNodes = character(0); removeReceivingNodes = character(0); morespacefor7 = 0; palatte =  c("yellow", "pink", "orange", "lavender", "dark green", "blue", "gold"); leg = T
source('loadCel.R')

Y = kmallu  # this file helps to order the plot.
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

mat= match(as.character(outneuron[-133]), as.character(inneuron))
m = sqrt(rowSums((u[-133,]-v[mat])^2))

# these two lines are to ensure zero degree nodes are not included in the kmeans step.

kmall = kmeans(rbind(u[-remrow,],v[-remcol,]), centers = 7, nstart = 100)
clustlab = kmall$cluster

# after removing the zero degree nodes, we must ensure that the rest of the code is unaffected.  So, put in a bunch of NA's for the 0-deg nodes. 
ucnew = rep(NA, nrow(u)); vcnew = rep(NA, nrow(v))
ucnew[-remrow] = clustlab[1:nrow(u[-remrow,])]
vcnew[-remcol] = clustlab[-(1:nrow(u[-remrow,]))]
clustlab = c(ucnew,vcnew)

# kmallu = clustlab[1:nrow(u)]
# newLabels = lm(Y[-remrow] ~ as.factor(kmallu[-remrow])-1)$coef
# newLabels = rank(newLabels)
# allzeros = c(remrow, nrow(u) + remcol)
# clustlab[-allzeros] = newLabels[clustlab[-allzeros]]

# tmp = clustlab
# tmp[clustlab==1] = 2
# tmp[clustlab== 2] = 1
# tmp[clustlab==3] = 5
# tmp[clustlab==4] = 7
# tmp[clustlab==5] = 3
# tmp[clustlab==6] = 4
# tmp[clustlab==7] = 6
# clustlab = tmp

kmallu = clustlab[1:nrow(u)]
kmallv = clustlab[-(1:nrow(u))]

kmallu[is.na(kmallu)] = 0 # add groups for NA.
kmallv[is.na(kmallv)] = 0



Zu = model.matrix(~ as.factor(kmallu) -1)
Zv = model.matrix(~ as.factor(kmallv) -1)

# find a different permutation.
adjm = table(kmallu[1:144], kmallv[1:144])
adjm = adjm[-1,]; adjm = adjm[,-1]  # remove NA's
g = graph.adjacency(adjm, mode = "directed")
ord = order(page.rank(g)$vec)

adjm[ord,ord]

Zu = Zu[,-1]; Zv = Zv[,-1]
Zu = Zu[,ord]; Zv = Zv[,ord]

# do the analysis on the unweighted graph.
Ab = A
Ab[A>0] = 1 


# B = solve(t(Zu)%*%Zu)%*%t(Zu)%*%Ab%*%Zv %*%solve(t(Zv)%*%Zv)
# # B = B[-1,]; B = B[,-1]  # remove the NA groups

# do the analysis on the weighted graph.
B = solve(t(Zu)%*%Zu)%*%t(Zu)%*%A%*%Zv %*%solve(t(Zv)%*%Zv)
# B = B[-1,]; B = B[,-1]  # remove the NA groups


pdf(file = "~/Dropbox/clusteringT+K/Convert to PNAS/celBmat.pdf", height =4 ,width = 4)
Btmp = B
# diag(Btmp) = 0
Btmp = t(Btmp[7:1,])
par(mar=c(5,4.5,4,7), las = 1)
image((Btmp), col = grey(seq(1,0,len=100)^.5), axes=F, main = quote(widehat("B")),
      xlab = "Receiving clusters", ylab = "Sending clusters")  
axis(1, at=seq(0,1,len = 7), labels=as.character(1:7), lty=0)
axis(2, at=seq(0,1,len = 7), labels=as.character(7:1), lty=0)
image.plot(t(Btmp), legend.only=T,col = grey(seq(1,0,len=100)^.5))
dev.off()



tab = table(kmallu[1:144], kmallv[1:144])
tab = tab[-1,]; tab = tab[,-1]
tab[ord,ord]
library(xtable)
xtable(tab[ord,ord])




adjm = table(kmallu[1:144], kmallv[1:144])
adjm = adjm[-1,]; adjm = adjm[,-1]
# adjm = as.matrix(B)
g = graph.adjacency(adjm, mode = "directed", weighted = T)
order(page.rank(g)$vec)

ord = order(page.rank(g)$vec)
tmp = adjm[ord,ord]
truescore = sum(tmp[lower.tri(tmp)])
truescore

score = rep(NA, 100000)
for(tick in 1:100000){
  atmp = diag(rep(0,7))
  for(i in 1:7) atmp[i,-i] = sample(adjm[i,-i])
  g = graph.adjacency(atmp, mode = "directed", weighted =T)
  ord = order(page.rank(g)$vec)
  tmp = atmp[ord,ord]
  score[tick] = sum(tmp[lower.tri(tmp)])
}
mean(score<=truescore)












# plot the node transitions from sending to receiving clusters. 

# par(mar=c(5,4.5,4,7), las = 1)
# image((B[7:1,]), col = grey(seq(1,0,len=100)^.5), axes=F, main = quote(widehat("B")),
#       xlab = "Receiving clusters", ylab = "Sending clusters")  
# axis(1, at=seq(0,1,len = 7), labels=as.character(1:7), lty=0)
# axis(2, at=seq(0,1,len = 7), labels=as.character(7:1), lty=0)
# image.plot(t(B[7:1,]), legend.only=T,col = grey(seq(1,0,len=100)^.5))




# pdf(file = "/Users/karlrohe/Dropbox/directedclustering/jrssb2/figures/celBmat.pdf", height =4 ,width = 8)
# par(mar=c(5,4.5,4,7), las = 1, mfrow = c(1,2))
# image(t(B[7:1,]), col = grey(seq(1,0,len=100)^.5), axes=F, main = quote(widehat("B")),
#       xlab = "Receiving clusters", ylab = "Sending clusters")  
# axis(1, at=seq(0,1,len = 7), labels=as.character(1:7), lty=0)
# axis(2, at=seq(0,1,len = 7), labels=as.character(7:1), lty=0)
# image.plot(t(B[7:1,]), legend.only=T,col = grey(seq(1,0,len=100)^.5))
# 
# image((t(B[7:1,] - t(B)[7:1,])), col = grey(seq(1,0,len=100)^.5), axes=F, main = quote(widehat("B") - widehat("B'")),
#       xlab = "Receiving clusters", ylab = "Sending clusters")  
# axis(1, at=seq(0,1,len = 7), labels=as.character(1:7), lty=0)
# axis(2, at=seq(0,1,len = 7), labels=as.character(7:1), lty=0)
# image.plot((t(B[7:1,] - t(B)[7:1,])), legend.only=T,col = grey(seq(1,0,len=100)^.5))
# dev.off()


tmp = B
colnames(tmp) = c()
rownames(tmp) = c()
round(tmp*100)


# # pdf(file = "/Users/karlrohe/Dropbox/directedclustering/jrssb2/figures/celBmat.pdf", height =4 ,width = 4)
# par(mar=c(5,4.5,4,7), las = 1)
# image((B[7:1,]), col = grey(seq(1,0,len=100)^.5), axes=F, main = quote(widehat("B")),
#       xlab = "Receiving clusters", ylab = "Sending clusters")  
# axis(1, at=seq(0,1,len = 7), labels=as.character(1:7), lty=0)
# axis(2, at=seq(0,1,len = 7), labels=as.character(7:1), lty=0)
# image.plot(t(B[7:1,]), legend.only=T,col = grey(seq(1,0,len=100)^.5))
# # dev.off()
