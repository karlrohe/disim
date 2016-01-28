rm(list = ls())
library(igraph)
library(xtable)


G = read.graph("data/polblogs.gml", format  = "gml")

# find largest connected component
clust = clusters(G)  
N = length(V(G))
bad = (1:N)[-which(clust$membership == which.max(clust$csize))]
G = delete.vertices(G, bad)

# these are the dem blogs
dem = (V(G)$value == 0)
dem = 1:N
dem = dem < 759  
dem = dem[-bad]


# compute the adjacency matrix from the igraph
A = get.adjacency(G)
A = as.matrix(A)
n = nrow(A)

####################################
######  Undirected Analysis   ######
####################################

As = A + t(A)  # A_ij is number of links between two blogs.
rss = rowSums(As)  # row sums symmetric
css = colSums(As)  # column sums symmetric
mds  = mean(rss)  # mean degree symmetric

Ls = diag(1/sqrt(rss + mds))%*%As%*%diag(1/sqrt(css + mds))
ss = svd(Ls)

us = ss$u[,1:2]  #  u symmetric
plot(us, pch = ".")
abline(0,0)
boxplot(us[,2]~dem)

usn = t(apply(us,1,function(x) return(x/sqrt(sum(x^2 + 10^(-10))))))
boxplot(usn[,2]~dem)

plot(usn)

# leverage scores:
lev = rowSums(us^2)
lowLev = which(log(lev)< -15)
points(usn[lowLev,], pch = 19, col = "red")
length(lowLev)


kms = kmeans(usn,2, nstart = 100)
points(kms$centers, pch = 19, col ="blue")

dem.est = rep(F, n)
Tab = table(kms$clust); demLabel = as.numeric(names(Tab[which.min(Tab)]))
dem.est[kms$clust == demLabel] = T  # XXX if there are strange label switch issues, check this line.
table(dem.est, dem)
correct = dem.est == dem
if(sum(correct) < n/2){
  flip = dem.est==T; 
  dem.est[flip] = F
  dem.est[!flip] = T
}
sum(dem.est == dem)
wrong = dem.est != dem

plot(density((log(lev))))
lines(density((log(lev[wrong]))), col = "red")  # misclustered nodes tend to have small leverage scores.





##################################
######  Directed Analysis   ######
##################################


rs = rowSums(A)
cs = colSums(A)
md  = mean(rs)

L = diag(1/sqrt(rs + md))%*%A%*%diag(1/sqrt(cs + md))

s= svd(L)
plot(s$d)
plot(s$u[,1:2])
plot(s$v[,1:2])
plot(s$u[,2], s$v[,2])

u = s$u[,1:2]; v = s$v[,1:2]
un = t(apply(u,1,function(x) return(x/sqrt(sum(x^2) + 10^(-10)))))

# "good" corresponds to nodes with degree greater than one.  
goodu = which(rowSums(A)!=0)  # here it is out-degree
vn = t(apply(v,1,function(x) return(x/sqrt(sum(x^2) + 10^(-10)))))
goodv = which(colSums(A)!=0)  # here it is in-degree
boxplot(vn[goodv,2]~dem[goodv])
boxplot(un[goodu,2]~dem[goodu])

good = intersect(goodu, goodv)  # "good" nodes must both send and receive at least one edge.
plot(un[good,2], vn[good,2])
all = rbind(un[goodu,],vn[goodv,])
km = kmeans(all,centers = 2, nstart = 1000)

ulabels = rep(NA,n)
ulabels[goodu] = km$clust[1:length(goodu)]

vlabels = rep(NA,n)
vlabels[goodv] = km$clust[-(1:length(goodu))]
table(ulabels, vlabels, dem)

table(ulabels)

# ensure that label 1 is dem (i.e. "political left") and label 2 is gop (i.e. "political right")
if(sum(diag(table(ulabels,dem))) > n/2){
  gops = ulabels==1; dems = ulabels ==2
  ulabels[gops] = 2; ulabels[dems] = 1
  gops = vlabels==1; dems = vlabels ==2
  vlabels[gops] = 2; vlabels[dems] = 1
}
table(ulabels)
table(vlabels)

dirlab= rep(NA, n)
demTF = function(lab){
  if(is.na(lab)) return("NA")
  if(lab == 1) return("dem")
  if(lab == 2) return("gop")
  
}
dirlabfun = function(u,v){
  return(paste(demTF(v),2,demTF(u), sep=""))
}
for(i in 1:n)  {
  dirlab[i] = dirlabfun(ulabels[i], vlabels[i])
  
}
table(ulabels)


table(dem.est, dirlab)

these = (1:n)[-grep("NA", dirlab)]  # removes the not-good elements.
table(dem.est[these], dirlab[these])

which(dirlab == "dem2dem" & dem.est == F)

table(dirlab[-grep("NA", dirlab)])

# now restrict to nodes with at least 3 in and 3 out edges.

hideg = (rs >2) & (cs>2)
these = intersect(these, which(hideg))
table(dem.est[hideg], dirlab[hideg])

names = get.vertex.attribute(G, "label")
names[dirlab == "dem2gop" & hideg]
names[dirlab == "gop2dem"  & hideg]
d2g = dirlab == "dem2gop" & hideg
g2d = dirlab == "gop2dem"  & hideg

x = which(d2g + g2d == 1)
to.dem = rowSums(A[x,dem])  
to.gop =rowSums(A[x,!dem])
from.dem =colSums(A[dem,x])
from.gop=colSums(A[!dem,x])

mixedNodes= data.frame(cbind(names[x], dirlab[x], from.dem, from.gop, to.dem, to.gop))



xtable(mixedNodes)



