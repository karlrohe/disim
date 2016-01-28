# install.packages("XLConnect")
# rm(list =ls())
library(XLConnect)
library(igraph)


# download data from:
# http://worms.aecom.yu.edu/PHP/chemicalmatrix.xlsx

wb<- loadWorkbook("data/chemicalmatrix.xlsx")
X <- readWorksheet(wb, sheet = "Sheet1")
a = X[,-(1:3)]
a = a[-(1:5),]
A = a[-1,]
A = A[,-1]
A = A[-171,]
A = A[,-237]
A =as.matrix(A)


colnames(A) = a[1,-c(1,237)]
# rownames(A) = a[-c(1,171),1]
A[is.na(A)] = 0
A = as.matrix(A)
class(A) <- "numeric"


A = rbind(A, matrix(0, nrow = ncol(A)-nrow(A), ncol = ncol(A)))
rownames(A) = colnames(A)


# some of the neurons neither send nor receive edges.  Remove them:
rem = which(rowSums(A)==0 & colSums(A) ==0)
A = A[-rem,]
A = A[,-rem]


outneuron = rownames(A)
inneuron = colnames(A)






# the data file lables 7 groups of neurons. 
# R reads them in from excel strange.  This unpacks it.
neuronType = as.matrix(X)[4,-(1:4)]
neuronType = neuronType[1:ncol(A)]

# there are multiple groups with the same character string, but they denote different
#  clusters because of the way they are ordered in the document.
#  clustNum keeps a unique labeling.
clustNum = rep(NA, ncol(A))  
clustNum[1] = 1
names(neuronType) = 1:length(neuronType)
lastType = neuronType[1]
for(i in 2:length(neuronType)){
  if(!is.na(neuronType[i])){
    lastType = neuronType[i]
    clustNum[i] = clustNum[i-1]+1
  }
  if(is.na(neuronType[i])){
    neuronType[i]= lastType
    clustNum[i] = clustNum[i-1]
  }
  
}




lab2 = read.csv("data/JarrellLabelsWithNotes.csv")
outClusterLabels = lab2[match(outneuron, lab2[,2]),3]
inClusterLabels = lab2[match(inneuron, lab2[,2]),3]


# char3 = substr(as.character(colnames(A)), 1, 3)
# 
# 
# 
# #character matching to color the motor groups.
# # groupA starts with aol,  dBWM,  pol, cdl, vBWM,  dgl, ail, pil
# clusterNum[clusterNum == 7] = "A"
# 
# # groupB starts with dsr, dsp, gec, grt, aob, pob, vsp, vsr, adp
# clustNum[which(char3 %in% c("dsr", "dsp", "gec", "grt", "aob", "pob", "vsp", "vsr", "adp"))]
# 
# 
# # groupL starts with vBWM, ail/pil
# groupL = c("vBW", "ail", "pil")
# clustNum[which(char3 %in% groupL)]
# 
# # groupD dgl, cdl/pol
# groupD = c("dgl", "cdl", "pol")
# clustNum[which(char3 %in% groupD)]
# 
# # groupE dBWM
# groupE = c("dBW")
# clustNum[which(char3 %in% groupE)]
# 
# # groupF sph, grt/gec/aob/pob, vsp/dsp/adp, int
# groupF = c("sph", "grt", "gec", "aob", "pob", "vsp", "dsp", "adp", "int")
# 







# ############  Run Di-Sim on c Elegans:





A = log(A+1)
n = nrow(A)

rs = rowSums(A); mrs = mean(rs[rs>0])
cs = colSums(A); mcs = mean(cs[cs>0])
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

# these points will not be used in the k-means step of di-sim.  This is equivalent to remove these rows/columns of the matrix
remrow = which(rowSums(A)==0)  
# remcol = which(colSums(A)==0)  all nodes have an incoming edge.

kmall = kmeans(rbind(u[-remrow,],v), centers = k, nstart = 100)
clustlab = kmall$cluster

# after removing the zero degree nodes, we must ensure that the rest of the code is unaffected.  So, put in a bunch of NA's for the 0-deg nodes. 
ucnew = rep(NA, n); 
ucnew[-remrow] = clustlab[1:nrow(u[-remrow,])]

vcnew = rep(NA, n)
vcnew = clustlab[-(1:nrow(u[-remrow,]))]

clustlab = c(ucnew,vcnew)


kmallu = clustlab[1:n]
kmallv = clustlab[-(1:n)]

kmallu[is.na(kmallu)] = 0 # add groups for NA.
# kmallv[is.na(kmallv)] = 0


# find a permutation of the group labels using pageRank
#   to see the feedforward structure.
adjm = table(kmallu, kmallv)
adjm = adjm[-1,]  # remove NA's
g = graph.adjacency(adjm, mode = "directed")
ord = order(page.rank(g)$vec)


tmp = clustlab
for(j in 1:k) tmp[clustlab==ord[j]] = j
clustlab = tmp

kmallu = clustlab[1:n];  kmallu[is.na(kmallu)] = 0 # add groups for NA.
kmallv = clustlab[-(1:n)]


# # ensure that the relabeling worked:
# adjm = table(kmallu, kmallv)
# adjm = adjm[-1,]  # remove NA's
# g = graph.adjacency(adjm, mode = "directed")
# ord = order(page.rank(g)$vec)
# ord  # this should read 1:7.


adjmFinal = adjm[ord,ord]











