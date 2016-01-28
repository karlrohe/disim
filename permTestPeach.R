rm(list = ls())
# install.packages("XLConnect")
rm(list =ls())
library(XLConnect)

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

library(igraph)


plotTitle =""; removeSendingNodes = character(0); removeReceivingNodes = character(0); morespacefor7 = 0; palatte =  c("yellow", "pink", "orange", "lavender", "dark green", "blue", "gold"); leg = T


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

kmall = kmeans(rbind(u[-remrow,],v), centers = 7, nstart = 100)
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
kmallv[is.na(kmallv)] = 0

Zu = model.matrix(~ as.factor(kmallu) -1)
Zv = model.matrix(~ as.factor(kmallv) -1)

# find a permutation of the group labels using pageRank
#   to see the feedforward structure.
adjm = table(kmallu, kmallv)
adjm = adjm[-1,]  # remove NA's
g = graph.adjacency(adjm, mode = "directed")
ord = order(page.rank(g)$vec)

adjm = adjm[ord,ord]



g = graph.adjacency(adjm)
ord = order(page.rank(g)$vec)
tmp = adjm[ord,ord]
truescore = sum(tmp[lower.tri(tmp)])
truescore

score = rep(NA, 1000000)
for(tick in 1:1000000){
  atmp = diag(rep(0,7))
  for(i in 1:7) atmp[i,-i] = sample(adjm[i,-i])
  g = graph.adjacency(atmp)
  ord = order(page.rank(g)$vec)
  tmp = atmp[ord,ord]
  score[tick] = sum(tmp[lower.tri(tmp)])
}
print(mean(score>truescore))
save.image(file = "permTest.RData")