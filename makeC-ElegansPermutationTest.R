
k = 7  # adjust as needed. 
source('makecElegans.R')

reps = 1000  # this is the number of times it will perform the resampling.

g = graph.adjacency(adjm)
ord = order(page.rank(g)$vec)
tmp = adjm[ord,ord]
truescore = sum(tmp[lower.tri(tmp)])
truescore

score = rep(NA, reps)  
for(tick in 1:reps){
  atmp = diag(rep(0,k))
  for(i in 1:k) atmp[i,-i] = sample(adjm[i,-i])
  g = graph.adjacency(atmp)
  ord = order(page.rank(g)$vec)
  tmp = atmp[ord,ord]
  score[tick] = sum(tmp[lower.tri(tmp)])
}
mean(score<=truescore)




# for the permutation test on B:

rownames(adjmFinal) = 1:k
colnames(adjmFinal) = 1:k
# library(xtable)
# xtable(adjmFinal)



Zu = model.matrix(~ as.factor(kmallu) -1)
Zv = model.matrix(~ as.factor(kmallv) -1)


# # do the analysis on the unweighted graph.
Ab = A
Ab[A>0] = 1 
B = solve(t(Zu)%*%Zu)%*%t(Zu)%*%Ab%*%Zv %*%solve(t(Zv)%*%Zv)
B = B[-1,]; 
# B = B[,-1]  # remove the NA groups

# # UNCOMMENT HERE TO do the analysis on the weighted graph.
# B = solve(t(Zu)%*%Zu)%*%t(Zu)%*%A%*%Zv %*%solve(t(Zv)%*%Zv)
# B = B[-1,] # remove the NA groups
# # B = B[,-1]  



g = graph.adjacency(B, weighted = T)
ord = order(page.rank(g)$vec)
tmp = B[ord,ord]
truescore = sum(tmp[lower.tri(tmp)])
truescore

score = rep(NA, reps)  
for(tick in 1:reps){
  atmp = diag(rep(0,k))
  for(i in 1:k) atmp[i,-i] = sample(B[i,-i])
  g = graph.adjacency(atmp, weighted = T)
  ord = order(page.rank(g)$vec)
  tmp = atmp[ord,ord]
  score[tick] = sum(tmp[lower.tri(tmp)])
}
mean(score<=truescore)