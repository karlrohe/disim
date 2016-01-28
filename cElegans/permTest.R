source('makecElegans.R')

reps = 1000  # this is the number of times it will perform the resampling.

g = graph.adjacency(adjm)
ord = order(page.rank(g)$vec)
tmp = adjm[ord,ord]
truescore = sum(tmp[lower.tri(tmp)])
truescore

score = rep(NA, reps)  
for(tick in 1:reps){
atmp = diag(rep(0,7))
for(i in 1:7) atmp[i,-i] = sample(adjm[i,-i])
g = graph.adjacency(atmp)
ord = order(page.rank(g)$vec)
tmp = atmp[ord,ord]
score[tick] = sum(tmp[lower.tri(tmp)])
}
mean(score>truescore)