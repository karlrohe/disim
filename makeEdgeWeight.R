# make figure of edge weights


source('makecElegans.R')
pdf(file = "figures/celEdgeWeight.pdf", height = 3, width = 6.5)
par(mfrow = c(1,3))
hist(G[G!=0], main ="Edge weights", xlab = "")
hist(sqrt(G[G!=0]), main = "SquareRoot(weight)", xlab = "")
hist(log(G[G!=0]+1), main = "Log(weight+1)", xlab = "")
dev.off()