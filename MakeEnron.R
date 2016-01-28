rm(list = ls())
M = read.table(file = "data/enron/messages.tsv", sep = "\t")
MM = as.matrix(M[,c(1,3,5)])
#* message id (1-21635)
#* filename
#* unix time (seconds since January 1, 1970)
#* message subject
#* sender's employee id (1-156)

R = read.table(file = "data/enron/recipients.tsv", sep = "\t")
#* message id (1-21635)
#* recipient number (1-57)
#* receiver's employee id (1-156)

E= read.table(file = "data/enron/employees.txt", sep = "\t")
dim(E)

# need to make the graph.
# time, from, to
n = length(R[,1])
G = matrix(0, nrow = n, ncol = 3)
MMindex = match(R[,1], MM[,1])
for(i in 1:n){
  if(complete.cases(MMindex[i]))  G[i,] = c(MM[MMindex[i],2:3], R[i,3])
}
G = G[complete.cases(MMindex),]

library(Matrix)

Gm = spMatrix(nrow = 156, ncol = 156, G[,2], G[,3], rep(1, nrow(G)))
Gm = as.matrix(Gm)
# these people are not connected to anybody in the graph
rem = which((colSums(Gm) ==0 )*(rowSums(Gm) ==0 ) ==1)
# remove them
Gm = Gm[-rem, -rem]
E = E[-rem,]




md = mean(rowSums(Gm))
Dr = diag(1/sqrt(rowSums(Gm) + md))
Dc = diag(1/sqrt(colSums(Gm) + md))
L = Dr%*%Gm%*%Dc

s = svd(L)


pdf(file = "figures/Enron.pdf", height = 3, width = 5)
par(mfrow = c(1,2))
par(las = 1)
plot(s$d[1:15], ylab = "singular values", main = "Eigengap after the\nsecond and fifth values.")

points(s$d[1:5], pch = 19, col = "grey")
points(s$d[1:2], pch = 19)

u = s$u[,1:2]; v = s$v[,1:2]
m2= sqrt(rowSums((u-v)^2))
u = s$u[,1:5]; v = s$v[,1:5]
m5 = sqrt(rowSums((u-v)^2))
n = length(m2)
xx = c(rep(2,n), rep(5,n))
xj = rnorm(xx,xx,.03)
plot(xj, c(m2,m5), xlim = c(1,6), xaxt="n", bty="n", yaxt = "n", ylab = "", xlab = "", main = "Asymmetry scores")
m2ind = which.max(m2); m5ind = which.max(m5)
points(xj[c(m2ind, m5ind+n)], c(m2,m5)[c(m2ind, m5ind+n)], col ="red", pch = 19)
axis(side = 1,at = c(2,5),labels = c("K=2","K=5"),tick = F)
axis(side = 2,at = seq(0,.75,by = .25), labels = seq(0,.75,by = .25),tick = T , las = 1)

dev.off()
