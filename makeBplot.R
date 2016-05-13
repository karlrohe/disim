rm(list = ls())

k = 7  # adjust as needed.  
source('makecElegans.R')
library(fields)  # for labeling figures
library(xtable)

rownames(adjmFinal) = 1:k
colnames(adjmFinal) = 1:k
xtable(adjmFinal)



Zu = model.matrix(~ as.factor(kmallu) -1)
Zv = model.matrix(~ as.factor(kmallv) -1)


# # do the analysis on the unweighted graph.
Ab = A
Ab[A>0] = 1 
B = solve(t(Zu)%*%Zu)%*%t(Zu)%*%Ab%*%Zv %*%solve(t(Zv)%*%Zv)
B = B[-1,]; 
# B = B[,-1]  # remove the NA groups

# # UNCOMMENT HERE TO do the analysis on the weighted graph, as in the paper.
# B = solve(t(Zu)%*%Zu)%*%t(Zu)%*%A%*%Zv %*%solve(t(Zv)%*%Zv)
# B = B[-1,] # remove the NA groups
# # B = B[,-1]  




# pdf(file = "figures/celBmat5.pdf", height =4 ,width = 4)
# pdf(file = "figures/celBmat.pdf", height =4 ,width = 4)
par(mar=c(5,4.5,4,7), las = 1)
image((B[k:1,]), col = grey(seq(1,0,len=100)^.5), axes=F, main = quote(widehat("B")),
      xlab = "Receiving clusters", ylab = "Sending clusters")  
axis(1, at=seq(0,1,len = k), labels=as.character(1:k), lty=0)
axis(2, at=seq(0,1,len = k), labels=as.character(k:1), lty=0)
image.plot(t(B[k:1,]), legend.only=T,col = grey(seq(1,0,len=100)^.5))
# dev.off()


# # UNCOMMENT TO inspect the asymmetry in B:
# B = B-t(B)
# par(mar=c(5,4.5,4,7), las = 1)
# image((B[k:1,]), col = grey(seq(1,0,len=100)^.5), axes=F, main = quote(widehat("B")),
#       xlab = "Receiving clusters", ylab = "Sending clusters")  
# axis(1, at=seq(0,1,len = k), labels=as.character(1:k), lty=0)
# axis(2, at=seq(0,1,len = k), labels=as.character(k:1), lty=0)
# image.plot(t(B[k:1,]), legend.only=T,col = grey(seq(1,0,len=100)^.5))
