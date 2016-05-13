# scree plot
rm(list = ls())
k=5  # this does not affect the code here.  just needed to load the data...
source('makecElegans.R')



pdf("figures/wormScree.pdf", height = 4, width =3)
n = 20
par(las = 1, bty="n", mar = c(0, 3,4,1))
plot(c(0,s$d[1:n]), col = c("white", rep("red", 5), rep("black", 2), rep("grey", n-7)), pch = 19, xlab ="", main= "top 20 singular values", ylab = "",xaxt="n", yaxp = c(.1,.6,5))
dev.off()
