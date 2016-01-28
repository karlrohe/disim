# setwd("/Users/karlrohe/Dropbox/DirectedClustering/jrssb2/code/celegans")
# setwd('')
x = read.csv("~/Dropbox/clusteringT+K/Convert to PNAS/cele/newCE.csv", na.strings = "XX")

A = x[3:nrow(x), 4:ncol(x)]
outneuron = x[-c(1:2),3]
inneuron = x[2,-c(1:3)]
# visual check to see that the neurons line up:
 # rbind(as.character(outneuron), (inneuron[1:152]))[,1:115]
n = 115
A = as.matrix(A)
bad = which(colSums(is.na(A)) > 150)
bad = c(bad, 116, 154, 177, 222, 223)  # these other 5 columns of A represent sums of edges from different subsets of nodes. 
inneuron[c(116, 154, 177, 222, 223)]
A = A[,-bad]
inneuron = inneuron[-bad]
for(i in 1:nrow(A)) for(j in 1:ncol(A)) if(A[i,j] =="") A[i,j] = 0
Anew = matrix(0, nrow = nrow(A), ncol = ncol(A))
for(i in 1:nrow(A)) for(j in 1:ncol(A)) Anew[i,j] = as.numeric(A[i,j])
A = Anew
A[which(is.na(A), arr.ind = T)] = 0

# these neurons don't have names and they don't have edges:
which(is.na(match(outneuron, unlist(inneuron))) )
# they are just rows of zeros in the adjacency matrix.  remove them.
bad = which(is.na(match(outneuron, unlist(inneuron))) )
A = A[-bad,]
outneuron = outneuron[-bad]
outneuron = unlist(outneuron)
inneuron = unlist(inneuron)

tmp = read.csv(file = "~/Dropbox/clusteringT+K/Convert to PNAS/cele/newCEwithLabels.csv", na.strings = "XX")
tmp = tmp[-c((1:2), bad+2),]
# this shows that the labels still match.
which(as.character(outneuron) != as.character(tmp[,3]))
cbind(as.character(outneuron) , as.character(tmp[,3]))
outClusterLabels = as.character(tmp[,4])
# outClusterLabels is a label for the out neurons. need to translate this to a label for in neurons 
which(as.character(outneuron) != as.character(inneuron[1:length(outneuron)]))  # names line up.
#cbind(as.character(outneuron), as.character(inneuron[1:length(outneuron)]))
inClusterLabels = rep("", length(inneuron))
inClusterLabels[1:length(outneuron)] = outClusterLabels

#character matching to color the motor groups.
# groupA starts with aol,  dBWM,  pol, cdl, vBWM,  dgl, ail, pil
# groupB starts with dsr, dsp, gec, grt, aob, pob, vsp, vsr, adp
# groupL starts with vBWM, ail/pil
groupL = c("vBW", "ail", "pil")
# groupD dgl, cdl/pol
groupD = c("dgl", "cdl", "pol")
# groupE dBWM
groupE = c("dBW")
# groupF sph, grt/gec/aob/pob, vsp/dsp/adp, int
groupF = c("sph", "grt", "gec", "aob", "pob", "vsp", "dsp", "adp", "int")

char3 = substr(as.character(inneuron), 1, 3)
inClusterLabels[which(char3 %in% groupL)] = "L"
inClusterLabels[which(char3 %in% groupD)] = "D"
inClusterLabels[which(char3 %in% groupE)] = "E"
inClusterLabels[which(char3 %in% groupF)] = "F"

# labelings = list(send = cbind(as.character(outneuron), outClusterLabels), rec = cbind(as.character(inneuron), inClusterLabels))
# labelings
# save(labelings, file = "JarrellNodeLabels.RData")

# load('originalOrder.RData')  # this file helps to order the rows of the figure.
