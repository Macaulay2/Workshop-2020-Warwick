library("TDA")
library("TDAstats")
library("Rfast")
library("proxy")

statsComp<-function(m){
  return(c(diag(m),c(m[1:5,6:10])))
}

genStats<-function(magn,numVar,r,sampleSize,numStats){
  randL<-runif(numVar*r*sampleSize,-magn,magn)
  sample<-array(randL,dim=c(r,numVar,sampleSize))
  statsM<-matrix(,nrow = sampleSize, ncol = numStats)
  for(i in 1:sampleSize){
    stats<-statsComp(t(sample[,,i])%*%sample[,,i])
    statsM[i,]<-stats
  }
  return(statsM)
}

r<-1
magn<- 1
numVar<-10
sampleSize<-1000
numStats<-35

stats10<-genStats(1,numVar,10,1000,numStats)
stats10test<-genStats(1,numVar,10,1000,numStats)
stats5<-genStats(1.2,numVar,5,1000,numStats)
stats4<-genStats(1.25,numVar,4,1000,numStats)
stats3<-genStats(1.3,numVar,3,1000,numStats)
stats2<-genStats(1.35,numVar,2,1000,numStats)
#stats1<-genStats(magn,numVar,1,sampleSize,numStats) -- need to figure out because 
#R automatically changes vectors between row and column

stats10.phom<-calculate_homology(stats10, dim=0)
joint10with5<-rbind(stats10,stats5)
joint10with10<-rbind(stats10,stats10test)
joint10with5.phom<-calculate_homology(joint10with5, dim=0)
joint10with10.phom<-calculate_homology(joint10with10, dim=0)
stats5.phom<-calculate_homology(stats5, dim=0)
concat10with5<-c(stats10.phom[,3],stats5.phom[,3])

hist(stats10.phom[,3])
hist(joint10with5.phom[,3])
hist(stats5.phom[,3])

hist(joint10with5.phom[,3])
hist(concat10with5)
hist(joint10with10.phom[,3])


stats4.phom<-calculate_homology(stats4, dim=0)
joint10with4<-rbind(stats10,stats4)
joint10with4.phom<-calculate_homology(joint10with4, dim=0)
concat10with4<-c(stats10.phom[,3],stats4.phom[,3])

hist(stats10.phom[,3])
hist(joint10with4.phom[,3])
hist(stats4.phom[,3])

hist(joint10with4.phom[,3])
hist(concat10with4)

stats3.phom<-calculate_homology(stats3, dim=0)
joint10with3<-rbind(stats10,stats3)
joint10with3.phom<-calculate_homology(joint10with3, dim=0)
concat10with3<-c(stats10.phom[,3],stats3.phom[,3])

hist(stats10.phom[,3])
hist(joint10with3.phom[,3])
hist(stats3.phom[,3])

hist(stats10.phom[,3])
hist(joint10with3.phom[,3])
hist(concat10with3)
hist(joint10with10.phom[,3])

stats2.phom<-calculate_homology(stats2, dim=0)
joint10with2<-rbind(stats10,stats2)
joint10with2.phom<-calculate_homology(joint10with2, dim=0)
concat10with2<-c(stats10.phom[,3],stats2.phom[,3])

hist(stats10.phom[,3])
hist(joint10with2.phom[,3])
hist(stats2.phom[,3])

hist(joint10with2.phom[,3])
hist(concat10with2)

maxDeath10with5=max(joint10with5.phom[,3])

maxDeath10=max(stats10.phom[,3])

#dist10to10=dist(stats10,stats10)
#diag(dist10to10)<-rep(10000,10)

dist10to5=dist(stats10,stats5)
#head(dist10to5,1)
minDist10to10=apply(dist10to10, 1, FUN = min)
minDist10to5=apply(dist10to5, 1, FUN = min)
which(minDist10to10>minDist10to5)
#any(minDist54>maxDeath5)
#head(minDist10to5,-1)-stats10.phom[,3]
which((head(minDist10to5,-1)-stats10.phom[,3])>0)
remaining10to5=which(minDist10to5>maxDeath10) 
length(remaining10to5)
#length 95 when 1000 vs 1000
#       60 when 1000 vs 2000
#       58 when 1000 vs 3000
#       46 when 1000 vs 5000
#       38 when 1000 vs 10000

dist10to4=dist(stats10,stats4)
minDist10to4=apply(dist10to4, 1, FUN = min)
#any(minDist54>maxDeath5)
remaining10to4=which(minDist10to4>maxDeath10) 
length(remaining10to4)

#length 275 when 1000 vs 1000
#length 228 when 1000 vs 2000
#       200 when 1000 vs 3000
#       185 when 1000 vs 5000
#       158 when 1000 vs 10000

stats10to4remaining.phom <- calculate_homology(stats10[remaining10to4,])
plot_persist(stats10to4remaining.phom)

##test with random matrices
M1000A<-matrix(runif(10000,-1,1),ncol=10)
M1000B<-matrix(runif(10000,-1,1),ncol=10)
M2000<-rbind(M1000A,M1000B)

M1000A.phom<-calculate_homology(M1000A, dim=0)
M2000.phom<-calculate_homology(M2000, dim=0)
M1000B.phom<-calculate_homology(M1000B, dim=0)

hist(M1000A.phom[,3])
hist(M2000.phom[,3])
hist(M1000B.phom[,3])

#PSDL5 <- t(read.csv("~/Dropbox/Math conferences/Workshop-2020-Warwick/AlgebraicStatistics/ColoredGraphicalModels/PSDL5.csv", header=FALSE))
#PSDL4 <- t(read.csv("~/Dropbox/Math conferences/Workshop-2020-Warwick/AlgebraicStatistics/ColoredGraphicalModels/PSDL4.csv", header=FALSE))
#PSDL3 <- t(read.csv("~/Dropbox/Math conferences/Workshop-2020-Warwick/AlgebraicStatistics/ColoredGraphicalModels/PSDL3.csv", header=FALSE))



#maxscale <- 5        # limit of the filtration
#maxdimension <- 1    # components and loops

#DiagRips <- ripsDiag(X = PSDL5, maxdimension, maxscale,library = c("GUDHI", "Dionysus"), location = TRUE, printProgress = FALSE)

PSDL5.phom <- calculate_homology(PSDL5, dim=0)
plot_barcode(PSDL5.phom)
plot_persist(PSDL5.phom)

PSDL4.phom <- calculate_homology(PSDL4)
plot_barcode(PSDL4.phom)
plot_persist(PSDL4.phom)

PSDL3.phom <- calculate_homology(PSDL3)
plot_barcode(PSDL3.phom)
plot_persist(PSDL3.phom)

#for (i in 1:100) {
#  PSDL5[i,]
#}
maxDeath5=max(PSDL5.phom[,3])


dist54=dist(PSDL5,PSDL4)
minDist54=apply(dist54, 1, FUN = min)
#any(minDist54>maxDeath5)
remaining54=which(minDist54>maxDeath5)
PSDL54remaining.phom <- calculate_homology(PSDL5[remaining54,])

dist53=dist(PSDL5,PSDL3)
minDist53=apply(dist53, 1, FUN = min)
remaining53=which(minDist53>maxDeath5)
PSDL53remaining.phom <- calculate_homology(PSDL5[remaining53,])
plot_barcode(PSDL53remaining.phom)
plot_persist(PSDL53remaining.phom)
