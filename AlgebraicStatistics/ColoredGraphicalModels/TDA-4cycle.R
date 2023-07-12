# reject null if p<0.05
library("TDAstats")
library("proxy")
library("Matrix")

genStats1<-function(magn,numVar,r,sampleSize,numStats){
  randL<-runif(numVar*r*sampleSize,-magn,magn)
  sample<-array(randL,dim=c(r,numVar,sampleSize))
  statsM<-matrix(,nrow = sampleSize, ncol = numStats)
  for(i in 1:sampleSize){
    stats<-statsComp((matrix(sample[,,i]))%*%t(matrix(sample[,,i])))
    statsM[i,]<-stats
  }
  return(statsM)
}

genStats<-function(magn,numVar,r,sampleSize,numStats){
  randL<-runif(numVar*r*sampleSize,-magn,magn)
  sample<-array(randL,dim=c(r,numVar,sampleSize))
  statsM<-matrix(,nrow = sampleSize, ncol = numStats)
  for(i in 1:sampleSize){
    stats<-statsComp(t(sample[,,i])%*%(sample[,,i]))
    statsM[i,]<-stats
  }
  return(statsM)
}

#magn=1
#numVar=4
#r=4
#sampleSize=1
#numStats=8
#pos=c(2)
#tol=100
genStatsNbourhood<-function(magn,numVar,r,sampleSize,numStats,pos, tol){
  sample<-array(,dim=c(r,numVar,sampleSize))
  for(i in 1:sampleSize){
    for (j in pos){
  sample[j,,i]<-runif(numVar,-magn,magn)
    }
    for (k in setdiff(1:r,pos))
  sample [k,,i]<-runif(numVar,-magn/tol,magn/tol)    
  }
  #randL<-runif(numVar*r*sampleSize,-magn,magn)
  #sample<-array(randL,dim=c(r,numVar,sampleSize))
  statsM<-matrix(,nrow = sampleSize, ncol = numStats)
  for(i in 1:sampleSize){
    stats<-statsComp(t(sample[,,i])%*%(sample[,,i]))
    statsM[i,]<-stats
  }
  return(statsM)
}

# uncolored 4-cycle
statVector<-c(1,5,6,10,11,13,15,16)
statsComp<-function(m){
  return(c(m)[statVector])
}

# magn=1
# numVar=4
# r=4
# sampleSize=1
# numStats=8
# m=(sample[,,1])%*%t(sample[,,1])
# (matrix(sample[,,1]))%*%t(matrix(sample[,,1]))
# stats<-statsComp(m)
# diag(m)
# sample[,,1]


#Generate sufficient stats for matrices of different rank
stats4<-genStats(0.1,4,4,10,length(statVector))
stats3<-genStats(0.1,4,3,10,length(statVector))
stats2<-genStats(0.1,4,2,10,length(statVector))
stats1<-genStats1(0.1,4,1,10,length(statVector))
jointAll<-rbind(stats4,stats3,stats2,stats1) 

median(stats4)
median(stats3)
median(stats2)
median(stats1)
median(jointAll)

dist4=dist(stats4,method="euclidean")
dist3=dist(stats3,method="euclidean")
dist2=dist(stats2,method="euclidean")
dist1=dist(stats1,method="euclidean")
distAll=dist(jointAll,method="euclidean")

median(dist4)
median(dist3)
median(dist2)
median(dist1)
median(distAll)

stats1<-genStats1(0.1,4,1,4000,length(statVector))
stats4N1pos1<-genStatsNbourhood(0.1,4,4,1000,length(statVector),c(1), 100)
stats4N1pos2<-genStatsNbourhood(0.1,4,4,1000,length(statVector),c(2), 100)
stats4N1pos3<-genStatsNbourhood(0.1,4,4,1000,length(statVector),c(3), 100)
stats4N1pos4<-genStatsNbourhood(0.1,4,4,1000,length(statVector),c(4), 100)
stats4N1<-rbind(stats4N1pos1,stats4N1pos2,stats4N1pos3,stats4N1pos4)
jointAll<-rbind(stats1,stats4N1)

c(median(stats1[,1]),median(stats1[,2]),median(stats1[,3]),median(stats1[,4]),
  median(stats1[,5]),median(stats1[,6]),median(stats1[,7]),median(stats1[,8]))
c(median(stats4N1[,1]),median(stats4N1[,2]),median(stats4N1[,3]),median(stats4N1[,4]),
  median(stats4N1[,5]),median(stats4N1[,6]),median(stats4N1[,7]),median(stats4N1[,8]))
c(median(jointAll[,1]),median(jointAll[,2]),median(jointAll[,3]),median(jointAll[,4]),
  median(jointAll[,5]),median(jointAll[,6]),median(jointAll[,7]),median(jointAll[,8]))

dist1=dist(stats1,method="euclidean")
dist4N1=dist(stats4N1,method="euclidean")
distAll=dist(jointAll,method="euclidean")

median(dist1)
median(dist4N1)
median(distAll)

write.table(stats4, "G18-4-01.csv", sep=",",  col.names=FALSE)
write.table(stats3, "G18-3-01.csv", sep=",",  col.names=FALSE)
write.table(stats2, "G18-2-01.csv", sep=",",  col.names=FALSE)
write.table(stats1, "G18-1-01.csv", sep=",",  col.names=FALSE)

#Compare basic descriptive properties of each statistic (pay attention to range!)
#sometimes rank 1 does not go to -1 for some stats
statId=3
summary(stats4[,statId])
summary(stats3[,statId])
summary(stats2[,statId])
summary(stats1[,statId])

wilcox.test(stats4,stats3) #reject the null
wilcox.test(stats4,stats2) #reject the null
wilcox.test(stats4,stats1) #reject the null
# Generate auxiliary rank 4 sufficient stats for mixing
# Generating separately to ensure independence 
stats4for4<-genStats(1,4,4,1000,length(statVector))
stats4for3<-genStats(1,4,4,1000,length(statVector))
stats4for2<-genStats(1,4,4,1000,length(statVector))
stats4for1<-genStats(1,4,4,1000,length(statVector))

# Generate auxiliary rank 3,2,1 sufficient stats for concatenating
# Generating separately to ensure independence 
stats3aux<-genStats(1,4,3,1000,length(statVector))
stats2aux<-genStats(1,4,2,1000,length(statVector))
stats1aux<-genStats(1,4,1,1000,length(statVector))


#Calculate homology of dim=0 for sufficients stats of different rank
stats4.phom<-calculate_homology(stats4, dim=0)
stats3.phom<-calculate_homology(stats3, dim=0)
stats2.phom<-calculate_homology(stats2, dim=0)
stats1.phom<-calculate_homology(stats1, dim=0)

stats4all.phom<-calculate_homology(stats4)
plot_barcode(stats4all.phom)
plot_persist(stats4all.phom)

joint2with1<-rbind(stats2,stats1) 
joint2with1.phom<-calculate_homology(joint2with1)
plot_barcode(joint2with1.phom)
plot_persist(joint2with1.phom)

stats1all.phom<-calculate_homology(stats1)
plot_persist(stats1all.phom)

stats3aux.phom<-calculate_homology(stats3aux, dim=0)
stats2aux.phom<-calculate_homology(stats2aux, dim=0)
stats1aux.phom<-calculate_homology(stats1aux, dim=0)

#null-hypothesis: X and Y have the same probability distribution
wilcox.test(stats4.phom[,3],stats3.phom[,3]) #reject
wilcox.test(stats4.phom[,3],stats2.phom[,3]) #reject
wilcox.test(stats4.phom[,3],stats1.phom[,3]) #reject

#Calculate homology of dim=0 for sufficients stats of union of different types of stats
#4 with 4 vs 4 with 3
joint4with3<-rbind(stats4for3,stats3) #need to use a different set of full rank stats to ensure independence
joint4with3.phom<-calculate_homology(joint4with3, dim=0)
joint4with4<-rbind(stats4for4,stats4) #need to use a different set of full rank stats to ensure independence
joint4with4.phom<-calculate_homology(joint4with4, dim=0)
wilcox.test(joint4with4.phom[,3],joint4with3.phom[,3]) #fail to reject the null

##plot histogram
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hg4with3=hist(joint4with3.phom[,3],plot=FALSE)
hg4=hist(joint4with4.phom[,3],plot=FALSE)
plot(hg4with3, col = c1)
plot(hg4, col = c2, add = TRUE)

#Calculate homology of dim=0 for sufficients stats of union of different types of stats
#4 with 4 vs 4 with 2
joint4with2<-rbind(stats4for2,stats2) #need to use a different set of full rank stats to ensure independence
joint4with2.phom<-calculate_homology(joint4with2, dim=0)
wilcox.test(joint4with4.phom[,3],joint4with2.phom[,3]) #fail to reject the null

##plot histogram
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hg4with2=hist(joint4with2.phom[,3],plot=FALSE)
hg4=hist(joint4with4.phom[,3],plot=FALSE)
plot(hg4with2, col = c1)
plot(hg4, col = c2, add = TRUE)

#Calculate homology of dim=0 for sufficients stats of union of different types of stats
#4 with 4 vs 4 with 1
joint4with1<-rbind(stats4for1,stats1) #need to use a different set of full rank stats to ensure independence
joint4with1.phom<-calculate_homology(joint4with1, dim=0)
wilcox.test(joint4with4.phom[,3],joint4with1.phom[,3]) #reject the null

##plot histogram
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hg4with1=hist(joint4with1.phom[,3],plot=FALSE)
hg4=hist(joint4with4.phom[,3],plot=FALSE)
plot(hg4with1, col = c1)
plot(hg4, col = c2, add = TRUE)



#Compare to the case when we concatenate results for different ranks vs truly mixing
concat4with3<-c(stats4.phom[,3],stats3aux.phom[,3])
wilcox.test(concat4with3,joint4with3.phom[,3]) #reject the null

concat4with2<-c(stats4.phom[,3],stats2aux.phom[,3])
wilcox.test(concat4with2,joint4with2.phom[,3]) #reject the null

concat4with1<-c(stats4.phom[,3],stats1aux.phom[,3])
wilcox.test(concat4with1,joint4with1.phom[,3]) #fail to reject the null

# Graph 17
numStats=7
statsComp<-function(m){
  return(c(m[1,1]+m[3,3],m[2,2],m[4,4],m[1,2],m[2,3],m[3,4],m[1,4]))
}

#Generate sufficient stats for matrices of different rank
stats4<-genStats(1,4,4,5000,numStats)
stats3<-genStats(1,4,3,5000,numStats)
stats2<-genStats(1.1,4,2,5000,numStats)
stats1<-genStats1(1.5,4,1,5000,numStats)

write.table(stats4, "G17-4.csv", sep=",",  col.names=FALSE)
write.table(stats3, "G17-3.csv", sep=",",  col.names=FALSE)
write.table(stats2, "G17-2.csv", sep=",",  col.names=FALSE)
write.table(stats1, "G17-1.csv", sep=",",  col.names=FALSE)

#Compare basic descriptive properties of each statistic (pay attention to range!)
#sometimes rank 1 does not go to -1 for some stats
statId=5
summary(stats4[,statId])
summary(stats3[,statId])
summary(stats2[,statId])
summary(stats1[,statId])

stats4all.phom<-calculate_homology(stats4)
plot_barcode(stats4all.phom)
plot_persist(stats4all.phom)

joint2with1<-rbind(stats2,stats1) 
joint2with1.phom<-calculate_homology(joint2with1)
plot_barcode(joint2with1.phom)
plot_persist(joint2with1.phom)

stats1all.phom<-calculate_homology(stats1)
plot_persist(stats1all.phom)

# Generate auxiliary rank 4 sufficient stats for mixing
# Generating separately to ensure independence 
stats4for4<-genStats(1,4,4,1000,numStats)
stats4for3<-genStats(1,4,4,1000,numStats)
stats4for2<-genStats(1,4,4,1000,numStats)
stats4for1<-genStats(1,4,4,1000,numStats)

# Generate auxiliary rank 3,2,1 sufficient stats for concatenating
# Generating separately to ensure independence 
stats3aux<-genStats(1,4,3,1000,numStats)
stats2aux<-genStats(1,4,2,1000,numStats)
stats1aux<-genStats(1,4,1,1000,numStats)


#Calculate homology of dim=0 for sufficients stats of different rank
stats4.phom<-calculate_homology(stats4, dim=0)
stats3.phom<-calculate_homology(stats3, dim=0)
stats2.phom<-calculate_homology(stats2, dim=0)
stats1.phom<-calculate_homology(stats1, dim=0)

stats3aux.phom<-calculate_homology(stats3aux, dim=0)
stats2aux.phom<-calculate_homology(stats2aux, dim=0)
stats1aux.phom<-calculate_homology(stats1aux, dim=0)

#null-hypothesis: X and Y have the same probability distribution
wilcox.test(stats4.phom[,3],stats3.phom[,3]) #fail to reject
wilcox.test(stats4.phom[,3],stats2.phom[,3]) #fail to reject
wilcox.test(stats4.phom[,3],stats1.phom[,3]) #reject

#Calculate homology of dim=0 for sufficients stats of union of different types of stats
#4 with 4 vs 4 with 3
joint4with3<-rbind(stats4for3,stats3) #need to use a different set of full rank stats to ensure independence
joint4with3.phom<-calculate_homology(joint4with3, dim=0)
joint4with4<-rbind(stats4for4,stats4) #need to use a different set of full rank stats to ensure independence
joint4with4.phom<-calculate_homology(joint4with4, dim=0)
wilcox.test(joint4with4.phom[,3],joint4with3.phom[,3]) #fail to reject the null

##plot histogram
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hg4with3=hist(joint4with3.phom[,3],plot=FALSE)
hg4=hist(joint4with4.phom[,3],plot=FALSE)
plot(hg4with3, col = c1)
plot(hg4, col = c2, add = TRUE)

#Calculate homology of dim=0 for sufficients stats of union of different types of stats
#4 with 4 vs 4 with 2
joint4with2<-rbind(stats4for2,stats2) #need to use a different set of full rank stats to ensure independence
joint4with2.phom<-calculate_homology(joint4with2, dim=0)
wilcox.test(joint4with4.phom[,3],joint4with2.phom[,3]) #fail to reject the null

##plot histogram
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hg4with2=hist(joint4with2.phom[,3],plot=FALSE)
hg4=hist(joint4with4.phom[,3],plot=FALSE)
plot(hg4with2, col = c1)
plot(hg4, col = c2, add = TRUE)

#Calculate homology of dim=0 for sufficients stats of union of different types of stats
#4 with 4 vs 4 with 1
joint4with1<-rbind(stats4for1,stats1) #need to use a different set of full rank stats to ensure independence
joint4with1.phom<-calculate_homology(joint4with1, dim=0)
wilcox.test(joint4with4.phom[,3],joint4with1.phom[,3]) #reject the null

##plot histogram
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hg4with1=hist(joint4with1.phom[,3],plot=FALSE)
hg4=hist(joint4with4.phom[,3],plot=FALSE)
plot(hg4with1, col = c1)
plot(hg4, col = c2, add = TRUE)



#Compare to the case when we concatenate results for different ranks vs truly mixing
concat4with3<-c(stats4.phom[,3],stats3aux.phom[,3])
wilcox.test(concat4with3,joint4with3.phom[,3]) #reject the null

concat4with2<-c(stats4.phom[,3],stats2aux.phom[,3])
wilcox.test(concat4with2,joint4with2.phom[,3]) #reject the null

concat4with1<-c(stats4.phom[,3],stats1aux.phom[,3])
wilcox.test(concat4with1,joint4with1.phom[,3]) #reject the null

# Graph 16
numStats=7
statsComp<-function(m){
  return(c(m[1,1]+m[2,2],m[3,3],m[4,4],m[1,2],m[2,3],m[3,4],m[1,4]))
}

tesPts=runif(4*1*1,-1,1)
testSample<-array(tesPts,dim=c(1,4,1))
testMatrix=matrix(testSample[,,1])%*%t(matrix(testSample[,,1]))
qr(testMatrix)$rank

#Generate sufficient stats for matrices of different rank
stats4<-genStats(1,4,4,5000,numStats)
stats3<-genStats(1.1,4,3,5000,numStats)
stats2<-genStats(1.3,4,2,5000,numStats)
stats1<-genStats1(1.7,4,1,5000,numStats)

# stats4<-genStats(1,4,4,5000,numStats)
# stats3<-genStats(1,4,3,5000,numStats)
# stats2<-genStats(1,4,2,5000,numStats)
# stats1<-genStats1(1,4,1,5000,numStats)

write.table(stats4, "G16-4.csv", sep=",",  col.names=FALSE)
write.table(stats3, "G16-3.csv", sep=",",  col.names=FALSE)
write.table(stats2, "G16-2.csv", sep=",",  col.names=FALSE)
write.table(stats1, "G16-1.csv", sep=",",  col.names=FALSE)

wilcox.test(stats4,stats3) #fail to reject
wilcox.test(stats4,stats2) #fail reject
wilcox.test(stats4,stats1) #fail reject

#Compare basic descriptive properties of each statistic (pay attention to range!)
#sometimes rank 1 does not go to -1 for some stats
statId=7
summary(stats4[,statId])
summary(stats3[,statId])
summary(stats2[,statId])
summary(stats1[,statId])

# Generate auxiliary rank 4 sufficient stats for mixing
# Generating separately to ensure independence 
stats4for4<-genStats(1,4,4,1000,numStats)
stats4for3<-genStats(1,4,4,1000,numStats)
stats4for2<-genStats(1,4,4,1000,numStats)
stats4for1<-genStats(1,4,4,1000,numStats)

# Generate auxiliary rank 3,2,1 sufficient stats for concatenating
# Generating separately to ensure independence 
stats3aux<-genStats(1,4,3,1000,numStats)
stats2aux<-genStats(1,4,2,1000,numStats)
stats1aux<-genStats(1,4,1,1000,numStats)


#Calculate homology of dim=0 for sufficients stats of different rank
stats4.phom<-calculate_homology(stats4, dim=0)
stats3.phom<-calculate_homology(stats3, dim=0)
stats2.phom<-calculate_homology(stats2, dim=0)
stats1.phom<-calculate_homology(stats1, dim=0)

stats3aux.phom<-calculate_homology(stats3aux, dim=0)
stats2aux.phom<-calculate_homology(stats2aux, dim=0)
stats1aux.phom<-calculate_homology(stats1aux, dim=0)

stats1all.phom<-calculate_homology(stats1, dim=7)
plot_persist(stats1all.phom)
plot_barcode(stats1all.phom)

#null-hypothesis: X and Y have the same probability distribution
wilcox.test(stats4.phom[,3],stats3.phom[,3]) #fail to reject
wilcox.test(stats4.phom[,3],stats2.phom[,3]) #fail reject
wilcox.test(stats4.phom[,3],stats1.phom[,3]) #fail reject

#Calculate homology of dim=0 for sufficients stats of union of different types of stats
#4 with 4 vs 4 with 3
joint4with3<-rbind(stats4for3,stats3) #need to use a different set of full rank stats to ensure independence
joint4with3.phom<-calculate_homology(joint4with3, dim=0)
joint4with4<-rbind(stats4for4,stats4) #need to use a different set of full rank stats to ensure independence
joint4with4.phom<-calculate_homology(joint4with4, dim=0)
wilcox.test(joint4with4.phom[,3],joint4with3.phom[,3]) #fail to reject the null

##plot histogram
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hg4with3=hist(joint4with3.phom[,3],plot=FALSE)
hg4=hist(joint4with4.phom[,3],plot=FALSE)
plot(hg4with3, col = c1)
plot(hg4, col = c2, add = TRUE)

#Calculate homology of dim=0 for sufficients stats of union of different types of stats
#4 with 4 vs 4 with 2
joint4with2<-rbind(stats4for2,stats2) #need to use a different set of full rank stats to ensure independence
joint4with2.phom<-calculate_homology(joint4with2, dim=0)
wilcox.test(joint4with4.phom[,3],joint4with2.phom[,3]) #fail to reject the null

##plot histogram
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hg4with2=hist(joint4with2.phom[,3],plot=FALSE)
hg4=hist(joint4with4.phom[,3],plot=FALSE)
plot(hg4with2, col = c1)
plot(hg4, col = c2, add = TRUE)

#Calculate homology of dim=0 for sufficients stats of union of different types of stats
#4 with 4 vs 4 with 1
joint4with1<-rbind(stats4for1,stats1) #need to use a different set of full rank stats to ensure independence
joint4with1.phom<-calculate_homology(joint4with1, dim=0)
wilcox.test(joint4with4.phom[,3],joint4with1.phom[,3]) #fail to reject

##plot histogram
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hg4with1=hist(joint4with1.phom[,3],plot=FALSE)
hg4=hist(joint4with4.phom[,3],plot=FALSE)
plot(hg4with1, col = c1)
plot(hg4, col = c2, add = TRUE)



#Compare to the case when we concatenate results for different ranks vs truly mixing
concat4with3<-c(stats4.phom[,3],stats3aux.phom[,3])
wilcox.test(concat4with3,joint4with3.phom[,3]) #reject the null

concat4with2<-c(stats4.phom[,3],stats2aux.phom[,3])
wilcox.test(concat4with2,joint4with2.phom[,3]) #reject the null

concat4with1<-c(stats4.phom[,3],stats1aux.phom[,3])
wilcox.test(concat4with1,joint4with1.phom[,3]) #reject the null



# Graph 15
numStats=6
statsComp<-function(m){
  return(c(m[1,1]+m[2,2]+m[3,3],m[4,4],m[1,2],m[2,3],m[3,4],m[1,4]))
}


#Generate sufficient stats for matrices of different rank
stats4<-genStats(1,4,4,5000,numStats)
stats3<-genStats(1.1,4,3,5000,numStats)
stats2<-genStats(1.3,4,2,5000,numStats)
stats1<-genStats1(1.7,4,1,5000,numStats)

write.table(stats4, "G15-4.csv", sep=",",  col.names=FALSE)
write.table(stats3, "G15-3.csv", sep=",",  col.names=FALSE)
write.table(stats2, "G15-2.csv", sep=",",  col.names=FALSE)
write.table(stats1, "G15-1.csv", sep=",",  col.names=FALSE)

#Compare basic descriptive properties of each statistic (pay attention to range!)
#sometimes rank 1 does not go to -1 for some stats
statId=2
summary(stats4[,statId])
summary(stats3[,statId])
summary(stats2[,statId])
summary(stats1[,statId])

# Generate auxiliary rank 4 sufficient stats for mixing
# Generating separately to ensure independence 
stats4for4<-genStats(1,4,4,1000,numStats)
stats4for3<-genStats(1,4,4,1000,numStats)
stats4for2<-genStats(1,4,4,1000,numStats)
stats4for1<-genStats(1,4,4,1000,numStats)

# Generate auxiliary rank 3,2,1 sufficient stats for concatenating
# Generating separately to ensure independence 
stats3aux<-genStats(0.95,4,3,1000,numStats)
stats2aux<-genStats(0.95,4,2,1000,numStats)
stats1aux<-genStats(0.95,4,1,1000,numStats)


#Calculate homology of dim=0 for sufficients stats of different rank
stats4.phom<-calculate_homology(stats4, dim=0)
stats3.phom<-calculate_homology(stats3, dim=0)
stats2.phom<-calculate_homology(stats2, dim=0)
stats1.phom<-calculate_homology(stats1, dim=0)

stats3aux.phom<-calculate_homology(stats3aux, dim=0)
stats2aux.phom<-calculate_homology(stats2aux, dim=0)
stats1aux.phom<-calculate_homology(stats1aux, dim=0)

#null-hypothesis: X and Y have the same probability distribution
wilcox.test(stats4.phom[,3],stats3.phom[,3]) #reject
wilcox.test(stats4.phom[,3],stats2.phom[,3]) #reject
wilcox.test(stats4.phom[,3],stats1.phom[,3]) #reject

#Calculate homology of dim=0 for sufficients stats of union of different types of stats
#4 with 4 vs 4 with 3
joint4with3<-rbind(stats4for3,stats3) #need to use a different set of full rank stats to ensure independence
joint4with3.phom<-calculate_homology(joint4with3, dim=0)
joint4with4<-rbind(stats4for4,stats4) #need to use a different set of full rank stats to ensure independence
joint4with4.phom<-calculate_homology(joint4with4, dim=0)
wilcox.test(joint4with4.phom[,3],joint4with3.phom[,3]) #fail to reject

##plot histogram
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hg4with3=hist(joint4with3.phom[,3],plot=FALSE)
hg4=hist(joint4with4.phom[,3],plot=FALSE)
plot(hg4with3, col = c1)
plot(hg4, col = c2, add = TRUE)

#Calculate homology of dim=0 for sufficient stats of union of different types of stats
#4 with 4 vs 4 with 2
joint4with2<-rbind(stats4for2,stats2) #need to use a different set of full rank stats to ensure independence
joint4with2.phom<-calculate_homology(joint4with2, dim=0)
wilcox.test(joint4with4.phom[,3],joint4with2.phom[,3]) #fail to reject

##plot histogram
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hg4with2=hist(joint4with2.phom[,3],plot=FALSE)
hg4=hist(joint4with4.phom[,3],plot=FALSE)
plot(hg4with2, col = c1)
plot(hg4, col = c2, add = TRUE)

#Calculate homology of dim=0 for sufficient stats of union of different types of stats
#4 with 4 vs 4 with 1
joint4with1<-rbind(stats4for1,stats1) #need to use a different set of full rank stats to ensure independence
joint4with1.phom<-calculate_homology(joint4with1, dim=0)
wilcox.test(joint4with4.phom[,3],joint4with1.phom[,3]) #reject

##plot histogram
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hg4with1=hist(joint4with1.phom[,3],plot=FALSE)
hg4=hist(joint4with4.phom[,3],plot=FALSE)
plot(hg4with1, col = c1)
plot(hg4, col = c2, add = TRUE)



#Compare to the case when we concatenate results for different ranks vs truly mixing
concat4with3<-c(stats4.phom[,3],stats3aux.phom[,3])
wilcox.test(concat4with3,joint4with3.phom[,3]) #fail to reject

concat4with2<-c(stats4.phom[,3],stats2aux.phom[,3])
wilcox.test(concat4with2,joint4with2.phom[,3]) #fail to reject

concat4with1<-c(stats4.phom[,3],stats1aux.phom[,3])
wilcox.test(concat4with1,joint4with1.phom[,3]) #fail to reject

#PLOTTING
#Compare stats
plot(stats4[,1:2],col="blue")
points(stats3[,1:2], col = "red")

plot(stats4[,c(1,3)],col="blue")
points(stats3[,c(1,3)], col = "red")

plot(stats4[,2:3],col="blue")
points(stats3[,2:3], col = "red")

plot(stats4[,1:2],col="blue")
points(stats1[,1:2], col = "red")

plot(stats4[,c(1,5)],col="blue")
points(stats1[,c(1,5)], col = "red")

plot(stats4[,2:3],col="blue")
points(stats1[,2:3], col = "red")

a1=rnorm(10)
a2=rnorm(10)
plot(a1)
points(a2,col=2)
plot(1,2,col="blue")
points(3,4, col = "red")
