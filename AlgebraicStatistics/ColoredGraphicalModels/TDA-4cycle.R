# reject null if p<0.05
library("TDAstats")

# uncolored 4-cycle
statVector<-c(1,5,6,10,11,13,15,16)
statsComp<-function(m){
  return(c(m)[statVector])
}

genStats<-function(magn,numVar,r,sampleSize,numStats){
  randL<-runif(numVar*r*sampleSize,-magn,magn)
  sample<-array(randL,dim=c(r,numVar,sampleSize))
  statsM<-matrix(,nrow = sampleSize, ncol = numStats)
  for(i in 1:sampleSize){
    stats<-statsComp((matrix(sample[,,i]))%*%t(matrix(sample[,,i])))
    statsM[i,]<-stats
  }
  return(statsM)
}


#Generate sufficient stats for matrices of different rank
stats4<-genStats(1,4,4,1000,length(statVector))
stats3<-genStats(1,4,3,1000,length(statVector))
stats2<-genStats(1,4,2,1000,length(statVector))
stats1<-genStats(1,4,1,1000,length(statVector))

#Compare basic descriptive properties of each statistic (pay attention to range!)
#sometimes rank 1 does not go to -1 for some stats
statId=8
summary(stats4[,statId])
summary(stats3[,statId])
summary(stats2[,statId])
summary(stats1[,statId])

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

stats3aux.phom<-calculate_homology(stats3aux, dim=0)
stats2aux.phom<-calculate_homology(stats2aux, dim=0)
stats1aux.phom<-calculate_homology(stats1aux, dim=0)

#null-hypothesis: X and Y have the same probability distribution
wilcox.test(stats4.phom[,3],stats3.phom[,3]) #fail to reject
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

genStats<-function(magn,numVar,r,sampleSize,numStats){
  randL<-runif(numVar*r*sampleSize,-magn,magn)
  sample<-array(randL,dim=c(r,numVar,sampleSize))
  statsM<-matrix(,nrow = sampleSize, ncol = numStats)
  for(i in 1:sampleSize){
    stats<-statsComp((matrix(sample[,,i]))%*%t(matrix(sample[,,i])))
    statsM[i,]<-stats
  }
  return(statsM)
}


#Generate sufficient stats for matrices of different rank
stats4<-genStats(1,4,4,1000,numStats)
stats3<-genStats(1,4,3,1000,numStats)
stats2<-genStats(1,4,2,1000,numStats)
stats1<-genStats(1,4,1,1000,numStats)

#Compare basic descriptive properties of each statistic (pay attention to range!)
#sometimes rank 1 does not go to -1 for some stats
#statId=8
#summary(stats4[,statId])
#summary(stats3[,statId])
#summary(stats2[,statId])
#summary(stats1[,statId])

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
wilcox.test(concat4with1,joint4with1.phom[,3]) #reject the null

# Graph 16
numStats=7
statsComp<-function(m){
  return(c(m[1,1]+m[2,2],m[3,3],m[4,4],m[1,2],m[2,3],m[3,4],m[1,4]))
}

genStats<-function(magn,numVar,r,sampleSize,numStats){
  randL<-runif(numVar*r*sampleSize,-magn,magn)
  sample<-array(randL,dim=c(r,numVar,sampleSize))
  statsM<-matrix(,nrow = sampleSize, ncol = numStats)
  for(i in 1:sampleSize){
    stats<-statsComp((matrix(sample[,,i]))%*%t(matrix(sample[,,i])))
    statsM[i,]<-stats
  }
  return(statsM)
}

tesPts=runif(4*1*1,-1,1)
testSample<-array(tesPts,dim=c(1,4,1))
testMatrix=matrix(testSample[,,1])%*%t(matrix(testSample[,,1]))
qr(testMatrix)$rank

#Generate sufficient stats for matrices of different rank
stats4<-genStats(1,4,4,1000,numStats)
stats3<-genStats(1,4,3,1000,numStats)
stats2<-genStats(1,4,2,1000,numStats)
stats1<-genStats(1,4,1,1000,numStats)

#Compare basic descriptive properties of each statistic (pay attention to range!)
#sometimes rank 1 does not go to -1 for some stats
#statId=8
#summary(stats4[,statId])
#summary(stats3[,statId])
#summary(stats2[,statId])
#summary(stats1[,statId])

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

genStats<-function(magn,numVar,r,sampleSize,numStats){
  randL<-runif(numVar*r*sampleSize,-magn,magn)
  sample<-array(randL,dim=c(r,numVar,sampleSize))
  statsM<-matrix(,nrow = sampleSize, ncol = numStats)
  for(i in 1:sampleSize){
    stats<-statsComp((matrix(sample[,,i]))%*%t(matrix(sample[,,i])))
    statsM[i,]<-stats
  }
  return(statsM)
}


#Generate sufficient stats for matrices of different rank
stats4<-genStats(1,4,4,1000,numStats)
stats3<-genStats(0.95,4,3,1000,numStats)
stats2<-genStats(0.95,4,2,1000,numStats)
stats1<-genStats(0.95,4,1,1000,numStats)

#Compare basic descriptive properties of each statistic (pay attention to range!)
#sometimes rank 1 does not go to -1 for some stats
#statId=8
#summary(stats4[,statId])
#summary(stats3[,statId])
#summary(stats2[,statId])
#summary(stats1[,statId])

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
