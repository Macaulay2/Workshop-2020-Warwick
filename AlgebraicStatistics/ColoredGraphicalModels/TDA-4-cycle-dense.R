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
  #return(sample)
}

genStatsNbourhood<-function(magn,numVar,r,sampleSize,numStats,pos, tol){
  sample<-array(,dim=c(r,numVar,sampleSize))
  for(i in 1:sampleSize){
    for (j in pos){
      sample[j,,i]<-runif(numVar,-magn,magn)
    }
    for (k in setdiff(1:r,pos))
      sample [k,,i]<-runif(numVar,-magn/tol,magn/tol)    
  }
  statsM<-matrix(,nrow = sampleSize, ncol = numStats)
  for(i in 1:sampleSize){
    stats<-statsComp(t(sample[,,i])%*%(sample[,,i]))
    statsM[i,]<-stats
  }
  return(statsM)
  #return(sample)
}

# uncolored 4-cycle
statVector<-c(1,5,6,10,11,13,15,16)
statsComp<-function(m){
  return(c(m)[statVector])
}

nVar<-4
r<-1
numR4Choice<-1000
numChoice<-choose(nVar,r)
numR4<-numChoice*numR4Choice
numR1<-numR4

stats1<-genStats1(0.1,nVar,r,numR1,length(statVector))
stats4N1pos1<-genStatsNbourhood(0.1,nVar,nVar,numR4Choice,length(statVector),c(1), 100)
stats4N1pos2<-genStatsNbourhood(0.1,nVar,nVar,numR4Choice,length(statVector),c(2), 100)
stats4N1pos3<-genStatsNbourhood(0.1,nVar,nVar,numR4Choice,length(statVector),c(3), 100)
stats4N1pos4<-genStatsNbourhood(0.1,nVar,nVar,numR4Choice,length(statVector),c(4), 100)
stats4N1<-rbind(stats4N1pos1,stats4N1pos2,stats4N1pos3,stats4N1pos4)
jointAll41<-rbind(stats1,stats4N1)

write.table(jointAll41, "G18-4-1.csv", sep=",",  col.names=FALSE)

distAll41=as.matrix(dist(jointAll41,method="euclidean"))+diag(10^6,numR4+numR1,numR4+numR1)
nearn41<-apply(distAll41, 1, FUN = min)
# median(nearn41)
# median(nearn41[1:numR1])
# median(nearn41[(numR1+1):(numR1+numR4)])

sum(nearn41[1:numR1]<=median(nearn41))/length(nearn41[1:numR1])
sum(nearn41[(numR1+1):(numR1+numR4)]<=median(nearn41))/length(nearn41[(numR1+1):(numR1+numR4)])

#hist(nearn[1:numR1])
#hist(nearn[(numR1+1):(numR1+numR4)])

nVar<-4
r<-2
numR4Choice<-1000
numChoice<-choose(nVar,r)
numR4<-numChoice*numR4Choice
numR2<-numR4

stats2<-genStats(0.1,nVar,r,numR2,length(statVector))
#rankMatrix(stats2[,,1])
#rankMatrix(t(stats2[,,1])%*%(stats2[,,1]))
stats4N2pos12<-genStatsNbourhood(0.1,nVar,nVar,numR4Choice,length(statVector),c(1,2), 100)
#rankMatrix(stats4N2pos12[,,1])
#rankMatrix(t(stats4N2pos12[,,1])%*%(stats4N2pos12[,,1]))
stats4N2pos13<-genStatsNbourhood(0.1,nVar,nVar,numR4Choice,length(statVector),c(1,3), 100)
stats4N2pos14<-genStatsNbourhood(0.1,nVar,nVar,numR4Choice,length(statVector),c(1,4), 100)
stats4N2pos23<-genStatsNbourhood(0.1,nVar,nVar,numR4Choice,length(statVector),c(2,3), 100)
stats4N2pos24<-genStatsNbourhood(0.1,nVar,nVar,numR4Choice,length(statVector),c(2,4), 100)
stats4N2pos34<-genStatsNbourhood(0.1,nVar,nVar,numR4Choice,length(statVector),c(3,4), 100)
stats4N2<-rbind(stats4N2pos12,stats4N2pos13,stats4N2pos14,stats4N2pos23,stats4N2pos24,stats4N2pos34)
jointAll42<-rbind(stats2,stats4N2)

write.table(jointAll42, "G18-4-2.csv", sep=",",  col.names=FALSE)

distAll42=as.matrix(dist(jointAll42,method="euclidean"))+diag(10^6,numR4+numR2,numR4+numR2)
nearn42<-apply(distAll42, 1, FUN = min)
#median(nearn)
#median(nearn[1:numR1])
#median(nearn[(numR1+1):(numR1+numR4)])

sum(nearn42[1:numR2]<=median(nearn42))/length(nearn42[1:numR2])
sum(nearn42[(numR2+1):(numR2+numR4)]<=median(nearn42))/length(nearn42[(numR2+1):(numR2+numR4)])

nVar<-4
r<-3
numR4Choice<-1000
numChoice<-choose(nVar,r)
numR4<-numChoice*numR4Choice
numR3<-numR4

stats3<-genStats(0.1,nVar,r,numR3,length(statVector))
#rankMatrix(stats2[,,1])
#rankMatrix(t(stats2[,,1])%*%(stats2[,,1]))
stats4N3pos123<-genStatsNbourhood(0.1,nVar,nVar,numR4Choice,length(statVector),c(1,2,3), 100)
#rankMatrix(stats4N2pos12[,,1])
#rankMatrix(t(stats4N2pos12[,,1])%*%(stats4N2pos12[,,1]))
stats4N3pos124<-genStatsNbourhood(0.1,nVar,nVar,numR4Choice,length(statVector),c(1,2,4), 100)
stats4N3pos134<-genStatsNbourhood(0.1,nVar,nVar,numR4Choice,length(statVector),c(1,3,4), 100)
stats4N3pos234<-genStatsNbourhood(0.1,nVar,nVar,numR4Choice,length(statVector),c(2,3,4), 100)
stats4N3<-rbind(stats4N3pos123,stats4N3pos124,stats4N3pos134,stats4N3pos234)
jointAll43<-rbind(stats3,stats4N3)

write.table(jointAll43, "G18-4-3.csv", sep=",",  col.names=FALSE)

distAll43=as.matrix(dist(jointAll43,method="euclidean"))+diag(10^6,numR4+numR3,numR4+numR3)
nearn43<-apply(distAll43, 1, FUN = min)
#median(nearn)
#median(nearn[1:numR1])
#median(nearn[(numR1+1):(numR1+numR4)])

sum(nearn43[1:numR3]<=median(nearn43))/length(nearn43[1:numR3])
sum(nearn43[(numR3+1):(numR3+numR4)]<=median(nearn43))/length(nearn43[(numR3+1):(numR3+numR4)])

