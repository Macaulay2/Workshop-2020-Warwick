
--GRAPH 9

restart
load "functions.m2"
R=QQ[l_1..l_6]
K=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
netList IG1_*

m=4 --rank of empirical covariance matrices
k=10 --number of points we want to try out
L=empiricalVanishingPolynomials(IG1,m,k,p,n,stats)
(L1,L2)=differentSign(IG1_1,m,k,p,n,stats);
netList L1 --initial point (suff statistics, evaluation)
netList L2 --if exists, point with different sign (suff statistics, evaluation, number of points we had to try)


v1=vector(flatten entries gens L1_0)
v2=vector(flatten entries gens L2_0)
P=interiorPoint(v1,v2,n,L_0,K)

--Sanity checks for a single point:
apply(flatten P,i->sub(i,RR))
sub(v1,RR)
sub(v2,RR)
sub(sub(sub(L_0,QQ[t_1..t_n]),matrix P),RR)


intP=P_0
tol=0.000001 
subs=rankCompletion(intP,p,rk,stats,S,Rtotal,tol)
apply(flatten entries gens minors(2,subs),i->sub(i,RR))
eigenvalues sub(subs,QQ)

        



--CAREFUL: the following data gave 2 points, one of which does not seem to be a convex
--combination. Can it be in the interior and not be a comb conv?
v1=vector {51241/7056, 14048465/1016064, 15/16, 115729/4704, 691/288, 129/28}
v2=vector {28813/1296, 223745/31752, 165701/6400, 24047/1512, 6681/280, 15719/2520}
apply(flatten P_0,i->sub(i,RR))
sub(v1,RR)
sub(v2,RR)
sub(sub(sub(L_0,QQ[t_1..t_n]),matrix{P_0}),RR)


apply(flatten P_1,i->sub(i,RR))
sub(v1,RR)
sub(v2,RR)
sub(sub(sub(L_0,QQ[t_1..t_n]),matrix{P_1}),RR)

intP=P_1
rk=1
rankCompletion(intP,p,rk,stats,S,Rtotal)
    

--FRED'S HEADS

restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_3,0,l_4},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_4,0,l_5,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
netList IG1_*

m=4 --rank of empirical covariance matrices
k=10 --number of points we want to try out
L=empiricalVanishingPolynomials(IG1,m,k,p,n,stats)
use ring IG1
(L1,L2)=differentSign(IG1_0,m,k,p,n,stats);
netList L1 --initial point (suff statistics, evaluation)
netList L2 --if exists, point with different sign (suff statistics, evaluation, number of points we had to try)


v1=vector(flatten entries gens L1_0)
v2=vector(flatten entries gens L2_0)
P=interiorPoint(v1,v2,n,L_0,K)


intP=P_0
tol=0.000001 
compl=rankCompletion(intP,p,rk,stats,S,Rtotal,tol)

checkRankCompletion(compl,rk,tol)




--GRAPH 13

restart
load "functions.m2"
R=QQ[l_1..l_6]
K=matrix{{l_1,l_3,0,l_6},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_6,0,l_5,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
netList IG1_*

-----test
Rt=QQ[support IG1]
f1=sub(IG1_0,Rt)
f2=sub(IG1_1,Rt)

f=f1^2+f2^2
m=4
k=10
(L1,L2)=differentSign(f,m,k,p,n,stats);
netList L1 --initial point (suff statistics, evaluation)
netList L2 --if exists, point with different sign (suff statistics, evaluation, number of points we had to try)
------

m=4 --rank of empirical covariance matrices
k=10 --number of points we want to try out
L=empiricalVanishingPolynomials(IG1,m,k,p,n,stats)
netList L
use ring IG1
(L1,L2)=differentSign(IG1_0,m,k,p,n,stats);
netList L1 --initial point (suff statistics, evaluation)
netList L2 --if exists, point with different sign (suff statistics, evaluation, number of points we had to try)


v1=vector(flatten entries gens L1_0)
v2=vector(flatten entries gens L2_0)

IG=sub(IG1,QQ[t_1..t_6])
sub(IG_0,transpose matrix v1)
sub(IG_0,transpose matrix v2)


sub(IG_1,transpose matrix v1)
sub(IG_1,transpose matrix v2)


P=interiorPoint(v1,v2,n,L_0,K)


intP=P_0
tol=0.000001 
compl=rankCompletion(intP,p,rk,stats,S,Rtotal,tol)

I=sub(IG1,QQ[support IG1])
dim I, degree I
