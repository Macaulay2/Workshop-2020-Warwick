------------------------------------------------------------
-- NON-COMPLETE GRAPHS ON 4 NODES
------------------------------------------------------------


-----------------------------------------------------------
--chordal graphs with 4 nodes
-----------------------------------------------------------
--closer to RCOP
restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_3,l_4,l_3},{l_3,l_2,l_5,0},{l_4,l_5,l_1,l_5},{l_3,0,l_5,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

--Elimination ideal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
netList IG1_*

IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2

--Intersection of the algebraic boundary with the interior
m=4 --rank of empirical covariance matrices
k=100 --number of points we want to try out
use ring IG1
L=empiricalVanishingPolynomials(IG1,m,k,p,n,stats)
(L1,L2)=differentSign(IG1_0,m,k,p,n,stats);
netList L1 --initial point (suff statistics, evaluation)
netList L2 --if exists, point with different sign (suff statistics, evaluation, number of points we had to try)

v1=vector(flatten entries gens L1_0)
v2=vector(flatten entries gens L2_0)
P=interiorPoint(v1,v2,n,L_0,K)


--Rank-1 completion
intP=P_0
tol=0.000001 
subs=rankCompletion(intP,p,rk,stats,S,Rtotal,tol)
--THIS IS TOO LONG!!!
apply(flatten entries gens minors(2,subs),i->sub(i,RR))
eigenvalues sub(subs,QQ)


--MLE for rk 1 matrices
empiricalMLEexistence(1,500,K) --score eq ideal does not have the right dimension
empiricalMLEexistence(2,500,K) --500


--Test ML degree for different ranks
restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_3,l_4,l_3},{l_3,l_2,l_5,0},{l_4,l_5,l_1,l_5},{l_3,0,l_5,l_2}}
p=4
n=1
X=random(QQ^p,QQ^n);              
S=(1/n)*X*transpose(X);
rank S
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
netList I_*
dim I, degree I
--(3,3)
J=saturate(I,det K);
dim J, degree J
--n=4 MLdeg=3
--n=3 MLdeg=3
--n=2 MLdeg=3
--n=1 MLdeg=mostly 2, sometimes 1, in rare occasions 0


--further from RCOP
restart
load "functions.m2"
R=QQ[l_1..l_6]
K=matrix{{l_1,l_3,l_4,l_3},{l_3,l_2,l_5,0},{l_4,l_5,l_1,l_6},{l_3,0,l_6,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

--Elimination ideal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
netList IG1_*

IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2

--Intersection of the algebraic boundary with the interior
m=4 --rank of empirical covariance matrices
k=100 --number of points we want to try out
use ring IG1
L=empiricalVanishingPolynomials(IG1,m,k,p,n,stats)
netList L
(L1,L2)=differentSign(IG1_0,m,k,p,n,stats);
netList L1 --initial point (suff statistics, evaluation)
netList L2 --if exists, point with different sign (suff statistics, evaluation, number of points we had to try)

v1=vector(flatten entries gens L1_0)
v2=vector(flatten entries gens L2_0)
P=interiorPoint(v1,v2,n,L_0,K)


--Rank-1 completion
intP=P_0
tol=0.000001 
subs=rankCompletion(intP,p,rk,stats,S,Rtotal,tol)
--THIS IS TOO LONG!!!
apply(flatten entries gens minors(2,subs),i->sub(i,RR))
eigenvalues sub(subs,QQ)


--MLE for rk 1 matrices
empiricalMLEexistence(1,500,K) --101,125
empiricalMLEexistence(2,500,K) --500



--4-cycle uncolored
restart
load "functions.m2"
R=QQ[l_1..l_8]
K=matrix{{l_1,l_5,0,l_6},{l_5,l_2,l_7,0},{0,l_7,l_3,l_8},{l_6,0,l_8,l_4}}
--Test ML degree for different ranks
p=4
n=1
X=random(QQ^p,QQ^n);              
S=(1/n)*X*transpose(X);
rank S
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
--netList (trim I)_*
dim I, degree I
--(5,10)
J=saturate(I,det K);
dim J, degree J
--n=4 MLdeg=5
--n=3 MLdeg=5
--n=2 MLdeg=mostly 4, rarely 3
--n=1 MLdeg=1

