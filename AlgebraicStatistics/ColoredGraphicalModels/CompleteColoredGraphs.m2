-----------------------------------------------------------
--3-cycle uncolored
-----------------------------------------------------------
--Test ML degree for different ranks
restart
load "functions.m2"
R=QQ[l_1..l_6]
K=matrix{{l_1,l_4,l_5},{l_4,l_2,l_6},{l_5,l_6,l_3}}
p=3
n=1
X=random(QQ^p,QQ^n);              
S=(1/n)*X*transpose(X);
rank S
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
netList I_*
dim I, degree I
--(3,4)
J=saturate(I,det K);
dim J, degree J
--n=3 MLdeg=1
--n=2 MLdeg=0
--n=1 MLdeg=0

-----------------------------------------------------------
--3-cycle with 3 vertices equal
-----------------------------------------------------------
restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_2,l_3},{l_2,l_1,l_4},{l_3,l_4,l_1}}
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
m=3 --rank of empirical covariance matrices
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
intP=P_1
tol=0.000000000001 
subs=rankCompletion(intP,p,rk,stats,S,Rtotal,tol)
use ring IG1
sub(sub(IG1_0,{t_1=>intP_0,t_2=>intP_1,t_3=>intP_2,t_4=>intP_3}),RR) 
apply(flatten entries gens minors(2,subs),i->sub(i,RR))
eigenvalues sub(subs,QQ)


--MLE for rk 1 matrices
empiricalMLEexistence(1,500,K) --333,351
empiricalMLEexistence(2,500,K) --500


--Test ML degree for different ranks
restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_2,l_3},{l_2,l_1,l_4},{l_3,l_4,l_1}}
p=3
n=1
X=random(QQ^p,QQ^n);              
S=(1/n)*X*transpose(X);
rank S
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
netList I_*
dim I, degree I
--(1,4)
J=saturate(I,det K);
dim J, degree J
--n=3 MLdeg=4
--n=2 MLdeg=4
--n=1 MLdeg=3

-----------------------------------------------------------
--3-cycle with 3 edges equal
-----------------------------------------------------------
restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_4,l_4},{l_4,l_2,l_4},{l_4,l_4,l_3}}
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
m=3 --rank of empirical covariance matrices
k=100 --number of points we want to try out
use ring IG1
L=empiricalVanishingPolynomials(IG1,m,k,p,n,stats)
(L1,L2)=differentSign(IG1_0,m,k,p,n,stats);
netList L1 --initial point (suff statistics, evaluation)
netList L2 --if exists, point with different sign (suff statistics, evaluation, number of points we had to try)

v1=vector(flatten entries gens L1_0)
v2=vector(flatten entries gens L2_0)
P=interiorPoint(v1,v2,n,L_0,K)

--Are the interior points in the singular locus of IG1?
intP=P_0
IG1=sub(IG1,QQ[t_1..t_4])
sing=trim (IG1+ideal flatten entries jacobian IG1)
netList sing_*
sub(sub(sing,{t_1=>intP_0,t_2=>intP_1,t_3=>intP_2,t_4=>intP_3}),RR) 
sub(sub(IG1,{t_1=>intP_0,t_2=>intP_1,t_3=>intP_2,t_4=>intP_3}),RR) 


--Rank-1 completion
intP=P_0
use ring IG1
sub(sub(IG1_0,{t_1=>intP_0,t_2=>intP_1,t_3=>intP_2,t_4=>intP_3}),RR) 
--Depending on the tolerance the rank completion exists or not     
tol=0.0000001 
subs=rankCompletion(intP,p,rk,stats,S,Rtotal,tol) --as we increase the tolerance, it stops working!!!
apply(flatten entries gens minors(2,subs),i->sub(i,RR))
eigenvalues sub(subs,QQ)

--MLE for rk 1 matrices
empiricalMLEexistence(1,500,K)  --0
empiricalMLENoPD(1,500,K)  --500
empiricalMLERealReg(1,500,K)  --499

--Test ML degree for different ranks
restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_4,l_4},{l_4,l_2,l_4},{l_4,l_4,l_3}}
p=3
n=1
X=random(QQ^p,QQ^n);              
S=(1/n)*X*transpose(X);
rank S
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
netList I_*
dim I, degree I
J=saturate(I,det K);
dim J, degree J
--n=3 MLdeg=4
--n=2 MLdeg=4
--n=1 MLdeg=3


criticalPoints=zeroDimSolve(J);
criticalMatrices=genListMatrix(criticalPoints,K);
netList criticalMatrices	
apply(criticalMatrices,i->eigenvalues i)

checkPD(criticalMatrices)
checkRealReg(criticalMatrices)

-- Multidegrees
RT=QQ[l_1..l_4]**QQ[t_1..t_4]
K=matrix{{l_1,l_4,l_4},{l_4,l_2,l_4},{l_4,l_4,l_3}}
tr=sum(toList apply(1..4,i->l_i*t_i))

I=ideal submatrix'(jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{tr}}),{4,5,6,7},)
netList I_*
dim I, degree I
multidegree I
J=time saturate(I,sub(ideal det K,RT));
netList J_*
dim J,degree J
multidegree J

IG1=sub(IG1,RT);
Irk1=I+IG1;
netList Irk1_*
betti (trim Irk1)
dim Irk1,degree Irk1
multidegree Irk1
Jrk1=time saturate(Irk1,sub(ideal det K,RT));
netList Jrk1_*
dim Jrk1,degree Jrk1
multidegree Jrk1

-----------------------------------------------------------
--3-cycle with 2 vertices equal
-----------------------------------------------------------
restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_2,l_3},{l_2,l_1,l_4},{l_3,l_4,l_5}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

--Elimination ideal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
netList IG1_*

--Intersection of the algebraic boundary with the interior
m=3 --rank of empirical covariance matrices
k=100 --number of points we want to try out
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
apply(flatten entries gens minors(2,subs),i->sub(i,RR))
eigenvalues sub(subs,QQ)


--MLE for rk 1 matrices
empiricalMLEexistence(1,1,K)

-----------------------------------------------------------
--3-cycle with 2 edges equal
-----------------------------------------------------------
restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_4,l_4},{l_4,l_2,l_5},{l_4,l_5,l_3}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

--Elimination ideal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
netList IG1_*
sing=trim ideal singularLocus IG1;
netList sing_*
betti sing
codim IG1,degree IG1
codim sing,degree sing

use R
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2

--Intersection of the algebraic boundary with the interior
m=3 --rank of empirical covariance matrices
k=100 --number of points we want to try out
use ring IG1
L=empiricalVanishingPolynomials(IG1,m,k,p,n,stats)
--NONE OF THE GENS CHANGES SIGNS!!!
(L1,L2)=differentSign(IG1_0,m,k,p,n,stats);
netList L1 --initial point (suff statistics, evaluation)
netList L2 --if exists, point with different sign (suff statistics, evaluation, number of points we had to try)

v1=vector(flatten entries gens L1_0)
v2=vector(flatten entries gens L2_0)
P=interiorPoint(v1,v2,n,L_0,K)


--Rank-1 completion
intP=P_0
tol=0.000001 
subs=rankCompletion(intP,p,rk,stats,S,Rtotal,tol) --Works!
apply(flatten entries gens minors(2,subs),i->sub(i,RR))
eigenvalues sub(subs,QQ)


--MLE for rk 1 matrices
empiricalMLEexistence(1,100,K) --score eq ideal doesn't have the right dim



X=random(QQ^3,QQ^1);              
S=(1/3)*X*transpose(X);
rank S
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
netList I_*
dim I, degree I

J=saturate(I,det K);
dim J, degree J
criticalPoints=zeroDimSolve(J);
criticalMatrices=genListMatrix(criticalPoints,K);
netList criticalMatrices	
apply(criticalMatrices,i->eigenvalues i)

checkPD(criticalMatrices)
checkRealReg(criticalMatrices)


--Test ML degree for different ranks
restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_4,l_4},{l_4,l_2,l_5},{l_4,l_5,l_3}}
p=3
n=2
X=random(QQ^p,QQ^n);              
S=(1/n)*X*transpose(X);
rank S
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
netList I_*
dim I, degree I
--(2,4)
J=saturate(I,det K);
dim J, degree J
--n=3 MLdeg=1
--n=2 MLdeg=1
--n=1 MLdeg=0
-----------------------------------------------------------
--3-cycle with 2 equal vertices with 2 equal edges including the one between the vertices
-----------------------------------------------------------
restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_2,l_3},{l_2,l_4,l_3},{l_3,l_3,l_4}}
empiricalMLEexistence(1,10,K)


-----------------------------------------------------------
--3-cycle with 2 equal vertices with 2 equal edges not including the one between the vertices
-----------------------------------------------------------
restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_2,l_2},{l_2,l_4,l_3},{l_2,l_3,l_4}}
empiricalMLEexistence(1,10,K)


-----------------------------------------------------------
--2-CLIQUES
-----------------------------------------------------------
restart
load "functions.m2"
R=QQ[l_1,l_2]
K=matrix{{l_1,l_2},{l_2,l_1}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

--Elimination ideal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
IG1
netList IG1_*


-----------------------------------------------------
--UNDERSTANDING RANG DEFFICIENT MATRICES IN K_G
-----------------------------------------------------

--all vertices same color
R=QQ[a,b,c]
K=matrix{{1/3,a,b},{a,1/3,c},{b,c,1/3}}
I1=minors(2,K)
dim I1,degree I1
sols=zeroDimSolve I1
netList sols

k=3
M=sub(K,{a=>(coordinates sols_k)_0,b=>(coordinates sols_k)_1,c=>(coordinates sols_k)_2})
eigenvalues M

I2=minors(3,K)
dim I2,degree I2


--all edges same color

R=QQ[a,b,c]
K=matrix{{a,c,c},{c,b,c},{c,c,1-a-b}}
I1=minors(2,K)
dim I1,degree I1
sols=zeroDimSolve I1
netList sols

k=3
M=sub(K,{a=>(coordinates sols_k)_0,b=>(coordinates sols_k)_1,c=>(coordinates sols_k)_2})
eigenvalues M


I2=minors(3,K)
dim I2,degree I2
