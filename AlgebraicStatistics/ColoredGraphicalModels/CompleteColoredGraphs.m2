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

(V,n,K2)=embeddedK(K);
algBoundary(V,n,K2)
boundaryComponents(K2,3,n)
boundaryComponents(K2,2,n)
boundaryComponents(K2,1,n)

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
empiricalMLEexistence(1,500,K) --333,351, 385 with real coefficients
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
eigenvalues S
minors(2,S)

I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
netList I_*
dim I, degree I
--(1,4)
J=saturate(I,det K);
dim J, degree J
--n=3 MLdeg=4
--n=2 MLdeg=4
--n=1 MLdeg=3


restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_2,l_3},{l_2,l_1,l_4},{l_3,l_4,l_1}}

(V,n,K2)=embeddedK(K);
algb=algBoundary(V,n,K2)
netList algb
boundaryComponents(K2,3,n)
boundaryComponents(K2,2,n)
boundaryComponents(K2,1,n) --I don't think it makes sense
minors(1,K)
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
empiricalMLEexistence(1,500,K)  --0, 92 with real coefficients
empiricalMLENoPD(1,500,K)  --500
empiricalMLERealReg(1,500,K)  --499

--Test ML degree for different ranks
restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_4,l_4},{l_4,l_2,l_4},{l_4,l_4,l_3}}
p=3
n=1
--X=random(QQ^p,QQ^n);
X=random(RR^p,RR^n);              

--X=transpose matrix {toList apply(0..p-1,i->lift(X_(i,0),QQ))}
X=transpose matrix {toList apply(0..p-1,i->promote(X_(i,0),QQ))}

S=(1/n)*X*transpose(X);
--S=matrix toList apply(0..2,i->toList apply(0..2,j->lift(S_(i,j),QQ)))
rank S
eigenvalues S
netList apply(flatten entries gens minors(2,S),i->sub(i,RR))
netList apply(flatten entries gens minors(3,S),i->sub(i,RR))

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

restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_4,l_4},{l_4,l_2,l_4},{l_4,l_4,l_3}}

(V,n,K2)=embeddedK(K);
algb=algBoundary(V,n,K2)
netList algb
boundaryComponents(K2,3,n)
boundaryComponents(K2,2,n)


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
L=empiricalVanishingPolynomials(IG1,m,k,p,n,stats) --ONLY 1 POLYNOMIAL CHANGES SIGNS
(L1,L2)=differentSign(IG1_0,m,k,p,n,stats);
netList L1 --initial point (suff statistics, evaluation)
netList L2 --if exists, point with different sign (suff statistics, evaluation, number of points we had to try)
v1=vector(flatten entries gens L1_0)
v2=vector(flatten entries gens L2_0)
P=interiorPoint(v1,v2,n,L_0,K)
intP=P_0
use ring IG1
sub(sub(IG1,{t_1=>intP_0,t_2=>intP_1,t_3=>intP_2,t_4=>intP_3}),RR) 


--Rank-1 completion: NOT POSSIBLE because interior points do not lie in IG1

--MLE for rk 1 matrices
empiricalMLEexistence(1,1,K) --WRONG DIMENSION


--Test ML degree for different ranks
restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_2,l_3},{l_2,l_1,l_4},{l_3,l_4,l_5}}
p=3
n=1
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

restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_2,l_3},{l_2,l_1,l_4},{l_3,l_4,l_5}}

(V,n,K2)=embeddedK(K);
algb=algBoundary(V,n,K2)
netList algb
boundaryComponents(K2,3,n)
boundaryComponents(K2,2,n)


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

(V,n,K2)=embeddedK(K);
algb=algBoundary(V,n,K2)
netList algb
boundaryComponents(K2,3,n)
boundaryComponents(K2,2,n)


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
--No interior points

--MLE for rk 1 matrices
empiricalMLEexistence(1,100,K) --score eq ideal doesn't have the right dim

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
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*


--Algebraic boundary
(V,n,K2)=embeddedK(K);
algb=algBoundary(V,n,K2)
netList algb
boundaryComponents(K2,3,n)
boundaryComponents(K2,2,n)


--Elimination ideal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
netList IG1_*

--Intersection of the algebraic boundary with the interior
m=3 --rank of empirical covariance matrices
k=100 --number of points we want to try out
use ring IG1
L=empiricalVanishingPolynomials(IG1,m,k,p,n,stats)

empiricalMLEexistence(1,500,K)--109

--Test ML degree for different ranks
restart
p=3
n=1
X=random(QQ^p,QQ^n);              
S=(1/n)*X*transpose(X);
rank S
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
netList I_*
dim I, degree I
--(1,5)
J=saturate(I,det K);
dim J, degree J
--n=3 MLdeg=3
--n=2 MLdeg=3
--n=1 MLdeg=2


-----------------------------------------------------------
--3-cycle with 2 equal vertices with 2 equal edges not including the one between the vertices
-----------------------------------------------------------
restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_2,l_2},{l_2,l_4,l_3},{l_2,l_3,l_4}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

--Algebraic boundary
(V,n,K2)=embeddedK(K);
algb=algBoundary(V,n,K2)
netList algb
boundaryComponents(K2,3,n)
boundaryComponents(K2,2,n)


--Elimination ideal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
netList IG1_*

IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2

--Intersection of the algebraic boundary with the interior
m=3 --rank of empirical covariance matrices
k=1000 --number of points we want to try out
use ring IG1
L=empiricalVanishingPolynomials(IG1,m,k,p,n,stats)

empiricalMLEexistence(1,500,K)


--Test ML degree for different ranks
p=3
n=3
X=random(QQ^p,QQ^n);              
S=(1/n)*X*transpose(X);
rank S
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
netList I_*
dim I, degree I
--(1,5)
J=saturate(I,det K);
dim J, degree J
--n=3 MLdeg=3
--n=2 MLdeg=3
--n=1 MLdeg=2


-----------------------------------------------------------
--3-cycle with 2 equal vertices and 3 equal edges 
-----------------------------------------------------------
restart
load "functions.m2"
R=QQ[l_1..l_3]
K=matrix{{l_1,l_3,l_3},{l_3,l_1,l_3},{l_3,l_3,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

--Algebraic boundary
(V,n,K2)=embeddedK(K);
algb=algBoundary(V,n,K2)
netList algb
boundaryComponents(K2,3,n)
boundaryComponents(K2,2,n)


--Elimination ideal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
IG1

IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2


--Test ML degree for different ranks
p=3
n=1
X=random(QQ^p,QQ^n);              
S=(1/n)*X*transpose(X);
rank S
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
netList I_*
dim I, degree I
--(1,2)
J=saturate(I,det K);
dim J, degree J
--n=3 MLdeg=2
--n=2 MLdeg=2
--n=1 MLdeg=2
 -----------------------------------------------------------
--3-cycle with 3 equal vertices and 2 equal edges 
-----------------------------------------------------------
restart
load "functions.m2"
R=QQ[l_1..l_3]
K=matrix{{l_1,l_2,l_3},{l_2,l_1,l_3},{l_3,l_3,l_1}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

--Algebraic boundary
(V,n,K2)=embeddedK(K);
algb=algBoundary(V,n,K2)
netList algb
boundaryComponents(K2,3,n)
boundaryComponents(K2,2,n)


--Elimination ideal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
IG1

IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2


--Test ML degree for different ranks
p=3
n=3
X=random(QQ^p,QQ^n);              
S=(1/n)*X*transpose(X);
rank S
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
netList I_*
dim I, degree I
--(1,2)
J=saturate(I,det K);
dim J, degree J
--n=3 MLdeg=2
--n=2 MLdeg=2
--n=1 MLdeg=2


 -----------------------------------------------------------
--3-cycle with 3 equal vertices and 3 equal edges 
-----------------------------------------------------------
restart
load "functions.m2"
R=QQ[l_1,l_2]
K=matrix{{l_1,l_2,l_2},{l_2,l_1,l_2},{l_2,l_2,l_1}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

--Algebraic boundary
(V,n,K2)=embeddedK(K);
algb=algBoundary(V,n,K2)
netList algb
boundaryComponents(K2,3,n)
boundaryComponents(K2,2,n)

--Elimination ideal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
IG1

IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2


--Test ML degree for different ranks
p=3
n=1
X=random(QQ^p,QQ^n);              
S=(1/n)*X*transpose(X);
rank S
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
netList I_*
dim I, degree I
--(1,1)
J=saturate(I,det K);
dim J, degree J
--n=3 MLdeg=1
--n=2 MLdeg=1
--n=1 MLdeg=1



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

-----------Test on data generation on the 3-cycle with all edges equal
restart
R=QQ[x,y,z]
U=matrix{{x,y,z},{y,z,x}}
X=transpose(U)*U
rank X
--principal 2-minors
X_(0,0)*X_(1,1)-X_(0,1)^2
X_(0,0)*X_(2,2)-X_(2,2)^2
X_(1,1)*X_(2,2)-X_(1,2)^2
loadPackage "SemidefiniteProgramming"
loadPackage "SOS"


restart
R=QQ[x_1..x_3,s_11,s_12,s_13,s_22,s_23,s_33]
U=matrix{{x_1,x_2,x_3}}
X=transpose(U)*U
rank X
S=genericSymmetricMatrix(R,s_11,3)

J=ideal{s_12-x_1*x_2,s_13-x_1*x_3,s_23-x_2*x_3,s_11+s_22+s_33-(x_1^2+x_2^2+x_3^2)}
netList J_*
codim J, degree J

I=eliminate({x_1,x_2,x_3},J)
netList I_*
codim I, degree I

isSubset()
