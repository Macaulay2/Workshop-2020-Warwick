
--GRAPH 14: 2+2 vertices equal (opposite) and no edges equal

restart
load "functions.m2"
R=QQ[l_1..l_6]
K=matrix{{l_1,l_3,0,l_6},{l_3,l_2,l_4,0},{0,l_4,l_1,l_5},{l_6,0,l_5,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*


-- Computation of the elimination ideals
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG3=sub(rankProjection(stats,3,S),ring(suffStat(K)))
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
netList IG1_*
toString IG1 -- ideal(t_6*t_4-t_5*t_3,t_6^2+t_5^2+t_4^2+t_3^2-4*t_2*t_1)

I1_G11=ideal(t_5^2+t_4^2-2*t_5*t_3+t_3^2-4*t_2*t_1)
isSubset(I1_G11,IG1)
isSubset(IG1,I1_G11)
g=(sub(I1_G11,{t_4=>t_4+t_6}))_0
g==IG1_1+2*IG1_0
--IG1_0 is a 2-minor, hence it vanishes in rank 1 (didn't happen before because t4=t4+t6)

-- Empirical analysis of sign of the generators of IG1
m=4 --rank of empirical covariance matrices
k=10000 --number of points we want to try out
L=empiricalVanishingPolynomials(IG1,m,k,p,n,stats) --IG1_0 changes signs
(L1,L2)=differentSign(IG1_1,m,k,p,n,stats);
netList L1 --initial point (suff statistics, evaluation)
netList L2 --if exists, point with different sign (suff statistics, evaluation, number of points we had to try)

m=1
U=random(QQ^4,QQ^m)
V=U*transpose(U)
rank V

scoord=matrix {flatten toList apply(0..p-1,i->toList apply(i..p-1,j->-V_(i,j)))}|matrix {toList apply(1..n,i->0)}
auxcoord=sub(stats,scoord)
tcoord=toList apply(1..n,i->t_i=>auxcoord_(i-1))
sub(f,tcoord)

--Infeasibility test (of the vanishing of the single generator of IG1)
RL=QQ[t_1..t_6,u11,u12,u13,u14,u22,u23,u24,u33,u34,u44]
LT=matrix{{u11,0,0,0},{u12,u22,0,0},{u13,u23,u33,0},{u14,u24,u34,u44}}
PSD=LT*transpose(LT)
f=sub(IG1_1,RL)
fu=sub(f,{t_1=>PSD_(0,0)+PSD_(2,2),t_2=>PSD_(1,1)+PSD_(3,3),t_3=>2*PSD_(0,1),t_4=>2*PSD_(1,2),t_5=>2*PSD_(2,3),t_6=>2*PSD_(0,3)})

loadPackage "SemidefiniteProgramming"
loadPackage "SumsOfSquares"
sol=solveSOS (-fu)
peek sol
s=sosPoly sol
peek s
coefs=s#coefficients
gene=s#generators
netList gene
prod=apply(coefs,gene,(i,j)->i*j^2)
sos=sum(prod)
-fu==sos --true: IG1_1 is a NSOS


--Computation of boundary components and algebraic boundary
restart
load "functions.m2"
R=QQ[l_1..l_6]
K=matrix{{l_1,l_3,0,l_6},{l_3,l_2,l_4,0},{0,l_4,l_1,l_5},{l_6,0,l_5,l_2}}

(V,n,K2)=embeddedK(K);
B3=time boundaryComponents(K2,4,n)
-- {}
B2=time boundaryComponents(K2,3,n)
netList B2                         --we should check that all of them are actually needed
                                   -- i.e. there's a real point in K_G (in the dual)
				   -- Is this equivalent to check that there's a real point here coming from a PSD matrix?
B1=time boundaryComponents(K2,2,n)
-- {}
algBoundary(V,n,K2)
netList oo 

IG1=ideal(t_6*t_4-t_5*t_3,t_6^2+t_5^2+t_4^2+t_3^2-4*t_2*t_1);

isSubset(B2_2,IG1) 
isSubset(B2_0,IG1) 
isSubset(B2_1,IG1) 

(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal);

--check change of sign for the other components of the algebraic boundary
m=4 --rank of empirical covariance matrices
k=10 --number of points we want to try out
--L=empiricalVanishingPolynomials(sub(B2_1,ring(suffStat(K))),m,k,p,n,stats)
--doesn't vanish (it's -SOS, see below)
L=empiricalVanishingPolynomials(sub(B2_0,ring(suffStat(K))),m,k,p,n,stats)
--does vanish (not SOS)

RL=QQ[t_1..t_6,u11,u12,u13,u14,u22,u23,u24,u33,u34,u44]
LT=matrix{{u11,0,0,0},{u12,u22,0,0},{u13,u23,u33,0},{u14,u24,u34,u44}}
PSD=LT*transpose(LT)
f=sub((B2_0)_0,RL)
f=sub((B2_1)_0,RL)
fu=sub(f,{t_1=>PSD_(0,0)+PSD_(2,2),t_2=>PSD_(1,1)+PSD_(3,3),t_3=>2*PSD_(0,1),t_4=>2*PSD_(1,2),t_5=>2*PSD_(2,3),t_6=>2*PSD_(0,3)})

loadPackage "SemidefiniteProgramming"
loadPackage "SumsOfSquares"
sol=solveSOS (-fu)
peek sol
s=sosPoly sol
peek s
coefs=s#coefficients
gene=s#generators
netList gene
prod=apply(coefs,gene,(i,j)->i*j^2)
sos=sum(prod)
-fu==sos --true


-----------------
restart
load "functions.m2"
R=QQ[l_1..l_6]
K=matrix{{l_1,l_3,0,l_6},{l_3,l_2,l_4,0},{0,l_4,l_1,l_5},{l_6,0,l_5,l_2}}

(V,n,K2)=embeddedK(K);

B2=time boundaryComponents(K2,3,n)
 -- used 0.671875 seconds
netList B2                         --we should check that all of them are actually needed

I2=minors(3,K2);
L2=minimalPrimes(I2);
netList L2
dualVariety(L2_0,n,l,t)
dualVariety(L2_1,n,l,t)
dualVariety(L2_2,n,l,t)
dualVariety(L2_3,n,l,t)
dualVariety(L2_4,n,l,t)

netList (L2_0)_* 
K20=matrix{{l_1,l_3,0,l_4},{l_3,l_2,l_4,0},{0,l_4,l_1,-l_3},{l_4,0,-l_3,l_2}}
trim(minors(3,K20)+ideal{l_1*l_2-l_5^2-l_6^2})
ex20=sub(K20,{l_1=>1,l_2=>0,l_3=>0,l_4=>0,l_5=>0,l_6=>0}) --rk 2 PSD matrix
netList unique detPrincipalMinors(K20)
factor det K20
--l1,l2>=0, l1l2=l3^2+l4^2

netList (L2_1)_* 
K21=matrix{{l_1,l_3,0,-l_4},{l_3,l_2,l_4,0},{0,l_4,l_1,l_3},{-l_4,0,l_3,l_2}}
trim(minors(3,K21)+ideal{l_1*l_2-l_3^2})
ex21=sub(K21,{l_1=>1,l_2=>0,l_3=>0,l_4=>0,l_5=>0,l_6=>0}) --rk 2 PSD matrix
netList unique detPrincipalMinors(K20)
factor det K20
--l1,l2>=0, l1l2=l3^2+l4^2

principalMinors(K21)
netList unique detPrincipalMinors(K21)
factor det K21 
--they have to satisfy that l1,l2>=0 and l1*l2>=l3^2, but l1*l2-l3^2=0
--hence l1,l2>=0 and l3^3=l1l2
eigenvalues sub(sub(K21,{l_1=>1,l_2=>1,l_3=>1}),QQ)

L2_2
K22=matrix{{0,l_3,0,l_6},{l_3,0,l_4,0},{0,l_4,0,l_5},{l_6,0,l_5,0}}
trim(minors(3,K22)+ideal{l_3*l_5-l_4*l_6})

-- compute list of ppal minors 
principalMinors=M->(
--M=sub(M,QQ[support M]);
mm:=numcols M;
ss:=drop(subsets mm,-1);
for s in ss list submatrix'(M,s,s)
);

-- compute list of determinants of ppal minors 
detPrincipalMinors=M->(
L:=principalMinors M;
for l in L list det l
);

principalMinors(K22)
unique detPrincipalMinors(K22)
factor det K22 
-- No PSD matrix belongs to this component

