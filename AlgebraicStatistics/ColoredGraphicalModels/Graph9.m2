-- Graph 9 Uhler
-- auxiliary functions (Olga's functions)

dualVariety=(IX,n,x,u)->(
    c:=codim IX;
    JacX:=diff(matrix{toList(x_1..x_n)}, transpose gens IX);
    AugJacX:=matrix{toList(u_1..u_n)}||JacX;
    SingX:=IX+minors(c,jacobian IX);
    conormalX:=saturate(IX+minors(c+1,AugJacX),SingX);
    dualX:=eliminate(toList(x_1..x_n),conormalX);
    dualX
    )

boundaryComponents=(K,p)->(
    I:=minors(p,K);
    minPrimes:=minimalPrimes I;
    m:= length minPrimes;
    allComponents:=for i to  m-1 list (dualComponent:=dualVariety(minPrimes_i,n,l,t),numgens trim dualComponent);
    boundaryComponents:= for i in allComponents list (if i_1==1 then i_0)
    )
algBoundary=(K)->(
    s=numgens target K;       
    delete (null, flatten (for p from 1 to s list boundaryComponents(K,p))) 
    )

---------------------------------------------------------------------------------------------------------------------------
--GRAPH 9
--STUDY THE SIGN OF THE IRREDUCIBLE FACTORS OF THE ALGEBRAIC BOUNDARY ON THE INTERIOR OF THE CONE OF SUFFICIENT STATISTICS
---------------------------------------------------------------------------------------------------------------------------

restart
n=6
R=QQ[l_1..l_n,t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
K=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}
numcols K
Sigma=genericSymmetricMatrix(R,s_11,4)

--conmputation in Mathematica shows that there are no points of rank 1 on the boundary
-- of the cone of concentration matrices (one can also directly see that if K has rank at most
--1, then it is necessary the 0 matrix)
boundaryComponents(K,2) -- not necessary!!!
bdry2=boundaryComponents(K,3) -- corresponds to I_G,1-- sum of dimenesions n-1 because we're working projectively
boundaryComponents(K,4)-- corresponds to I_G,2=0
algBoundary(K)

I2=minors(3,K)
minPrimes2=minimalPrimes I2
isSubset(minPrimes2_1,(minPrimes2_0+minPrimes2_3))
for i to  (length minPrimes2)-1 list dualVariety(minPrimes2_i,n,l,t)

I_G1=ideal(4*t_2*t_3-t_5^2-t_6^2,t_3*t_4^2-t_1*t_5^2-2*t_1*t_5*t_6-t_1*t_6^2) --irreducible
isSubset(bdry2_1,(bdry2_0+bdry2_3))
---test with new functions
load "functions.m2"
(V,n,K2)=embeddedK(K)
n
gens ring K2
gens R

algBoundary(V,n,K2)

--empiricalVanishingPolynomials
restart
naiveMatrixProduct=(M1,M2)->(
    m=numcols M1-1;
    matrix toList apply(0..m, i -> toList apply(0..m, j -> M1_(i,j)*M2_(i,j))))

R=QQ[l_1..l_6]
K=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}

p=numcols K
n=#support K
w=flatten toList apply(1..p, i -> toList apply(i..p, j -> (i,j)));
v= apply (w, ij -> s_ij);
Rtotal=QQ[v,reverse toList apply(1..n, i->t_i)]
--gens Rtotal
--toString gens Rtotal
S=genericSymmetricMatrix(Rtotal,s_(1,1),p)
--support S
M=mutableMatrix {{apply(1..n,i->0)}}
stat=ideal{}
for i from 0 to n-1 do (
    M_(0,i)=1;
    aux=matrix M;
    stat=stat+(t_(i+1)-sum (flatten entries naiveMatrixProduct(sub(K,aux),S)));
    M_(0,i)=0;
)
netList stat_*

adjacencyMat=(A)->(
    matrix toList apply(0..(p-1),i->toList apply(0..(p-1),j-> if A_(i,j)!=0 then 1 else 0))
)


M=mutableMatrix {{apply(1..n,i->0)}}
for i from 1 to n do (
    M_(0,i-1)=1;
    aux=matrix M;
    A=adjacencyMat naiveMatrixProduct(sub(K,aux),S);
    print(A,sum(flatten entries A));
    M_(0,i-1)=0;
)
coord


restart
load "functions.m2"
R=QQ[l_1..l_6]
K=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}
(p,n,Rtotal,S)=coloredData(K);
stats=sub(suffStat(K),Rtotal);
ring(suffStat(K))
netList stats_*

IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
netList IG1_*

m=4
k=10
L=empiricalVanishingPolynomials(IG1,m,k,p,n,stats)
(L1,L2)=differentSign(IG1_1,m,k,p,n,stats)


empiricalVanishingPolynomials=(f,m,k,p,n,stats)->(
     X:=random(QQ^p,QQ^m);    
     SE:=(1/m)*X*transpose(X);
     aux:=matrix {join(flatten toList apply(0..(p-1),i->toList apply(i..(p-1),j->-SE_(i,j))),
	toList apply(1..n,i->0))};
     suffE:=sub(stats,aux);
     subs:=toList apply(1..n,i->t_i=>suffE_(i-1));
     eval:=sub(f,subs);
     --initialize vars
     sign:=1;
     suffEnew:=suffE;
     evalnew:=eval;
     i:=0;
     while (i<k and sign==1) list (
          i=i+1;
          X=random(QQ^p,QQ^m);
          SE=(1/m)*X*transpose(X);
          aux=matrix {join(flatten toList apply(0..(p-1),i->toList apply(i..(p-1),j->-SE_(i,j))),
	   toList apply(1..n,i->0))};
          suffEnew=sub(stats,aux);
          subs=toList apply(1..n,i->t_i=>suffEnew_(i-1));
          evalnew=sub(f,subs);
          if sub(evalnew*eval,QQ)<0 then sign=-1;
          if sign==1 then continue; (suffEnew,evalnew,i));
     return(suffE,eval) 
     )


m=4
r=1
p=numcols K
X=random(QQ^p,QQ^m)
SE=(1/m)*X*transpose(X)
toList apply(1..n,i->0)
reverse toList apply(1..n, i->t_i)
aux=matrix {join(flatten toList apply(0..(p-1),i->toList apply(i..(p-1),j->-SE_(i,j))),
	toList apply(1..n,i->0))}
suffE=sub(stats,aux)
subs=toList apply(1..n,i->t_i=>suffE_(i-1))
eval=sub(IG1_r,subs)

sign=1

--while p list x do z
--while i < 10 list (i = i+1; if odd i then continue; i^2)

k=100
i=0
while (i<k and sign==1) list (
    i=i+1;
    X=random(QQ^p,QQ^m);
    SE=(1/m)*X*transpose(X);
    toList apply(1..n,i->0);
    reverse toList apply(1..n, i->t_i);
    aux=matrix {join(flatten toList apply(0..(p-1),i->toList apply(i..(p-1),j->-SE_(i,j))),
	toList apply(1..n,i->0))};
    suffEnew=sub(stats,aux);
    subs=toList apply(1..n,i->t_i=>suffEnew_(i-1));
    evalnew=sub(IG1_r,subs);
    if sub(evalnew*eval,QQ)<0 then sign=-1;
    if sign==1 then continue; (suffEnew,evalnew,i))

k=10
i=0
netList (while (i<k) list (
    i=i+1;
    X=random(QQ^p,QQ^m);
    SE=(1/m)*X*transpose(X);
    toList apply(1..n,i->0);
    reverse toList apply(1..n, i->t_i);
    aux=matrix {join(flatten toList apply(0..(p-1),i->toList apply(i..(p-1),j->-SE_(i,j))),
	toList apply(1..n,i->0))};
    suffE=sub(stats,aux);
    subs=toList apply(1..n,i->t_i=>suffE_(i-1));
    evalnew=sub(IG1_r,subs);
    if sub(evalnew*eval,QQ)<0 then sign=-1;
    if sign==1 then continue; (suffE,evalnew,sign)))



-- algebraic boundary 
HG=algBoundary(K);
netList HG
isPrime(HG_2)
toString (HG_2)_0

gens ring(K)
boundaryComponents(K,2)
boundaryComponents(K,3)
boundaryComponents(K,1)

I=HG_2+HG_3+HG_4
netList I_*

I==ideal{HG_2,HG_4}
I==ideal{HG_3,HG_4}
netList (flatten entries gens gb I)

-- EXTENSION THEOREM ------------------------------------------------------------------------
RE=QQ[drop(gens R,6),t_8,t_7]
--sufficient statistics
varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
Istat=ideal(t_1-s_11,t_2-s_22-s_44,t_3-s_33,t_4-2*(s_12+s_14),t_5-2*s_23,t_6-2*s_34,t_7-2*s_13,t_8-2*s_24)
Istatorig=ideal(t_1-s_11,t_2-s_22-s_44,t_3-s_33,t_4-2*(s_12+s_14),t_5-2*s_23,t_6-2*s_34)
S=sub(Sigma,RE)
IG1=eliminate(varList,minors(2,S)+Istat)
IG1orig=eliminate(varList,minors(2,S)+Istatorig)
IG1orig==sub(I,RE) --true

eliminate({t_7,t_8},I_G1)==I_G1orig
netList I_G1_*
netList I_G1orig_*
flatten entries gens I_G1orig
idealHG=sub(product HG,RE)
isSubset(idealHG,I_G1orig)  --- V(I_G1orig) is contained in the algebraic boundary
(I_G1orig)_0 % gb idealHG
(I_G1orig)_1 % gb idealHG

f1=t_6^2+t_5^2-4*t_3*t_2
f2=t_4^2*t_3-2*t_6*t_5*t_1-4*t_3*t_2*t_1

J=ideal{f1,f2}
I==J
------------------------------------------------------------------------------------------
--Chose your favourite rank n and the number of elements k you want to try
--This procedure prints out the the values of the sufficient statistics for k random matrices of rank n
n=4
k=5
for i from 1 to k do (
--    X=random(QQ^4,QQ^n);              
    X=sub(random(ZZ^4,ZZ^n),QQ);
    S=(1/n)*X*transpose(X);
    aux={t_1=>S_(0,0),t_2=>S_(1,1)+S_(3,3),t_3=>S_(2,2),t_4=>2*(S_(0,1)+S_(0,3)),t_5=>2*S_(1,2),t_6=>2*S_(2,3)};
    print S;
    print toString aux;
--    print apply(0..4,i->sub(HG_i,aux)); 
--    print (sub(HG_2,aux))_0;
--    print(apply(flatten entries gens IG1orig, i->sub(i,aux)))
      print(apply(flatten entries gens J, i->sub(i,aux)))
    )

--Below we have two sample covariance matrices of rank 4 
-- whose sufficient statistics take different signs when 
-- evaluated at f3

--S=matrix{{251/4,131/4,20,23},{131/4,101/4,15,13},{20,15,41/4,4},{23,13,4,20}}
    
| 251/4 131/4 20   23 |
| 131/4 101/4 15   13 |
| 20    15    41/4 4  |
| 23    13    4    20 |

P=matrix{{251/4,131/4,20, 23},{131/4,101/4,15,13},{20, 15, 41/4, 4},{23,13,4,20}}
rank P

P1={t_1 => 251/4, t_2 => 181/4, t_3 => 41/4, t_4 => 223/2, t_5 => 30, t_6 => 8}
for i to 4 do print (sub(HG_i,P1))_0
--for i from 0 to 1 do print sub(J_i,P1)
v1=vector{251/4,181/4,41/4,223/2,30,8}

--S=matrix{{27/2,22,22,43/4},{22,175/4,37,99/4},{22,37,73/2,75/4},{43/4,99/4,75/4,61/4}}

| 27/2 22    22   43/4 |
| 22   175/4 37   99/4 |
| 22   37    73/2 75/4 |
| 43/4 99/4  75/4 61/4 |

P2={t_1 => 27/2, t_2 => 59, t_3 => 73/2, t_4 => 131/2, t_5 => 74, t_6 => 75/2}
for i to 4 do print (sub(HG_i,P2))_0
--for i from 0 to 1 do print sub(J_i,P2)
v2=vector{27/2,59,73/2,131/2,74,75/2}

--TRY IT OUT WITH A DIFFERENT CHOICE OF POINTS
| 61/4 57/4 12    45/2  |
| 57/4 37/2 57/4  26    |
| 12   57/4 63/2  135/4 |
| 45/2 26   135/4 109/2 |
P1={t_1 => 61/4, t_2 => 73, t_3 => 63/2, t_4 => 147/2, t_5 => 57/2, t_6 => 135/2}
for i to 4 do print (sub(HG_i,P1))_0
--for i from 0 to 1 do print sub(J_i,P1)
v1=vector{61/4,73,63/2,147/2,57/2,135/2}

| 87/2  15   34    169/4 |
| 15    19   33/2  45/2  |
| 34    33/2 115/4 143/4 |
| 169/4 45/2 143/4 185/4 |
P2={t_1 => 87/2, t_2 => 261/4, t_3 => 115/4, t_4 => 229/2, t_5 => 33, t_6 => 143/2}
for i to 4 do print (sub(HG_i,P2))_0
--for i from 0 to 1 do print sub(J_i,P2)
v2=vector{87/2,261/4,115/4,229/2,33,143/2}

--MORE CHOICHES OF POINTS
P1={t_1 => 57/4, t_2 => 79/2, t_3 => 35/2, t_4 => 103/2, t_5 => 21, t_6 => 69/2}
for i to 4 do print (sub(HG_i,P1))_0
--for i from 0 to 1 do print sub(J_i,P1)

v1=vector{57/4,79/2,35/2,103/2,21,69/2}
P2={t_1 => 35, t_2 => 209/4, t_3 => 81/2, t_4 => 213/2, t_5 => 56, t_6 => 60}
for i to 4 do print (sub(HG_i,P2))_0
--for i from 0 to 1 do print sub(J_i,P1)
v2=vector{35,209/4,81/2,213/2,56,60}


-- POINTS USEFUL FOR EXTENSION THEOREM
P1={t_1 => 55, t_2 => 329/4, t_3 => 35, t_4 => 365/2, t_5 => 68, t_6 => 109/2}
for i from 0 to 1 do print sub(J_i,P1)
v1=vector{55,329/4,35,365/2,68,109/2}
P2={t_1 => 25, t_2 => 277/4, t_3 => 131/4, t_4 => 95, t_5 => 111/2, t_6 => 58}
for i from 0 to 1 do print sub(J_i,P2)
v2=vector{25,277/4,131/4,95,111/2,58}
-- Point in the intersection of the segment between P1 and P2 and relevant
test={44.4548, 77.6804, 34.2091, 151.743, 63.6062, 55.7303}

--Define line through 2 points (with auxiliary variable)
restart
Rt=QQ[t_1..t_6,m]
v=sub(vector delete(m,gens Rt),Rt)
--run line of code for v1 and v2 above
P=sub(v1,Rt)
Q=sub(v2,Rt)
paramr=ideal flatten entries (v-P-m*(Q-P))
r=eliminate(m,paramr)

--Line through 2 points (in ring with variables t1...t6)
Rtt=QQ[t_1..t_6]
r=sub(r,Rtt)
betti trim r
dim r
sub(r,{t_1 => v1_0, t_2 => v1_1, t_3 => v1_2, t_4 => v1_3, t_5 => v1_4, t_6 => v1_5})
sub(r,{t_1 => v2_0, t_2 => v2_1, t_3 => v2_2, t_4 => v2_3, t_5 => v2_4, t_6 => v2_5})


-- Intersection of the line and the factor of the algebraic boundary (f3)
--f=t_3*t_4^2-t_1*t_5^2-2*t_1*t_5*t_6-t_1*t_6^2 --usual setting
f=t_4^2*t_3-2*t_6*t_5*t_1-4*t_3*t_2*t_1   --extension setting

I=ideal(f,r)
betti trim I
dim I
degree I

loadPackage "EigenSolver"
sols=zeroDimSolve I;
netList sols  -- 3 real solutions with first choice of points
              -- 1 real sol, 2 complex sols for second choice of points

coordinates sols_0
test=apply(coordinates sols_0,j->lift(j,QQ))

-- How close to vanishing at the line and the cubic are actually these solutions?
aux={}
k=0
for i from 1 to 6 do aux=append(aux,t_i=>sub((coordinates sols_k)_(i-1),CC));
for i from 1 to 6 do aux=append(aux,t_i=>(coordinates sols_k)_(i-1));
for i from 1 to 6 do aux=append(aux,t_i=>lift((coordinates sols_k)_(i-1),QQ));
for i from 1 to 6 do aux=append(aux,t_i=>sub(test_(i-1),CC));

netList aux
sub(r,aux)
sub(f,aux)

ROlga=QQ[t_8,t_7,t_6,t_5,t_4,t_3,t_2,t_1]
IG1orig=ideal{t_6^2+t_5^2-4*t_3*t_2,t_4^2*t_3-2*t_6*t_5*t_1-4*t_3*t_2*t_1}
IG1=ideal{t_4^2-4*t_8*t_1-4*t_2*t_1,
    t_7*t_4-2*t_6*t_1-2*t_5*t_1,t_6*t_5-2*t_8*t_3,t_6^2+t_5^2-4*t_3*t_2,
    t_7*t_6+t_7*t_5-2*t_4*t_3,t_7^2-4*t_3*t_1,2*t_8*t_7-t_6*t_4-t_5*t_4+2*t_7*t_2}

test=apply(test,i->sub(i,ROlga))
partialSol=sub(IG1orig,{t_1 =>test_0, t_2 =>test_1, t_3 => test_2, t_4 => test_3, t_5 => test_4, t_6 => test_5})
sub(IG1,{t_1 =>(coordinates sols_k)_0, t_2 =>(coordinates sols_k)_1, t_3 => (coordinates sols_k)_2, t_4 => (coordinates sols_k)_3, t_5 => (coordinates sols_k)_4, t_6 => (coordinates sols_k)_5})


------------------
-- How can we check whether these points are in the interior 
-- of the cone of algebraic statistics?
--------------------
-- Option 1: check if one of them is a convex combination of the other two
Pr=apply(flatten entries P,i->sub(i,RR))
Qr=apply(flatten entries Q,i->sub(i,RR))
netList sols
--sols_k could be a convex combination of P and Q because in all coordinates 
-- the value of sols_k is in between the corresponding value of P and Q

k=0
A=sub(vector apply(coordinates sols_k,i->lift(i,QQ)),Rt)
conv= ideal flatten entries (A-m*Q-(1-m)*P)
dim conv
-- doesn't seem to be the case... COULD IT BE A NUMERICAL ERROR? Because otherwise how would
-- it be possible that a point is in the line through P and Q, its coordinates have 
-- values in between but still is not a convex combination

-- Manual check: find m for the first equation and substitute
toString conv
--ideal((197/4)*m-2959052059/200230340,-(55/4)*m+2764866769/670123772,-(105/4)*m+120116117/15249496,46*m-719969693/52160278,-44*m+262528697/19884177,-(59/2)*m+133422799/15072717)
--(197/4)*m-2959052059/200230340 => m = 2959052059/9861344245
aux=sub(conv,m=>(113712381*4)/(7694582*197))
apply(0..5,i->sub(aux_i,RR))


--Option 2: find 4-rank pd matrix with this suff statistics via MLE computation.
coord=apply(coordinates sols_k,i->lift(i,QQ)) --lift to QQ
--Step 1: symmetric matrix with same suff stat
RPD=QQ[a,b]
S=matrix{{coord_0,coord_3/4,a,coord_3/4},
         {coord_3/4,coord_1/2,coord_4/2,b},
	 {a,coord_4/2,coord_2,coord_5/2},
	 {coord_3/4,b,coord_5/2,coord_1/2}}
{S_(0,0),S_(1,1)+S_(3,3),S_(2,2),2*(S_(0,1)+S_(0,3)),2*S_(1,2),2*S_(2,3)}==coord
apply({0,1,2,3},apply({0,1,2,3},))

S_(0,0)
S_(0,0)*S_(1,1)-S_(0,1)^2
det submatrix'(S,{3},{3})
I1=ideal {det submatrix'(S,{3},{3}),b}
dim I1,degree I1
solsMinor3=zeroDimSolve I1;
netList solsMinor3 
-- for -6.984 < a < 24.5852, det submatrix'(S,{3},{3})>0 
det sub(S,{a=>0})
I2=ideal{det sub(S,{a=>0}),a}
dim I2,degree I2
solsMinor4=zeroDimSolve I2;
netList solsMinor4 
-- no real solutions for a=0
Sa=sub(S,{a=>0,b=>0})
eigenvalues Sa

RMLE=QQ[l_1..l_6];
K=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}
rank Sa

trace(Sa*K)
Raux=QQ[l_1..l_6,a,b]
S
sub(S,Raux)*sub(K,Raux)
tr=trace(sub(S,Raux)*sub(K,Raux))

f=sub(tr,RMLE)

I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{f}})};
J=saturate(I,det K);
dim J,degree J -- ML-degree in table

debug loadPackage "GraphicalModelsMLE"
needsPackage "EigenSolver"
solutions=zeroDimSolve(J);
solutions	


M=genListMatrix(solutions,K)
L=checkPD(M)
(maxPt,E)=maxMLE(L,Sa)

SS=inverse E
eigenvalues SS
SS
apply({S_(0,0),S_(1,1)+S_(3,3),S_(2,2),2*(S_(0,1)+S_(0,3)),2*S_(1,2),2*S_(2,3)},i->sub(i,RR))
{SS_(0,0),SS_(1,1)+SS_(3,3),SS_(2,2),2*(SS_(0,1)+SS_(0,3)),2*SS_(1,2),2*SS_(2,3)}
apply(coord,i->sub(i,RR))
-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------
-- SAME COMPUTATIONS WITH COMPLETION OF GRAPH 9
restart
n=8
R=QQ[l_1..l_n,t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
K=matrix{{l_1,l_4,l_7,l_4},{l_4,l_2,l_5,l_8},{l_7,l_5,l_3,l_6},{l_4,l_8,l_6,l_2}}
Sigma=genericSymmetricMatrix(R,s_11,4)


-- algebraic boundary 
HG=algBoundary(K);
Istat=ideal(t_1-s_11,t_2-s_22-s_44,t_3-s_33,t_4-2*(s_12+s_14),t_5-2*s_23,t_6-2*s_34,t_7-2*s_13,t_8-2*s_24)
IG3=eliminate(support Sigma,minors(4,Sigma)+Istat)
IG2=eliminate(support Sigma,minors(3,Sigma)+Istat)
netList IG2_*
eliminate({t_7,t_8},IG2)
IG1=eliminate(support Sigma,minors(2,Sigma)+Istat)
netList IG1_*
netList (trim IG1)_*
eliminate({t_7,t_8},IG1)

n=4
k=100
for i from 1 to k do (
--    X=random(QQ^4,QQ^n);              
    X=sub(random(ZZ^4,ZZ^n),QQ);
    S=(1/n)*X*transpose(X);
    aux={t_1=>S_(0,0),t_2=>S_(1,1)+S_(3,3),t_3=>S_(2,2),t_4=>2*(S_(0,1)+S_(0,3)),t_5=>2*S_(1,2),t_6=>2*S_(2,3),t_7=>2*S_(0,2),t_8=>2*S_(1,3)};
--    print S;
--    print toString aux;
--    print apply(0..4,i->sub(HG_i,aux)); 
--    print (sub(HG_2,aux))_0;
--    print(apply(flatten entries gens IG1orig, i->sub(i,aux)))
      print(apply(flatten entries gens IG1, i->sub(i,aux)))
    )

-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------
-- Compute elimination ideals
Istat=ideal{t_1-s_11,t_2-s_22-s_44,t_3-s_33,t_4-2*(s_12+s_14),t_5-2*s_23,t_6-2*s_34};
netList Istat_*

I_G3=eliminate(support Sigma,minors(4,Sigma)+Istat)
I_G2=eliminate(support Sigma,minors(3,Sigma)+Istat)
I_G1=eliminate(support Sigma,minors(2,Sigma)+Istat)

netList (I_G1)_*
isPrime I_G1

RChol=QQ[l11,l12,l13,l14,l22,l23,l24,l33,l34,l44]

Rtotal=QQ[s_11..s_14,s_22..s_24,s_33..s_34,s_44,gens RChol,t_1..t_n]
S=sub(Sigma,Rtotal)
L=matrix{{l11,0,0,0},{l12,l22,0,0},{l13,l23,l33,0},{l14,l24,l34,l44}}
SPD=L*transpose(L)
IChol=trim ideal flatten entries (S-SPD);
eliminate(support S,minors(4,S)+IChol)

Istat=ideal{t_1-s_11,t_2-s_22-s_44,t_3-s_33,t_4-2*(s_12+s_14),t_5-2*s_23,t_6-2*s_34};
netList Istat_*
test=eliminate(support S,Istat+IChol)
netList test_*
netList (trim test)_*
minors(4,SPD)

IChol3=eliminate(support SPD,minors(4,SPD)+test)
IChol2=eliminate(support SPD,minors(3,SPD)+test)
IChol3=eliminate(support SPD,minors(2,SPD)+test)


---try out SOS
loadPackage "SemidefiniteProgramming"
loadPackage "SumsOfSquares"
sol=solveSOS (I_G1)_1
peek sol
s=sosPoly sol
peek s
coefs=s#coefficients
gene=s#generators

prod=apply(coefs,gene,(i,j)->i*j^2)
sos=sum(prod)
F_3==sos --true


-- TEST: Elimination ideal with rk constraints

--irr factor of the algebraic boundary that changes signs in the interior
f=t_3*t_4^2-t_1*t_5^2-2*t_1*t_5*t_6-t_1*t_6^2
-- t_3*t_2^2-t_1*(t_5+t_6)^2
codim ideal f
degree ideal f
isSubset(HG_3,I_G1) --true, I_G1=HG_2+HG_4
isSubset(HG_2,HG_3+HG_4) --true
isSubset(HG_4,HG_2+HG_3) --false

--Compute MLE for different sizes of samples
restart
m=4;
L=flatten for i from 1 to binomial(m+1,2) list l_i;
Raux=QQ[L];

K=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}
R=QQ[support K];
K=sub(K,R);

--n=1;
--X=random(QQ^m,QQ^n);              
--S=(1/m)*X*transpose(X);
--S=matrix{{251/4,131/4,20,23},{131/4,101/4,15,13},{20,15,41/4,4},{23,13,4,20}}
--S=matrix{{27/2,22,22,43/4},{22,175/4,37,99/4},{22,37,73/2,75/4},{43/4,99/4,75/4,61/4}}
S=matrix{{239859/5000, 20, 5, 57697/2000},{ 20, 20, 432029/20000, 3},{ 5, 432029/20000, 181267/10000, 168519/20000},{ 57697/2000, 3, 168519/20000, 293759/10000}}

rank S
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
J=saturate(I,det K);
degree J -- ML-degree in table

debug loadPackage "GraphicalModelsMLE"
needsPackage "EigenSolver"
sols=zeroDimSolve(J);
sols	
M=genListMatrix(sols,K)
L=checkPD(M)
(maxPt,E)=maxMLE(L,S)

SS=inverse E
apply({S_(0,0),S_(1,1)+S_(3,3),S_(2,2),2*(S_(0,1)+S_(0,3)),2*S_(1,2),2*S_(2,3)},i->sub(i,RR))
{SS_(0,0),SS_(1,1)+SS_(3,3),SS_(2,2),2*(SS_(0,1)+SS_(0,3)),2*SS_(1,2),2*SS_(2,3)}


--Can we understand whether it has a solution or not based on the values of the 4
-- factors of the algebraic boundary that are not part of the elimination ideal
-- in rank 1?
t1=S_(0,0)+S_(1,1);
t2=S_(2,2)+S_(3,3);
t3=2*S_(0,1);
t4=2*(S_(1,2)+S_(0,3));
t5=2*S_(2,3);

aux={t1+t3,t1-t3,t2+t5,t2-t5}


n=1
count={};
for i from 1 to 10000 do(
X=random(QQ^m,QQ^n);              
S=(1/m)*X*transpose(X);
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
J=saturate(I,det K);
sols=zeroDimSolve(J);
M=genListMatrix(sols,K);
L=checkPD(M);
count=append(count,{degree J,#L});
    )

summary=tally count
--Tally{{1, 0} => 10100}

------------------------------------------------------------------------------------
--How do suff stat look like when they cannot be completed into a full rank matrix?
------------------------------------------------------------------------------------
restart
m=4;
L=flatten for i from 1 to binomial(m+1,2) list l_i;
Raux=QQ[L];

K=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}
R=QQ[support K];
K=sub(K,R);

n=1;
--X=random(QQ^m,QQ^n);              
--S=(1/m)*X*transpose(X);

S=matrix {{9/4, 3/16, 3/4, 3/8}, {3/16, 1/64, 1/16, 1/32}, {3/4, 1/16, 1/4, 1/8}, {3/8, 1/32, 1/8, 1/16}}
rank S
eigenvalues S

{S_(0,0),S_(1,1)+S_(3,3),S_(2,2),2*(S_(0,1)+S_(0,3)),2*S_(1,2),2*S_(2,3)}
M=matrix{{9/4,1/4,0,5/16},
        {1/4,2/64,1/16,0},
	{0,1/16,1/4,1/8},
	{5/16,0,1/8,3/64}}
{M_(0,0),M_(1,1)+M_(3,3),M_(2,2),2*(M_(0,1)+M_(0,3)),2*M_(1,2),2*M_(2,3)}
rank M
eigenvalues M

M=matrix{{9/4,1/4,3,5/16},
        {1/4,2/64,1/16,11},
	{3,1/16,1/4,1/8},
	{5/16,11,1/8,3/64}}
{M_(0,0),M_(1,1)+M_(3,3),M_(2,2),2*(M_(0,1)+M_(0,3)),2*M_(1,2),2*M_(2,3)}
rank M
eigenvalues M

S=M

I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
J=saturate(I,det K);
degree J -- ML-degree in table

debug loadPackage "GraphicalModelsMLE"
needsPackage "EigenSolver"
sols=zeroDimSolve(J);
sols	
M=genListMatrix(sols,K)
L=checkPD(M)
(maxPt,E)=maxMLE(L,S)


---------------------RANK 1 EXTENSION

R=QQ[s_11..s_14,s_22..s_24,s_33..s_34,s_44]
Sigma=genericSymmetricMatrix(R,s_11,4)
IdMin=minors(2,Sigma)
dim IdMin, degree IdMin
coord={1214817414/27327025, 807645352/10397027, 469043467/13711069, 918171032/6050823, 482712667/7589086, 1074513611/19280609}
Istat=ideal{1214817414/27327025-s_11,807645352/10397027-s_22-s_44,469043467/13711069-s_33,
    918171032/6050823-2*(s_12+s_14),482712667/7589086-2*s_23,1074513611/19280609-2*s_34};
dim Istat, degree Istat
I=(trim IdMin)+Istat
dim I,degree Istat
netList I_*

ideal(s_34^2-s_33*s_44,s_24*s_34-s_23*s_44,s_14*s_34-s_13*s_44,s_24*s_33-s_23*s_34,s_14*s_33-s_13*s_34,s_24^2-s_22*s_44,s_23*s_24-s_22*s_34,s_14*s_24-s_12*s_44,s_13*s_24-s_12*s_34,s_
      23^2-s_22*s_33,s_14*s_23-s_12*s_34,s_13*s_23-s_12*s_33,s_14*s_22-s_12*s_24,s_13*s_22-s_12*s_23,s_14^2-s_11*s_44,s_13*s_14-s_11*s_34,s_12*s_14-s_11*s_24,s_13^2-s_11*s_33,s_12*s_13-s_11
      *s_23,s_12^2-s_11*s_22,-s_11+1214817414/27327025,-s_22-s_44+807645352/10397027,
      -s_33+469043467/13711069,-2*s_12-2*s_14+918171032/6050823,-2*s_23+482712667/7589086,-2*s_34+1074513611/
      19280609)

aux=ideal apply(flatten entries gens I, i->sub(i,{s_11=>1214817414/27327025,
	    s_33=>469043467/13711069,s_23=>482712667/(2*7589086),
	    s_34=>1074513611/(2*19280609)}))
netList aux
toString aux_0
-(469043467/13711069)*s_44+1154579500224259321/1486967533643524
toString aux_3
(469043467/13711069)*s_24-518681330893610537/585288799333496

aux2=ideal apply(flatten entries gens aux, i->sub(i,{s_44=>1154579500224259321*13711069/(469043467*1486967533643524),
s_24=>518681330893610537*13711069/(585288799333496*469043467)}))

aux3=ideal apply(flatten entries gens aux, i->sub(i,{s_22=>398703859514842939224267097468193/7251431509877722661419244634116}))
sub(aux3_9,RR)
