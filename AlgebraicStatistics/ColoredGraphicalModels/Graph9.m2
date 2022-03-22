-- Graph 9 Uhler
-- auxiliary functions (Olga's functions)
dualVariety=(IX,n,x,u)->(
    c:=codim IX;
    JacX:=submatrix(transpose jacobian IX,toList(0..n-1));
    AugJacX:=matrix{toList(u_1..u_n)}||JacX;
    SingX:=minors(c,JacX);
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
Sigma=genericSymmetricMatrix(R,s_11,4)

-- algebraic boundary 
HG=algBoundary(K);
netList HG
isPrime(HG_2)
toString (HG_2)_0

--Chose your favourite rank n and the number of elements k you want to try
--This procedure prints out the the values of the sufficient statistics for k random matrices of rank n
n=4
k=3
for i from 1 to k do (
--    X=random(QQ^4,QQ^n);              
    X=sub(random(ZZ^4,ZZ^n),QQ);
    S=(1/n)*X*transpose(X);
    aux={t_1=>S_(0,0),t_2=>S_(1,1)+S_(3,3),t_3=>S_(2,2),t_4=>2*(S_(0,1)+S_(0,3)),t_5=>2*S_(1,2),t_6=>2*S_(2,3)};
    print S;
    print toString aux;
--    print apply(0..4,i->sub(HG_i,aux)); 
    print (sub(HG_2,aux))_0;
    )

--Below we have two sample covariance matrices of rank 4 
-- whose sufficient statistics take different signs when 
-- evaluated at f3

| 251/4 131/4 20   23 |
| 131/4 101/4 15   13 |
| 20    15    41/4 4  |
| 23    13    4    20 |

P1={t_1 => 251/4, t_2 => 181/4, t_3 => 41/4, t_4 => 223/2, t_5 => 30, t_6 => 8}
for i to 5 do print (sub(HG_i,P1))_0

| 27/2 22    22   43/4 |
| 22   175/4 37   99/4 |
| 22   37    73/2 75/4 |
| 43/4 99/4  75/4 61/4 |

P2={t_1 => 27/2, t_2 => 59, t_3 => 73/2, t_4 => 131/2, t_5 => 74, t_6 => 75/2}
for i to 5 do print (sub(HG_i,P2))_0

--Define line through 2 points (with auxiliary variable)
restart
Rt=QQ[t_1..t_6,m]
v=sub(vector delete(m,gens Rt),Rt)
P=sub(vector{251/4,181/4,41/4,223/2,30,8},Rt)
Q=sub(vector{27/2,59,73/2,131/2,74,75/2},Rt)
paramr=ideal flatten entries (v-P-m*(Q-P))
r=eliminate(m,paramr)

--Line through 2 points (in ring with variables t1...t6)
Rtt=QQ[t_1..t_6]
r=sub(r,Rtt)
betti trim r
dim r
sub(r,{t_1 => 251/4, t_2 => 181/4, t_3 => 41/4, t_4 => 223/2, t_5 => 30, t_6 => 8})
sub(r,{t_1 => 27/2, t_2 => 59, t_3 => 73/2, t_4 => 131/2, t_5 => 74, t_6 => 75/2})


-- Intersection of the line and the factor of the algebraic boundary (f3)
f=t_3*t_4^2-t_1*t_5^2-2*t_1*t_5*t_6-t_1*t_6^2
I=ideal(f,r)
betti trim I
dim I
degree I

loadPackage "EigenSolver"
sols=zeroDimSolve I;
netList sols  -- 3 real solutions

-- How close to vanishing at the line and the cubic are actually these solutions?
aux={}
k=2
for i from 1 to 6 do aux=append(aux,t_i=>sub((coordinates sols_k)_(i-1),CC));
netList aux
sub(r,aux)
sub(f,aux)

-- How can we check whether these points are in the interior 
-- of the cone of algebraic statistics?
-- Option 1: check if one of them is a convex combination of the other two
Pr=apply(flatten entries P,i->sub(i,RR))
Qr=apply(flatten entries Q,i->sub(i,RR))
netList sols
--sols_2 could be a convex combination of P and Q because in all coordinates 
-- the value of sols_2 is in between the corresponding value of P and Q

k=2
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
aux=sub(conv,m=>2959052059/9861344245)
apply(0..5,i->sub(aux_i,RR))
-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------
-- Compute elimination ideals
Istat=ideal{t_1-s_11,t_2-s_22-s_44,t_3-s_33,t_4-2*(s_12+s_14),t_5-2*s_23,t_6-2*s_34};
netList Istat_*

I_G3=eliminate(support Sigma,minors(4,Sigma)+Istat)
I_G2=eliminate(support Sigma,minors(3,Sigma)+Istat)
I_G1=eliminate(support Sigma,minors(2,Sigma)+Istat)
netList (I_G1)_*

--irr factor of the algebraic boundary that changes signs in the interior
f=t_3*t_4^2-t_1*t_5^2-2*t_1*t_5*t_6-t_1*t_6^2
-- t_3*t_2^2-t_1*(t_5+t_6)^2
codim ideal f
degree ideal f
isSubset(HG_3,I_G1) --true, I_G1=HG_2+HG_4
isSubset(HG_2,HG_3+HG_4) --true
isSubset(HG_4,HG_2+HG_3) --false

loadPackage "Divisor"
isSmooth HG_2 --false
isSmooth (HG_2,IsGraded=>true) --false
-- Not smooth, also without the origin

isSmooth ideal{t_1^2-t_2^3} --false
singularLocus ideal{t_1^2-t_2^3} 

isSmooth ideal{t_1^2+t_2^2-1} --true
singularLocus ideal{t_1^2+t_2^2-1} --don't get the answer: shouldn't the singular locus be empty if it's smooth?

isSmooth ideal{t_1^2+t_2^2-t_1^3} --false
isSmooth ideal{t_1+t_2^2-t_1^3} --true

--Cartan's umbrella
isSmooth ideal{t_3*(t_1^2+t_2^2)-t_1^3}
isSmooth (ideal{t_3*(t_1^2+t_2^2)-t_1^3},IsGraded=>true) --not smooth, even without zero

--Compute MLE for different sizes of samples
restart
m=4;
L=flatten for i from 1 to binomial(m+1,2) list l_i;
Raux=QQ[L];

K=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}
R=QQ[support K];
K=sub(K,R);

n=1;
X=random(QQ^m,QQ^n);              
S=(1/m)*X*transpose(X);
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
