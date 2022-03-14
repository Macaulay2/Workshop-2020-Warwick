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

n=6
R=QQ[l_1..l_n,t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
K=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}
Sigma=genericSymmetricMatrix(R,s_11,4)

-- algebraic boundary 
HG=algBoundary(K);
netList HG

aux={t_1=>S_(0,0),t_2=>S_(1,1)+S_(3,3),t_3=>S_(2,2),t_4=>2*(S_(0,1)+S_(0,3)),t_5=>2*S_(1,2),t_6=>2*S_(2,3)}

sub((HG)_0,aux)

sub((HG)_1,aux)

sub((HG)_2,aux)

sub((HG)_3,aux)

sub((HG)_4,aux)

--Chose your favourite rank n and the number of elements k you want to try
--This procedure prints out the the values of the sufficient statistics for k random matrices of rank n
n=2
k=10
for i from 1 to k do (
    X=random(QQ^4,QQ^n);              
    S=(1/4)*X*transpose(X);
    aux={t_1=>S_(0,0),t_2=>S_(1,1)+S_(3,3),t_3=>S_(2,2),t_4=>2*(S_(0,1)+S_(0,3)),t_5=>2*S_(1,2),t_6=>2*S_(2,3)};
    print apply(0..4,i->sub(HG_i,aux)); 
    )


-- Compute elimination ideals
Istat=ideal(t_1-s_11,t_2-s_22-s_44,t_3-s_33,t_4-2*(s_12+s_14),t_5-2*s_23,t_6-2*s_34);
netList Istat_*

I_G3=eliminate(support Sigma,minors(4,Sigma)+Istat)
I_G2=eliminate(support Sigma,minors(3,Sigma)+Istat)
I_G1=eliminate(support Sigma,minors(2,Sigma)+Istat)
netList (I_G1)_*


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
