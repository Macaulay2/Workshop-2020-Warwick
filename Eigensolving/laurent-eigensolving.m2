--------------------------------------
-- solve a univariate polynomial
--------------------------------------
restart
R=RR[x_0,x_1]

-- input polynomial; x_0=0 is the point at infinity
f=(x_1-2*x_0)*(x_1-x_0)^2
f=1; for i from 1 to 11 do f=f*(x_1-i*x_0); 
f=random(50,R)
f=(x_1-x_0)^10
f=1; for i from 0 to 10 do f=f*(x_1-i*x_0); -- need generalized eignevalues

-- solving
d=(degree f)_0; mb=basis(d,R)
(C,M)=coefficients(f,Monomials=>mb); -- M is the coefficient matrix
K=syz transpose M;
D0=K^{0..(d-1)}; -- the companion matrix 
D1=K^{1..d}; -- the identity
D0=sub(D0,RR);D1=sub(D1,RR); rank D0 == d -- sub in RR
listSol=eigenvalues (inverse(D0)*D1) -- not optimal numerically
-- (SVD(D0))_0


-- Verification
for j from 0 to #listSol - 1 do print sub(f,{x_0=>1,x_1=>listSol_j})

-- Remark: K could have been computed as an approximate kernel of a linear system

--------------------------------------
-- Solving two points in a P1xP1
--------------------------------------
restart
R=QQ[x_0,x_1,y_0,y_1,Degrees=>{{1,0},{1,0},{0,1},{0,1}}]

-- set up an example
P=ideal((x_0-2*x_1),(y_0-y_1)) -- point (1:1/2)x(1:1)
Q=ideal((3*x_0-x_1),(5*y_0-y_1)) -- point (1:3)x(1:5)
I=intersect(P,Q) 
JJ=super basis({2,2},I)
J=JJ*(random(QQ^(numcols JJ),QQ^3))
J= ideal J
transpose gens J
B=ideal(x_0,x_1)*ideal(y_0,y_1)
JB=saturate(J,B) 
primaryDecomposition JB -- gives back the input points


-- Elimination matrix
nu={3,4}; hilbertFunction(nu,J)-hilbertFunction(nu,JB)  -- (3,4) or (4,3) is the saturation index
satind={3,4}
psi=(res J).dd_1
elimMat=matrix basis(satind,psi)
-- transpose basis(satind,R)
K=syz transpose elimMat

-- solve in the yi's
D0=K^{0,1} -- indexed by the basis 1,y_1
D1=K^{1,2} -- indexed by the basis y_1,y_1^2 
eigenvalues (sub(inverse(D0)*D1,RR)) -- ,y_0=1, y_1 values

-- solve in the xi's
D0=K^{0,5} -- indexed by the basis 1,x_1
D1=K^{5,10} -- indexed by the basis x_1,x_1^2
eigenvalues (sub(inverse(D0)*D1,RR)) -- ,x_0=1, x_1 values

-- solve both at the same time to avoid a matching step
D0=K^{0,1} -- indexed by the basis 1,y_1
D1=K^{5,6} -- indexed by the basis x_1,x_1*y_1
(EVal,EVec)=eigenvectors (sub(inverse(D0)*D1,RR))
{EVal_0 , EVec_(1,0)/EVec_(0,0) }
{EVal_1 , EVec_(1,1)/EVec_(0,1) }


------------------------------------------------
-- Solve the intersection of two curves in a P2
------------------------------------------------

restart
R=QQ[x_0,x_1,x_2]
d1=2; f1=random(d1,R)
d2=2; f2=random(d2,R)
I=ideal(f1,f2)
psi=(res I).dd_1
nu=d1+d2-2 -- HF(R/I;nu) coincide with HP(R/I;nu)
elimMat=matrix basis(nu,psi); -- Matrix
basis(nu,R)
numrows(elimMat) - rank(elimMat) -- number of roots
-- need to be completed