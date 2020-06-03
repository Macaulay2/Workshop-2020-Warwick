-- study of bi-homogenization
QQ[x,y]
F = {x^2+y^2-1, x^3-y^2} -- defines 6 points in A^2
dim ideal F
degree ideal F
degree radical ideal F
QQ[x]
roots (x^2+x^3-1)

R = QQ[x0,x1,y0,y1,Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
hF = {x1^2*y0^2+y1^2*x0^2-x0^2*y0^2, y0^2*x1^3-x0^3*y1^2} -- defines ... A^1 x A^1
-*
ass radical ideal hF
primaryDecomposition ideal hF
*-

-*
QQ[x0,x1,y0,y1,Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
F = ideal {x1^2+y1^2-1, x1^3-y1^2}
homogenize(F,x0) -- can't homogenize if multi-graded!!!
*-

--------------------------------------
-- Solving points in a P1xP1
--------------------------------------
restart
loadPackage "NumericalAlgebraicGeometry"
R = QQ[x0,x1,y0,y1,Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
-- this example below requires generalized eigenvalues
-- hF = {x1^2*y0^2+y1^2*x0^2-x0^2*y0^2, y0^2*x1^3-x0^3*y1^2} -- defines ... A^1 x A^1
-- J = ideal hF
J=ideal(random({2,30},R),random({2,30},R));

--solve with NAG for comparison
S=CC[x1,y1]
rho=map(S,R,{1,x1,1,y1})
JS=rho(J)
NAGsol=solveSystem {JS_0,JS_1}
netList NAGsol
use R

-- start the method based on eigencomputations
-- matrix reverse for i to 7 list for j to 7 list hilbertFunction({j,i},J)
-- satind={4,3} 
-- Here is an automatic estimation assuming a complete intersection : 2 equations
deg=degrees J; satind={(deg_0)_0+(deg_1)_0-1,(deg_0)_1+(deg_1)_1-1}
psi=map(R^{{0,0}}, R^(-(J_*/degree)), gens J);
elimMat=matrix basis(satind,psi);
-- mb=transpose basis(satind,R) --basis used to index the rows of elimMat
-- K = syz transpose elimMat -- "exact" basis of coker of linear map
-- For a better stability, compute the numerical kernel
elimMatTr=transpose sub(elimMat,RR);
nc=rank source elimMatTr; numrk=numericalRank elimMatTr
(S,U,Vt)=SVD(elimMatTr); K=transpose Vt^{numrk..nc-1}; -- numerical kernel
numroots=rank source K -- expected number of roots

D0=K^{0..numroots-1}; -- indexed by the basis 1,y_1,..
--sub(mb^{0..numroots-1},{x0=>1,y0=>1})
D1=K^{satind_1+1..satind_1+numroots}; -- indexed by the basis x_1*1,x_1*y_1,... 
-- sub(mb^{satind_1+1..satind_1+numroots},{x0=>1,y0=>1})
-- SVD D0
(EVal,EVec)=eigenvectors (inverse(D0)*D1); EVect=D0*EVec;
listSol=apply(#EVal, i -> {EVal_i,EVect_(1,i)/EVect_(0,i)}) -- roots (x_1,y_1) assuming x0=y0=1

-- Comparison
netList listSol
netList NAGsol
