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
-- Solving two points in a P1xP1
--------------------------------------
J = ideal hF
satind={4,4}
basis(satind,R)
basis(satind,R^(-(J_*/degree)))
psi=map(R^{{0,0}}, R^(-(J_*/degree)), gens J)
elimMat=matrix basis(satind,psi)
K = syz transpose elimMat -- basis of coker of linear map
hilbertFunction(satind,R/J)

!!!!!!!!!!!!!1

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

