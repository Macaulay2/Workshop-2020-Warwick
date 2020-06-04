restart;

load "GraphicalModelsMLE.m2"

R=QQ[x,y];
FR=frac R;
F=1/(x^2+y^2);
M=entries JacobianMatrixOfRationalFunction(F)
N=transpose {{-2*x/(x^2 + y^2)^2,-2*y/(x^2 + y^2)^2 }}
assert(M === N)

needsPackage("Graphs");
needsPackage("GraphicalModels");
G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
R=gaussianRing(G)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
J=scoreEquationsFromCovarianceMatrix(R,U);
I=ideal(20*p_(3,4)+39,50*p_(4,4)-271,440104*p_(3,3)-742363,230*p_(2,2)-203,16*p_(1,1)-115,5*l_(3,4)+2,110026*l_(2,3)-2575,55013*l_(1,3)-600,115*l_(1,2)+26);
assert(J===I)

-- NumericalAlgebraicGeometry version 1.11
needsPackage "NumericalAlgebraicGeometry"
s= solveSystem J_*


--Alternative
needsPackage "MonodromySolver"
s2=sparseMonodromySolve(polySystem J_*)
assert(s==s2)

--Comparison 
elapsedTime solveSystem J_*
elapsedTime sparseMonodromySolve(polySystem J_*)


matRtoR2 = (M,F) -> (
    E:=entries(M);    
    return matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)))    
);



--Example 2.1.13 OWL
restart
needsPackage "GraphicalModelsMLE"
needsPackage "GraphicalModels"
needsPackage "NumericalAlgebraicGeometry"
G=graph{{1,2},{2,3},{3,4},{1,4}}
R=gaussianRing(G)
K=undirectedEdgesMatrix R
d=4
numS=lift(d*(d+1)/2,ZZ);
lpRvar=apply(numgens(R)-numS,i->(gens(R))_i);
R2=coefficientRing(R)[lpRvar]
-- R2 = QQ[a,b,c,e,f,g,h,l]

R2map=apply(numgens(R),i-> if i<= numgens(R)-numS-1 then (gens(R2))_i else 0)
F=map(R2,R,R2map)
K2=F(K)
X=random(ZZ^4,ZZ^4)
S=X*transpose(X)
I=ideal{jacobian ideal{determinant(K2)}-determinant(K2)*jacobian(ideal{trace(K2*S)})}
J=saturate(I,ideal{determinant(K2)})
dim J, degree J
options solveSystem
s=solveSystem J_*
solveSystem(J_*,Precision=>infinity)
sparseMonodromySolve(polySystem J_*)

J2 = trim J
solveSystem(J_*,Software=>BERTINI)

--| 1 5 5 0 |
--| 7 2 7 4 |
--| 9 4 0 5 |
--| 7 4 7 8 |

U={random(ZZ^1,ZZ^4),random(ZZ^1,ZZ^4),random(ZZ^1,ZZ^4),random(ZZ^1,ZZ^4)}
scoreEquationsFromCovarianceMatrix(R,U)
-- expects R to be created from a MixedGraph!!!


