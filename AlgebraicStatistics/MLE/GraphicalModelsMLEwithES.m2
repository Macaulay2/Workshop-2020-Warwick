restart;

load "GraphicalModelsMLE.m2"
load "EigenSolver.m2"
--needsPackage("EigenSolver",FileName => "~/Workshop-2020-Warwick/Eigensolving/EigenSolver.m2")

Example 1 (from tests)
G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
R=gaussianRing(G)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
J=scoreEquationsFromCovarianceMatrix(R,U);
dim J,degree J


--With brand new eigensolver
s=zeroDimSolve J


-- NumericalAlgebraicGeometry version 1.11
needsPackage "NumericalAlgebraicGeometry"
s= solveSystem J_*


--Alternatives
needsPackage "MonodromySolver"
s2=sparseMonodromySolve(polySystem J_*)
assert(s==s2)

--Comparison 
elapsedTime solveSystem J_*
elapsedTime sparseMonodromySolve(polySystem J_*)


--Example 2.1.13 OWL: undirected graph
restart
needsPackage "GraphicalModelsMLE"
needsPackage "EigenSolver"
--needsPackage("EigenSolver",FileName => "~/Workshop-2020-Warwick/Eigensolving/EigenSolver.m2")
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
--X={{1,5,5,0},{7,2,7,4},{9,4,0,5},{7,4,7,8}}
S=X*transpose(X)
I=ideal{jacobian ideal{determinant(K2)}-determinant(K2)*jacobian(ideal{trace(K2*S)})}
J=saturate(I,ideal{determinant(K2)})
dim J, degree J

sols=zeroDimSolve J
#sols
realPoints sols


--options solveSystem
needsPackage "NumericalAlgebraicGeometry"
s=solveSystem J_*
solveSystem(J_*,Precision=>infinity)
sparseMonodromySolve(polySystem J_*)

J2 = trim J
solveSystem(J_*,Software=>BERTINI)

--| 1 5 5 0 |
--| 7 2 7 4 |
--| 9 4 0 5 |
--| 7 4 7 8 |





PDcheck = L -> 
( 
    mat = {};
    for l in L do
    (
    	flag = 0;
    	L1 = eigenvalues l;
    	for t in L1 do 
    	(	 
	    if 0 >= t then flag = 1;
     	);
        if flag == 0 then mat = mat | {l} ;
    );
    if mat == {} then print("none of the matrices are pd");
    return mat;
    
    -- input - list of matrices 
    -- output - list of positive definite matrices
);

