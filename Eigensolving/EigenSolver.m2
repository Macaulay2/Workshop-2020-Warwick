newPackage(
    "EigenSolver",
    Version => "0.1", 
    Date => "June 2020",
    Authors => {{Name => "", 
	    Email => "", 
	    HomePage => ""}},
    Headline => "polynomial system solver relying on eigen-computations",
    PackageExports => {"NumericalAlgebraicGeometry"},
    AuxiliaryFiles => false,
    DebuggingMode => false
    )
    
export {
    "zeroDimSolve",
    "Basis",
    "Multiplier"
}

zeroDimSolve = method(Options => {
    symbol Basis => null,
    symbol Multiplier => 0,
    Strategy => "Stickelberger"})
zeroDimSolve Ideal := List => opts -> I -> ( -- Assume I is radical
    if opts.Strategy == "Stickelberger" then return eigSolve1(I, opts)
    else if opts.Strategy == "elimMatrixP1P1" then return eigSolveP1P1(I, opts)
    else if opts.Strategy == "elimMatrixP2" then return eigSolveP2(I, opts)
    else error "Unknown strategy";
)

reduceCoeffs = (A, I, m) -> transpose sub(last coefficients(A % I, Monomials => m), coefficientRing ring I)

eigSolve1 = method(Options => options zeroDimSolve) -- based on Tomas' function SolvePolyEqByEigV
eigSolve1 Ideal := List => opts -> I -> (
    R := ring I;
    B := if opts.Basis =!= null then opts.Basis else sub(basis(R/I), R);
    f := if opts.Multiplier =!= 0 then opts.Multiplier else random(1, R);
    M := reduceCoeffs(f*B, I, B);
    (E, V) := eigenvectors M;
    V = V*inverse diagonalMatrix V^{0};
    (entries transpose(reduceCoeffs(vars R, I, B)*V))/(p -> point{p})
)

eigSolveP1P1 = method(Options => options zeroDimSolve) -- based on Laurent's code
eigSolveP1P1 Ideal := List => opts -> J -> ( -- currently assumes dim R = 2 et 2 equations, homogenizes to P^1 x P^1
    R := ring J;
    FF := coefficientRing R;
    Lvar := vars R;
    (xx0, yy0) := (symbol xx0, symbol yy0);
    nR := FF[xx0,Lvar_(0,0),yy0,Lvar_(0,1)];
    E0 := homogenize(homogenize(sub(J_0,nR),xx0,{1,1,0,0}),yy0,{0,0,1,1});
    E1 := homogenize(homogenize(sub(J_1,nR),xx0,{1,1,0,0}),yy0,{0,0,1,1});
    JnR := ideal(E0,E1); -- homogenized equations
    nnR := FF[xx0,Lvar_(0,0),yy0,Lvar_(0,1),Degrees=>{{1,0},{1,0},{0,1},{0,1}}];
    JnnR := sub(JnR,nnR); -- homogeneized equations in the bigraded ring
    deg := degrees JnnR;
    satind := {(deg_0)_0+(deg_1)_0-1,(deg_0)_1+(deg_1)_1-1};
    psi := map(nnR^{{0,0}}, nnR^(-(JnnR_*/degree)), gens JnnR);
    elimMat := matrix basis(satind,psi);
    elimMatTr := transpose sub(elimMat,CC);
    numrk := numericalRank elimMatTr;
    K := numericalKernel(elimMatTr,getDefault(Tolerance));
    numroots := rank source K; -- expected number of roots
    D0 := K^{0..numroots-1}; -- indexed by the basis 1,y_1,..
    D1 := K^{satind_1+1..satind_1+numroots}; -- indexed by the basis x_1*1,x_1*y_1,... 
    (EVal,EVec) := eigenvectors (inverse(D0)*D1); -- NEED GENERALIZED EIGENVALUES
    EVect := D0*EVec;
    apply(#EVal, i -> point{{EVal_i,EVect_(1,i)/EVect_(0,i)}}) -- roots (x_1,y_1) assuming x0=y0=1
)

eigSolveP2 = method(Options => options zeroDimSolve) -- based on Laurent's code
eigSolveP2 Ideal := List => opts -> J -> ( -- currently assumes dim R = 2 et 2 (nonlinear) equations, homogenizes to P^2
    R := ring J;
    FF := coefficientRing R;
    Lvar := vars R;
    xx0 := symbol xx0;
    nR := FF[xx0,gens R];
    E0 := homogenize(sub(J_0,nR),xx0);
    E1 := homogenize(sub(J_1,nR),xx0);
    JnR := ideal(E0,E1); -- homogenized equations
    deg := degrees JnR;
    satind := (deg_0)_0+(deg_1)_0-1;
    psi := map(nR^1, nR^(-(JnR_*/degree)), gens JnR);
    B := flatten entries basis(satind, nR);
    H := hashTable apply(#B, i -> B#i => i);
    elimMat := matrix basis(satind,psi);
    elimMatTr := transpose sub(elimMat,CC);
    K := numericalKernel(elimMatTr,getDefault(Tolerance));
    numroots := rank source K; -- expected number of roots
    B0 := select(numroots, B, b -> member(nR_1*b // nR_0, B));
    D0 := K^(apply(B0, b -> H#b));
    B1 := apply(B0, b -> (b*nR_1)//nR_0);
    D1 := K^(apply(B1, b -> H#b));
    (EVal,EVec) := eigenvectors (inverse(D0)*D1); -- NEED GENERALIZED EIGENVALUES (provisionally assume no solutions at infinity, separated by x_1)
    EVect := D0*EVec;
    apply(#EVal, i -> point{{EVal_i,EVect_(2,i)/EVect_(0,i)}}) -- roots (x_1,y_1) assuming x0=y0=1    
)

-- needs "laurent-eigensolving.m2"
-- needs "documentation.m2"


TEST ///
-- Problem 2.11 in https://math.berkeley.edu/~bernd/cbms.pdf
n = 3
R = QQ[a_1..a_(n-1),b_1..b_(n-1)]
S = R[x,y]/ideal(x*y-1)
f = sum(n-1, i -> R_i*x^(i+1) + R_(n-1+i)*y^(i+1)) + x^n + y^n
I = sub(ideal apply(toList(2..2*n-1), i -> last coefficients(f^i, Monomials => {1_S})), R)
sols = zeroDimSolve I
realPts = realPoints sols
assert(#sols == 66 and #realPts == 6)
sort apply(sols, p -> norm evaluate(gens I, p))
///

TEST ///
needsPackage "ExampleSystems"
cyc6 = ideal cyclic(6, QQ)
elapsedTime sols = zeroDimSolve cyc6;
sort apply(sols, p -> norm evaluate(gens cyc6, p))
///

TEST ///
A = QQ[x0,x1,y0,y1,Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
I = ideal(random({2,5},A),random({2,5},A))
R = QQ[x1,y1]
spe = map(R,A,matrix{{1,x1,1,y1}})
J=spe I
listSol = zeroDimSolve(J, Strategy => "elimMatrixP1P1")
assert(#listSol == 20)

B = QQ[x0,x1,x2]
I = ideal(random(2,B), random(2,B))
R = QQ[x1,x2]
spe = map(R,B,matrix{{1,x1,x2}})
J = spe I
elapsedTime listSol = zeroDimSolve(J, Strategy => "elimMatrixP2")
sort apply(listSol, p -> norm evaluate(gens J, p))
assert(#listSol == product(I_*/degree/first))
///

TEST /// -- Example from Roser and Olga
needsPackage "GraphicalModels"
G=graph{{1,2},{2,3},{3,4},{1,4}}
R=gaussianRing(G)
K=undirectedEdgesMatrix R
d=4
numS=lift(d*(d+1)/2,ZZ);
lpRvar=apply(numgens(R)-numS,i->(gens(R))_i);
R2=coefficientRing(R)[lpRvar]
R2map=apply(numgens(R),i-> if i<= numgens(R)-numS-1 then (gens(R2))_i else 0)
F=map(R2,R,R2map)
K2=F(K)
X=random(ZZ^4,ZZ^4)
S=X*transpose(X)
I=ideal{jacobian ideal{determinant(K2)}-determinant(K2)*jacobian(ideal{trace(K2*S)})}
J=saturate(I,ideal{determinant(K2)})
sols = zeroDimSolve J
assert(#sols == 5)
///

end--

restart
needsPackage "EigenSolver"
loadPackage("EigenSolver", Reload => true)
installPackage("EigenSolver", RemakeAllDocumentation => true)
check "EigenSolver"
