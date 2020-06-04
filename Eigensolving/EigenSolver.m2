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
zeroDimSolve Ideal := List => opts -> I -> (
    if opts.Strategy == "Stickelberger" then return eigSolve1(I, opts)
    else if opts.Strategy == "syzygy" then return eigSolve2(I, opts)
    else error "Unknown strategy";
)

reduceCoeffs = (A, I, m) -> transpose sub(last coefficients(A % I, Monomials => m), coefficientRing ring I)

eigSolve1 = method(Options => options zeroDimSolve) -- based on Tomas' function SolvePolyEqByEigV
eigSolve1 Ideal := List => opts -> I -> (
    R := ring I;
    B := if opts.Basis != null then opts.Basis else sub(basis(R/I), R);
    f := if opts.Multiplier != 0 then opts.Multiplier else random(1, R);
    M := reduceCoeffs(f*B, I, B);
    (E, V) := eigenvectors M;
    V = V*inverse diagonalMatrix V^{0};
    entries transpose(reduceCoeffs(vars R, I, B)*V)
)

eigSolve2 = method(Options => options zeroDimSolve) -- based on Laurent's code
<<<<<<< HEAD
eigSolve2 Ideal := List => opts -> J -> ( -- currently assumes R = P^1 x P^1, J bihomogeneous complete intersection
=======
eigSolve2 Ideal := List => opts -> J -> ( -- currently assumes dim R = 2 et 2 equations
>>>>>>> fa905878fd83d7e0722be12c64d89f672a25934d
    R := ring J;
    FF := coefficientRing R;
    Lvar := vars R; 
    nR=FF[xx0,Lvar_(0,0),yy0,Lvar_(0,1)];
    E0:=homogenize(homogenize(sub(J_0,nR),xx0,{1,1,0,0}),yy0,{0,0,1,1});
    E1:=homogenize(homogenize(sub(J_1,nR),xx0,{1,1,0,0}),yy0,{0,0,1,1});
    JnR := ideal(E0,E1); -- homogenized equations
    nnR := FF[xx0,Lvar_(0,0),yy0,Lvar_(0,1),Degrees=>{{1,0},{1,0},{0,1},{0,1}}];
    JnnR := sub(JnR,nnR); -- homogeneized equations in the bigraded ring
    deg := degrees JnnR;
    satind := {(deg_0)_0+(deg_1)_0-1,(deg_0)_1+(deg_1)_1-1};
<<<<<<< HEAD
    psi := map(R^{{0,0}}, R^(-(J_*/degree)), gens J);
    elimMat := transpose sub(matrix basis(satind,psi),CC);
    K := numericalKernel(elimMat,getDefault(Tolerance));
=======
    psi := map(nnR^{{0,0}}, nnR^(-(JnnR_*/degree)), gens JnnR);
    elimMat := matrix basis(satind,psi);
    elimMatTr := transpose sub(elimMat,CC);
    numrk := numericalRank elimMatTr;
    K := numericalKernel(elimMatTr,getDefault(Tolerance));
>>>>>>> fa905878fd83d7e0722be12c64d89f672a25934d
    numroots := rank source K; -- expected number of roots
    D0 := K^{0..numroots-1}; -- indexed by the basis 1,y_1,..
    D1 := K^{satind_1+1..satind_1+numroots}; -- indexed by the basis x_1*1,x_1*y_1,... 
    (EVal,EVec) := eigenvectors (inverse(D0)*D1); -- NEED GENERALIZED EIGENVALUES
    EVect := D0*EVec;
    apply(#EVal, i -> point{{EVal_i,EVect_(1,i)/EVect_(0,i)}}) -- roots (x_1,y_1) assuming x0=y0=1
)

needs "laurent-eigensolving.m2"
needs "documentation.m2"


TEST ///
-- Problem 2.11 in https://math.berkeley.edu/~bernd/cbms.pdf
n = 3
R = QQ[a_1..a_(n-1),b_1..b_(n-1)]
S = R[x,y]/ideal(x*y-1)
f = sum(n-1, i -> R_i*x^(i+1) + R_(n-1+i)*y^(i+1)) + x^n + y^n
I = sub(ideal apply(toList(2..2*n-1), i -> last coefficients(f^i, Monomials => {1_S})), R)
s = zeroDimSolve I
realPts = realPoints(s/(v -> point matrix{v}))
assert(#realPts == 6)
G = gens sub(I, CC[gens R])
-- tally apply(50, i -> (
    -- s = zeroDimSolve I;
    -- realPts = realPoints(s/(v -> point {v}));
    -- all(realPts, p -> clean(1e-8, evaluate(G, p)) == 0)
-- ))
///

<<<<<<< HEAD
TEST /// -- Complete intersection on P^1 x P^1
R = QQ[x0,x1,y0,y1,Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
J = ideal(random({2,5},R),random({2,5},R))
=======
TEST ///
A = QQ[x0,x1,y0,y1,Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
I = ideal(random({2,5},A),random({2,5},A))
R = QQ[x1,y1]
spe = map(R,A,matrix{{1,x1,1,y1}})
J=spe I
>>>>>>> fa905878fd83d7e0722be12c64d89f672a25934d
listSol = zeroDimSolve(J, Strategy => "syzygy")
assert(#listSol == 20)
///

end--

restart
needsPackage "EigenSolver"
loadPackage("EigenSolver", Reload => true)
installPackage("EigenSolver", RemakeAllDocumentation => true)
check "EigenSolver"
