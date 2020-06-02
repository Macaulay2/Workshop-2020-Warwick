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
    else if opts.Strategy == "syzygy" then return null
    else error "Unknown strategy";
)

reduceCoeffs = (A, I, m) -> transpose sub(last coefficients(A % I, Monomials => m), ultimate(coefficientRing, ring I))

eigSolve1 = method(Options => options zeroDimSolve) -- based on Tomas' function SolvePolyEqByEigV
eigSolve1 Ideal := Matrix => opts -> I -> (
    R := ring I;
    B := if opts.Basis != null then opts.Basis else sub(basis(R/I), R);
    f := if opts.Multiplier != 0 then opts.Multiplier else random(1, R);
    M := reduceCoeffs(f*B, I, B);
    (E, V) := eigenvectors M;
    V = V*inverse diagonalMatrix V^{0};
    entries transpose(reduceCoeffs(vars R, I, B)*V)
)

-- needs "tomas-eigensolving.m2" -- errors when loading
needs "laurent-eigensolving.m2"
needs "documentation.m2"


TEST ///
-- Problem 2.11 in https://math.berkeley.edu/~bernd/cbms.pdf
n = 3
-- R = QQ[a_1..a_(n-1),b_1..b_(n-1), MonomialOrder => Eliminate (2*n-3)] -- doesn't work with elimination order
R = QQ[a_1..a_(n-1),b_1..b_(n-1)]
S = R[x,y]/ideal(x*y-1)
f = sum(n-1, i -> R_i*x^(i+1) + R_(n-1+i)*y^(i+1)) + x^n + y^n
I = sub(ideal apply(toList(2..2*n-1), i -> last coefficients(f^i, Monomials => {1_S})), R)
s = zeroDimSolve I
realPts = realPoints(s/(v -> point matrix{v}))
assert(#realPts == 6)
-- G = gens sub(I, CC[gens R])
-- tally apply(50, i -> (
    -- s = zeroDimSolve I;
    -- realPts = realPoints(s/(v -> point matrix{v}));
    -- all(realPts, p -> clean(1e-8, evaluate(G, p)) == 0)
-- )) -- seems around 6% fails
///

end--

restart
needsPackage "EigenSolver"
loadPackage("EigenSolver", Reload => true)
installPackage("EigenSolver", RemakeAllDocumentation => true)
check "EigenSolver"
