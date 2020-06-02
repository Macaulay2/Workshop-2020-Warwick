newPackage(
    "EigenSolver",
    Version => "0.1", 
    Date => "June 2020",
    Authors => {{Name => "", 
	    Email => "", 
	    HomePage => ""}},
    Headline => "polynomial system solver relying on eigen-computations",
    DebuggingMode => false
    )
    
export {}

needs "tomas-eigensolving.m2"
needs "laurent-eigensolving.m2"    
needs "documentation.m2"

end--


-- Problem 2.11 in https://math.berkeley.edu/~bernd/cbms.pdf
n = 3
-- R = QQ[a_1..a_(n-1),b_1..b_(n-1), MonomialOrder => Eliminate (2*n-3)] -- doesn't work with elimination order
R = QQ[a_1..a_(n-1),b_1..b_(n-1)]
S = R[x,y]/ideal(x*y-1)
f = sum(n-1, i -> R_i*x^(i+1) + R_(n-1+i)*y^(i+1)) + x^n + y^n
I = sub(ideal apply(toList(2..2*n-1), i -> last coefficients(f^i, Monomials => {1_S})), R)
first SolvePolyEqByEigV I
