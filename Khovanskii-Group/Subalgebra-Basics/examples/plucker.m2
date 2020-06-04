-- Plucker algebra example following:
-- "Toric Degenerations of Grassmannians from Matching Fields" - Mohammadi & Shaw

restart
needsPackage "SubalgebraBases"

n = 4;
k = 2;

D = apply((1,1) .. (k,n), (i,j) -> n*k - (n-j)*(i-1))
--Diagonal matching field weight matrix

R = QQ[x_(1,1) .. x_(k,n), MonomialOrder => {Weights => D}, Global => false];
X = transpose genericMatrix(R, n, k);
M = matrix { for J in subsets(n, k) list det X_J };
A = subring M; -- Plucker algebra

setWeight(A, toList D)
S = QQ[apply(toSequence\subsets(n, k), J -> p_J) ]
S = presentationRing(A, S)
peek A.cache
presentation A -- Plucker relations
leadTerm(1, presentation A)

-*
=====
Notes
=====

In order use setWeight successfully you must first call
setWeight before the presentation or presentationRing

*-




-- The above does the same as the following
inducedWeights = for p in flatten entries M list (
    E := (exponents leadTerm p)#0;
    sum apply(E, D, (i, j) -> i*j)
    )

newPresRing = newRing(presentationRing A, MonomialOrder => {Weights => inducedWeights})
pluckerI = sub(presentation A, newPresRing)
I = ideal pluckerI
leadTerm(1, I)




