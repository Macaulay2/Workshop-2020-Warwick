-- Showcase File

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

S = QQ[apply(toSequence\subsets(n, k), J -> p_J) ]
presentationRing(A, S)
peek A.cache
presentation A -- Plucker relations

inducedWeights = for p in flatten entries M list (
    E := (exponents leadTerm p)#0;
    sum apply(E, D, (i, j) -> i*j)
    )

newPresRing = newRing(presentationRing A, MonomialOrder => {Weights => inducedWeights})
pluckerI = sub(presentation A, newPresRing)
I = ideal pluckerI
leadTerm(1, I)

sB = subalgebraBasis(A, Limit => 10);


-- invariants example showing %
--
restart
newDoc = 2020
needs "SubalgebraBases"

R = QQ[x1, x2, x3];
S = QQ[e1, e2, e3, y];
A = subring {
    x1 + x2 + x3, 
    x1*x2 + x1*x3 + x2*x3,
    x1*x2*x3,
    (x1 - x2)*(x1 - x3)*(x2 - x3)
    }
d=1
f = x1^(d+1)*x2^d+x2^(d+1)*x3^d+x3^(d+1)*x1^d
sub(f//A, presentationRing A)
presentation A
f%A
