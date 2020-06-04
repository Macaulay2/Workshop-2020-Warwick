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

---

--Example: Maximal minors of a generic 2x4 matrix.

--produces ideal of generic minors
genericMinors = (k, m, n) -> (
    -- k by k minors of a generic m by n matrix
    R = ZZ/101[x_(1,1) .. x_(m,n), MonomialOrder=>Lex];
		X = transpose genericMatrix(R, n, m);
    gens minors(k, X));

I = ideal genericMinors(2,2,4);
S = subring I_* 
subalgebraBasis S

peek S
peek S.cache
peek S.cache.SubalgComputations

toricSyz = S.cache.SubalgComputations.SyzygyIdeal
describe ring toricSyz
G = gens gb S.cache.SubalgComputations.SyzygyIdeal
rels = selectInSubring(1,G) --does elimination order on toric ideal to get only algebraic relations

liftedRel = (S.cache.SubalgComputations.Substitution(rels))_(0,0) -- substitutes original f's back into relation
P = ring(liftedRel)
phi = map(R,P,matrix{gens R | I_*})

phi(liftedRel) - phi(rels_(0,0))

factor phi(rels)_(0,0)

leadTerm phi(liftedRel)<leadTerm phi(leadTerm rels_(0,0))

---
--Example: 3-minors of a 3x5 matrix

I = ideal genericMinors(3,3,5); --maximal minors of a 2x4 matrix
S = subring I_* --make a subring out of gens I
subalgebraBasis S --compute SAGBI basis. Now information in S has been changed.

S.cache.SubalgComputations.SyzygyIdeal --to get toric ideal
G = gens gb S.cache.SubalgComputations.SyzygyIdeal --get GB of toric ideal
rels = sort (selectInSubring(1,G)) --does elimination order on toric ideal to get only algebraic relations

liftedRels = sort (S.cache.SubalgComputations.Substitution(rels)) -- substitutes original f's back into relation

P = ring(liftedRels)
phi = map(R,P,matrix{gens R | I_*})
factor phi(liftedRels_(0,1))

phi(liftedRels_(0,0))-phi(rels_(0,0))

leadTerm phi(liftedRels_(0,0))<leadTerm phi(leadTerm rels_(0,0))

