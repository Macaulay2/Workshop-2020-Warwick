-- First, let's compute the number of lines on a cubic surface in P^3.
restart
needsPackage "Schubert2"
-- Grassmannian Gr(2,4):
G = flagBundle({2,2}, VariableNames => {b,c}) 
(S,Q) = bundles G
rank S
rank Q
chern S
trim ideal intersectionRing G
A = intersectionRing G
chern S
chern Q
F = symmetricPower(3, Q)
chern F
chern_4 F
integral chern_4 F
integral chern F

-- Let's create an abstract variety X
-- and then do Riemann-Roch for a divisor on X
restart
needsPackage "Schubert2"
  -- Consider a 3-fold, with canonical class K, and chern classes c1=-K, c2, c3
  -- and with a divisor D.
  A = QQ[a, Degrees=>{0}][K, c_2, c_3, Degrees=>{1,2,3}, Join=>false][D, Join=>false]
  X = abstractVariety(3, A)
  X.TangentBundle = abstractSheaf(X, Rank=>3, ChernClass=>1-K+c_2+c_3)
  chern tangentBundle X
  OO(D)
  chern oo
  chi OO(D)
  chi OO(0_A)
  -- uses Grothendieck-Riemann-Roch:
  symmetricPower(a, tangentBundle X)
  chern oo
  symmetricPower(3, tangentBundle X)
  chern oo

-- toric variety
restart
needsPackage "NormalToricVarieties"
needsPackage "Schubert2"
X = smoothFanoToricVariety(4, 10)
isSmooth X
picardGroup X
pt = base(a,b,c)
Xa = abstractVariety(X, pt)
IX = intersectionRing Xa
describe IX
S = ring X
basis(1, IX)
G = OO(a*t_4 + b*t_5 + c*t_6)
chi G
sub(chi G, {a => 10, b => -9, c => -3})
sub(chi G, {a => 0, b => 1, c => 0})
D = X_5
for i from 0 to 4 list HH^i(X, OO(D))

Y = sectionZeroLocus OO(sum gens IX)
IY = intersectionRing Y
use IY
describe IY
chern exteriorPower(dim Y, cotangentBundle Y) -- chern classes of KY are all 0
F = OO(a*t_4 + b*t_5 + c*t_6)
chi F
chi OO(t_5)

for i from 0 to 3 list HH^i(X, OO(X_2))
for i from 0 to 3 list HH^i(X, OO(X_2 + toricDivisor X))
chi OO(t_2)  -- OO(X_2) is OO_X(D_2) which is OO(t_2), as chern class
             -- of OO(D_2) is t_2, computed from Hirzebruch-Riemann-Roch.

restart
needsPackage "Schubert2"
methods abstractVariety
help oo
A = QQ[a,b,c,d,Degrees=>{0,1,1,2}]
X = abstractVariety(2, A)
X.TangentBundle = abstractSheaf(X, Rank=>2, ChernClass => 1 + c + d)
integral(b*c)
integral d
integral(c^2)
integral(b^2)
TX = tangentBundle X
chern TX
ch TX
OOB = OO(b)
chern OOB
chi OOB


-- Now let's compute the number of lines on a smooth hypersurface X
-- of degree d in P^n
-- When d+1 = 2(n-1) = dim Gr(2, n+1), we expect a finite number
-- of lines on X.  What is this number?
-- Now let's try several such
n = 5 -- so d = 2n-3
X = flagBundle({n-1,2}, VariableNames => {,c})
dim X
(S,Q) = bundles X
A = intersectionRing X
rank S
rank Q
F = symmetricPower(2*n-3, Q)
chern F
integral chern F 

hashTable for n from 3 to 12 list (
    X = flagBundle({n-1,2});
    (S,Q) = bundles X;
    (2*n-3, n) => integral chern symmetricPower(2*n-3, Q)
    )


