---------------------------------------------
restart
printWidth = 45;
-*
---------------------------------------------
A TORIC VARIETY is an integral scheme such 
that an algebraic torus forms a Zariski open 
subscheme and the natural action of this 
torus on itself extends to the entire scheme.
---------------------------------------------
Michel Demazure (1970) provided the formal 
definition while studying the Cremona group.
---------------------------------------------
Miles Reid coined the phrase "toric variety" 
in his translation of Vladimir Danilov's 
influential 1978 survey.
---------------------------------------------
Toric varieties have been / are an important 
testing ground in the development of
  (*)  resolution of singularities,
  (*)  birational geometry,
  (*)  polyhedral geometry, and
  (*)  mirror symmetry.
---------------------------------------------
Because toric varieties are parametrized by 
monomials, they almost invariably arise in 
applications such as
  (*)  integer programming
  (*)  geometric modeling
  (*)  coding theory
  (*)  algebraic statistics
---------------------------------------------
Normal toric varieties correspond to strongly 
convex rational polyhedral fans. This makes 
the theory of normal toric varieties very 
explicit and computable.
---------------------------------------------
*-
needsPackage "NormalToricVarieties";
-- The projective plane is a toric variety
P = toricProjectiveSpace 2
dim P
rays P  -- primitive lattice points on rays
max P   -- maximal cones in fan
ring P  -- Cox ring
---------------------------------------------
-- Toric varieties include more exotic spaces
X = normalToricVariety(
    {{1,0,0},{0,1,0},{0,0,1},{0,-1,2},{0,0,-1},{-1,1,-1},{-1,0,-1},{-1,-1,0}},
    {{0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,1,5},{1,2,7},{2,3,7},{3,4,7},{4,5,6},{4,6,7},{5,6,7},{1,5,7}});
# rays X
isComplete X            -- compact
isProjective X          -- but not projective
isSmooth X
---------------------------------------------
-- Every normal toric variety has a 
-- resolution of singularities given by 
-- another normal toric variety.
Y = weightedProjectiveSpace {1,2,3,4};
(dim Y, isSmooth Y, isSimplicial Y)
Z = makeSmooth Y;     -- creates a resolution
isSmooth Z
# (set rays Z - set rays Y)
phi = map(Y, Z, 1)   -- birational surjection
---------------------------------------------
-- We can work with torus-invariant divisors. 
Y_0   -- irreducible divisor given by 0th ray
isEffective Y_0
isCartier Y_0
isQQCartier Y_0
isCartier (60*Y_0)
isAmple (60*Y_0)
pullback(phi, 60*Y_0)
---------------------------------------------
-- We can compute the cohomology of sheaves
P2 = P ** P   -- product of projective planes
dim P2
OO_P2
HH^0(P2, OO_P2)
OO_P2(-3,0)
HH^2(P2, OO_P2(-3,0))
HH^2(P2, OO_P2(0,-3))
HH^4(P2, OO_P2(-3,-3))
prune exteriorPower(4, cotangentSheaf P2) 
OO toricDivisor P2
---------------------------------------------