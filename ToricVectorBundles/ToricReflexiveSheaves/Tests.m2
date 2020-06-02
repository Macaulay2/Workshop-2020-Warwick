-- -*- coding: utf-8 -*-
------------------------------------------------------------------------------
-- Copyright 2017-20 Gregory G. Smith
--
-- This program is free software: you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by the Free
-- Software Foundation, either version 3 of the License, or (at your option)
-- any later version.
--
-- This program is distributed in the hope that it will be useful, but WITHOUT
-- ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
-- FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
-- more details.
--
-- You should have received a copy of the GNU General Public License along
-- with this program.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------------
-- TESTS
------------------------------------------------------------------------------
TEST ///
  -- TEST 0: verify everything works for the tangent bundle on projective 2-space
  (X,R) = (toricProjectiveSpace 2, QQ[e_1,e_2]);
  W = {{(-e_1-e_2,1),(e_1,0)},{(e_1,1),(e_2,0)},{(e_2,1),(e_1,0)}};
  T = toricReflexiveSheaf(W, X);
  assert isWellDefined T
  assert not isArithmeticallyFree T
  assert isLocallyFree T
  assChar = associatedCharacters T;
  assert(set assChar#0 === set ({0,-1},{1,-1}))
  assert(set assChar#1 === set ({ -1,0},{ -1,1}))
  assert(set assChar#2 === set ({1,0},{0,1}))  
  assert isGloballyGenerated T
  assert(rank T == dim X)
  assert((R,{1}) === ambient T)
  assert(X === variety T)
  assert(subspace(0,infinity,T) == ideal(0_R))
  assert(subspace(0,1,T) == ideal(R_0 + R_1))
  assert(subspace(0,0,T) == ideal(R_0, R_1))
  assert(subspace(0,-infinity,T) == ideal(R_0, R_1))    
  assert(subspace(1,infinity,T) == ideal(0_R))  
  assert(subspace(1,1,T) == ideal(R_0))
  assert(subspace(1,0,T) == ideal(R_0, R_1))
  assert(subspace(1,-infinity,T) == ideal(R_0, R_1))      
  assert(subspace(2,infinity,T) == ideal(0_R))    
  assert(subspace(2,1,T) == ideal(R_1))
  assert(subspace(2,0,T) == ideal(R_0, R_1))
  assert(subspace(2,-infinity,T) == ideal(R_0, R_1))  
  assert({R_0+R_1,R_0,R_1} == groundSet T)
  F = cover T;
  assert(rank F == # groundSet T)
  assert isWellDefined F
  assert isArithmeticallyFree F
  assert isLocallyFree F
  assert(gens F == 1)
  f = gens T;
  assert(isWellDefined f)
  assert(source f === F)
  assert(target f === T)
  assert(matrix f == matrix {{1_QQ, 1, 0}, {1, 0, 1}})  
  (R,d) = ambient F;
  assert(numgens R == 3)
  assert(d == {1})
  assert(X === variety F)	 
  describe F  
  assert(subspace(0,infinity,F) == ideal(0_R))
  assert(subspace(0,1,F) == ideal(R_0))
  assert(subspace(0,0,F) == ideal(R_0, R_1, R_2))
  assert(subspace(0,-infinity,F) == ideal(R_0, R_1, R_2))    
  assert(subspace(1,infinity,F) == ideal(0_R))  
  assert(subspace(1,1,F) == ideal(R_1))
  assert(subspace(1,0,F) == ideal(R_0, R_1, R_2))
  assert(subspace(1,-infinity,F) == ideal(R_0, R_1, R_2))      
  assert(subspace(2,infinity,F) == ideal(0_R))    
  assert(subspace(2,1,F) == ideal(R_2))
  assert(subspace(2,0,F) == ideal(R_0, R_1, R_2))
  assert(subspace(2,-infinity,F) == ideal(R_0, R_1, R_2))    
  assert(apply(# rays X, i -> X_i) === apply(components F, E -> toricDivisor E))
  assert(isWellDefined directSum apply(# rays X, i -> toricReflexiveSheaf X_i))
  assert(isWellDefined(T ** (X_0+X_1+X_2)))
///

TEST ///
  -- TEST 1: verify everything for direct sums of line bundles
  X = smoothFanoToricVariety(3,15)
  L = toricReflexiveSheaf (- toricDivisor X)
  assert(variety L === X)
  assert(rank L === 1)
  assert({1} === (ambient L)#1)
  assert(coefficientRing (ambient L)#0 === QQ)
  assert(numgens (ambient L)#0 === 1)
  assert isWellDefined L
  assert isArithmeticallyFree L
  assert isLocallyFree L
  assChar = associatedCharacters L;
  assert all(#max X, 
    i -> all(entries (matrix ((rays X)_((max X)#i) ) * (transpose matrix assChar#i)), 
      e -> e === {1}))
  assert isGloballyGenerated L
  assert(toricDivisor L == - toricDivisor X)
  assert(# groundSet L === 1)
  describe det L
  describe L
  assert(det L === L)
  assert(toricDivisor symmetricPower(0,L) == 0*X_0)
  assert(toricDivisor symmetricPower(2,L) == -2 * toricDivisor X)
  assert(toricDivisor symmetricPower(3,L) == -3 * toricDivisor X)
  assert(isWellDefined gens L)
  assert(matrix gens L == id_(QQ^1))
  D = X_0 + 2*X_1 + 3*X_2 + 4*X_3 + 2*X_4 + 2*X_5 + X_6
  L' = L ** D
  assert(isWellDefined L')
  assert(toricDivisor L' == D - toricDivisor X)  
  L'' = toricReflexiveSheaf( 2*X_0 + X_1 + 2*X_2 + X_3 + 2*X_4 + X_5 + 2*X_6)
  assert(isWellDefined L')    
  E = L'' ++ L
  assert(isWellDefined E)
  assert(rank E === 2)
  assert(variety E === X)
  assert(isArithmeticallyFree E)
  assert(set components E - set(L'', L) === set {})
  f = map(E,E,3)
  assert(isWellDefined f)
  assert(matrix f == 3*id_(QQ^2))
  assert(rank ker f === 0)
  assert(toricDivisor det E == toricDivisor L'' + toricDivisor L)
///

TEST ///
  -- TEST 2: a rank three toric vector bundle that is ample but not globally generated
  (X,R) = (toricProjectiveSpace 2, QQ[e_1,e_2,e_3]);
  W = {{(e_1-e_2,5), (e_2-e_3,2), (e_1,-3)}, {(e_1,8), (e_2,0), (e_3,-1)}, 
    {(e_3,5), (e_2,0), (e_1,-4)}};
  E = toricReflexiveSheaf(W, X);
  assert isWellDefined E
  assert not isArithmeticallyFree E
  assert isLocallyFree E
  assChar = associatedCharacters E;
  assert(set assChar#0 === set ({8,-5},{0,-5},{ -1,-1}))
  assert(set assChar#1 === set ({ -1,-4},{ -2,0},{ -2,5}))
  assert(set assChar#2 === set ({8,-4},{0,0},{ -1,5}))  
  assert not isGloballyGenerated E
  describe E
  assert(rank E === 3)
  assert(groundSet E === {e_1, -e_1+e_2, e_2, -e_2+e_3, e_3})
  assert(apply(components cover E, L -> toricDivisor L) === { 
      -3*X_0+8*X_1-4*X_2, 5*X_0-4*X_2, -3*X_0, 2*X_0-X_1, -3*X_0-X_1+5*X_2})
  Eneg = symmetricPower(-1,E);
  assert isWellDefined Eneg
  rank Eneg === 0
  E0 = symmetricPower(0,E);
  assert(isWellDefined E0)
  assert(rank E0 === 1)
  E1 = symmetricPower(1,E)
  assert(E1 === E)
  E2 = symmetricPower(2,E);
  assert(isWellDefined E2)
  assert(rank E2 === binomial(rank E + 2 - 1, 2))
  ER = ring (groundSet E2)#0
  --assert(groundSet E2 == {ER_0^2-2*ER_0*ER_1+ER_1^2, ER_0^2,
  --    -ER_0^2+ER_0*ER_1, ER_0*ER_1-ER_1^2-ER_0*ER_2+ER_1*ER_2, ER_0*ER_1,
  --    -ER_0*ER_1+ER_0*ER_2,-ER_0*ER_2+ER_1^2, ER_1^2-2*ER_1*ER_2+ER_2^2,
  --    -ER_1^2+ER_1*ER_2, ER_0*ER_2, -ER_0*ER_2+ER_1*ER_2, ER_1*ER_2,
  --    -ER_1*ER_2+ER_2^2, ER_2^2})
  F2 = cover E2;
  assert(rank F2 == # groundSet E2)
  assert(isWellDefined F2)
  assert(isArithmeticallyFree F2)
  f = gens E2;
  assert(isWellDefined f)
  assert(source f === F2)
  assert(target f === E2)
///

TEST ///
  -- TEST 3 : Euler sequence on projective 3-space
  (X,R) = (toricProjectiveSpace 3, QQ[e_1..e_3]);
  W = {{(-e_1-e_2-e_3,1),(e_1,0),(e_2,0)},{(e_1,1),(e_2,0),(e_3,0)},{(e_2,1),(e_1,0),(e_3,0)},{(e_3,1),(e_1,0),(e_2,0)}};
  T = toricReflexiveSheaf(W,X);
  assert(isWellDefined T)
  assert(not isArithmeticallyFree T)
  F0 = cover T
  assert(isWellDefined F0)
  assert(isArithmeticallyFree F0)
  assert(rank F0 === 4)
  d0 = gens T
  assert(isWellDefined d0)
  K0 = ker d0
  assert(isWellDefined K0)
  assert(rank K0 === 1)
  iota = map(F0, K0, 1)
  assert(isWellDefined iota)
  F1 = cover K0
  assert(isWellDefined F1)
  assert(rank F1 == rank K0)
  f1 = gens K0
  assert(isWellDefined f1)
  d1 = iota * f1
  assert(isWellDefined d1)
  assert(d0 * d1 == 0)
///  

TEST ///
  -- TEST 4 : verifies routine on the zero sheaf
  X = toricProjectiveSpace 2;
  Z = toricReflexiveSheaf X;
  assert isWellDefined Z
  assert isArithmeticallyFree Z
  assert isLocallyFree Z
  assert(associatedCharacters Z === apply(max X, s -> {}))
  assert isGloballyGenerated Z
  assert(coefficientRing ((ambient Z)#0) === QQ and
      numgens ((ambient Z)#0) === 0 and {1} === (ambient Z)#1)
  assert(rank Z === 0)
  assert(variety Z === X)
  assert(components Z === {Z})
  assert(Z ++ Z === Z)
  assert(directSum(11:Z) === Z)
  assert(subspace(0,0,Z) == 0)
  assert(subspace(0,-4,Z) == 0)
  assert(subspace(0,infinity,Z) == 0)
  assert(subspace(0,-infinity,Z) == 0)
  --assert isWellDefined map(Z,Z,1)
  assert(groundSet Z === {})
  assert isArithmeticallyFree Z
  assert(cover Z === Z)
  cover Z
  assert(gens Z == 0)
  assert(rank exteriorPower(-1,Z) === 0)
  assert(exteriorPower(1,Z) ==  0)  
  --assert(exteriorPower(1,Z) ===  Z)
  assert(rank exteriorPower(2,Z) === 0) 
  assert(rank symmetricPower(-1,Z) === 0)     
  --assert(rank symmetricPower(2,Z) === 0)       
///

TEST ///
  -- TEST 5 : verifies toricTangentBundle routine
  T = toricTangentBundle toricProjectiveSpace 2;
  assert isWellDefined T
  assert isLocallyFree T
  (R,d) = ambient T;
  for i from 0 to # rays variety T - 1 do (
    assert(subspace(i,0,T) == ideal basis(d,R));
    assert(subspace(i,1,T) == ideal(vars R * transpose matrix{(rays variety T)#i}));
    assert(subspace(i,2,T) == 0))
  T = toricTangentBundle toricProjectiveSpace 4;
  assert isWellDefined T
  (R,d) = ambient T;
  for i from 0 to # rays variety T - 1 do (
    assert(subspace(i,0,T) == ideal basis(d,R));
    assert(subspace(i,1,T) == ideal(vars R * transpose matrix{(rays variety T)#i}));
    assert(subspace(i,2,T) == 0))
  T = toricTangentBundle smoothFanoToricVariety(3,15);
  assert isWellDefined T  
  (R,d) = ambient T;
  for i from 0 to # rays variety T - 1 do (
    assert(subspace(i,0,T) == ideal basis(d,R));
    assert(subspace(i,1,T) == ideal(vars R * transpose matrix{(rays variety T)#i}));
    assert(subspace(i,2,T) == 0))
///

TEST ///
  -- TEST 6: testing isGloballyGenerated
  X = toricProjectiveSpace 2;
  T = toricTangentBundle X;
  assert isGloballyGenerated T
  Om1 = toricReflexiveSheaf(-X_0);
  assert not isGloballyGenerated Om1
  assert not isGloballyGenerated(T ++ Om1)
///

------------------------------------------------------------------------------
end---------------------------------------------------------------------------
------------------------------------------------------------------------------


------------------------------------------------------------------------------
-- methods to be added to package

dual ToricReflexiveSheaf := ToricReflexiveSheaf => E -> (  
ToricReflexiveSheaf ** ToricReflexiveSheaf := ToricReflexiveSheaf => (E,F) -> (
sheafHom (ToricReflexiveSheaf, ToricReflexiveSheaf) := ToricReflexiveSheaf => (E,F) -> (
euler ToricReflexiveSheaf := ZZ => E -> (
toricEuler ToricReflexiveSheaf := RingElement => E -> (

isAmple ToricReflexiveSheaf := Boolean => E -> (
isVeryAmple ToricReflexiveSheaf := Boolean => E -> (
isBasePointFree ToricReflexiveSheaf := Boolean => E -> (
cohomology ToricReflexiveSheaf := 

cokernel ToricReflexiveSheafMap
image ToricReflexiveSheafMap
coimage ToricReflexiveSheafMap
ToricReflexiveSheaf _ Array := ToricReflexiveSheafMap => (E,v) -> (
ToricReflexiveSheaf ^ Array := ToricReflexiveSheafMap => (E,v) -> (  
  
chowRing, chern classes, etc.

------------------------------------------------------------------------------
uninstallPackage "Parliaments"
restart
installPackage "Parliaments"
check "Parliaments"


-- Example 3.8
X = toricProjectiveSpace(1) ** toricProjectiveSpace (1);
S = QQ[e_1..e_3];
W = {{(e_1+e_2,1),(e_2,0),(e_3,-1)},{(e_1+e_2,0),(e_2,2),(e_3,1)},
  {(e_1+e_3,1),(e_2,-1),(e_3,0)},{(e_1+e_3,0),(e_2,2),(e_3,1)}}
E = toricReflexiveSheaf(W,X);
describe E
subspacePoset E
groundSet E
A = arrangement groundSet E
flats A
ideal X
rays X
length 
(R,d) = ambient E


-- Example: direct sum of two line bundles
X = toricProjectiveSpace(2);
S = QQ[e_1..e_2];
W = {{(e_1,1),(e_2,2)},{(e_1,0),(e_2,1)},{(e_1,0),(e_2,1)},}
E = toricReflexiveSheaf(W,X)
isWellDefined E
describe E
subspacePoset E
groundSet E

-- Example: tangent bundle ++ a line bundle
X = toricProjectiveSpace(2);
S = QQ[e_1..e_3];
W = {{(-e_1-e_2,3),(e_3,1),(e_1,-1)},{(e_1,3),(e_3,1),(e_2,-1)},
  {(e_2,4),(e_3,1),(e_1,0)}}
E = toricReflexiveSheaf(W,X)
isWellDefined E
subspacePoset E
groundSet E
needsPackage "HyperplaneArrangements";
A = arrangement groundSet E
flats A

X = toricProjectiveSpace 2;
E = mukaiLazarsfeldBundle (2*X_1)
exteriorPower(2,E)
isWellDefined E
E ** (2*X_0 - 2*X_1)
mukaiLazarsfeldBundle (2*X_0)
E#0

toricReflexiveSheaf({{},{},{}},X)


------------------------------------------------------------------
-- BERNT IVAR's EXAMPLE
restart
installPackage "Parliaments"
X = normalToricVariety({{1,0},{1,1},{1,2},{0,1},{ -1,0},{0,-1}},{{0,1},{1,2},{2,3},{3,4},{4,5},{0,5}})
isWellDefined X
isSmooth X
S = QQ[e_1,e_2];
L = {{(e_2,130),(e_1,339)}, {(e_2,-246),(e_1,324)}, {(e_2,-375),(e_1,442)}, {(e_2,61),(e_1,232)}, 
  {(e_1,79),(e_2,535)}, {(e_2,91),(e_1+e_2,414)}}
E = toricReflexiveSheaf(L,X)  
isWellDefined E
groundSet E

describe E

E#0
subspace(E,0,
  

needsPackage "Polyhedra"

------------------------------------------------------------------------------
XXX
------------------------------------------------------------------------------
uninstallPackage "Parliaments"
restart
path = prepend("~/MyDocuments/Research/ToricVectorBundles/", path)
installPackage "Parliaments"
check "Parliaments"

restart
needsPackage "Parliaments"

X = toricProjectiveSpace 2
E = toricTangentBundle X
outerNormals = matrix rays X;
assChar = associatedCharacters E
G = groundSet E
(R,d) = ambient E
(M,C) = coefficients( matrix {G}, Monomials => basis(d,R))
divisors = apply(components cover E, L -> toricDivisor L)
Sigma = max X

i = 0
j = 0
u = assChar#i#j
(max X)#i
charVector = transpose matrix {u}
(first transpose entries(outerNormals * charVector - matrix vector divisors#j))_((max X)#i)

coordinates = submatrix(outerNormals * charVector - matrix vector divisors#j,(max X)#i,)
first transpose entries coordinates

(X,R) = (toricProjectiveSpace 2, QQ[e_1,e_2,e_3]);
W = {{(e_1-e_2,5), (e_2-e_3,2), (e_1,-3)}, {(e_1,8), (e_2,0), (e_3,-1)}, 
  {(e_3,5), (e_2,0), (e_1,-4)}};
E = toricReflexiveSheaf(W, X);
outerNormals = matrix rays X
assChar = associatedCharacters E
G = groundSet E
(R,d) = ambient E
(M,C) = coefficients( matrix {G}, Monomials => basis(d,R))
divisors = apply(components cover E, L -> toricDivisor L)
i = 2
max X
u = assChar#i#1
      potentialBases := for u in assChar#i list (
      	charVector = transpose matrix {u}
      	select(#G, j -> (
	    j = 1
	    coordinates = outerNormals*charVector - matrix vector divisors#j
	    coordinates_(1,0)
	    member(1, (max X)#i)
	    
	    << coordinates << endl;
	    localCoordinates := transpose submatrix(coordinates,(max X)#i,);
	    all(first entries coordinates, c -> c == 0))
	  )
      any(boxProduct potentialBases, s -> det C_s =!= 0))));  




isGloballyGenerated = method()
isGloballyGenerated ToricReflexiveSheaf := Boolean => E -> (
  if E == 0 then return true;
  X := variety E;
  outerNormals := matrix rays X;
  n := # rays X;
  assChar := associatedCharacters E;
  G := groundSet E;
  (R,d) := ambient E;
  (M,C) := coefficients( matrix {G}, Monomials => basis(d,R));
  divisors := apply(components cover E, L -> toricDivisor L);
  all(#assChar, i -> (
      potentialBases := for u in assChar#i list (
      	charVector := transpose matrix {u};
      	select(#G, j -> (
	    coordinates := outerNormals*charVector - matrix vector divisors#j;
	    all(n, k -> if member(k, (max X)#i) then coordinates_(k,0) == 0 
	      else coordinates_(k,0) <= 0))));
      any(boxProduct potentialBases, s -> det C_s =!= 0))));  




X = smoothFanoToricVariety(3,15)
L = toricReflexiveSheaf (- toricDivisor X)
OO toricDivisor L
describe L
sigma = (max X)#0
U = normalToricVariety((rays X)_sigma, {toList(0..#sigma - 1)} )
isWellDefined U
LU = toricReflexiveSheaf(for i in sigma list pairs L#i,U);
isWellDefined LU
G = groundSet LU
reverse sort unique values LU#0
G#0 % subspace(0, 1, LU)
C = transpose matrix for g in G list (
    for i to # rays U - 1 list (
      weights := reverse sort unique values LU#i;
      j := 0;
      while g % subspace(i, weights#j, LU) != 0 do j = j+1;
      j));
   <<

    entries transpose ((inverse matrix rays U) * C)


transpose inverse matrix rays U

associatedCharacters T
max X
rays X
G = groundSet TU
C = for g in G list (
  for i from 0 to # rays U - 1 list (
    weights := reverse sort unique values TU#i;
    j := 0;
    while g % subspace(i, weights#j, TU) != 0 do j = j+1;
    j))
entries transpose ((inverse matrix rays U)*matrix C)




W = {{(-e_1-e_2,1),(e_1,0)},{(e_1,1),(e_2,0)},{(e_2,1),(e_1,0)}};
T = toricReflexiveSheaf(W,X);

(X,R) = (toricProjectiveSpace 2, QQ[e_1, e_2])
D = X_0
ML =  D
rank ML
isWellDefined ML
F0 = cover ML
rank F0
epsilon = gens ML
K0 = ker epsilon
isWellDefined K0
F1 = cover K0
rank F1 
gens K0
d1 = map(F0, K0, 1) * (gens K0)
K1 = ker d1
K1 == 0 
epsilon * d1
C = complex {matrix epsilon, matrix d1}
degree C_0
(C, Weights => {-1})
help betti


C0 = (target matrix epsilon) ** QQ[]
C1 = (source matrix epsilon) ** QQ[]
C2 = (source matrix d1) ** QQ[]
degrees C1
degrees matrix d1
help (map, Module, Module, Matrix)
C = chainComplex (map(C0,C1,matrix epsilon),map(C1,C2, matrix d1))
C0 = target matrix epsilon
C1 = (source matrix epsilon) ** QQ[]
C2 = (source matrix d1) ** QQ[]
map(C0,C1,(matrix epsilon) ** QQ[])
betti C
break

options betti


needsPackage "Complexes"

D = 2*X_0
ML = mukaiLazarsfeldBundle D
rank ML
ambient ML
isWellDefined ML
F0 = cover ML
rank F0
epsilon = gens ML
K0 = ker epsilon
isWellDefined K0
F1 = cover K0
rank F1 
gens K0
d1 = map(F0, K0, 1) * (gens K0)
K1 = ker d1
K1 == 0 
epsilon * d1

betti 
chainComplex(matrix epsilon, matrix d1)

break


------------------------------------------------------------------------------
-- ZZZ
restart
needsPackage "Parliaments";
X = toricProjectiveSpace 2;
D = 2*X_0;
ML = mukaiLazarsfeldBundle D
E = exteriorPower(2, ML)
res E
res (E ** D)
res (ML ** D)

C = time res E;
C
display C

D = {{0, rank target maps#0}} | apply(#maps, 
  i -> {i+1, rank source maps#i} | apply(components source maps#i, L -> toricDivisor L))
transpose 
transpose netList
Table transpose D
help table

display
net 


D


------------------------------------------------------------------------------
-- SYZYGIES OF VERONESE EMBEDDINGS?
--

restart
needsPackage "Parliaments"

display = C -> (
  D := netList( transpose {{ -1, rank target C#0}}, 
    Boxes => false, Alignment => Center, VerticalSpace => 1);
  for i to #C - 1 do (
    D = D | "  " | (netList ( {i, rank source C#i} | apply(components source C#i, 
	  L -> HH^i(variety L, OO toricDivisor L)), 
  	Boxes => false, Alignment => Center, VerticalSpace => 1)));
  D);

-- PP^2;  OO(3)
X = toricProjectiveSpace 2;
D = 3 * X_0;
ML = mukaiLazarsfeldBundle D;
C = res (ML ** D);
display C
time display res (exteriorPower(2,ML) ** D)
time display res (exteriorPower(3,ML) ** D)
time display res (exteriorPower(4,ML) ** D)
time display res (exteriorPower(5,ML) ** D)
time display res (exteriorPower(6,ML) ** D)
time display res (exteriorPower(7,ML) ** D)
time display res (exteriorPower(8,ML) ** D)

-- PP^2;  OO(4)
X = toricProjectiveSpace 2;
D = 4 * X_0;
ML = mukaiLazarsfeldBundle D;
C = res (ML ** D)
display C
time display res (exteriorPower(2,ML) ** D)
time display res (exteriorPower(3,ML) ** D)
time display res (exteriorPower(4,ML) ** D)
time display res (exteriorPower(5,ML) ** D)
time display res (exteriorPower(6,ML) ** D)
time display res (exteriorPower(7,ML) ** D)
time display res (exteriorPower(8,ML) ** D)
time display res (exteriorPower(9,ML) ** D)
time display res (exteriorPower(10,ML) ** D)
time display res (exteriorPower(11,ML) ** D)

-- PP^2;  OO(5)
X = toricProjectiveSpace 2;
D = 5 * X_0;
ML = mukaiLazarsfeldBundle D;
C = res (ML ** D);
display C
time display res (exteriorPower(2,ML) ** D)
time display res (exteriorPower(3,ML) ** D)
time display res (exteriorPower(4,ML) ** D)
time display res (exteriorPower(5,ML) ** D)
time display res (exteriorPower(6,ML) ** D)
time display res (exteriorPower(7,ML) ** D)
time display res (exteriorPower(8,ML) ** D)
time display res (exteriorPower(9,ML) ** D)
time display res (exteriorPower(10,ML) ** D)
time display res (exteriorPower(11,ML) ** D)
time display res (exteriorPower(12,ML) ** D)
time display res (exteriorPower(13,ML) ** D)


-- PP^3;  OO(3)
X = toricProjectiveSpace 3;
D = 3 * X_0;
ML = mukaiLazarsfeldBundle D;
C = res (ML ** D);
display C
time display res (exteriorPower(2,ML) ** D)
time display res (exteriorPower(3,ML) ** D)
time display res (exteriorPower(4,ML) ** D)
time display res (exteriorPower(5,ML) ** D)
time display res (exteriorPower(6,ML) ** D)
time display res (exteriorPower(7,ML) ** D)
time display res (exteriorPower(8,ML) ** D)

-- PP^3;  OO(4)
X = toricProjectiveSpace 3;
D = 4 * X_0;
ML = mukaiLazarsfeldBundle D;
C = res (ML ** D);
display C
time display res (exteriorPower(2,ML) ** D)
time display res (exteriorPower(3,ML) ** D)
time display res (exteriorPower(4,ML) ** D)
time display res (exteriorPower(5,ML) ** D)
time display res (exteriorPower(6,ML) ** D)
time display res (exteriorPower(7,ML) ** D)
time display res (exteriorPower(8,ML) ** D)
time display res (exteriorPower(9,ML) ** D)
time display res (exteriorPower(10,ML) ** D)
time display res (exteriorPower(11,ML) ** D)

-- PP^3;  OO(5)
X = toricProjectiveSpace 3;
D = 5 * X_0;
ML = mukaiLazarsfeldBundle D;
C = res (ML ** D);
display C
time display res (exteriorPower(2,ML) ** D)
time display res (exteriorPower(3,ML) ** D)
time display res (exteriorPower(4,ML) ** D)
time display res (exteriorPower(5,ML) ** D)
time display res (exteriorPower(6,ML) ** D)
time display res (exteriorPower(7,ML) ** D)
time display res (exteriorPower(8,ML) ** D)
time display res (exteriorPower(9,ML) ** D)
time display res (exteriorPower(10,ML) ** D)
time display res (exteriorPower(11,ML) ** D)
time display res (exteriorPower(12,ML) ** D)
time display res (exteriorPower(13,ML) ** D)



-----------------------------------------------------------------
-- Andreas Hochenegger question
restart
needsPackage "Parliaments";
X = toricProjectiveSpace 2;
Omega = toricCotangentBundle X;
E = symmetricPower(2,Omega) ** (X_0+X_1+X_2)
parliament = apply(components cover E, L -> toricDivisor L)
#parliament
for D in parliament list HH^0(X, OO D)



------------------------------------------------------------------------------
-- Hiro example
restart
needsPackage "Parliaments";
(X,R) = (hirzebruchSurface 1, QQ[e_1,e_2]);
W = {{(e_1,4),(e_2,-2)}, {(e_1,3),(e_2,2)}, {(e_2,5),(e_1,0)}, {(e_1+e_2,3),(e_1,-1)}}
E = toricReflexiveSheaf(W,X);
describe E
resolution E
characters E
orbits(X, 0)
associatedCharacters E
toricDivisor exteriorPower(2,E)
intersectionRing X



(X,R) = (toricProjectiveSpace 2, QQ[y_1,y_2]);
E = toricTangentBundle X
describe E
A = ZZ/32003[gens ring X | gens (ambient E)#0];
B = ZZ/32003[z_0..z_8];
f = map(A,B, {A_0*A_1*A_2*A_3, A_0^2*A_2*A_3, A_0*A_2^2*A_3,
	 A_0*A_1^2*A_4, A_0^2*A_1*A_4, A_0*A_1*A_2*A_4,
	 A_1^2*A_2*(A_3+A_4), A_0*A_1*A_2*(A_3+A_4), A_1*A_2^2*(A_3+A_4)})
I = saturate(ideal mingens ker f, ideal gens B)
netlist =  I_*
codim I
hilbertPolynomial(I, Projective => false)

C = ZZ/32003[gens ring X | gens B];
gens C
M = matrix{
    {C_0*C_1*C_2, C_0^2*C_2, C_0*C_2^2, 0_C, 0_C, 0_C, C_1^2*C_2, C_0*C_1*C_2, C_1*C_2^2},
    {0_C, 0_C, 0_C, C_0*C_1^2, C_0^2*C_1, C_0*C_1*C_2, C_1^2*C_2, C_0*C_1*C_2, C_1*C_2^2},
    {C_3,C_4,C_5, C_6, C_7, C_8, C_9, C_10, C_11}}
saturate(minors(3,M), ideal(C_0,C_1,C_2))

M = matrix{{ 

(X,R) = (toricProjectiveSpace 2, QQ[y_1,y_2]);
E = toricTangentBundle X
describe E
A = ZZ/32003[gens ring X | gens (ambient E)#0];
B = ZZ/32003[z_0..z_8];
f = map(A,B, {A_0*A_1*A_2*A_3, A_0^2*A_2*A_3, A_0*A_2^2*A_3,
	 A_0*A_1^2*A_4, A_0^2*A_1*A_4, A_0*A_1*A_2*A_4,
	 A_1^2*A_2*(A_3+A_4), A_0*A_1*A_2*(A_3+A_4), A_1*A_2^2*(A_3+A_4)})
I = saturate(ideal mingens ker f, ideal gens B)
codim I
hilbertPolynomial(I, Projective => false)

