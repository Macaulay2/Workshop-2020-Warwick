
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
-- CODE
------------------------------------------------------------------------------
--- kludge to access parts of the 'Core'
hasAttribute = value Core#"private dictionary"#"hasAttribute";
getAttribute = value Core#"private dictionary"#"getAttribute";
ReverseDictionary = value Core#"private dictionary"#"ReverseDictionary";

------------------------------------------------------------------------------
-- Basic features of the toric reflexive sheaf datatype
------------------------------------------------------------------------------
ToricReflexiveSheaf = new Type of HashTable; 
ToricReflexiveSheaf.synonym = "toric reflexive sheaf"
ToricReflexiveSheaf.GlobalAssignHook = globalAssignFunction
ToricReflexiveSheaf.GlobalReleaseHook = globalReleaseFunction
net ToricReflexiveSheaf := net @@ expression
expression ToricReflexiveSheaf := E -> if hasAttribute(E, ReverseDictionary) 
    then expression getAttribute(E, ReverseDictionary) else 
    Describe (expression  toricReflexiveSheaf) (expression "...", expression variety E)
describe ToricReflexiveSheaf := E -> (
    local klyachkoData;
    if rank E === 0 then return net expression 0
    else (
    	n := # rays variety E;
    	for i to n-1 do (
      	    flag := {{toString(i) | ": "}, {null}};
      	    P := sort apply(pairs E#i, p -> {p#1,p#0});
      	    for p in P do (
      		flag = {flag#0 | {p#1}, flag#1 | {p#0}});
      	    if i === 0 then klyachkoData = flag
      	    else klyachkoData = klyachkoData | flag
	    );
    	net Table klyachkoData
	)
    )
ToricReflexiveSheaf#{Standard,AfterPrint} = 
ToricReflexiveSheaf#{Standard,AfterNoPrint} = E -> (
    << endl;
    << concatenate(interpreterDepth:"o") << lineNumber;
    << " : ToricReflexiveSheaf on " << variety E << endl;
    )  

ambient ToricReflexiveSheaf := Sequence => E -> E.ambient
variety ToricReflexiveSheaf := NormalToricVariety => E -> E.variety
normalToricVariety ToricReflexiveSheaf := NormalToricVariety => opts -> E -> E.variety
rank ToricReflexiveSheaf := ZZ => E -> E.rank

isWellDefined ToricReflexiveSheaf := Boolean => E -> (
    n := # rays variety E;  
    -- check keys
    K := keys E;
    expectedKeys := set (toList(0..n-1) | 
	{symbol ambient, symbol rank, symbol variety, symbol cache}
	);
    if set K =!= expectedKeys then (
    	if debugLevel > 0 then (
      	    added := toList(K - expectedKeys);
      	    missing := toList(expectedKeys - K);
      	    if #added > 0 then << "-- unexpected key(s): " << toString added << endl;
      	    if #missing > 0 then << "-- missing key(s): " << toString missing << endl
	    );
    	return false
	);
    -- check types
    for i from 0 to n-1 do (
    	if not instance(E#i, HashTable) then (
      	    if debugLevel > 0 then 
      	        << "-- expected 'E#" << i << "' to be a HashTable" << endl;
      	    return false
	    );
    	if not all(keys E#i, k -> instance(k, RingElement)) then (
	    if debugLevel > 0 then 
	        << "-- expected the keys of 'E#" << i << "' to be ring elements" << endl;
	    return false
	    );
    	if not all(values E#i, v -> instance(v, ZZ)) then (
	    if debugLevel > 0 then 
	        << "-- expected the values of 'E#" << i << "' to be integers" << endl;
	    return false
	    );
	);
    if not instance(E.ambient, Sequence) or length E.ambient =!= 2 then (
	if debugLevel > 0 then 
  	    << "-- expected 'E.ambient' to be a Sequence with two entries" << endl;
	return false
	);
    if not instance(E.ambient#0, PolynomialRing) then (
    	if debugLevel > 0 then 
      	    << "-- expected 'E.ambient#0' to be a PolynomialRing" << endl;
	return false
	);
    if not instance(E.ambient#1, List) then (
    	if debugLevel > 0 then 
	    << "-- expected 'E.ambient#1' to be a List" << endl;  
    	return false
	);
    if not all(length E.ambient#1, j -> instance(j, ZZ)) then (
      	if debugLevel > 0 then 
	    << "-- expected 'E.ambient#1' to be a list of integers" << endl;  
    	return false
	);
    if length E.ambient#1 =!= degreeLength E.ambient#0 then (
      	if debugLevel > 0 then 
	    << "-- expected 'E.ambient#1' to be a degree in the ambient ring" << endl;
	return false
	);
    if not instance(E.rank, ZZ) then (
    	if debugLevel > 0 then 
      	    << "-- expected 'E.rank' to be a ZZ" << endl; 
	return false
	);  	
    if not instance(E.variety, NormalToricVariety) then (
    	if debugLevel > 0 then 
      	    << "-- expected 'E.variety' to be a NormalToricVariety" << endl;
    	return false
	);
    if not instance(E.cache, CacheTable) then (
    	if debugLevel > 0 then 
	    << "-- expected `E.cache' to be a CacheTable" << endl;    
    	return false
	);
    -- check flags
    (R, d) := ambient E;
    for i from 0 to n-1 do (
    	if not all(keys E#i, k -> ring k === R) then (
      	    if debugLevel > 0 then (
		<< "-- expected the keys of 'E#"; 
		<< i << "' to be elements in the ambient ring" << endl
		);
      	    return false
	    );
    	if not all(keys E#i, k -> isHomogeneous k) then (
      	    if debugLevel > 0 then (
		<< "-- expected the keys of 'E#";
		<< i << "' to be homogeneous" << endl
		);
      	    return false
	    );
    	if not all(keys E#i, k -> degree k === d) then (
      	    if debugLevel > 0 then (
		<< "-- expected the keys of 'E#";
		<< i << "' to have degree equal to the ambient degree" << endl
		);
      	    return false
	    );    
    	if rank E =!= rank source mingens ideal keys E#i then (
      	    if debugLevel > 0 then (
		<< "-- expected the keys of 'E#" << i;
		<< "' to form a basis for a vector space of dimension " << rank E << endl
		);
      	    return false
	    )
	);
    I := ideal keys E#0;
    if not all(n, i -> I == ideal keys E#i) then (
    	if debugLevel > 0 then (
      	    << "-- expected the ring elements on each ray ";
	    << "to generate the same ideal" << endl
	    );
    	return false
	);    
    true
    )

ToricReflexiveSheaf == ZZ := Boolean => (E, m) -> (
    if m =!= 0 then 
	error "attempted to compare a sheaf to a nonzero integer";
    rank E === 0
    )
ZZ == ToricReflexiveSheaf := Boolean => (m, E) -> E == m


toricReflexiveSheaf = method()
toricReflexiveSheaf (List, NormalToricVariety) := ToricReflexiveSheaf => (L, X) -> (
    n := # L;
    if n =!= # rays X then   
      	error "expected a filtration for each torus-invariant divisor";
    r := #L#0;
    if r === 0 then error "expected a nonempty filtration";
    R := ring L#0#0#0;
    d := degree L#0#0#0;
    new ToricReflexiveSheaf from hashTable (
      	apply(n, i -> i => hashTable L#i) | {
      	    symbol ambient => (R, d),
      	    symbol rank => r,
      	    symbol variety => X,
      	    symbol cache => new CacheTable}
	)
    )
-- constructor for the zero sheaf
toricReflexiveSheaf NormalToricVariety := ToricReflexiveSheaf => X -> (
    n := # rays X;
    (R, d) := (QQ[], {1});
    new ToricReflexiveSheaf from hashTable (
      	apply(n, i -> i => hashTable {}) | {
      	    symbol ambient => (R, d),
      	    symbol rank => 0,
      	    symbol variety => X,
      	    symbol cache => new CacheTable}
	)
    )
toricReflexiveSheaf ToricDivisor := ToricReflexiveSheaf => (
    cacheValue symbol ToricReflexiveSheaf) (D -> (
    	c := entries vector D;
    	e := symbol e;
    	R := QQ(monoid[e]);
    	n := # rays variety D;
    	W := for i to n-1 list {{R_0,c#i}};
    	X := variety D;  
    	toricReflexiveSheaf(W, X) 
	)
    )

toricDivisor ToricReflexiveSheaf := ToricDivisor => opts -> (
    (cacheValue symbol toricDivisor) (E -> (
    	    r := rank E;
    	    if r =!= 1 then error "expect a rank-one sheaf";
    	    X := variety E;
    	    n := # rays X;
    	    g := first groundSet E;
    	    D := sum for i from 0 to n-1 list (
      		weights := reverse sort unique values E#i;
      		j := 0;
      		while g % subspace(i,weights#j, E) != 0 do j = j + 1;
      		weights#j * X_i
		);
    	    D 
	    )
	)
    )

toricTangentBundle = method();
toricTangentBundle NormalToricVariety := ToricReflexiveSheaf => X -> (
    d := dim X;
    e := symbol e;
    R := QQ(monoid[e_1..e_d]);
    raysX := rays X;
    n := # rays X;
    W := for i from 0 to n-1 list (
    	f := (vars R * transpose matrix {raysX#i})_(0,0);
    	basisSet := {(f, 1)};
    	for g in select(gens R, h -> h // f == 0) do basisSet = basisSet | {(g, 0)};
    	basisSet
	);
    toricReflexiveSheaf(W, X)
    )

toricCotangentBundle = method();
toricCotangentBundle NormalToricVariety := ToricReflexiveSheaf => X -> (
    d := dim X;
    e := symbol e;
    R := QQ(monoid[e_1..e_d]);
    raysX := rays X;
    n := # rays X;
    W := for i from 0 to n-1 list (
    	M := matrix { raysX#i };
    	basisSet := { ( (vars R * transpose M)_(0,0), -1 ) };
    	for g in first entries (vars R * gens ker M) do (
      	    basisSet = basisSet | { (g,0) } );
    	basisSet
	);
    toricReflexiveSheaf(W, X)
    )

-*
TODO: weights are correct but vectors are not; need to construct complements
dual ToricReflexiveSheaf := ToricReflexiveSheaf => {} >> o -> E -> (
  X := variety E;
  n := # rays X;
  W := {};
  for i from 0 to n-1 do W = W | {apply(pairs E#i, p -> (p#0, - p#1))};
  toricReflexiveSheaf(W,X))
*-

trim ToricReflexiveSheaf := ToricReflexiveSheaf => opts -> E -> (
    if rank E === 0 then return E;
    chosenBasis := gens gb subspace(0, -infinity, E);
    r := rank E;
    KK := coefficientRing (ambient E)#0;
    e := symbol e;
    newR := KK(monoid[e_1..e_r]);
    X := variety E;
    n := # rays X;
    W := for i from 0 to n-1 list (
    	apply(pairs E#i, 
      	    p -> ( (vars newR * sub(p#0 // chosenBasis, KK))_(0,0), p#1 ) ));
    F := toricReflexiveSheaf(W,X);
    (R, d) := ambient E;
    BE := basis(d,R);
    F.cache.trim = map(E, F, 
    	sub( (coefficients( chosenBasis, Monomials => BE) )#1, KK));
    F
    )

ToricReflexiveSheaf.directSum = args -> (
    sheaves := select(args, E -> rank E =!= 0);
    m := #sheaves;
    if m === 0 then return args#0;
    X := variety sheaves#0;  
    if not all(sheaves, E -> variety E === X) then 
    	error "expected all sheaves to be over the same variety";
    p := apply(toList sheaves, E -> numgens (ambient E)#0);
    s := apply(m, k -> sum(k, j -> p#j));
    e := symbol e;
    R := QQ(monoid[e_0..e_(sum(p)-1)]);
    phi := new HashTable from apply(m, k -> 
      	    k => map(R, (ambient sheaves#k)#0, apply(toList(s#k .. s#k+p#k-1), j -> R_j))
	);
    n := # rays X;
    W := for i from 0 to n-1 list flatten apply(m, 
    	k -> apply(pairs sheaves#k#i, (r,w) -> (phi#k(r),w)));
    E := toricReflexiveSheaf(W, X);
    E.cache.components = toList sheaves;
    E
    )
ToricReflexiveSheaf ++ ToricReflexiveSheaf := ToricReflexiveSheaf => (E,F) -> directSum(E,F)
directSum ToricReflexiveSheaf := E -> directSum (1 : E)
components ToricReflexiveSheaf := E -> if E.cache.?components then 
    E.cache.components else {E}

-*
TODO: 
ToricReflexiveSheaf ** ToricReflexiveSheaf := ToricReflexiveSheaf => (E, F) -> (
    ??
    )
*-

tensor (ToricReflexiveSheaf, ToricReflexiveSheaf) := ToricReflexiveSheaf => opts ->(E1, E2) -> (
    -- The logic of this code can compute tensor product of arbitrary
    -- number of arguments. To interface with M2 tensor, it only allows two
    -- arguments.

    -- sheaves :=  select(args, E -> rank E =!= 0);
    sheaves := (E1, E2);
    numSheaves := #sheaves;
    ambients := apply(sheaves, E -> ambient E);
    -- if numSheaves === 0 then return args#0;
    X := variety sheaves#0;
    if not all(sheaves, G -> variety G === X) then
        error "expected all sheaves to be over the same variety";
    fiberDims := apply(sheaves,
	               G -> rank source basis ((ambient G)#1, (ambient G)#0));
    -- Construct a new ambient ring.
    R := tensor(apply(ambients, i -> i#0));

    startIndices := apply(numSheaves, k -> sum(k, j -> fiberDims#j));
    maps := apply(numSheaves,
	          i -> map(R, ambients#i#0, apply(toList ((startIndices#i)..(startIndices#i + fiberDims#i-1)),
			  j -> R_j)));

    -- Compute tensor products of basis elements and their weights.
    W := apply(# rays X, i -> (
	basisPairs := apply(numSheaves, j -> apply( pairs sheaves#j#i, p -> (maps_j(p_0), p_1)));
	makeTensors := (l1, l2) ->
	    flatten apply (l2, p -> apply(l1, q -> (p_0*q_0, p_1+q_1)));
	fold(makeTensors, basisPairs_0, basisPairs_{1..numSheaves-1})
	)
    );
    E := toricReflexiveSheaf(W, X);
    E.cache.components = toList sheaves;
    E
    )

ToricReflexiveSheaf ** ToricReflexiveSheaf := ToricReflexiveSheaf => (E,F) -> tensor(E,F)


ToricReflexiveSheaf ** ToricDivisor := ToricReflexiveSheaf => (E,D) -> (
    c := entries vector D;
    X := variety D;
    n := # rays X;
    W := for i from 0 to n-1 list for p in pairs E#i list (p#0,p#1+c#i);
    toricReflexiveSheaf(W, X)
    )
ToricDivisor ** ToricReflexiveSheaf := ToricReflexiveSheaf => (D,E) -> E ** D

exteriorPower (ZZ, ToricReflexiveSheaf) := ToricReflexiveSheaf => opts -> (p, E') -> (  
    E := trim E';
    X := variety E;
    r := rank E;
    (R, d) := ambient E;
    if p < 0 or p > r then return toricReflexiveSheaf X;
    if p === 0 then return toricReflexiveSheaf (0 * X_0);
    if p === 1 then return E';
    KK := coefficientRing (ambient E)#0;
    e := symbol e;
    newR := KK(monoid[e_1..e_r, SkewCommutative => true]);
    t := symbol t;
    degRing := KK(monoid[t]);
    n := # rays X;
    W := {};
    for i to n-1 do (
	psi := sub( 
	    (coefficients(matrix{keys E#i}, Monomials => vars R))#1, degRing);
	psi = exteriorPower(p, map(degRing^(-values E#i), degRing^r, psi));
	weights := degrees target psi;
	psi = basis(p, newR) * sub(psi, KK);
	W = W | { 
	    for k to binomial(r,p) -1 list ( psi_(0,k), weights#k#0)} );
    toricReflexiveSheaf(W, X) 
    ) 

determinant ToricReflexiveSheaf := ToricReflexiveSheaf => opts -> E -> 
    exteriorPower(rank E, E)

symmetricPower (ZZ, ToricReflexiveSheaf) := ToricReflexiveSheaf => (p, E') -> (
    -- should think harder about the matroid labels in the general case
    if (ambient E')#1 =!= {1} then E := trim E' else E = E';
    X := variety E;
    r := rank E;
    (R,d) := ambient E;
    if p < 0 then return toricReflexiveSheaf X;
    if p === 0 then return toricReflexiveSheaf(0 * X_0);
    if p === 1 then return E';
    if rank E === 0 then return toricReflexiveSheaf(R, p*d, X);
    W := join for i to # rays X - 1 list (
	psi := map(R^1, R^(-values E#i), 
	    matrix{apply(keys E#i, g -> matrix{{g}})});
	psi = symmetricPower(p, psi);
	weights := degrees source psi;
	apply(binomial(r+p-1,p), k -> (psi_(0,k), weights#k#0)));
    toricReflexiveSheaf(W,X) )

subspace = method();
subspace (ZZ, ZZ, ToricReflexiveSheaf) := 
subspace (ZZ, InfiniteNumber, ToricReflexiveSheaf) := Ideal => (i,j,E) -> (
  subspaceGens := select(keys E#i, f -> E#i#f >= j);
  (R,d) := ambient E;  
  if #subspaceGens === 0 then ideal(0_R)
  else ideal subspaceGens)

-- local routine; the input L is expected to be a list of lists
boxProduct = L -> (
  m := #L;
  if m === 0 then {}
  else if m === 1 then apply(L#0, l -> {l})
  else if m === 2 then flatten table(L#0, L#1, (k,l) -> {k} | {l})
  else flatten table(first L, boxProduct drop(L,1),  (k,l) -> {k} | l))

subspacePoset = method();
subspacePoset ToricReflexiveSheaf := HashTable => E -> (  
    n := #rays variety E;
    (R, d) := ambient E;
    L := boxProduct apply(n, i -> sort unique values E#i);
    KK := coefficientRing R;
    BE := basis(d,R);
    V := keys set for l in L list gens intersect apply(n, 
	i -> image sub(
	    (coefficients(gens subspace(i,l#i,E), Monomials => BE))#1, KK
	    )
	);
    H := new MutableHashTable;
    for v in V do (
    	m := rank source v;
    	if H#?m then H#m = H#m | {v}
    	else if m === 0 then continue
    	else H#m = {v});
    hashTable apply(keys H, k -> k => H#k)
    )

groundSet = method();
groundSet ToricReflexiveSheaf := List => (cacheValue symbol groundSet) (E -> (  
  	if rank E === 0 then return {};
  	P := subspacePoset E;
  	if P#?1 then gSet := P#1 else (
    	    k := 2;
    	    while not P#?k do k = k + 1;
    	    gSet = {(first P#k)_{0}}
	    );
  	for k from 2 to rank E do (
    	    if P#?k then for V in P#k do (
      		spots := positions(gSet, f -> f % V == 0);
      		for j from 0 to rank source V - 1 do (
		    if spots === {} then M := 0 * V_{j}
		    else M = matrix{gSet_spots};
      		    g := V_{j} % M;
      		    if g != 0 then (
    	  		gSet = gSet | {g};
    	  		spots = spots | {#gSet-1} 
			) 
		    )
		)
	    );
  	(R, d) := ambient E;
  	BE := basis(d,R);
  	gSet = apply(gSet, g -> (basis(d,R) * g)_(0,0) );
  	reverse sort gSet)
    )

isArithmeticallyFree = method();
isArithmeticallyFree ToricReflexiveSheaf := Boolean => E -> 
    rank E == # groundSet E

isLocallyFree = method();
isLocallyFree ToricReflexiveSheaf := Boolean => E -> (
    if E == 0 then return true;
    X := variety E;
    all(max X, sigma -> (
      	    U := normalToricVariety( (rays X)_sigma, {sigma});
      	    EU := toricReflexiveSheaf( apply(sigma, i -> pairs E#i), U);
      	    isArithmeticallyFree EU
	    )
	)
    )

associatedCharacters = method();
associatedCharacters ToricReflexiveSheaf := List => (cacheValue symbol associatedCharacters) (E -> (
    	X := variety E;
    	if E == 0 then return apply(max X, sigma -> {});
    	for sigma in max X list (
    	    U := normalToricVariety((rays X)_sigma, {toList(0..#sigma-1)});
    	    EU := toricReflexiveSheaf(for i in sigma list pairs E#i,U);
    	    G := groundSet EU;
    	    C := transpose matrix for g in G list (
      		for i to # rays U - 1 list (
		    weights := reverse sort unique values EU#i;
		    j := 0;
		    while g % subspace(i, weights#j, EU) != 0 do j = j+1;
		    weights#j
		    )
		);
    	    if rank source C =!= rank EU then 
		error "expected sheaf to be locally free";
    	    if not isSmooth U then 
		error "expected sheaf on a smooth toric variety";
    	    entries transpose ((inverse matrix rays U) * C)
	    )
	)
    )

isGloballyGenerated = method()
isGloballyGenerated ToricReflexiveSheaf := Boolean => E -> (
    if E == 0 then return true;
    X := variety E;
    n := # rays X;  
    outerNormals := matrix rays X;
    assChar := associatedCharacters E;
    G := groundSet E;
    (R, d) := ambient E;
    (M, C) := coefficients( matrix {G}, Monomials => basis(d,R));
    divisors := apply(components cover E, L -> toricDivisor L);
    all(#assChar, i -> (
      	    potentialBases := for u in assChar#i list (
      		charVector := transpose matrix {u};
      		select(#G, j -> (
	    		coordinates := outerNormals*charVector - matrix vector divisors#j;
	    		all(n, k -> if member(k, (max X)#i) then 
				coordinates_(k,0) == 0
	      		    else coordinates_(k,0) <= 0
			    )
			)
		    )
		);
      	    any(boxProduct potentialBases, s -> det C_s =!= 0)
	    )
	)
    );  

separatesJets = method()
separatesJets ToricReflexiveSheaf := ZZ => E -> (
 testing := false;
 if testing then << "METHOD: separatesJets" << endl;
-- exclude trivial case
 if E == 0 then return 0;
-- setup groundSet, parliament of polytopes, character
 X := variety E;
 n := dim X;
 maxConeE := apply(max X, s -> (rays X)_s);
 gsE := groundSet E;
 parE := hashTable pack(2, mingle(gsE, apply(components cover E, D -> polytope toricDivisor D) ));
 onefaces := applyValues(parE, p -> facesAsPolyhedra(n-1,p));
-- note: change of sign!
 assCharE := hashTable pack(2, mingle(maxConeE, apply(associatedCharacters E, l->apply(l, u -> -transpose matrix {u})) ));
-- exclude another trivial case
 if min apply(values parE, dim) < 0 then (
  if testing then << "There are empty polytopes, so not even 0-jets are separated." << endl;
  return -1;
 );
-- working hypothesis: separates all jets (expectation will be lowered during computation)
 lJets := infinity;
-- if one of the polytopes is not full dimensional, lower expectation immediately
 if min apply(values parE, dim) < n then (
  if testing then << "There are polytopes which are not of full dimension, so at most separates 0-jets." << endl;
  lJets = 0;
 );
-- this list keeps track whether we hit all the polytopes
 indexedPolytopes := {};
 for sigma in maxConeE do ( 
  if testing then << "Test on the maximal cone:" << endl << sigma << endl;
  admPolys := new MutableHashTable from parE;
  for u in assCharE#sigma do (
   if testing then << "for the component u of the character:" << endl << u << endl;
   uVertexOfPolys := hashTable select(pairs admPolys, p -> isFace(convexHull u, p_1));
   if #uVertexOfPolys == 0 then (
    if testing then << "is not vertex of an available polytope." << endl;
    return -1;
   )
   else (
    if testing then << "is vertex of polytopes: " << endl << applyValues(uVertexOfPolys, vertices) << endl;
   );
   maxLengthRealisedAt := 0;
   if lJets > 0 then (
    uMaxEdgeLength := -infinity;
    for p in pairs uVertexOfPolys do (
     edgesContainU := select(onefaces#(p_0), f -> contains(f, u));
     if testing then << "Edges of polytope containing u: " << endl << apply(edgesContainU, vertices) << endl;
     coneAtU := dualCone coneFromVData fold(apply(edgesContainU, q -> vertices(q + convexHull (-u) )), (k,l) -> k|l);
     isSameCone := coneAtU == coneFromVData transpose matrix sigma;
     if testing then << "The dual cone at this vertex is:" << endl << rays coneAtU << endl
                     << "this is the same as the maximal cone: " << isSameCone << endl;
     localEdgeLength := if isSameCone then (
      edgesVert := apply(edgesContainU,vertices);
      edgesLengths := apply(edgesVert, e -> gcd flatten entries (e_1 - e_0)); 
      if testing then << "The edges have (lattice) lengths: " << edgesLengths << endl;
      min edgesLengths
     )
     else
      0;
     if localEdgeLength > uMaxEdgeLength then (
      uMaxEdgeLength = localEdgeLength;
      maxLengthRealisedAt = p;
     )
    );
    if testing then << "Maximal edge length for u is " << uMaxEdgeLength << " realised at: " << endl << vertices maxLengthRealisedAt_1 << endl;
    lJets = min(lJets, uMaxEdgeLength);
   )
   -- case lJets == 0;
   else (
    maxLengthRealisedAt = first pairs uVertexOfPolys;
   );
   remove(admPolys, maxLengthRealisedAt_0);
   indexedPolytopes = append(indexedPolytopes, maxLengthRealisedAt_0);
  );
 );
 if set indexedPolytopes === set keys parE then
  lift(lJets,ZZ)
 else
  -1
)


cover ToricReflexiveSheaf := ToricReflexiveSheaf => (cacheValue symbol cover) (E -> (
	if E == 0 then return E;
    	X := variety E;
    	n := # rays X;
    	(RE, dE) := ambient E;    
    	gSet := groundSet E;
    	if gSet === {} then (
      	    Z := toricReflexiveSheaf(RE,dE,X);
      	    E.cache.generators = map(E, Z, 0);
      	    return Z);
    	divisors := for g in gSet list (
      	    sum for i from 0 to n-1 list (
    		weights := reverse sort unique values E#i;
    		j := 0;
    		while g % subspace(i,weights#j, E) != 0 do j = j + 1;
    		weights#j *X_i
		)
	    );
    	F := directSum apply(divisors, D -> toricReflexiveSheaf D);
    	(M, C) := coefficients(matrix{gSet}, Monomials => basis(dE, RE));
    	KK := coefficientRing RE;
    	E.cache.generators = map(E, F, entries sub(C,KK));
    	F
	)
    )

generators ToricReflexiveSheaf := ToricReflexiveSheafMap => opts -> E -> (
  if not E.cache.?generators then cover E;
  E.cache.generators)


------------------------------------------------------------------------------
-- Basic features of the toric reflexive sheaf map datatype
------------------------------------------------------------------------------
ToricReflexiveSheafMap = new Type of HashTable
ToricReflexiveSheafMap.synonym = "equivariant map of toric reflexive sheaves"
source ToricReflexiveSheafMap := ToricReflexiveSheaf => f -> f.source
target ToricReflexiveSheafMap := ToricReflexiveSheaf =>  f -> f.target
matrix ToricReflexiveSheafMap := Matrix => opts -> f -> f.matrix
variety ToricReflexiveSheafMap := NormalToricVariety => f -> variety source f

net ToricReflexiveSheafMap := f -> net matrix f
ToricReflexiveSheafMap#{Standard,AfterPrint} = 
ToricReflexiveSheafMap#{Standard,AfterNoPrint} = f -> (
    << endl;				  -- double space
    << concatenate(interpreterDepth:"o") << lineNumber << " : ToricReflexiveSheafMap ";
    << net target f << " <--- " << net source f << endl;)

map(ToricReflexiveSheaf, ToricReflexiveSheaf, List) := ToricReflexiveSheafMap => opts -> (E, F, L) -> (
    (RF, dF) := ambient F;  
    (RE, dE) := ambient E;  
    rF := rank source basis(dF, RF);
    rE := rank source basis(dE, RE);    
    KK := coefficientRing RF;
    phi := map(KK^rE, KK^rF, L);
    new ToricReflexiveSheafMap from {
      	symbol source => F,
      	symbol target => E,
      	symbol matrix => phi,
      	symbol cache => new CacheTable}
    )
map(ToricReflexiveSheaf, ToricReflexiveSheaf, Matrix) := ToricReflexiveSheafMap => 
    opts -> (E, F, f) -> map(E, F, entries f)
map(ToricReflexiveSheaf, ToricReflexiveSheaf, Function) := ToricReflexiveSheafMap => opts -> (E, F, f) -> (
    (RF, dF) := ambient F;  
    (RE, dE) := ambient E;  
    rF := rank source basis(dF, RF);
    rE := rank source basis(dE, RE);    
    map(E, F, table(rE, rF, f))
    )
map(ToricReflexiveSheaf, ToricReflexiveSheaf, Number) :=
map(ToricReflexiveSheaf, ToricReflexiveSheaf, RingElement) := ToricReflexiveSheafMap => opts -> (E, F, g) -> (
    (RF, dF) := ambient F;  
    (RE, dE) := ambient E;  
    rF := rank source basis(dF, RF);
    rE := rank source basis(dE, RE);  
    KK := coefficientRing RF; 
    if g == 0 then map(E, F, table(rE, rF, (i,j) -> 0_KK))
    else if rF === rE then map(E, F, g * id_(KK^rE))
    else error "expected 0, or ambient vector spaces of the same dimension"
    )
  
ToricReflexiveSheafMap * ToricReflexiveSheafMap := ToricReflexiveSheafMap => (g,f) -> (
    if source g =!= target f then error "maps are not composable";
    map(target g, source f, matrix g * matrix f)
    )
  
ToricReflexiveSheafMap == ZZ := Boolean => (f,m) -> matrix f == m
ZZ == ToricReflexiveSheafMap := Boolean => (m,f) -> f == m
    
isWellDefined ToricReflexiveSheafMap := Boolean => f -> (
    -- check keys
    K := keys f;
    expectedKeys := set{ symbol source, symbol target, symbol matrix, symbol cache};
    if set K =!= expectedKeys then (
    	if debugLevel > 0 then (
      	    added := toList(K - expectedKeys);
      	    missing := toList(expectedKeys - K);
      	    if #added > 0 then 
            << "-- unexpected key(s): " << toString added << endl;
      	    if #missing > 0 then 
            << "-- missing keys(): " << toString missing << endl);
    	return false
	);
    -- check types
    if not instance(f.source, ToricReflexiveSheaf) then (
    	if debugLevel > 0 then 
      	    << "-- expect 'f.source' to be a ToricReflexiveSheaf" << endl;
    	return false
	);
    if not instance(f.target, ToricReflexiveSheaf) then (
     	if debugLevel > 0 then 
	    << "-- expect 'f.target' to be a ToricReflexiveSheaf" << endl;
     	return false
	);
    if not instance(f.matrix, Matrix) then (
    	if debugLevel > 0 then 
      	    << "-- expect 'f.matrix' to be a Matrix " << endl;
    	return false
	); 
    F := source f;
    E := target f;
    (RF, dF) := ambient F;  
    (RE, dE) := ambient E;
    KK := coefficientRing RF; 
    if KK =!= coefficientRing RE then (
    	if debugLevel > 0 then 
      	    << "-- expected the coefficient rings of the ambient rings to be equal" << endl;
    	return false
	);      
    if not all(flatten entries matrix f, r -> instance(r, KK)) then (
    	if debugLevel > 0 then 
      	    << "-- expected entries belong to coefficient ring of the ambient ring" << endl;
    	return false
	);
    if not instance(E.cache, CacheTable) then (
    	if debugLevel > 0 then 
      	    << "-- expected `E.cache' to be a CacheTable" << endl;
    	return false
	);    
    -- check underlying varieties
    if variety F =!= variety E then (
    	if debugLevel > 0 then 
      	    << "-- expected sheaves over the same variety" << endl;
    	return false);
    -- check compatiblity of flags
    n := # rays variety F;
    BF := basis(dF, RF);  
    BE := basis(dE, RE);
    for i from 0 to n-1 do (
    	weights := sort unique (values F#i | values E#i);
    	for w in weights do (
      	    IF := subspace(i,w,F);
      	    (M, C) := coefficients(mingens IF, Monomials => BF);
      	    IE := ideal (BE * (matrix f) * sub(C, KK));
      	    if not isSubset(IE, subspace(i,w,E)) then (
    		if debugLevel > 0 then (
      	  	    << "-- expected " << toString w << "-th part of the filtrations";
      	  	    << " on the " << toString i << "-th ray to be compatible" << endl);
    		return false
		)
	    )
	);    
    true)

kernel ToricReflexiveSheafMap := ToricReflexiveSheaf => opts -> (
    cacheValue symbol kernel) (f -> (
    	E := source f;
    	(R, d) := ambient E;
    	I := ideal( basis(d, R) * gens ker matrix f );
    	X := variety E;    
    	if I == 0 then return toricReflexiveSheaf X;
    	n := # rays X;
    	W := for i from 0 to n-1 list (
      	    basisSet := {};
      	    weights := reverse sort unique values E#i;
      	    for w in weights do (
		V := ideal select( (ideal mingens intersect(I, subspace(i,w,E)))_*, 
	  	    g -> degree g === d);
		if V === ideal() then continue else (
	  	    spots := positions(basisSet, f -> f#0 % V == 0);
	  	    C := unique flatten entries (gens V % (matrix{apply(basisSet_spots, 
		  		    p -> p#0)} ** ring V));
	  	    C = for c in C list if c == 0 then continue else (c,w); 
	  	    if C == {} then continue 
	  	    else basisSet = basisSet | C));
      	    basisSet
	    );
    	toricReflexiveSheaf(W, X)
	)
    )

mukaiLazarsfeldBundle = method();
mukaiLazarsfeldBundle ToricDivisor := ToricReflexiveSheaf => D -> (
    coeffs := entries vector D;
    X := variety D;
    n := # rays X;
    S := ring X;
    monos := first entries basis(sum apply(n, i -> coeffs#i * degree S_i), S);
    E := directSum apply(monos, g -> (  
      	    coeffs = first exponents g;
      	    toricReflexiveSheaf sum apply(n, i -> - coeffs#i * X_i)
	    )
	);
    f := map(toricReflexiveSheaf D, E, matrix{toList(rank E : 1_QQ)});
    ker f)


-*
Complex = new Type of List;
Complex.synonym = "chain complex of toric reflexive sheaves"  
net Complex := C ->  (
  D := netList( transpose {{ -1, rank target C#0}}, 
    Boxes => false, Alignment => Center, VerticalSpace => 1);
  for i to #C - 1 do (
    D = D | "  " | (netList ( {i, rank source C#i} | apply(components source C#i, L -> toricDivisor L), 
  	Boxes => false, Alignment => Center, VerticalSpace => 1)));
  D);
Complex#{Standard,AfterPrint} = 
Complex#{Standard,AfterNoPrint} = C -> (
  << endl;
  << concatenate(interpreterDepth:"o") << lineNumber << " : Chain complex of toric reflexive sheaves over ";
  << variety C#0 << endl;)  

resolution ToricReflexiveSheaf := Complex => opts -> E -> (
  (R,d) := ambient E;
  nabla := map(toricReflexiveSheaf(R,d, variety E), E, 0);
  maps := {};
  while true do (
    F := source nabla;
    K := ker nabla;
    if K == 0 then break;
    G := cover K;
    nabla = map(F,K,1) * gens K;
    maps = maps | {nabla});
  new Complex from maps);
*-
