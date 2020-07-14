export {
    "Subring",
    "subring",
    "liftedPresentation",
    "liftedPresentationRing",
    "presentationRing",
    "getWeight",
    "setWeight"
    }


-- Organization
-- 1) Subring type (and associated methods)
-- 2) experimental valuation type

-- 1) Subring type (and associated basic methods)

-- todo: eventually, we might want Subring to inherit from Ring
-- but inheriting from "Ring" is not straightforward, so HashTable for now
Subring = new Type of HashTable

-- constructor
subring = method()
subring Matrix := M -> (
    R := ring M;
    cTable := new CacheTable from{
	SubalgComputations => new MutableHashTable from {},
	SagbiGens => matrix(R, {{}}),
	SagbiDegrees => {},
	SagbiDone => false
	}; 
    new Subring from {
    	"AmbientRing" => R,
    	"Generators" => M,
	cache => cTable
	}    
    )
subring List := L -> subring matrix{L}

-- overloaded "getters"
gens Subring := o -> A -> A#"Generators"
numgens Subring := A -> numcols gens A
ambient Subring := A -> A#"AmbientRing"
net Subring := A -> "subring of " | toString(ambient A)

presentationRing = method()
presentationRing Subring := A -> (
    B := ambient A;
    k := coefficientRing B;
    e := symbol e;
    nA := numgens A;
    return presentationRing(A, k[apply(nA, i -> e_i)]);
    )

presentationRing (Subring, Ring) := (A, newPresRing) -> (
    if not A.cache#?"PresentationRing" then (
	B := ambient A;
	nA := numgens A;
	assert(nA == numgens newPresRing);
	if A.cache#?"AmbientWeight" then (
	    D := A.cache#"AmbientWeight";
	    weightRing := newRing(B, MonomialOrder => { Weights => D});
	    inducedWeights := for p in flatten entries gens A list (
	    	p' := sub(p, weightRing);
	    	E := (exponents leadTerm p')#0;
	    	sum apply(E, D, (i,j) -> i*j)
	    	);
	    A.cache#"PresentationRing" = newRing(newPresRing, MonomialOrder => {Weights => inducedWeights});
	    ) else (
	    A.cache#"PresentationRing" = newPresRing;
	    );
	);
    A.cache#"PresentationRing"
    )

presentationRing (List, Matrix) := (pVars, M) -> (
    assert (#pVars === numcols M);
    newRing (ring M, Variables => pVars)
    )

presentationRing Matrix := M -> (
    n := numcols M;
    pVars := toList vars (0..n-1);
    presentationRing (pVars, M)
    )

-*
    Input:  - matrix with n columns over polynomial ring R
            - a list of n new variables
            - optional: a prescribed term order
    Output: a polynomial ring with with variables coming from R and the list
    *-

liftedPresentationRing = method (
    TypicalValue => PolynomialRing,
    Options => {MonomialOrder => null})
liftedPresentationRing (List, Matrix) := PolynomialRing => o -> (pVars, M) -> (
    assert (#pVars === numcols M);
    X := (ring M)_*;
    newVars := X | pVars;
    liftedPR := if o.MonomialOrder === null then
    newRing (ring M, Variables => newVars) else
    newRing (ring M, Variables => newVars, MonomialOrder => o.MonomialOrder);
    liftedPR
    )
liftedPresentationRing Matrix := PolynomialRing => o -> M -> (
    pVars := (presentationRing M)_*;
    liftedPresentationRing (pVars, M, MonomialOrder => o.MonomialOrder))


-- computes the presentation of the subring in the presentation ring
liftedPresentation = method()
liftedPresentation Subring := (cacheValue "LiftedPresentation")(A -> (
    B := ambient A;
    P := presentationRing A;
    G := gens A;
    k := coefficientRing B;
    (nB, nA) := (numgens B, numgens A);
    -- introduce nA "tag variables" w/ monomial order that eliminates non-tag variables

    -- e := symbol e;
    -- C := k[gens B | apply(nA, i -> e_i), MonomialOrder => append(getMO B, Eliminate nB)];
    C := k[gens B | gens P, MonomialOrder => append(getMonomialOrder B, Eliminate nB)];
    B2C := map(C,B,(vars C)_{0..nB-1});
    ideal(B2C G - (vars C)_{nB..numgens C-1})
    ))


-- computes relations of presentation using gb
presentation Subring := A -> (
    if not A.cache#?"LiftedGB" then (
	A.cache#"LiftedGB" = gb liftedPresentation A;
	);
    presentationGens := selectInSubring(1, gens A.cache#"LiftedGB");
    P := presentationRing A;
    sub(presentationGens, P)
    )


-- quotient ring given by a presentation
ring Subring := A -> (
    I := ideal presentation A;
    P := ring I;
    P/I
    )
options Subring := A -> A.cache#"Options"

-- output: r in ambient of A such that f = a + r w/ a in A, r "minimal"
-- See proposition 3.6.1 in "Computational Commutative Algebra" book 1 by Kreuzer and Robbiano.
RingElement % Subring := (f, A) -> (
    
    if hasComputedSagbi A then(
	return subduction(A, f);
	);
    
    pA := presentationRing A;
    tagVars := take(gens pA, numgens A);
    tagSubs := flatten entries sub(gens A, ring liftedPresentation A);
    subTable := apply(tagVars, tagSubs, (ei, fi) -> ei => fi);
    nonRemainder := sub(sub(sub(f // A, pA), subTable), ambient A);
    use ring f;
    f - nonRemainder
    )

-- output: for A=k[g1,..,gk], p(e1,...,ek) st f = p(g1,..,gk) + r
RingElement // Subring := (f, A) -> (
    I := liftedPresentation A;
    T := ring I;
    sub(f, T) % I
    )

-- probably needs to change! (why? :P)
member (RingElement, Subring) := (f, A) -> (
    r := f%A;
    R := ambient A;
    r == 0_R
    )

-- 2) experimental valuation type

Valuation = new Type of HashTable
-- what should a valuation need to know?
source Valuation := v -> v#source
target Valuation := v -> v#target
net Valuation := v -> "valuation: " | toString target v | " <-- " | toString source v
Valuation RingElement := (v, f) -> (
    assert(v#?evaluate and ring f === source v);
    assert instance(v#evaluate, Function);
    v#evaluate f
    )

MonomialValuation = new Type of Valuation
monomialValuation = method()
monomialValuation Ring := R -> new MonomialValuation from {
    source => R,
    target => ZZ^(numgens R),
    evaluate => (f -> matrix exponents leadTerm f)
    }
leadTerm (MonomialValuation, RingElement) := (v, f) -> (
    assert(ring f === source v);
    leadTerm f
    )


-- set a weight on the ambient ring that induces a weight on the presentation ring
setWeight = method()
setWeight (Subring, List) := (A, W) -> (
    B := ambient A;
    assert(numgens B == length W);
    if not A.cache#?"AmbientWeight" then (
	A.cache#"AmbientWeight" = W;
	) else (
	print "Weight has already been set"
	);
    )

getWeight = method()
getWeight Subring := A -> (
    if not A.cache#?"AmbientWeight" then (
	return "None set"
	);
    A.cache#"AmbientWeight"
    )


-*
R=QQ[x,y]
f = x+y^2
v = monomialValuation R
source v
target v
v f
leadTerm(v,f) < y
*-

end--

