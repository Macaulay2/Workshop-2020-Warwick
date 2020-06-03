export {
    "Subring",
    "subring",
    "liftedPresentation",
    "presentationRing"
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
    new Subring from {
    	"AmbientRing" => R,
    	"Generators" => M,
	cache => new CacheTable from {}
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
    if not A.cache#?"PresentationRing" then (
	B := ambient A;
	k := coefficientRing B;
	e := symbol e;
	nA := numgens A;
	A.cache#"PresentationRing" = k[apply(nA, i-> e_i)];
	);
    A.cache#"PresentationRing"
    )


-- lifted presentation using the cacheValue function
liftedPresentation = method()
liftedPresentation Subring := (cacheValue "LiftedPresentation")(A -> (
    B := ambient A;
    P := presentationRing A;
    G := gens A;
    k := coefficientRing B;
    (nB, nA) := (numgens B, numgens A);
    -- introduce nA "tag variables" w/ monomial order that eliminates non-tag variables
    e := symbol e;
    -- C := k[gens B | apply(nA, i -> e_i), MonomialOrder => append(getMO B, Eliminate nB)];
    C := k[gens B | gens P, MonomialOrder => append(getMO B, Eliminate nB)];
    B2C := map(C,B,(vars C)_{0..nB-1});
    ideal(B2C G - (vars C)_{nB..numgens C-1})
    ))

-- leaving the old liftedPresentation for comparison
-- delete after cacheValue discussion (20200603)
-*
liftedPresentation Subring := A -> (
    if not A.cache#?"LiftedPresentation" then (
    	B := ambient A;
	P := presentationRing A;
	G := gens A;
    	k := coefficientRing B;
	(nB, nA) := (numgens B, numgens A);
	-- introduce nA "tag variables" w/ monomial order that eliminates non-tag variables
        -- C := k[gens B | gens P, MonomialOrder => append(getMO B, Eliminate nB)];
	C := k[gens B | gens P, MonomialOrder => {Eliminate nB}];
	B2C := map(C,B,(vars C)_{0..nB-1});
    	A.cache#"LiftedPresentation" = ideal(B2C G - (vars C)_{nB..numgens C-1});
	);
    A.cache#"LiftedPresentation"
    )
 *-

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
-- these need to be implemented

-- output: r in ambient of A such that f = a + r w/ a in A, r "minimal"
RingElement % Subring := (f, A) -> (
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

-- probably needs to change!
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

