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


liftedPresentation = method()
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

end--

TEST ///
R = QQ[x1, x2, x3];
S = QQ[e1, e2, e3, y];
f = map(R, S, {x1 + x2 + x3, x1*x2 + x1*x3 + x2*x3, x1*x2*x3,
(x1 - x2)*(x1 - x3)*(x2 - x3)});
A = subring matrix f;
T = ring liftedPresentation A;
assert (presentation A == matrix {{T_3^2*T_4^2-4*T_3^3*T_5-4*T_4^3+18*T_3*T_4*T_5-27*T_5^2-T_6^2}})
///