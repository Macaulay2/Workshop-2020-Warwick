export {
    "Subring",
    "subring",
    "liftedPresentation"
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

liftedPresentation = method()
liftedPresentation Subring := A -> (    
    if not A.cache#?"LiftedPresentation" then (
    	B := ambient A;
	G := gens A;
    	k := coefficientRing B;
	(nB, nA) := (numgens B, numgens A);
	-- introduce nA "tag variables" w/ monomial order that eliminates non-tag variables
	e := symbol e;	       
    	C := k[gens B | apply(nA, i -> e_i), MonomialOrder => append(getMO B, Eliminate nB)];
	B2C := map(C,B,(vars C)_{0..nB-1});
    	A.cache#"LiftedPresentation" = ideal(B2C G - (vars C)_{nB..numgens C-1});
	);
    A.cache"LiftedPresentation"
    )

-- computes an ideal of relations
presentation Subring := A -> selectInSubring(1, liftedPresentation A)

-- quotient ring given by a presentation
ring Subring := A -> (
    I := ideal presentation A;
    R := ring I;
    R/I
    )
options Subring := A -> A.cache#"Options"
-- these need to be implemented

-- output: r in ambient of A such that f = a + r w/ a in A, r "minimal"  
RingElement % Subring := (f, A) -> (
    ret := f;
    assert(ring ret == ambient A);
    ret
    )

-- output: for A=k[g1,..,gk], p(e1,...,ek) st f = p(g1,..,gk) + r
RingElement // Subring := (f, A) -> (
    I := presentation A;
    f % I
    )

-- probably needs to change!
member (RingElement, Subring) := (f, A) -> (
    r := f%A;
    R := ambient A;
    r == 0_R
    )