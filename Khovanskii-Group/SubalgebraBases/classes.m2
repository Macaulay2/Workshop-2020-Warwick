export {
    "Subring",
    "subring"
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
-- computes an ideal of relations
presentation (Matrix, Subring) := (G, A) -> (
    if not A.cache.?presentation then (
    	B := ambient A;
    	k := coefficientRing B;
	(nB, nA) := (numgens B, numgens A);
	-- introduce nA "tag variables" w/ monomial order that eliminates non-tag variables
    	Cmonoid := monoid [Variables => nB + nA,
	                   MonomialOrder => (options B).MonomialOrder];
    	C := k Cmonoid;
	B2C := map(C,B,(vars C)_{0..nB-1});
    	I := ideal(B2C G - (vars C)_{nB..numgens C-1});
    	A.cache.presentation = selectInSubring(1, gens gb I);
	);
    A.cache.presentation
    )
-- quotient ring given by a presentation
ring Subring := A -> (
    I := ideal presentation A;
    R := ring I;
    R/I
    )
options Subring := A -> A.cache#"Options"
-- these need to be implemented
RingElement % Subring := (f, A) -> f
RingElement // Subring := (f, A) -> f
member (RingElement, Subring) := (f, A) -> (
    r := f%A;
    R := ambient A;
    r == 0_R
    )