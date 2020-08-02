export {
    "Subring",
    "subring",
    "PresRing",
    "makePresRing",
    "getWeight",
    "setWeight",
    "presentationRing"
    }

-- todo: eventually, we might want Subring to inherit from Ring
-- but inheriting from "Ring" is not straightforward, so HashTable for now
Subring = new Type of HashTable


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
	-- The PresRing of a Subring instance is immutable because the generators are immutable.
	"PresRing" => makePresRing(R, M),
	"isSagbi" => false,
	"isPartialSagbi" => false,
	"partialDegree" => 0,
	cache => cTable
	}    
    )
subring List := L -> subring matrix{L}

-- overloaded "getters"
gens Subring := o -> A -> A#"Generators"
numgens Subring := A -> numcols gens A
ambient Subring := A -> A#"AmbientRing"
net Subring := A -> "subring of " | toString(ambient A)

-- TODO: Write a Subring equality operator which calculates the mathematical equality of subalgebras 
-- TODO: Write good tests for % and //.

-- This type is compatible with internal maps that are generated in the Sagbi algorithm.
-- Originally, this was stored directly in the cache of an instance of Subring. 
-- The problem with that solution is there is a need to use these maps outside of the Sagbi algorithm computations.
PresRing = new Type of HashTable

net PresRing := pres -> (    
    tense := pres#"TensorRing";
    A := numcols vars tense;
    B := numcols selectInSubring(1, vars tense);
    "PresRing instance ("|toString(B)|" generators in "|toString(A-B)|" variables)"
    )

-- gensR are elements of R generating some subalgebra.
-- R is a polynomial ring.
makePresRing = method(TypicalValue => PresRing)  
makePresRing(Ring, Matrix) := (R, gensR) ->( 
    if(R =!= ring(gensR)) then(
	error "The generators of the subalgebra must be in the ring R.";
	);
    makePresRing(R, first entries gensR)
    );
 
makePresRing(Ring, List) := (R, gensR) ->( 
    
    if(ring(matrix({gensR})) =!= R) then(
	error "The generators of the subalgebra must be in the ring R.";
	);
    
    ambR := R;
    nBaseGens := numgens ambR;
    nSubalgGens := length gensR;
    
    -- Create a ring with combined generators of base and subalgebra.  
    MonoidAmbient := monoid ambR;
    CoeffField := coefficientRing ambR;
    
    -- Construct the monoid of a ring with variables corresponding to generators of the ambient ring and the subalgebra.
    -- Has an elimination order that eliminates the generators of the ambient ring.
    -- The degrees of generators are set so that the SyzygyIdeal is homogeneous.
    newOrder := prepend(Eliminate nBaseGens, MonoidAmbient.Options.MonomialOrder);


    -- TODO: It would be nice if the variables had better names    
    --sym1 := getSymbol "x";
    --sym2 := getSymbol "y";
    NewVariables := monoid[
        Variables=>nBaseGens+nSubalgGens,
	--Variables => (sym1_1..sym1_nBaseGens)|(sym2_1..sym2_nSubalgGens),
	Degrees=>join(degrees source vars ambR, degrees source matrix({gensR})),
        MonomialOrder => newOrder];
        
    TensorRing := CoeffField NewVariables;	    
    
    presVars := monoid[
        Variables=>nSubalgGens,
        Degrees=> degrees source matrix({gensR})
        ];
    
    (pres, InclusionPres) := selectVariables(toList (nBaseGens..(nBaseGens+nSubalgGens-1)),TensorRing);
    
    ProjectionPres := map(pres, TensorRing,
        matrix {toList (nBaseGens:(0_pres))} | (vars pres));
    
    ProjectionInclusion := map(TensorRing, TensorRing,
        (matrix {toList(nBaseGens:0_(TensorRing))}) |
	(vars TensorRing)_{nBaseGens .. nBaseGens+nSubalgGens-1});
    
    ProjectionBase := map(ambR, TensorRing,
        (vars ambR) | matrix {toList(nSubalgGens:0_(ambR))});
    
    InclusionBase := map(TensorRing, ambR,
        (vars TensorRing)_{0..nBaseGens-1});
    
    SubstitutionPres := map(ambR, pres, matrix({gensR}));

    Substitution := map(TensorRing, TensorRing,
        (vars TensorRing)_{0..nBaseGens-1} | InclusionBase(matrix({gensR})));
    
    SyzygyIdeal := ideal(
        (vars TensorRing)_{nBaseGens..nBaseGens+nSubalgGens-1}-
	InclusionBase(leadTerm matrix({gensR})));
    
    submap := Substitution;
    genVars := (vars TensorRing)_{numgens ambient R..numgens TensorRing-1};
    liftedPres := ideal(submap(genVars) - genVars);
    FullSub := ProjectionBase*Substitution;
     
    ht := new HashTable from {
	"TensorRing" => TensorRing,
	"ProjectionInclusion" => ProjectionInclusion,
	"ProjectionBase" => ProjectionBase,
	"InclusionBase" => InclusionBase,
	"Substitution" => Substitution,
	"FullSub" => FullSub,
	"PresentRing" => pres,
	"ProjectionPres"=> ProjectionPres,
	"InclusionPres" => InclusionPres,
	"SubstitutionPres" => SubstitutionPres,
	"SyzygyIdeal" => SyzygyIdeal,
	"LiftedPres" => liftedPres
	};
    
    new PresRing from ht
    );

makePresRing(Subring) := subR -> (
    subR#"PresRing"
    );

-- computes relations of presentation using gb
presentation Subring := A -> (
    if not A.cache#?"LiftedPresGB" then (
	pres := A#"PresRing";
	A.cache#"LiftedPresGB" = gb (pres#"LiftedPres");
	);
    selectInSubring(1, gens A.cache#"LiftedGB");
    );

-- quotient ring given by a presentation
ring Subring := A -> (
    I := ideal presentation A;
    P := ring I;
    P/I
    );
-- f % Subring is never going to be an element of the subalgebra, hence the ouput
-- is in the lower variables of TensorRing.
-- input: f in ambient A or TensorRing of A. 
-- output: r in TensorRing of A such that f = a + r w/ a in A, r "minimal"
RingElement % Subring := (f, A) -> (
    pres := A#"PresRing";
    if ring f === ambient A then(
	f = (pres#"InclusionBase")(f);
	) else if ring f =!= pres#"TensorRing" then(
	error "The RingElement f must be in either TensorRing or ambient A.";
	);
    ans := (subduction(A, f));
    ans    
    );

-- f // Subring is always going to be inside of the subalgebra, hence the output
-- should be in the upper variables of TensorRing.
-- input: f in ambient A or TensorRing of A. 
-- output: a in TensorRing of A such that f = a + r w/ a in A, r "minimal."
RingElement // Subring := (f, A) -> (
    pres := A#"PresRing";
    tense := pres#"TensorRing";
    if ring f === ambient A then(
	f = (pres#"InclusionBase")(f);
	) else if ring f =!= tense then(
	error "The RingElement f must be in either the TensorRing or ambient ring of A.";
	);
    result := f - (f % A);
    I := pres#"LiftedPres";
    result % I
    );

Matrix % Subring := (M, A) -> (
    pres := A#"PresRing";
    ents := for i from 0 to numrows M - 1 list(
	for j from 0 to numcols M - 1 list(M_(i,j) % A)
	);
    matrix(pres#"TensorRing", ents)
    );
Matrix // Subring := (M, A) -> (
    pres := A#"PresRing";
    ents := for i from 0 to numrows M - 1 list(
	for j from 0 to numcols M - 1 list(M_(i,j) // A)
	);
    matrix(pres#"TensorRing", ents)
    );

-- Perhaps it is a bug that this will sometimes throw an error when it should return false.
member (RingElement, Subring) := (f, A) -> (
    r := f%A;
    r == 0_(A#"PresRing"#"TensorRing")
    );

-----------------------------------------------------------------
-- experimental valuation type
-----------------------------------------------------------------

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

