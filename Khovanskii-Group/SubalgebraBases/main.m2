debug Core -- gets rid of "raw" error during installation. probably a better way...

export {
    "subalgebraBasis",
    "subduction",
    "sagbi" => "subalgebraBasis", 
    "PrintLevel",
    "SagbiDegrees",
    "SubalgComputations",
    "SagbiGens",
    "SagbiDone"
    }

-- A wrapper around rawSubduction. Intended for internal use only.
   -- M should be a 1-row matrix with entries in the tensor ring of subR.
subduct = method(TypicalValue => Matrix)
subduct(Subring, Matrix) := (subR, M) -> (   
    subalgComp := subR.cache.SubalgComputations;
    ambR := ambient subR;
    
    if entries M === {{}} then (
	return M;
	);
    
    -- Substitution encodes the set used as the subalgebra basis.
    F := subalgComp#"Substitution";
    J := subalgComp#"sagbiGB";
    numblocks := rawMonoidNumberOfBlocks raw monoid ambR;
    rawSubduction(numblocks, raw M, raw F, raw J)
    );

-- If S is a SAGBI basis, the result is zero IFF f is an element of the subalgebra generated by S. 
   -- S is a set of RingElements (as in Robbiano and Sweedler.)
   -- f is a RingElement to perform subduction on, relative to S.
subduction = method(TypicalValue => RingElement)
subduction(Matrix, RingElement) := (S, f) -> (
    subR := subring S;
    
    -- This writes data we need to subR.cache.SubalgComputations.
    appendToBasis(subR, S);
    subalgComp := subR.cache.SubalgComputations;  
    F := subalgComp#"Substitution";
    J := gb(subalgComp#"SyzygyIdeal");
    numblocks := rawMonoidNumberOfBlocks raw monoid ambient subR;
    fMat := matrix({{subalgComp#"InclusionBase"(f)}});    
    result := rawSubduction(numblocks, raw fMat, raw F, raw J);
    result = promote(result_(0,0), source (subalgComp#"ProjectionBase"));
    result = subalgComp#"ProjectionBase"(result);
    
    result
    );

-- Main function for computing a subalgebraBasis.
-- For now, subduction is performed in the engine (in older version, Strategy toggles between engine and top-level implementation)
subalgebraBasis = method(Options => {
    Strategy => null,
    Limit => 100,
    PrintLevel => 0})

subalgebraBasis List := o -> L -> (
    subalgebraBasis(o, subring L)
    );

subalgebraBasis Subring := o -> R -> (
    -- baseRing is the ring of the input matrix
    -- semiRing is the free semigroup ring formed from the SAGBI generators
    -- tensorRing is the ring formed from the tensor of the base ring and semigroup ring
    -- Pending is a list of lists, sorting elements of the algebra by degree.

    -- ProjectionInclusion is the projection from the tensor ring to the semiRing by sending the base ring generators to 0.
    -- projectionToBase is the projection from the tensor ring to the baseRing by sending the SAGBI generators to 0.
    -- inclusionOfBase is the inclusion of the baseRing into the tensorRing
        
    R.cache.SubalgComputations = new MutableHashTable;
    subalgComp := R.cache.SubalgComputations;

    R.cache.SagbiDegrees = {};
    
    subalgComp#"TensorRing" = null;  -- RS
    subalgComp#"SyzygyIdeal" = null; -- J
    subalgComp#"sagbiGB" = null;
    
    subalgComp#"ProjectionInclusion" = null; -- RStoS
    subalgComp#"ProjectionBase" = null;      -- RStoR
    subalgComp#"InclusionBase" = null;       -- RtoRS
    subalgComp#"Substitution" = null;        -- Gmap
    
    currDegree := null;     -- d
    nLoops := null;         -- nloops
    R.cache.SagbiDone = false;
    syzygyPairs := null;
    newElems := null;
    
    subalgComp#"Pending" = new MutableList from toList(o.Limit+1:{});
    R.cache.SagbiGens = matrix(ambient R,{{}});

    -- Get the maximum degree of the generators. This is used as a stopping condition.
    maxGensDeg := (max degrees source gens R)_0;

    -- Only look at generators below degree limit.  Add those generators to the SubalgebraGenerators
    reducedGens := compress submatBelowDegree(gens R, o.Limit+1);
    insertPending(R, reducedGens, o.Limit);
    -- Remove elements of coefficient ring
    (subalgComp#"Pending")#0 = {};
    processPending(R, o.Limit);
    currDegree = subalgComp#"CurrentLowest" + 1;
        
    while currDegree <= o.Limit and not R.cache.SagbiDone do (  	
        -- Construct a Groebner basis to eliminiate the base elements generators from the SyzygyIdeal.
	-- SyzygyIdeal is an ideal consisting of (variables repesenting subalgebra generators) minus (their leading term). 
	subalgComp#"sagbiGB" = gb(subalgComp#"SyzygyIdeal", DegreeLimit => currDegree);
	
	-- This will select the entries of sagbiGB that do not involve any of (leadTerms subalgComp#"SyzygyIdeal") and also
	-- have degree equal to currDegree. So, they will be exclusively in the higher block of variables of TensorRing.
	zeroGens := submatByDegree(mingens ideal selectInSubring(1, gens (subalgComp#"sagbiGB")), currDegree);
	    		
	-- Plug the generators into the degree currDegree polynomials that eliminate the lead terms (I.e. zeroGens.) 
	-- This changes it from a polynomial in the generators to a polynomial in the variables of the ambient ring.
       	syzygyPairs = subalgComp#"Substitution"(zeroGens);

	-- Have we previously found any syzygies of degree currDegree?
        if subalgComp#"Pending"#currDegree != {} then (
            syzygyPairs = syzygyPairs | subalgComp#"InclusionBase"(matrix{subalgComp#"Pending"#currDegree});
            subalgComp#"Pending"#currDegree = {};
            );
	
       	subd := subduct(R, syzygyPairs);
	
       	if entries subd != {{}} then (
	    -- converts back to the variables of the ambient ring.
	    subducted := (subalgComp#"ProjectionBase")(map(subalgComp#"TensorRing",subd));
	    newElems = compress subducted;
            ) else (
	    newElems = subd;
	    );
	
	if numcols newElems > 0 then (
	    -- Put newElems in pending and update subalgComp. 
            insertPending(R, newElems, o.Limit);
    	    processPending(R, o.Limit);
	    currDegree = subalgComp#"CurrentLowest";
            ) else (
	    -- "rawStatus1 raw (subalgComp#"sagbiGB") == 6" means that the GB is a complete GB (as if DegreeLimit was not specified.)
	    if sum toList apply(subalgComp#"Pending", i -> #i) == 0 and rawStatus1 raw (subalgComp#"sagbiGB") == 6 and currDegree > maxGensDeg then (
                R.cache.SagbiDone = true;
                if (o.PrintLevel > 0) then (
		    print("Finite SAGBI basis was found.");
		    );
            	);
            );
	currDegree = currDegree + 1;
    	);
    if currDegree > o.Limit and o.PrintLevel > 0 then (
	print("Limit was reached before a finite SAGBI basis was found");
    	);
    R.cache.SagbiGens
)

subalgebraBasis Matrix := o -> gensMatrix -> (
    R := subring gensMatrix;
    subalgebraBasis(R,o)
    )
