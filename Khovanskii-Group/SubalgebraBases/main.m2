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

-- A wrapper around rawSubduction.
-- Note: See "Using SAGBI bases to compute invariants" by Stillman and Tsai
-- for a description of how this subduction function is supposed to work.
subduction = method(TypicalValue => Matrix)
subduction(Subring, Matrix) := (subR, M) -> (      
    
    subalgComp := subR.cache.SubalgComputations;
    Mtensor := M;
    ambR := ambient subR;
        
    if ring M === subalgComp#"TensorRing" then (
    	-- It is already in the tensor ring, nothing to do.
	) else if ring M === ambR then(
	Mtensor = subalgComp#"InclusionBase"(M);
	) else(
	error "M should have entries from ambient(subR) (or possibly subR.cache.SubalgComputations.TensorRing.)";
	);
    
    -- Check that subalgComp.sagbiGB was computed in a degree that is high enough to handle M.
    degreesM := apply(first entries Mtensor, mono -> first degree mono);
    degM := 0;
    if degreesM === {} then (
	-- M was an empty matrix or otherwise didn't make sense.
	return M; 
	) else(
	degM = max degreesM;
	);
    degGB := {};
    
    if subalgComp#?"sagbiGB" == true then (
	degGB = (subalgComp#"sagbiGB")#"stopping options".DegreeLimit;
	) else (
	subalgComp#"sagbiGB" = gb(subalgComp#"SyzygyIdeal", DegreeLimit => degM);
	degGB = degM;
	);
    
    if degGB === {} then (
	-- This means DegreeLimit wasn't specified and it computed the full GB. Nothing to do.
	)else if degM > degGB then(
	-- degGB is not large enough to handle degM.
	-- This assumes that subalgComp.SyzygyIdeal is up to date.
       	subalgComp#"sagbiGB" = gb(subalgComp#"SyzygyIdeal", DegreeLimit => degM);
	);
    
    F := subalgComp#"Substitution";
    C := subalgComp#"sagbiGB";
    numblocks := rawMonoidNumberOfBlocks raw monoid ambR;
    rawSubduction(numblocks, raw Mtensor, raw F, raw C)
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

-- caches computation results inside of R
subalgebraBasis Subring := o -> R -> (
-- Declaration of variables
    -- baseRing is the ring of the input matrix
    -- sagbiGens is a list/matrix of all the generators
    -- semiRing is the free semigroup ring formed from the SAGBI generators
        -- I don't think that this is ever used.  Should be deleted!
    -- tensorRing is the ring formed from the tensor of the base ring and semigroup ring
    -- Pending is a list of lists, sorting elements of the algebra by degree

    -- ProjectionInclusion is the projection from the tensor ring to the semiRing by sending the base ring generators to 0.
    -- projectionToBase is the projection from the tensor ring to the baseRing by sending the SAGBI generators to 0.
    -- inclusionOfBase is the inclusion of the baseRing into the tensorRing

    -- currDegree is a number representing the current degree of interest
    -- nLoops is the number of loops of the algorithm
    -- maxDegree is the maximum degree of interest by the algorithm
    -- nNewGenerators is the number of new generators found in a loop
    
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
    sagbiGB := null;
    syzygyPairs := null;
    newElems := null;

    subalgComp#"Pending" = new MutableList from toList(o.Limit+1:{});

    -- Create an empty matrix of generators.
    R.cache.SagbiGens = matrix(ambient R,{{}});

    -- Get the maximum degree of the generators.
    -- This is used as a stopping condition.
    maxGensDeg := (max degrees source gens R)_0;

    -- Only look at generators below degree limit.  Add those generators to the SubalgebraGenerators
    reducedGens := compress submatBelowDegree(gens R, o.Limit+1);
    insertPending(R, reducedGens, o.Limit);

    -- Remove elements of coefficient ring
    (subalgComp#"Pending")#0 = {};
    
    -- Compute the initial value of subalgComp.
    processPending(R, o.Limit);
    currDegree = subalgComp#"CurrentLowest" + 1;
        
    while currDegree <= o.Limit and not R.cache.SagbiDone do (  	
        -- Construct a Groebner basis to eliminiate the base elements generators from the SyzygyIdeal.
	
    	-- At this point, the entries of subalgComp.SyzygyIdeal look like this:
	-- p_i - [one of the generators of SagbiGens]
    	-- where none of the p_i are involved in any of the [one of the generators of SagbiGens] terms
	subalgComp#"sagbiGB" = gb(subalgComp#"SyzygyIdeal", DegreeLimit => currDegree);
	zeroGens := submatByDegrees(mingens ideal selectInSubring(1, gens (subalgComp#"sagbiGB")), currDegree);	        
       	syzygyPairs = subalgComp#"Substitution"(zeroGens);
	
        if subalgComp#"Pending"#currDegree != {} then (
            syzygyPairs = syzygyPairs | subalgComp#"InclusionBase"(matrix{subalgComp#"Pending"#currDegree});
            subalgComp#"Pending"#currDegree = {};
            );
	
       	subd := subduction(R, syzygyPairs);
	if entries subd != {{}} then (
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
-- old way: intermediate results aren't cached
subalgebraBasis Matrix := o -> gensMatrix -> (
    R := subring gensMatrix;
    subalgebraBasis(R,o)
    )