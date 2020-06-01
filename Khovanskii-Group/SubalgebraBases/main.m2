debug Core -- gets rid of "raw" error during installation. probably a better way...

export {
    "subalgebraBasis",
    "sagbi" => "subalgebraBasis", 
    "PrintLevel",
    -- things that get cached in the computation: do we really want to export all of these?
    "ProjectionBase",
    "SyzygyIdeal",
    "Pending",
    "Substitution",
    "SagbiDegrees",
    "TensorRing",
    "SubalgComputations",
    "InclusionBase",
    "ProjectionInclusion",
    "SagbiGens",
    "SagbiDone"
    }

-- Main function for computing a subalgebraBasis.
-- For now, subduction is performed in the engine (in older version, Strategy toggles between engine and top-level implementation)
subalgebraBasis = method(Options => {
    Strategy => null,
    Limit => 100,
    PrintLevel => 0})
-- caches computation results inside of R
subalgebraBasis Subring := o -> R -> (
-- Declaration of variables
    -- baseRing is the ring of the input matrix
    -- sagbiGens is a list/matrix of all the generators
    -- semiRing is the free semigroup ring formed from the SAGBI generators
        -- I don't think that this is ever used.  Should be deleted!
    -- tensorRing is the ring formed from the tensor of the base ring and semigroup ring
    -- Pending is a list of lists, sorting elements of the algebra by degree

    -- projectionInclusion is the projection from the tensor ring to the semiRing by sending the base ring generators to 0.
    -- projectionToBase is the projection from the tensor ring to the baseRing by sending the SAGBI generators to 0.
    -- inclusionOfBase is the inclusion of the baseRing into the tensorRing

    -- currDegree is a number representing the current degree of interest
    -- nLoops is the number of loops of the algorithm
    -- maxDegree is the maximum degree of interest by the algorithm
    -- nNewGenerators is the number of new generators found in a loop

    -- CHECK IF THESE EXIST to avoid recomputing???

    R.cache.SubalgComputations = new MutableHashTable;

    subalgComp := R.cache.SubalgComputations;

    R.cache.SagbiDegrees = {};
    subalgComp.TensorRing = null;  -- RS
    subalgComp.SyzygyIdeal = null; -- J

    subalgComp.ProjectionInclusion = null; -- RStoS
    subalgComp.ProjectionBase = null;      -- RStoR
    subalgComp.InclusionBase = null;       -- RtoRS
    subalgComp.Substitution = null;        -- Gmap

    currDegree := null;     -- d
    nLoops := null;         -- nloops
    R.cache.SagbiDone = false;
    sagbiGB := null;
    syzygyPairs := null;
    newElems := null;

    subalgComp.Pending = new MutableList from toList(o.Limit+1:{}); -- Pending

    -- Create an empty matrix of generators.
    R.cache.SagbiGens = matrix(ambient R,{{}});

    -- Get the maximum degree of the generators.
        -- This is used as a stopping condition.
    maxGensDeg := (max degrees source gens R)_0;

    -- Only look at generators below degree limit.  Add those generators to the SubalgebraGenerators
    reducedGens := compress submatrixBelowDegree(gens R, o.Limit+1);
    insertPending(R, reducedGens, o.Limit);

    -- Remove elements of coefficient ring
    subalgComp.Pending#0 = {};

    -- Get the lowest degree of the pending list.  Add 1 and initialize to number of loops
    currDegree = grabLowestDegree(R, o.Limit) + 1;

    nLoops = currDegree;
    local subducted;
    ambR := ambient R;
    -- While the number of loops is within the limit and the isDone flag is false, continue to process
    while nLoops <= o.Limit and not R.cache.SagbiDone do (
        nLoops = nLoops + 1;
    
        -- Construct a Groebner basis to eliminiate the base elements generators from the SyzygyIdeal.
        sagbiGB = gb(subalgComp.SyzygyIdeal, DegreeLimit=>currDegree);
        syzygyPairs = subalgComp.Substitution(submatrixByDegrees(mingens ideal selectInSubring(1, gens sagbiGB), currDegree));
        if subalgComp.Pending#currDegree != {} then (
            syzygyPairs = syzygyPairs | subalgComp.InclusionBase(matrix{subalgComp.Pending#currDegree});
            subalgComp.Pending#currDegree = {};
        );
        
        subducted = subalgComp.ProjectionBase(map(subalgComp.TensorRing,rawSubduction(rawMonoidNumberOfBlocks raw monoid ambR, raw syzygyPairs, raw subalgComp.Substitution, raw sagbiGB)));
        newElems = compress subducted;
        if numcols newElems > 0 then (
            insertPending(R, newElems, o.Limit);
            currDegree = grabLowestDegree(R, o.Limit);
        )
        else (
            if sum toList apply(subalgComp.Pending, i -> #i) === 0 and rawStatus1 raw sagbiGB == 6 and currDegree > maxGensDeg then (
                R.cache.SagbiDone = true;
                if (o.PrintLevel > 0) then << "SAGBI basis is FINITE!" << endl;
            	)
            );
            currDegree = currDegree + 1;
    	);
    R.cache.SagbiGens
)
-- old way: intermediate results aren't cached
subalgebraBasis Matrix := o -> gensMatrix -> (
    R := subring gensMatrix;
    subalgebraBasis(R,o)
)
