-- return the monomial order stashed inside of a ring
getMonomialOrder = R -> (options R).MonomialOrder

-- Sorts and adds the elements of the matrix "candidates" to the pending list of R
    -- R is a subalgebra
    -- candidates is a matrix of elements of the subalgebra.
    -- Algorithm makes a pass through the elements in the first row of "candidates" and places them in the correct sublist of subalgComp#"Pending".
insertPending = (R, candidates, maxDegree) -> (
    subalgComp := R.cache.SubalgComputations;
    
    if subalgComp#?"Pending" == false then(
	subalgComp#"Pending" = new MutableList from (maxDegree:{});
	);
    
    for i from 0 to (numcols candidates)-1 do(
        -- get the entry of the column and its degree
        candidate := candidates_(0,i);
        level := (degree candidate)_0;
        -- if the degree isn't too large, then add f to the correct list
        if level <= maxDegree then (
	    (subalgComp#"Pending")#level = append((subalgComp#"Pending")#level, candidate);
            );
    	);
    )

-- Finds the lowest nonempty list in Pending
    -- R is a subalgebra
    -- Algorithm makes a pass through the lists of Pending until it finds something nonempty
lowestDegree = (R, maxDegree) -> (
    subalgComp := R.cache.SubalgComputations;
    -- i steps through the lists of Pending
    i := 0;
    while i <= maxDegree and (subalgComp#"Pending")#i === {} do i=i+1;
    i
    )

-- Adds newGens to R.cache.SagbiGens. Updates the appropriate rings/maps in R.cache.SubalgComputations.
    -- R is of Type Subring
    -- newGens is a matrix of generators to be added
appendToBasis = (R, newGens) -> (
    ambR := ambient R;
    subalgComp := R.cache.SubalgComputations;
    
    R.cache.SagbiGens = R.cache.SagbiGens | newGens;
    R.cache.SagbiDegrees = R.cache.SagbiDegrees | flatten degrees source newGens;
        
    -- Find the number of generators of the ambient ring and the current list of subalgebra generators
    nBaseGens := numgens ambR;
    nSubalgGens := numcols R.cache.SagbiGens;
    
    -- A monoid with an elimination order that can be used to eliminate the generators of the base.
    NewVariables := monoid[
        Variables=>nBaseGens+nSubalgGens,
        Degrees=>join(degrees source vars ambR, degrees source R.cache.SagbiGens),
        MonomialOrder => Eliminate nBaseGens
	];
    
    CoeffField := coefficientRing ambR;
    subalgComp#"TensorRing" = CoeffField NewVariables;
        
    -- ProjectionInclusion sets the variables corresponding to the base equal to 0.  The result is in the tensor ring.
    subalgComp#"ProjectionInclusion" = map(subalgComp#"TensorRing", subalgComp#"TensorRing",
        matrix {toList(nBaseGens:0_(subalgComp#"TensorRing"))} |
	(vars subalgComp#"TensorRing")_{nBaseGens .. nBaseGens+nSubalgGens-1});
    
    -- ProjectionBase sets the variables corresponding to the subalgebra generators equal to 0 and maps into the ambient ring.
    subalgComp#"ProjectionBase" = map(ambR, subalgComp#"TensorRing",
        (vars ambR) | matrix {toList(nSubalgGens:0_(ambR))});
    
    -- InclusionBase is the inclusion map from the base ring to the tensor ring.  The variables are mapped to themselves
    subalgComp#"InclusionBase" = map(subalgComp#"TensorRing", ambR,
        (vars subalgComp#"TensorRing")_{0..nBaseGens-1});
    
    -- Replaces elements of the tensor ring with their formulas in terms of the base ring.
    subalgComp#"Substitution" = map(subalgComp#"TensorRing", subalgComp#"TensorRing",
        (vars subalgComp#"TensorRing")_{0..nBaseGens-1} | subalgComp#"InclusionBase"(R.cache.SagbiGens));
    
    -- An ideal consisting of variables repesenting subalgebra generators minus their leading term
    subalgComp#"SyzygyIdeal" = ideal(
        (vars subalgComp#"TensorRing")_{nBaseGens..nBaseGens+nSubalgGens-1}-
	subalgComp#"InclusionBase"(leadTerm R.cache.SagbiGens));
    )

--Accepts a 1-row matrix inputMatrix and returns a matrix of columns of inputMatrix whose entries all have total degree less than maxDegree
submatBelowDegree = (inputMatrix,maxDegree) -> (
    selectedCols := positions(0..numcols inputMatrix - 1,
        i -> (degrees source inputMatrix)_i < {maxDegree});
    inputMatrix_selectedCols
    )

--Accepts a 1-row matrix inputMatrix and returns a matrix of columns of inputMatrix where the highest degree entry has total degree equal to currDegree
submatByDegree = (inputMatrix, currDegree) -> (
    selectedCols := positions(0..numcols inputMatrix - 1,
        i -> (degrees source inputMatrix)_i === {currDegree});
    inputMatrix_selectedCols
    )


-- Reduces the lowest degree in subalgComp#"Pending", updating subalgComp#"Pending" and subalgComp#"sagbiGB".
-- The various maps, tensor ring, and syzygy ideal are updated to reflect this change.
-- !!!Assumes that the pending list has been subducted!!!
   -- R is the subalgebra.
   -- maxDegree is the degree limit.
processPending = (R, maxDegree) -> (

    subalgComp := R.cache.SubalgComputations;
    currentLowest := lowestDegree(R, maxDegree);
    
    if currentLowest <= maxDegree then (
	-- remove redundant elements of the lowest degree in subalgComp#"Pending".
	reducedGenerators := gens gb(matrix{(subalgComp#"Pending")#currentLowest}, DegreeLimit=>currentLowest);
    	(subalgComp#"Pending")#currentLowest = {};
    	insertPending(R, reducedGenerators, maxDegree);
    	-- Find the lowest degree elements after reduction.
    	currentLowest = lowestDegree(R, maxDegree);
    	-- Add new generators to the basis
    	if currentLowest <= maxDegree then (
            appendToBasis(R, matrix{(subalgComp#"Pending")#currentLowest});
            (subalgComp#"Pending")#currentLowest = {};
	    );
    	);
    subalgComp#"CurrentLowest" = currentLowest;
    )