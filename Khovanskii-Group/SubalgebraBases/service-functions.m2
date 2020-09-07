-- return the monomial order stashed inside of a ring
getMonomialOrder = R -> (options R).MonomialOrder

-- Makes a pass through the elements in the first row of "candidates" and places them in the correct sublist of subalgComp#"Pending".
    -- R is a subalgebra
    -- candidates is a matrix of elements of the subalgebra.
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
-- Makes a pass through the lists of Pending until it finds something nonempty
    -- R is a subalgebra
    -- maxDegree is an integer
lowestDegree = (R, maxDegree) -> (
    subalgComp := R.cache.SubalgComputations;
    i := 0;
    while i <= maxDegree and (subalgComp#"Pending")#i === {} do i=i+1;
    i
    )

-- Adds newGens to R.cache.SagbiGens. Updates the appropriate rings/maps in R.cache.SubalgComputations.
    -- R is of Type Subring
    -- newGens is a 1-row matrix of generators to be added
appendToBasis = (R, newGens) -> (
    subalgComp := R.cache.SubalgComputations;
    R.cache.SagbiGens = R.cache.SagbiGens | newGens;
    R.cache.SagbiDegrees = R.cache.SagbiDegrees | flatten degrees source newGens;
    subalgComp#"PartialSagbi" = subring(R.cache.SagbiGens);
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
	lowest := matrix({(subalgComp#"Pending")#currentLowest});
	reducedGenerators := gens gb(lowest, DegreeLimit=>currentLowest);
	
	if reducedGenerators_(0,0) == 1 then(
	    error "Could not process the lowest degree S-polynomials. Perhaps autosubduction was turned off?"
	    );
	
	(subalgComp#"Pending")#currentLowest = {};
    	insertPending(R, reducedGenerators, maxDegree);
    	currentLowest = lowestDegree(R, maxDegree);
	
    	if currentLowest <= maxDegree then (	    
            appendToBasis(R, matrix{(subalgComp#"Pending")#currentLowest});
            (subalgComp#"Pending")#currentLowest = {};
	    );
    	);
    subalgComp#"CurrentLowest" = currentLowest;
    )
