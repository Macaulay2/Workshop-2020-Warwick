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

-- Adds generators to the set of generators and computes the syzygies of the generators.  Also defines the appropriate ring maps for future use.
    -- R is of Type Subring
    -- newGens is a matrix of generators to be added
appendToBasis = (R, newGens) -> (
    -- unpack immutable fields
    ambR := ambient R;
    -- this gets modified
    subalgComp := R.cache.SubalgComputations;
    
    -- Add the new generators to the subalgebra generators
    R.cache.SagbiGens = R.cache.SagbiGens | newGens;
    R.cache.SagbiDegrees = R.cache.SagbiDegrees | flatten degrees source newGens;
        
    -- Find the number of generators of the ambient ring and the current list of subalgebra generators
    nBaseGens := numgens ambR;
    nSubalgGens := numcols R.cache.SagbiGens;
    
    -- Create a ring with combined generators of base and subalgebra.  Monoid is needed for constructing a monomial order and the coefficient ring is used to construct the new ring.
    MonoidAmbient := monoid ambR;
    CoeffField := coefficientRing ambR;
    
    -- Add on an elimination order that eliminates the generators of the base.
    -- Create a monoid with variables for both nBaseGens and nSubalgGens.
    -- Degrees of generators are set so that the SyzygyIdeal is homogeneous.
    newOrder := append(MonoidAmbient.Options.MonomialOrder, Weights=>nBaseGens:1);
    
    NewVariables := monoid[
        Variables=>nBaseGens+nSubalgGens,
        Degrees=>join(degrees source vars ambR, degrees source R.cache.SagbiGens),
        MonomialOrder => newOrder];
    
    -- Construct the free monoid ring with coefficients in CoeffField and and variables for both the base and subalgebra.
    subalgComp#"TensorRing" = CoeffField NewVariables;
    
    -- Construct maps between our rings to allow us to move polynomials around
    -- ProjectionInclusion sets the variables corresponding to the base equal to 0.  The result is in the tensor ring.
    -- ProjectionBase sets the variables corresponding to the subalgebra generators equal to 0 and maps into the ambient ring.
    -- InclusionBase is the inclusion map from the base ring to the tensor ring.  The variables are mapped to themselves
    -- Substitution replaces elements of the tensor ring with their formulas in terms of the base ring.
    subalgComp#"ProjectionInclusion" = map(subalgComp#"TensorRing", subalgComp#"TensorRing",
        matrix {toList(nBaseGens:0_(subalgComp#"TensorRing"))} |
	(vars subalgComp#"TensorRing")_{nBaseGens .. nBaseGens+nSubalgGens-1});
    
    subalgComp#"ProjectionBase" = map(ambR, subalgComp#"TensorRing",
        (vars ambR) | matrix {toList(nSubalgGens:0_(ambR))});
    
    subalgComp#"InclusionBase" = map(subalgComp#"TensorRing", ambR,
        (vars subalgComp#"TensorRing")_{0..nBaseGens-1});
    
    subalgComp#"Substitution" = map(subalgComp#"TensorRing", subalgComp#"TensorRing",
        (vars subalgComp#"TensorRing")_{0..nBaseGens-1} | subalgComp#"InclusionBase"(R.cache.SagbiGens));
    
    -- Construct an ideal consisting of variables repesenting subalgebra generators minus their leading term
    subalgComp#"SyzygyIdeal" = ideal(
        (vars subalgComp#"TensorRing")_{nBaseGens..nBaseGens+nSubalgGens-1}-
	subalgComp#"InclusionBase"(leadTerm R.cache.SagbiGens));
    )    

rowReduce = (elems, d) -> (
    return gens gb(elems, DegreeLimit=>d);
    )

--Accepts a 1-row matrix inputMatrix and returns a matrix of columns of inputMatrix whose entries all have total degree less than maxDegree
submatBelowDegree = (inputMatrix,maxDegree) -> (
    selectedCols := positions(0..numcols inputMatrix - 1,
        i -> (degrees source inputMatrix)_i < {maxDegree});
    inputMatrix_selectedCols
    )

--Accepts a 1-row matrix inputMatrix and returns a matrix of columns of inputMatrix where the highest degree entry has total degree equal to currDegree
submatByDegrees = (inputMatrix,currDegree) -> (
    selectedCols := positions(0..numcols inputMatrix - 1,
        i -> (degrees source inputMatrix)_i === {currDegree});
    inputMatrix_selectedCols
    )

-- Reduces the lowest degree list in the pending list.  Adds the results to Pending.  The new lowest degree list in pending is added to the subalgebra basis.  Returns the number of elements added.
    -- !!!Assumes that the pending list has been subducted!!!
    -- R is the subalgebra
processPending = (R, maxDegree) -> (

    subalgComp := R.cache.SubalgComputations;
    -- Finds the current lowest degree of the pending list.
    currentLowest := lowestDegree(R, maxDegree);
    
    -- If there are elements in the pending list, then work on them.
    local reducedGenerators;
    if currentLowest <= maxDegree then (    
	
	reducedGenerators = gens gb(matrix{(subalgComp#"Pending")#currentLowest}, DegreeLimit=>currentLowest);
    	(subalgComp#"Pending")#currentLowest = {};
    	insertPending(R, reducedGenerators, maxDegree);
    	-- Find the lowest degree elements after reduction.
    	currentLowest = lowestDegree(R, maxDegree);
    	-- Count the number of new generators and add them to the basis
    	numNewGenerators := 0;
    	if currentLowest <= maxDegree then (
            numNewGenerators = #(subalgComp#"Pending")#currentLowest;
            appendToBasis(R, matrix{(subalgComp#"Pending")#currentLowest});
            (subalgComp#"Pending")#currentLowest = {};
	    );
    	-- If number of new generators is zero, then nothing was added because pending was empty.  There is no way for pending to be empty unless currentLowest is maxDegree + 1.
    	);
    subalgComp#"CurrentLowest" = currentLowest;
    currentLowest
    )
