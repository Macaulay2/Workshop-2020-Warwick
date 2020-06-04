--Example: Maximal minors of a generic 2x4 matrix.

--produces ideal of generic minors
genericMinors = (k, m, n) -> (
    -- k by k minors of a generic m by n matrix
    R = ZZ/101[x_(0,0) .. x_(m-1,n-1)];
		X = transpose genericMatrix(R, n, m);
    gens minors(k, X));

I = ideal genericMinors(2,2,4); --maximal minors of a 2x4 matrix

S = subring I_* --make a subring out of gens I

subalgebraBasis S --compute SAGBI basis. Now information in S has been changed.

peek S.cache.SubalgComputations --to see information from computation

peek S.cache.SubalgComputations.SyzygyIdeal --to get toric ideal

G = gens gb S.cache.SubalgComputation.SyzygyIdeal --get GB of toric ideal

rels = selectInSubring(1,G) --does elimination order on toric ideal to get only algebraic relations

S.cache.SubalgComputations.Substitution(rels) -- substitutes original f's back into relation
--what is this output?
--Can one subduct this output with the original generators to get the Plucker relation?