--Input: (I) an ideal
--Output: (I) the ideal of algebraic relations among initial terms of gens I
--Description: Finds toric syzygies among initial terms of I

toricSyz = method()
toricSyz(Ideal) := Ideal => (I) -> (
    gensI = flatten entries (gens I);
    inI = {};
    for i from 0 to length gensI - 1 do (
    	inI = inI|{leadTerm (gensI)_i};
    	);
    dummyVars = ZZ/101[T_0..T_((length inI)-1)];
    F = map((ring I), dummyVars, inI);
    L = kernel F
    )

-----
--Example: Maximal minors of a generic 2x4 matrix.
--We expect to get 2 terms of the standard Plucker relation

genericMinors = (k, m, n) -> (
    -- k by k minors of a generic m by n matrix
    R = ZZ/101[x_(0,0) .. x_(m-1,n-1)];
		X = transpose genericMatrix(R, n, m);
    gens minors(k, X));

I = ideal genericMinors(2,2,4);

toricSyz(I)

--Next step: How to lift this algebraic relation to the Plucker relation
--using the subduction algorithm?