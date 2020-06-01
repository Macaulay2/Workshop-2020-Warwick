-*
Examples from the book Groebner Bases and Convex Polytopes by Bernd Sturmfels.
Chapter 11 "Canonical Subalgebra Bases"
*-

-- example 11.2
--complete?

-- example 11.7
restart
needsPackage "SubalgebraBases"
R= QQ[t_(1,1)..t_(3,3),MonomialOrder=>Lex]
T=transpose genericMatrix(R,3,3)
S32=sort subsets(3,2)
kf = matrix{flatten apply(S32, S1 -> (
	    apply(S32, S2 -> (
		    (i1, j1) := (first S1, last S1);
		    (i2, j2) := (first S2, last S2);
		    det T_{i1,j1}^{i2,j2}
		    )
    		)
	    )
	)}
A = subring kf
sbA = subalgebraBasis A
peek A.cache

-- no finite subalg basis
restart
R=QQ[x,y,MonomialOrder=>Lex]
(options R).MonomialOrder
needsPackage "SubalgebraBases"
methods subring
A = subring {x+y,x*y,x*y^2}
elapsedTime sB = subalgebraBasis(A, Limit=>30);
elapsedTime sB = subalgebraBasis(A, Limit=>30); -- this should take less time