needsPackage "NumericalHilbert"
needsPackage "NumericalAlgebraicGeometry"

multMat = method()
multMat(Ideal, RingElement) := (I, f) -> (
    if dim I =!= 0 then error("ideal must be zero-dimensional");
    R := ring I;
    monomBasis := basis(R/I);
    M = module ideal monomBasis;
    sub(matrix map(M,M,f), QQ)
    )

eigenSpace = method()
eigenSpace(Matrix, Number) := (M, e) -> (
    (eValues, eVectors) := eigenvectors M;
    ePositions := positions(eValues, i -> i == e);
    eVectors_ePositions
    )

checkIntersect = method()
checkIntersect(List) := L -> (
    M := last accumulate((i,j) -> i||j, L);
    if isFullNumericalRank M then false
    else true
    )
		

eigenSpaceList = method()
eigenSpaceList(Matrix) := M -> (
    (eValues, eVectors) := eigenvectors M;
    eValues = toList unique eValues;
    apply(eValues, i -> (i, eigenSpace(M, i)))
    )

stickelbergerSolve = method()
stickelbergerSolve Ideal := I -> (
    if dim I =!= 0 then error("ideal must be zero-dimensional");
    R := ring I;
    J := ideal gens gb I;
    gbList := flatten entries gens J;
    eList := apply(gens R, i -> eigenSpaceList multMat(J, i)); 
    outList = {};
    tmpList = {{}};
    for x in eList do (
    	outList = tmpList;
    	tmpList = {{}};
    	for y in x do (
	    for z in outList do (
	    	newitem = {z};
	    	if newitem == {{}} then newitem = {y}
	    	else newitem = {newitem#0, y};
	    	tmpList = append(tmpList,flatten newitem);
	    	);
	    );
    	outList = tmpList;
    	);
    eigList := delete( ,apply(outList, i -> 
	    if length i == length eList then i else null));
    matrixList := apply(eigList, i -> 
	apply(i, j -> transpose last j));
    delete( ,apply(length matrixList, i -> if checkIntersect matrixList#i then (
	    apply(eigList#i, j -> first j)
	    )
	else null
	))
    )
  

end

restart
needs("Stickelberger.m2")
needsPackage "EigenSolver"

R = QQ[x,y]
I = ideal (x^2-1,y^3-1)
stickelbergerSolve I

B = QQ[x0,x1,x2]
I = ideal(random(3,B), random(3,B))
R = QQ[x1,x2]
spe = map(R,B,matrix{{1,x1,x2}})
J = spe I



-- takes too long time
n = 3
R = QQ[a_1..a_(n-1),b_1..b_(n-1)]
S = R[x,y]/ideal(x*y-1)
f = sum(n-1, i -> R_i*x^(i+1) + R_(n-1+i)*y^(i+1)) + x^n + y^n
I = sub(ideal apply(toList(2..2*n-1), i -> last coefficients(f^i, Monomials => {1_S})), R)


-- multiple root example
R = QQ[x,y,z]
I = ideal(x^3 - 3*x^2 *y + 3*x*y^2 - y^3 - z^2, z^3 -3*z^2*x + 3*z*x^2 -x^3-y^2, y^3 - 3*y^2 * z + 3*y*z^2 - z^3 -x^2)
