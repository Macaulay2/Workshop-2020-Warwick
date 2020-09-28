needsPackage "NumericalAlgebraicGeometry"

barebonesP2 = method() -- compute 'satind' from homogenization and compute inputs for barebonesSolve using satind
barebonesP2 Ideal := J -> ( -- currently assumes dim R = 2 et 2 (nonlinear) equations, homogenizes to P^2
    I := ideal gens gb J;
    R := ring J;
    D := max flatten apply(flatten entries gens I, i -> degree i);
    quotientBasis := flatten entries sub(basis(R/I),R);
    degDmonoms := flatten entries matrix{for d to D list basis(d,R)};
    S := degDmonoms;
    maxIdeg := max (J_*/degree/sum);
    lowerDegGens := select(flatten entries gens I, g -> sum degree g < maxIdeg);
    outputVariables := gens R;
    x := R_0;
    signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
    barebonesSolve(J,S,x,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis)
    )     




barebonesP1P1 = method() 
barebonesP1P1 Ideal := J -> ( -- currently assumes dim R = 2 et 2 equations, homogenizes to P^1 x P^1
    I := ideal gens gb J;
    R := ring J;
    D := max flatten apply(flatten entries gens I, i -> degree i);
    quotientBasis := flatten entries sub(basis(R/I),R);
    degDmonoms := flatten entries matrix{for d to D list basis(d,R)};
    S := degDmonoms;
    maxIdeg := max (J_*/degree/sum);
    lowerDegGens := select(flatten entries gens I, g -> sum degree g < maxIdeg);
    outputVariables := gens R;
    x := R_0;
    signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
    barebonesSolve(J,S,x,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis)
    )



    
    
    

barebonesSolve = method(Options=>{OutputVariables=>null,QuotientBasis=>null}) 
-* 
In:  I, ideal of a polynomial ring R
     S, a list of polynomials whose cosets span R/I as a vector space
     x, the multiplier for which companion matrix is to be constructed
     signature, the "signature" of the method --- a list of pairs (m,f) where m in R and f in I_* 
     outputVariables, a subset of gens R
     quotientBasis, a subset of S that forms a basis of R/I 
Out: a list of vectors of size dim(R/I), 
     each vector contains the value of outputVariables for a point in Spec(R/I)
Assumptions:
     I is radical
     S is linearly independent in R
     x*quotientBasis is contained in S
Description:
     Compute a matrix C whose columns span 
       $$ (<S> \cap <m*f|(m,f)\in signature>)^\perp $$
     and whose rows correspond to elements of $S$
     Solve a generalized eigenvalue problem for appropriately chosen 
     maximal square submatrices A and xA of C.            
*-
barebonesSolve(Ideal,List,RingElement,List) := List => 
opts -> (I,S,x,signature) -> (
    R := ring I;
    CR := coefficientRing R;
    outputVariables := opts.OutputVariables; 
    if outputVariables===null then error "TO DO: set the default";
    quotientBasis := apply(opts.QuotientBasis, g -> lift(g,R));
    if quotientBasis===null then error "TO DO: set the default";
    assert(S_0 === 1_R);
    monomialsAndCoefficients := coefficients matrix {S | apply(signature, mf->first mf * last mf)};
    monoms := first monomialsAndCoefficients;
    coeffs := sub(last monomialsAndCoefficients, CR);
    Icoeffs := coeffs_{#S..(numcols coeffs - 1)};
    Scoeffs := coeffs_{0..#S-1};
    L := image Scoeffs;
    IcapL := intersect(L, image Icoeffs); -- this should be done numerically
    projL := (last SVD sub(transpose mingens L,CC))^{0..numgens L-1}; -- orthogonal projection to L
    projC := (last SVD(transpose(projL*sub(mingens IcapL,CC))))^{numgens IcapL ..numgens L-1}; -- orthogonal projection to (I \cap L)^\perp
    assert(numrows projC == #quotientBasis); --TO DO: error message
    ind := new HashTable from apply(#S,i->S_i=>i);
    projection := projC*projL; -- this is an ofthodonal projection (but doesn't need to be)  
    A := projection*Scoeffs_(apply(quotientBasis,g->ind#g));
    assert(det A!=0); --TO DO: numerical check instead
    xA := projection*Scoeffs_(apply(x*quotientBasis,g->ind#g));
    -- solve the eigen-problem: xA - a*A drops rank 
    (E,V) := eigenvectors(inverse A * xA); --TO DO: generalized e-solver
    assert isSubset(outputVariables,quotientBasis); --TO DO: in fact, this is not necessary
    apply(numcols V,i->apply(outputVariables, y->(
    		    yA := projection*Scoeffs_(apply(y*quotientBasis,g->ind#g));
    	    	    (solve(A*V_{i},yA*V_{i},ClosestFit=>true))_(0,0)
	    	    )) 
	)	
)
end


restart
load "projbarebone.m2"
B = QQ[x0,x1,x2]
I = ideal(random(4,B), random(4,B))
R = QQ[x1,x2]
spe = map(R,B,matrix{{1,x1,x2}})
J = spe I
barebonesP2 J



restart
load "projbarebone.m2"
A = QQ[x0,x,y0,y,Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
I = ideal(random({2,5},A),random({2,5},A))
R = QQ[x,y]
spe = map(R,A,matrix{{1,x,1,y}})
J = spe I
barebonesP1P1 J

