restart
debugLevel = 1
Precision => 1000
needsPackage "EigenSolver"
needsPackage "ExampleSystems"
needsPackage "NumericalCertification"
load "barebonesSolve.m2"


projOnto = method()
projOnto(Matrix, Matrix) := (u, v) -> (
    dotProdu := first flatten entries (conjugate transpose u * u);
    dotProduv := first flatten entries (conjugate transpose u * v);
    (dotProduv/dotProdu) * u
    )

gramSchmidt = method()
gramSchmidt(Matrix) := V -> (
    k := {};
    u1 := V_{0};
    u := u1;
    mat := 1/(norm(2,u)) * u;
    uMat := u;
    last apply(numcols V - 1, i -> (
	    k = apply(i+1, j -> projOnto(uMat_{j}, V_{i+1}));
	    u = V_{i+1} - sum k;
	    uMat = uMat|u;
	    e = 1/(norm(2,u)) * u;
	    mat = mat|e
	    )
	)
    )


schurDecomposition = method()
schurDecomposition(Matrix) := M -> (
    dummyM := M;
    Ulist := apply(numcols M - 1, i -> (
	    (E,V) := eigenvectors dummyM;
--	    print dummyM;
	    U = gramSchmidt V;
	    dummyM = submatrix(transpose(conjugate(U)) * dummyM * U, {1..numcols M - 1 - i}, {1..numcols M - 1 - i});
	    (id_(CC^(i))|0*random(CC^(i),CC^(numcols M - i)))||(transpose(0*random(CC^(i),CC^(numcols M - i)))|U)
	    )
	);
    Q := product Ulist;
    (Q, (conjugate transpose Q) * M * Q)
    )

generalizedEigenvalueProblem = method()
generalizedEigenvalueProblem(List) := L -> (
    len := #L;
    coeffList := apply(len-1, i -> random RR);
    coeffList = append(coeffList, 1-sum coeffList);
    linComb := sum apply(len, i -> coeffList#i * L#i);
    (Q,T) := schurDecomposition linComb;
    apply(L, i -> (conjugate transpose Q) * i * Q)
    )
    
    


barebonesSolveRevised = method(Options=>{OutputVariables=>null,QuotientBasis=>null}) 
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
barebonesSolveRevised(Ideal,List,List) := List => 
opts -> (I,S,signature) -> (
    outputVariables := opts.OutputVariables; 
    if outputVariables===null then error "TO DO: set the default";
    quotientBasis := opts.QuotientBasis;
    if quotientBasis===null then error "TO DO: set the default";
    R := ring I;
    CR := coefficientRing R;
    assert(S_0 === 1_R);
    monomialsAndCoefficients := coefficients matrix {S | apply(signature, mf->first mf * last mf)};
    monoms := first monomialsAndCoefficients;
    coeffs := sub(last monomialsAndCoefficients, CR);
    Icoeffs := coeffs_{#S..(numcols coeffs - 1)};
    Scoeffs := coeffs_{0..#S-1};
    L := image Scoeffs;
    if debugLevel > 0 then print (numgens L, numcols monoms);
    projL := (last SVD sub(transpose mingens L,CC))^{0..numgens L-1}; -- orthogonal projection to L
    if debugLevel > 0 then print numcols projL;
    IcapL := intersect(L, image Icoeffs); -- this should be done numerically
    if debugLevel > 0 then print numgens IcapL;
    projC := (last SVD(transpose(projL*sub(mingens IcapL,CC))))^{numgens IcapL..numgens L-1}; -- orthogonal projection to (I \cap L)^\perp
    assert(numrows projC == #quotientBasis); --TO DO: error message
    ind := new HashTable from apply(#S,i->S_i=>i);
    projection := projC*projL; -- this is an ofthodonal projection (but doesn't need to be)  
    A := projection*Scoeffs_(apply(quotientBasis,g->ind#g));
    invA := inverse A;
    BList := apply(outputVariables, y -> invA * projection*Scoeffs_(apply(y*quotientBasis, g-> ind#g)));
    eSolve := generalizedEigenvalueProblem BList;
    apply(numcols first eSolve, i -> apply(eSolve, j -> j_(i,i)))
    )


end


TEST///
-- Schur decomposition test
load "testExampleSystems.m2"
A = random(CC^5, CC^5);
eps = 1e-10;
(Q,T) = schurDecomposition A;
assert(clean_eps((inverse Q) - (conjugate transpose Q))==0)
assert(clean_eps(A - Q* T *(conjugate transpose Q))==0)
T = clean_eps T;
assert all(flatten apply(numcols T, i -> apply(i, j -> T_(i,j) == 0)), identity)
///


-- This script tests various examples in ExampleSystems.



-- bellido 
-- There were 40 solutions found in 0.546875 seconds (with a Bezout bound of 64).
load "testExampleSystems.m2"
J = ideal bellido(QQ)
maxIdeg = max (J_*/degree/sum); -- max degree of gens of I : 2
elapsedTime sols = zeroDimSolve J -- 0.0377411 seconds
#sols
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.3385
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)



I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max flatten apply(flatten entries gens I, i -> degree i); -- max degree of gb of I : 5
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
elapsedTime barebonesSolve(J,S,x,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- took too long time on lin alg.



-- boon
-- There were 8 solutions found in 10.1037 seconds (with a Bezout bound of 1024).
load "testExampleSystems.m2"
J = ideal boon(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of I : 4
elapsedTime sols = zeroDimSolve J -- 0.0177399 seconds
#sols
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.06555
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)
    

I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum)  --max flatten apply(flatten entries gens I, i -> degree i); -- max degree of gb of I : 4
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
--x = R_0;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));

elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 2.89916
sols = apply(oo, i -> point {i}); -- 8 solutions
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.05665





-- butcher
-- There were 368 solutions found in 82.589 seconds (with a Bezout bound of 4608).
-- dim 3
load "testExampleSystems.m2"
J = ideal butcher(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of I : 4
elapsedTime sols = zeroDimSolve J -- fail since J has dim 3.
#sols
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))





-- caprasse
-- There were 64 solutions found in 0.525506 seconds (with a Bezout bound of 144).
load "testExampleSystems.m2"
J = ideal caprasse(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of I : 4
elapsedTime sols = zeroDimSolve J; -- 0.0310881 seconds
#sols -- 56 not 64.
sort apply(sols, p -> norm evaluate(gens J, p)) -- 8 of them are incorrect solutions
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.277
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)
#((values oo)#2) -- number of certified solutions : 24


I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum)  -- 7
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg+4);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
--signature = flatten apply(degDmonoms,m->(apply(I_*,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 93.5149
sols = apply(oo, i -> point {i}); -- 56 solutions
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.24913
-- nothing certified


-- cassou
-- There were 2 solutions found in 39.0405 seconds (with a Bezout bound of 1344).
load "testExampleSystems.m2"
J = ideal cassou(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of I : 8
elapsedTime sols = zeroDimSolve J; -- 0.0164434 seconds
#sols -- 16 not 2.
sort apply(sols, p -> norm evaluate(gens J, p)) 
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.1129
-- nothing certified
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols) -- all of them are not approximate solutions.


theirs = solveSystem polySystem J
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), theirs) -- not certified. 
-- Why?
computeConstants(polySystem sub(J, CC[gens ring J]), first theirs) -- small beta, large gamma.
clean_1e-5 evaluate(jacobian J, first theirs) -- Jacobian with large entries.


I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum);  -- 8
--D = max flatten apply(flatten entries gens I, i -> degree i); -- max degree of gb of I : 2
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- returns large residual solutions
-- 159.448 seconds
sols = apply(oo, i -> point {i}); -- 16 solutions
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.09471
-- nothing certified
elapsedTime alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)
-- nothing certified.


-- cohn3
-- There were 72 solutions found in 13.8985 seconds (with a Bezout bound of 1080).
-- dim J = 1
load "testExampleSystems.m2"
J = ideal cohn3(QQ)
R = ring J;
maxIdeg = max (J_*/degree/sum) -- max degree of gens of I : 6
elapsedTime sols = zeroDimSolve J;



-- cyclic 5
-- There were 70 solutions found in 2.737 seconds with 5 variables.
load "testExampleSystems.m2"
J = ideal cyclic(5,QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of I : 5
elapsedTime sols = zeroDimSolve J; -- 0.0394186 seconds
#sols -- 70
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.336583
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)

I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum);  -- 8
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- not working, it might take too long time on linear alg.



-- cyclic 3
load "testExampleSystems.m2"
J = ideal cyclic(3,QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 3
elapsedTime sols = zeroDimSolve J; -- 0.00886598 seconds
#sols -- 6
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.03941
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)

I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum);  -- 4
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(I_*,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 0.0683171
sols = apply(oo, i -> point {i}); -- 6 solutions
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.0402972
elapsedTime alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)





-- eco8
-- There were 64 solutions found in 9.197 seconds (with a Bezout bound of 1458).
load "testExampleSystems.m2"
J = ideal eco8(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 3
elapsedTime sols = zeroDimSolve J; -- 0.503384 seconds
#sols -- 64
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.474047
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)

I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum);  -- 5
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(I_*,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
barebonesSolve(I,S,x,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- not working, it might take too long time on linear alg.





-- geneig
-- There were 10 solutions found in 1.447 seconds (with a Bezout bound of 243).
load "testExampleSystems.m2"
J = ideal geneig(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 3
elapsedTime sols = zeroDimSolve J; -- 0.0299402 seconds
#sols -- 10
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.0880489
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)

I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum)  -- 3
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 15.9242
sols = apply(oo, i -> point {i}); -- 10 solutions
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.0880489
elapsedTime alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)
-- all solutions are certified. 



-- katsura
-- There were 512 solutions found in 2.804 seconds with 10 variables.
load "testExampleSystems.m2"
J = ideal katsura(10,QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 2
elapsedTime sols = zeroDimSolve J; -- took too long time to compute gb.


-- katsura 8
load "testExampleSystems.m2"
J = ideal katsura(8,QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 2
elapsedTime sols = zeroDimSolve J; -- 3.37624
#sols -- 128
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.0880489
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)

I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum)  -- 8
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- took too long time on lin alg.



-- katsura 5
load "testExampleSystems.m2"
J = ideal katsura(5,QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 2
elapsedTime sols = zeroDimSolve J; -- 0.0226005
#sols -- 16
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.03941
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)

I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum)  -- 5
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 20.2684
sols = apply(oo, i -> point {i});
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.107208
elapsedTime alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)
-- all certified.




-- ku10
-- There were 4 solutions found in 2.266 seconds (with a Bezout bound of 1024).
load "testExampleSystems.m2"
J = ideal ku10(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 2
elapsedTime sols = zeroDimSolve J; -- 0.0318752
#sols -- 2 not 4
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.0287836
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)

I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum)  -- 2
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 0.72186
sols = apply(oo, i -> point {i}); -- 2 solutions
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.0282631
elapsedTime alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)
-- all certified.



-- lorentz
-- There were 12 solutions found in 0.0562 seconds (with a Bezout bound of 16).
load "testExampleSystems.m2"
J = ideal lorentz(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 2
elapsedTime sols = zeroDimSolve J; -- 0.0119325
#sols -- 11 not 12
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.0776831
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)

I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum)  -- 4
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 0.498564
sols = apply(oo, i -> point {i}); -- 11 solutions
-- 0.953831 seconds == RMK : using lowerDegGens makes solving faster ==
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.0776831
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)
-- all solutions are certified.

theirs = solveSystem polySystem J;
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), theirs)





-- rabmo
-- There were 58 solutions found in 1642.07 seconds (with a Bezout bound of 36000).
load "testExampleSystems.m2"
J = ideal rabmo(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 5
elapsedTime sols = zeroDimSolve J; -- 0.433235
#sols -- 16 not 58
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.134484
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)

I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum)  -- 5
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- out of memory







-- randomSystem
load "testExampleSystems.m2"
J = ideal randomSystem(3,5,QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 5
elapsedTime sols = zeroDimSolve J; -- 0.754872
#sols -- 125
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 1.14966
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols) -- random number of solutions are certified

elapsedTime theirs = solveSystem polySystem J; -- 125, 0.182215 seconds.
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), theirs) -- 1.124
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), theirs) -- all certified


I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum)  -- 13
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- took too long time on lin alg.








-- reimer5
-- There were 146 solutions found in 6.540 seconds (with a Bezout bound of 720).
load "testExampleSystems.m2"
J = ideal reimer5(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 6
elapsedTime sols = zeroDimSolve J; -- 1.5318
#sols -- 144
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.783844
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)


I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum)  -- 10
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- took too long time on lin alg.



-- rose
-- There were 136 solutions found in 4.954 seconds (with a Bezout bound of 216).
load "testExampleSystems.m2"
J = ideal rose(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 9
elapsedTime sols = zeroDimSolve J; -- 0.0676541
#sols -- 136
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.687595
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols) -- 128 solutions are certified


I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum)  -- 11
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- took too long time on lin alg.




-- sendra
-- There were 46 solutions found in 0.111 seconds (with a Bezout bound of 49).
load "testExampleSystems.m2"
J = ideal sendra(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 7
elapsedTime sols = zeroDimSolve J; -- 0.0200045
#sols -- 46
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.228829
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)


I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum)  -- 11
ㅐdegDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
sols = apply(oo, i -> point {i}); -- 46 solutions
-- 11.8806 seconds 
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.236832
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)
-- random number of solutions are certified.




-- trinks
-- There were 10 solutions found in 0.203 seconds (with a Bezout bound of 24).
load "testExampleSystems.m2"
J = ideal trinks(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 3
elapsedTime sols = zeroDimSolve J; -- 0.0154688
#sols -- 10
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.0776831
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)
#((values oo)#2) -- number of certified solutions : 10

I = ideal gens gb J;
R = ring I;
quotientBasis = apply(flatten entries basis(R/I), i -> sub(i, R));
D = max flatten apply(flatten entries gens I, i -> degree i) -- max degree of gb of I : 3
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
sols = apply(oo, i -> point {i}); -- 10 solutions
-- 3.25115 seconds 
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.0776831
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)
#((values oo)#2) -- number of certified solutions : 10





-- virasoro
-- There were 256 solutions found in 0.356 seconds (with a Bezout bound of 256).
load "testExampleSystems.m2"
J = ideal virasoro(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 2
elapsedTime sols = zeroDimSolve J; -- 7.41958
#sols -- 256
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 1.95857
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)
#((values oo)#2) -- random number of solutions certified 

I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum)  -- 9
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- took too long time on lin alg.



-- wright
-- There were 32 solutions found in 0.0944 seconds (with a Bezout bound of 32).
load "testExampleSystems.m2"
J = ideal wright(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 2
elapsedTime sols = zeroDimSolve J; -- 0.0204144
#sols -- 32
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.166867
alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)
#((values oo)#2) -- number of certified solutions : 32

I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum)  -- 6
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 23.7229 seconds
elapsedTime certifyRegularSolution(polySystem sub(J, CC[gens ring J]), sols) -- 0.0776831 
elapsedTime alphaTheoryCertification(polySystem sub(J, CC[gens ring J]), sols)
-- number of certified solutions : 32




