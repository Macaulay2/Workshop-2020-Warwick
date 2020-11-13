restart
debugLevel = 1
Precision => 1000
needsPackage "EigenSolver"
needsPackage "ExampleSystems"


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
    

    
ingredients = method()
-- Function returning inputs for barebonseSolveRevised from a given input ideal.
ingredients(Ideal) := J -> (
    maxIdeg := max (J_*/degree/sum); -- maximal degree of a given generator of J
    I := ideal gens gb J; 
    R := ring I;
    outputVariables := gens R;
    quotientBasis := flatten entries sub(basis(R/I), R);
    D := max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum); -- choosing a proper degree for monomials from quotientBasis, gb J, and J.
    S := flatten entries matrix{for d to D list basis(d,R)};
    lowerDegGens := select(flatten entries gens I, g -> sum degree g < maxIdeg); -- among gb, choose elements whose degree is less than maxIdeg
    signature := flatten apply(S, m->(apply(J_*|lowerDegGens,f->(m,f))));
    (quotientBasis, S, signature, outputVariables)
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





-- boon
-- There were 8 solutions found in 10.1037 seconds (with a Bezout bound of 1024).
load "testExampleSystems.m2"
J = ideal boon(QQ)
elapsedTime sols = zeroDimSolve J -- 0.0177399 seconds
#sols
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))
    
(quotientBasis, S, signature, outputVariables) = ingredients J;
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 2.89916
sols = apply(oo, i -> point {i}); -- 8 solutions
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))








-- cyclic 3
load "testExampleSystems.m2"
J = ideal cyclic(3,QQ)
elapsedTime sols = zeroDimSolve J; -- 0.00886598 seconds
#sols -- 6
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))

(quotientBasis, S, signature, outputVariables) = ingredients J;
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 0.0683171
sols = apply(oo, i -> point {i}); -- 6 solutions
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))







-- katsura 5
load "testExampleSystems.m2"
J = ideal katsura(5,QQ)
elapsedTime sols = zeroDimSolve J; -- 0.0226005
#sols -- 16
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))

(quotientBasis, S, signature, outputVariables) = ingredients J;
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 20.2684
sols = apply(oo, i -> point {i});
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))




-- ku10
-- There were 4 solutions found in 2.266 seconds (with a Bezout bound of 1024).
load "testExampleSystems.m2"
J = ideal ku10(QQ)
elapsedTime sols = zeroDimSolve J; -- 0.0318752
#sols -- 2 not 4
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))

(quotientBasis, S, signature, outputVariables) = ingredients J;
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 0.72186
sols = apply(oo, i -> point {i}); -- 2 solutions
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))



-- lorentz
-- There were 12 solutions found in 0.0562 seconds (with a Bezout bound of 16).
load "testExampleSystems.m2"
J = ideal lorentz(QQ)
elapsedTime sols = zeroDimSolve J; -- 0.0119325
#sols -- 11 not 12
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))

(quotientBasis, S, signature, outputVariables) = ingredients J;
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 0.498564
sols = apply(oo, i -> point {i}); -- 11 solutions
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))







-- sendra
-- There were 46 solutions found in 0.111 seconds (with a Bezout bound of 49).
load "testExampleSystems.m2"
J = ideal sendra(QQ)
elapsedTime sols = zeroDimSolve J; -- 0.0200045
#sols -- 46
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))


(quotientBasis, S, signature, outputVariables) = ingredients J;
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
sols = apply(oo, i -> point {i}); -- 46 solutions
-- 11.8806 seconds 
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))

-- refine 'sols' using Newton's method.
sols = apply(sols, i -> newton(polySystem J, i))
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))




-- trinks
-- There were 10 solutions found in 0.203 seconds (with a Bezout bound of 24).
load "testExampleSystems.m2"
J = ideal trinks(QQ)
elapsedTime sols = zeroDimSolve J; -- 0.0154688
#sols -- 10
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))

(quotientBasis, S, signature, outputVariables) = ingredients J;
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
sols = apply(oo, i -> point {i}); -- 10 solutions
-- 3.25115 seconds 
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))

-- refine 'sols' using Newton's method.
sols = apply(sols, i -> newton(polySystem J, i))
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))










-- wright
-- There were 32 solutions found in 0.0944 seconds (with a Bezout bound of 32).
load "testExampleSystems.m2"
J = ideal wright(QQ)
elapsedTime sols = zeroDimSolve J; -- 0.0204144
#sols -- 32
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))

(quotientBasis, S, signature, outputVariables) = ingredients J;
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 23.7229 seconds
sols = apply(oo, i -> point {i}); -- 10 solutions
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))



-- geneig
-- There were 10 solutions found in 1.447 seconds (with a Bezout bound of 243).
load "testExampleSystems.m2"
J = ideal geneig(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 3
elapsedTime sols = zeroDimSolve J; -- 0.0299402 seconds
#sols -- 10
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))

(quotientBasis, S, signature, outputVariables) = ingredients J;
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 15.9242
sols = apply(oo, i -> point {i}); -- 10 solutions
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))



-- redcyc5
load "testExampleSystems.m2"
R = QQ[z0,y1,y2,y3,y4]
J = ideal(1 + y1 + y2 + y3 + y4,
 y1 + y1*y2 + y2*y3 + y3*y4 + y4,
 y1*y2 + y1*y2*y3 + y2*y3*y4 + y3*y4 + y4*y1,
 y1*y2*y3 + y1*y2*y3*y4 + y2*y3*y4 + y3*y4*y1 + y4*y1*y2,
 z0**5*y1*y2*y3*y4  - 1)
elapsedTime sols = zeroDimSolve J; -- 0.0277107
#sols -- 14
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))

(quotientBasis, S, signature, outputVariables) = ingredients J;
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
sols2 = apply(oo, i -> point {i}); -- 14 solutions
-- 91.9528 seconds
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))



-----------------------
-- examples that are too large to compute or having positive dimension.
-----------------------


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

-- eco8
-- There were 64 solutions found in 9.197 seconds (with a Bezout bound of 1458).
load "testExampleSystems.m2"
J = ideal eco8(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 3
elapsedTime sols = zeroDimSolve J; -- 0.503384 seconds
#sols -- 64
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))

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

(quotientBasis, S, signature, outputVariables) = ingredients J;
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- took too long time on lin alg.

-- rabmo
-- There were 58 solutions found in 1642.07 seconds (with a Bezout bound of 36000).
load "testExampleSystems.m2"
J = ideal rabmo(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 5
elapsedTime sols = zeroDimSolve J; -- 0.433235
#sols -- 16 not 58
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))

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

-- virasoro
-- There were 256 solutions found in 0.356 seconds (with a Bezout bound of 256).
load "testExampleSystems.m2"
J = ideal virasoro(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 2
elapsedTime sols = zeroDimSolve J; -- 7.41958
#sols -- 256
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))

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



-- rose
-- There were 136 solutions found in 4.954 seconds (with a Bezout bound of 216).
load "testExampleSystems.m2"
J = ideal rose(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of J : 9
elapsedTime sols = zeroDimSolve J; -- 0.0676541
#sols -- 136
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))


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



-- bellido 
-- There were 40 solutions found in 0.546875 seconds (with a Bezout bound of 64).
load "testExampleSystems.m2"
J = ideal bellido(QQ)
maxIdeg = max (J_*/degree/sum); -- max degree of gens of I : 2
elapsedTime sols = zeroDimSolve J -- 0.0377411 seconds
#sols
sort apply(sols, p -> norm evaluate(gens J, p))
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))



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
-- took too long time on lin alg.



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



-----------------------
-- non radical examples
-----------------------

-- caprasse
-- There were 64 solutions found in 0.525506 seconds (with a Bezout bound of 144).
load "testExampleSystems.m2"
J = ideal caprasse(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of I : 4
elapsedTime sols = zeroDimSolve J; -- 0.0310881 seconds
#sols -- 56 not 64.
sort apply(sols, p -> norm evaluate(gens J, p)) -- 8 of them are incorrect solutions
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))


I = ideal gens gb J;
R = ring I;
quotientBasis = flatten entries sub(basis(R/I), R);
D = max ({1 + max(quotientBasis/degree/sum)}|(I_*|J_*)/degree/sum)  -- 7
degDmonoms = flatten entries matrix{for d to D list basis(d,R)};
S = degDmonoms;
lowerDegGens = select(flatten entries gens I, g -> sum degree g < maxIdeg);
outputVariables = gens R;
signature = flatten apply(degDmonoms,m->(apply(J_*|lowerDegGens,f->(m,f))));
--signature = flatten apply(degDmonoms,m->(apply(I_*,f->(m,f))));
elapsedTime sols = barebonesSolveRevised(J,S,signature,OutputVariables=>outputVariables,QuotientBasis=>quotientBasis) 
-- 93.5149
sols = apply(oo, i -> point {i}); -- 56 solutions


-- cassou
-- There were 2 solutions found in 39.0405 seconds (with a Bezout bound of 1344).
load "testExampleSystems.m2"
J = ideal cassou(QQ)
maxIdeg = max (J_*/degree/sum) -- max degree of gens of I : 8
elapsedTime sols = zeroDimSolve J; -- 0.0164434 seconds
#sols -- 16 not 2.
sort apply(sols, p -> norm evaluate(gens J, p)) 
apply(sols, p -> areEqual(0,norm evaluate(gens J, p)))


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

