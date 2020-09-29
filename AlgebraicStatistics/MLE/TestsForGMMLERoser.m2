restart;

installPackage("GraphicalModelsMLE")
check GraphicalModelsMLE
help GraphicalModelsMLE

----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------

--EXAMPLE 2.1.13 OWL
restart
needsPackage("GraphicalModelsMLE")
G=graph{{1,2},{2,3},{3,4},{1,4}}
R=gaussianRing(G)
-- U=random(ZZ^4,ZZ^4)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
--Generate score equations for an undirected graph on 4 nodes  with random observations
J=time scoreEquationsFromCovarianceMatrixUndir(R,U)
-- used 0.125 seconds
dim J, degree J
--Find the optimal covariance matrix
MLEsolver(J,R)

--Step by step
    --solve system with eigensolver
    sols=zeroDimSolve(J);
    --compute concentration matrix
    K=undirectedEdgesMatrix(R)
    --evaluate concentration matrix on solutions
    M=evalConcentrationMatrix(sols,R);
    --consider only PD matrices
    L=PDcheck M
    E=inverse L_0
       

-- test for undirected as  mixedGraph
restart
needsPackage("GraphicalModelsMLE")
G=graph{{1,2},{2,3},{3,4},{1,4}}
g=mixedGraph G
R2=gaussianRing(g)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
J2=time scoreEquationsFromCovarianceMatrix(R2,U)
--with usual saturation
-- used 11.1563 seconds
--iterating saturation (denoms)
-- used 0.421875 seconds
J2=time scoreEquationsFromCovarianceMatrix(R2,U,doSaturate =>false)
-- used 0.3125 seconds
dim J2, degree J2
test=ideal{58*k_(4,4)+47*k_(1,4)-21*k_(3,4)-8, 27*k_(3,3)+14*k_(2,3)-42*k_(3,4)-16,
     10*k_(2,2)-13*k_(1,2)+7*k_(2,3)-8, 115*k_(1,1)-26*k_(1,2)+94*k_(1,4)-16,
     49587264*k_(3,4)^2+159578097*k_(1,2)-401770875*k_(1,4)+86063425*k_(2,3)-241038399*k_(3,4)-
     9279488, 289984*k_(2,3)*k_(3,4)-2996077*k_(1,2)+3236687*k_(1,4)-1267133*k_(2,3)+2001475*k_(3
     ,4), 289984*k_(1,4)*k_(3,4)+572663*k_(1,2)-1267133*k_(1,4)+247207*k_(2,3)-713625*k_(3,4),
     289984*k_(1,2)*k_(3,4)-1267133*k_(1,2)+1588223*k_(1,4)-634637*k_(2,3)+786099*k_(3,4),
     12469312*k_(2,3)^2+159578097*k_(1,2)-401770875*k_(1,4)+94182977*k_(2,3)-216679743*k_(3,4)-
     9279488, 289984*k_(1,4)*k_(2,3)-1428105*k_(1,2)+2001475*k_(1,4)-713625*k_(2,3)+983079*k_(3,4
     ), 289984*k_(1,2)*k_(2,3)+2001475*k_(1,2)-3960705*k_(1,4)+786099*k_(2,3)-2523789*k_(3,4),
     163260992*k_(1,4)^2+159578097*k_(1,2)-347253883*k_(1,4)+86063425*k_(2,3)-216679743*k_(3,4)-
     9279488, 289984*k_(1,2)*k_(1,4)-713625*k_(1,2)+786099*k_(1,4)-302505*k_(2,3)+482391*k_(3,4),
     58866752*k_(1,2)^2+144498929*k_(1,2)-401770875*k_(1,4)+86063425*k_(2,3)-216679743*k_(3,4)-
     9279488};

J2==test
--Find the optimal covariance matrix
MLEsolver(J2,R2)

----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
    
--EXAMPLE 3.2.1 OWL
restart
needsPackage("GraphicalModelsMLE")
 G=graph{{1,2},{2,4},{3,4},{2,3}}
 R=gaussianRing(G)
-- U=random(ZZ^4,ZZ^4)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
--Generate score equations for an undirected graph on 4 nodes  with random observations
J=time scoreEquationsFromCovarianceMatrixUndir(R,U)
dim J, degree J
aux=toExternalString J
test=ideal{8*k_(3,4)-7, 8*k_(2,3)-7, 32*k_(2,4)-75, 203*k_(1,2)-52, 32*k_(4,4)-43, 2*k_(3,3)-3, 32480*k_(2,2)-184381, 203*k_(1,1)-40}

-- test for undirected as  mixedGraph
restart
needsPackage("GraphicalModelsMLE")
G=graph{{1,2},{2,4},{3,4},{2,3}}
g=mixedGraph G
R2=gaussianRing(g)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
J2=time scoreEquationsFromCovarianceMatrix(R2,U)
-- used 0.328125 seconds
-- used 0.234375 seconds
J2=time scoreEquationsFromCovarianceMatrix(R2,U,doSaturate =>false)
-- used 0.1875 seconds
dim J2, degree J2
test=ideal{8*k_(3,4)-7, 8*k_(2,3)-7, 32*k_(2,4)-75, 203*k_(1,2)-52, 32*k_(4,4)-43, 2*k_(3,3)-3, 32480*k_(2,2)-184381, 203*k_(1,1)-40}
J2==test
--Find the optimal covariance matrix
MLEsolver(J2,R2)


----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------


--EXAMPLE 3.2.8 OWL: chain graph from figure 3.2.5(a) turned into an undirected graph
restart
needsPackage("GraphicalModelsMLE")
G=graph{{2,3},{3,4},{1,3}}
R=gaussianRing(G)
undirectedEdgesMatrix(R)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
--U=random(ZZ^4,ZZ^4)
J=time scoreEquationsFromCovarianceMatrixUndir(R,U)
 -- used 0.03125 seconds
dim J, degree J
aux=toExternalString J
test=ideal{57*k_(3,4)-28, 43*k_(2,3)+28, 283*k_(1,3)-58, 19*k_(4,4)-6, 6242697*k_(3,3)-11951602, 43*k_(2,2)-54, 283*k_(1,1)-54}
J==test

MLEsolver(J,R)

-- test for undirected as  mixedGraph
restart
needsPackage("GraphicalModelsMLE")
G=graph{{2,3},{3,4},{1,3}}
g=mixedGraph(G)
R2=gaussianRing(g)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
--U=random(ZZ^4,ZZ^4)
J2=time scoreEquationsFromCovarianceMatrix(R2,U)
-- used 0.078125 seconds
dim J2, degree J2
test=ideal{57*k_(3,4)-28, 43*k_(2,3)+28, 283*k_(1,3)-58, 19*k_(4,4)-6, 6242697*k_(3,3)-11951602, 43*k_(2,2)-54, 283*k_(1,1)-54}
J2==test


--EXAMPLE 3.2.11 OWL: chain graph from figure 3.2.5(b) turned into an undirected graph
restart
needsPackage("GraphicalModelsMLE")
G=graph{{1,2},{2,5},{5,6},{2,4},{4,5},{3,4}}
R=gaussianRing(G)
undirectedEdgesMatrix(R)
--U=random(ZZ^6,ZZ^6)
--auxU=toExternalString U
U=matrix{{1, 2, 9, 6, 0, 0}, {2, 7, 7, 3, 2, 2}, {6, 3, 4, 1, 5, 5}, {5, 5, 8, 8, 7, 6}, {3, 2, 3, 8, 7, 5}, {8, 0, 5, 3, 8, 5}}
J=time scoreEquationsFromCovarianceMatrixUndir(R,U)
-- used 0.09375 seconds
dim J, degree J
aux=toExternalString J
test=ideal{452*k_(5,6)+627, 28733*k_(4,5)+639, 1703*k_(3,4)+72, 28733*k_(2,4)+309,
       28733*k_(2,5)-1781, 524*k_(1,2)-51, 452*k_(6,6)-915, 3961131380*k_(5,5)-4311459839, 12575600843*k_(4,4)-1906007136, 3406*k_(3,3)-771, 557075404*k_(2,2)-148839607, 524*k_(1,1)-111}
J==test

MLEsolver(J,R)

-- test for undirected as  mixedGraph
restart
needsPackage("GraphicalModelsMLE")
G=graph{{1,2},{2,5},{5,6},{2,4},{4,5},{3,4}}
g=mixedGraph G
R2=gaussianRing(g)
U=matrix{{1, 2, 9, 6, 0, 0}, {2, 7, 7, 3, 2, 2}, {6, 3, 4, 1, 5, 5}, {5, 5, 8, 8, 7, 6}, {3, 2, 3, 8, 7, 5}, {8, 0, 5, 3, 8, 5}}
J2=time scoreEquationsFromCovarianceMatrix(R2,U,doSaturate =>false)

J2=time scoreEquationsFromCovarianceMatrix(R2,U)
-- COMPUTATION DIDN'T FINISH IN A SHORT AMOUNT OF TIME
-- used 4.53125 seconds

dim J2, degree J2
test=ideal{452*k_(5,6)+627, 28733*k_(4,5)+639, 1703*k_(3,4)+72, 28733*k_(2,4)+309,
       28733*k_(2,5)-1781, 524*k_(1,2)-51, 452*k_(6,6)-915, 3961131380*k_(5,5)-4311459839, 12575600843*k_(4,4)-1906007136, 3406*k_(3,3)-771, 557075404*k_(2,2)-148839607, 524*k_(1,1)-111}
J2==test

JnoSat=time scoreEquationsFromCovarianceMatrix(R2,U,doSaturate=>false);

restart
debug loadPackage("GraphicalModelsMLE")

R=R2
--   R := gaussianRing(G);  
    ----------------------------------------------------
    -- Extract information about the graph
    ---------------------------------------------------- 
    -- Lambda
    L = directedEdgesMatrix R
    -- K 
    K = undirectedEdgesMatrix R
    -- Psi
    P = bidirectedEdgesMatrix R
    
    ----------------------------------------------------
    -- Create an auxiliary ring and its fraction field
    -- which do not have the s variables
    ----------------------------------------------------
    -- d is equal to the number of vertices in G
    d = numRows L
    -- create a new ring, lpR, which does not have the s variables
    (F,lpR)=changeRing(d,R)
    -- create its fraction field
    FR = frac(lpR)
    
    -----------------------------------------------------
    -- Construct Omega
    -----------------------------------------------------
    -- Kinv
    K=substitute(K, FR)
    Kinv=inverse K
    P=substitute(P,FR)
       
     --Omega
    if K==0 then W:=P else (if P==0 then W=Kinv else W = directSum(Kinv,P));
    --W:= directSum(Kinv,P);
    
    -- move to FR, the fraction field of lpR
    L= substitute (L,FR)
    
    -- Sigma
    if L==0 then S:=W else (
	IdL := inverse (id_(FR^d)-L);
    	S = (transpose IdL) * W * IdL
	);
    if S == Kinv then Sinv:= K else Sinv = inverse S; 
    
    -- Sample covariance matrix
    V = sampleCovarianceMatrix(U);
     
    -- Compute ideal J   
    C1 = trace(Sinv * V)
    C1derivative = JacobianMatrixOfRationalFunction(C1)
    LL=JacobianMatrixOfRationalFunction (det Sinv)*matrix{{1/det(Sinv)}} - (C1derivative)
--  LL=JacobianMatrixOfRationalFunction (det S)*matrix{{-1/det(S)}} - (C1derivative)
    LL=flatten entries(LL)
    denoms = apply(#LL, i -> lift(denominator(LL_i), lpR))
--    prod = product(denoms);
    J=ideal apply(#LL, i -> lift(numerator(LL_i),lpR))
    for i from 0 to (#denoms-1) do (J=time saturate(J,denoms_i)); 
    
    netList J_*


    LL2=JacobianMatrixOfRationalFunction(det Sinv)-det(Sinv)*C1derivative;
    LL2=flatten entries(LL2);
    J=ideal{LL2};
    
    
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------

--EXAMPLE 3.3.7 OWL
restart
needsPackage("GraphicalModelsMLE")
G=graph{{1,2},{1,3},{1,4}}
R=gaussianRing(G)
undirectedEdgesMatrix(R)
--U=random(ZZ^4,ZZ^4)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
J=scoreEquationsFromCovarianceMatrixUndir(R,U)
dim J, degree J
aux=toExternalString J
test=ideal{563*k_(1,4)+188, 283*k_(1,3)-58, 203*k_(1,2)-52, 563*k_(4,4)-230, 283*k_(3,3)-230,
     203*k_(2,2)-230, 3719535505*k_(1,1)-1940386226}
J==test
MLEsolver(J,R)

-- test for undirected as  mixedGraph
restart
needsPackage("GraphicalModelsMLE")
G=graph{{1,2},{1,3},{1,4}}
g=mixedGraph G
R2=gaussianRing(g)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
J2=scoreEquationsFromCovarianceMatrix(R2,U)
dim J2, degree J2
test=ideal{563*k_(1,4)+188, 283*k_(1,3)-58, 203*k_(1,2)-52, 563*k_(4,4)-230, 283*k_(3,3)-230,
     203*k_(2,2)-230, 3719535505*k_(1,1)-1940386226}
J2==test


----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------


--7.4 OWL conjecture on maximum likelihood of degree of Gaussian cycles
restart
needsPackage("GraphicalModelsMLE")
m=6;
L={{1,m}};
for l from 1 to m-1 do(L=L|{{l,l+1}});
G=graph(L)
R=gaussianRing(G)
undirectedEdgesMatrix(R)
U=random(ZZ^m,ZZ^m)
J=scoreEquationsFromCovarianceMatrixUndir(R,U)
dim J, degree J
MLEsolver(J,R)

conj=(m-3)*2^(m-2)+1
degree J===conj

--Could be computed up to m=7. For m=8 even in the MPI server the computation is killed.