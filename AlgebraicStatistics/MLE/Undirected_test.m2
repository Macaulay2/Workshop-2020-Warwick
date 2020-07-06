restart;

--Generate score equations for an undirected graph on 4 nodes  with random observations
load "GraphicalModelsMLE.m2"

-- Run internal functions for explicit checking
matZZtoQQ = (M) -> (
    E:=entries(M);    
    return matrix apply(#E, i -> apply(#(E_i), j -> (1/1)*E_i_j))    
);

    matRtolpR = (M,F) -> (
    E:=entries(M);    
    return matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)))    
);

--Examples and tests

 G=graph{{1,2},{2,3},{3,4},{1,4}}
 R=gaussianRing(G)
 U=random(ZZ^4,ZZ^4)
 ring(U)
  n := numRows U;
   -- converting it to list of matrix; rows of matrix correponds to the elements of the list
   X = for i to n-1 list U^{i};
   X
   
J=scoreEquationsFromCovarianceMatrixUndir(R,U)
J1 = scoreEquationsFromCovarianceMatrixUndir(R,X)

L=MLEsolver(J,R)
L1=MLEsolver(J1,R)
dim J
degree J

-- Test with rational matrix
 G=graph{{1,2},{2,3},{3,4},{1,4}}
 R=gaussianRing(G)
 U=random(QQ^4,QQ^4)
 ring(U)
J=scoreEquationsFromCovarianceMatrixUndir(R,U)

dim J
degree J

--Find the optimal concentration matrix
L=MLEsolver(J,R)

--Test MLE max
S=U*transpose(U);
n = #U;
MLEmax(R,L,S,n)


--Test which matrices in a list are PD
L={matrix{{1,0},{0,1}},matrix{{-2,0},{0,1}}};				
PDcheck(L)

--Do we need matRtolpR?
-- run lines 7-15 first!
G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
R=gaussianRing(G)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
V := sampleCovarianceMatrix(U);
 L := directedEdgesMatrix R;
 d := numRows L;
 -- Omega
    W := bidirectedEdgesMatrix R;
    -- move to a new ring, lpR, which does not have the s variables
    numSvars:=lift(d*(d+1)/2,ZZ);
    --lp rings is the ring without the s variables
    lpRvarlist:=apply(numgens(R)-numSvars,i->(gens(R))_i);
    KK:=coefficientRing(R);
    lpR:=KK[lpRvarlist];
    lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
    F:=map(lpR,R,lpRTarget);
    L1 = matRtolpR(L,F);
    L2=F(L);
    ring(L1)
    ring(L2)
ring(L1)===ring(L2)


-- Check why we need a new ring
 G=graph{{1,2},{2,3},{3,4},{1,4}}
 R=gaussianRing(G)
 U=random(ZZ^4,ZZ^4)
 -- run lines 7-10 first!
 if ring(U)===ZZ then U=matZZtoQQ(U);
 S=U*transpose(U);
 K=undirectedEdgesMatrix R;
 I=ideal{jacobian ideal{determinant(K)}-determinant(K)*jacobian(ideal{trace(K*S)})};
 J=saturate(I,ideal{determinant(K)});
 dim J
 degree J

-- Test changeRing on Mixed
 G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
 R = gaussianRing(G)
 U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
 J= scoreEquationsFromCovarianceMatrix(R,U)

 I=ideal(20*p_(3,4)+39,50*p_(4,4)-271,440104*p_(3,3)-742363,230*p_(2,2)-203,16*p_(1,1)-115,5*l_(3,4)+2,110026*l_(2,3)-2575,55013*l_(1,3)-600,115*l_(1,2)+26);
 J===I 
