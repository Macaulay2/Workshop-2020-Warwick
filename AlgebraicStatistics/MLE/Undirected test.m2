restart;

load "GraphicalModelsMLE.m2"
 G=graph{{1,2},{2,3},{3,4},{1,4}}
 R=gaussianRing(G)
 U=random(ZZ^4,ZZ^4)
J=scoreEquationsFromCovarianceMatrixUndir(R,U)
dim J
degree J

    S:=U*transpose(U)
    --S=matZZtoQQ(V)
    --V := sampleCovarianceMatrix(U);
    -- Concentration matrix K
    K=undirectedEdgesMatrix R;
    -- d is equal to the number of vertices in G
    d := numRows K;
    -- move to a new ring, lpR, which does not have the s variables
    numSvars:=lift(d*(d+1)/2,ZZ);
     --lp ring is the ring without the s variables
    lpRvarlist:=apply(numgens(R)-numSvars,i->(gens(R))_i);
    KK:=coefficientRing(R);
    lpR:=KK[lpRvarlist];
    lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
    F:=map(lpR,R,lpRTarget);
    K=F(K)
    I=ideal{jacobian ideal{determinant(K)}-determinant(K)*jacobian(ideal{trace(K*S)})}
    J=saturate(I,ideal{determinant(K)})

degree J
