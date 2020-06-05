restart;

load "GraphicalModelsMLE.m2"
 G=graph{{1,2},{2,3},{3,4},{1,4}}
 R=gaussianRing(G)
 U=random(ZZ^4,ZZ^4)
J=scoreEquationsFromCovarianceMatrixUndir(R,U)
dim J
degree J

MLEsolver(J,R)



