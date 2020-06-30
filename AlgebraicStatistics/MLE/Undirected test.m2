restart;

--Generate score equations for an undirected graph on 4 nodes  with random observations
load "GraphicalModelsMLE.m2"
 G=graph{{1,2},{2,3},{3,4},{1,4}}
 R=gaussianRing(G)
 U=random(ZZ^4,ZZ^4)
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



