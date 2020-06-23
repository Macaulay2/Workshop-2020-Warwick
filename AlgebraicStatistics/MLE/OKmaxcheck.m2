restart;
load "GraphicalModelsMLE.m2";

G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
R = gaussianRing(G)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
J=scoreEquationsFromCovarianceMatrix(R,U);
MLEsolver(J,R);
S=U*transpose(U);
n = #U;
MLEmax(R,L,S,n)
