restart;

installPackage("GraphicalModelsMLERoser")
check GraphicalModelsMLERoser
help GraphicalModelsMLE


--EXAMPLE 2.1.13 OWL
restart
needsPackage("GraphicalModelsMLERoser")
 G=graph{{1,2},{2,3},{3,4},{1,4}}
 R=gaussianRing(G)
 U=random(ZZ^4,ZZ^4)
--Generate score equations for an undirected graph on 4 nodes  with random observations
J=scoreEquationsFromCovarianceMatrixUndir(R,U)
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
       
    
--EXAMPLE 3.2.1 OWL
restart
needsPackage("GraphicalModelsMLERoser")
 G=graph{{1,2},{2,4},{3,4},{2,3}}
 R=gaussianRing(G)
 undirectedEdgesMatrix(R)
 U=random(ZZ^4,ZZ^4)
--Generate score equations for an undirected graph on 4 nodes  with random observations
J=scoreEquationsFromCovarianceMatrixUndir(R,U)
dim J, degree J

--Find the optimal covariance matrix
MLEsolver(J,R)

--EXAMPLE 3.2.8 OWL: chain graph from figure 3.2.5(a) turned into an undirected graph
restart
needsPackage("GraphicalModelsMLERoser")
G=graph{{2,3},{3,4},{1,3}}
R=gaussianRing(G)
undirectedEdgesMatrix(R)
U=random(ZZ^4,ZZ^4)
J=scoreEquationsFromCovarianceMatrixUndir(R,U)
dim J, degree J
MLEsolver(J,R)

--EXAMPLE 3.2.11 OWL: chain graph from figure 3.2.5(b) turned into an undirected graph
restart
needsPackage("GraphicalModelsMLERoser")
G=graph{{1,2},{2,5},{5,6},{2,4},{4,5},{3,4}}
R=gaussianRing(G)
undirectedEdgesMatrix(R)
U=random(ZZ^6,ZZ^6)
J=scoreEquationsFromCovarianceMatrixUndir(R,U)
dim J, degree J
MLEsolver(J,R)


--EXAMPLE 3.3.7 OWL
restart
needsPackage("GraphicalModelsMLERoser")
G=graph{{1,2},{1,3},{1,4}}
R=gaussianRing(G)
undirectedEdgesMatrix(R)
U=random(ZZ^4,ZZ^4)
J=scoreEquationsFromCovarianceMatrixUndir(R,U)
dim J, degree J
MLEsolver(J,R)


--7.4 OWL conjecture on maximum likelihood of degree of Gaussian cycles
restart
needsPackage("GraphicalModelsMLERoser")
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

