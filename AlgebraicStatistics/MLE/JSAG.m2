--Example 4-cycle
restart
loadPackage "GraphicalModelsMLE";    
G=graph{{1,2},{2,3},{3,4},{4,1}};
--U=random(RR^4,RR^4);
--S=sampleCovarianceMatrix(U)
--toString S
S=matrix {{.105409, -.0745495, -.0186132, .0621907}, {-.0745495, .0783734, -.00844503,
     -.0872842}, {-.0186132, -.00844503, .128307, .0230245}, {.0621907, -.0872842, .0230245,
     .109849}};
solverMLE(G,S,SampleData=>false)

--Example with all types of multiple edges
restart
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}},graph {{1,2}});
S =  matrix {{7/20, 13/50, -3/50, -19/100}, {13/50, 73/100, -7/100, -9/100},
      {-3/50, -7/100, 2/5, 3/50},{-19/100, -9/100, 3/50, 59/100}};
solverMLE(G,S,SampleData=>false,ConcentrationMatrix=>true)


-- Example with multiple local maxima I: SUR (seemingly unrelated regressions) Drton 2004
restart
debug loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}});
S=matrix {{.68366, .0716539, .819476, .48852}, {.0716539, .224382, .309413, .45075}, {.819476, .309413, 1.23168, .971008}, {.48852, .45075, .971008, 1.12175}}
(J,Sinv)=scoreEquationsInternal(gaussianRing G,S,SampleData=>false);
sols=zeroDimSolve(J)
netList sols
--3 real solutions
M=genListMatrix(sols,Sinv)
L=checkPD(M)
--3 PD matrices

--Build Hessian matrix to determine whether they are local minima, maxima or saddle points
V=roundMatrix(53,S)
JAC=jacobianMatrixOfRationalFunction(det Sinv)*matrix{{1/det(Sinv)}}-jacobianMatrixOfRationalFunction(trace(V*Sinv))
H=for f in flatten entries JAC list flatten entries jacobianMatrixOfRationalFunction(f)
H=matrix H
entries H==entries transpose H
netList sols
MM=genListMatrix({sols_0,sols_1,sols_2},H)
for m in MM list eigenvalues m
-- 2 local maxima (all eigenvalues are negative) and 1 saddle point (both negative and positive) 



-- Example with multiple local maxima II: NO LUCK YET
restart
debug loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph{{1,2},{2,3},{2,4}},bigraph{{1,3},{2,4}});
U=random(RR^4,RR^4);
S=sampleCovarianceMatrix(U)
(J,Sinv)=scoreEquationsInternal(gaussianRing G,S,SampleData=>false);
dim J, degree J
-- (0,2)
sols=zeroDimSolve(J);
netList sols
--2 real solutions (all solutions are real)
M=genListMatrix(sols,Sinv);
L=checkPD(M)
--2 PD matrices
maxMLE(L,S)

--Build Hessian matrix to determine whether they are local minima, maxima or saddle points
V=roundMatrix(53,S);
JAC=jacobianMatrixOfRationalFunction(det Sinv)*matrix{{1/det(Sinv)}}-jacobianMatrixOfRationalFunction(trace(V*Sinv));
H=for f in flatten entries JAC list flatten entries jacobianMatrixOfRationalFunction(f);
H=matrix H;
entries H==entries transpose H;
netList sols;
MM=genListMatrix(sols,H);
for m in MM list eigenvalues m
-- 1 local maxima (all eigenvalues are negative) and 1 saddle point (both negative and positive) 

