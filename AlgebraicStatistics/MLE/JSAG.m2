uninstallPackage "Graphs"
uninstallPackage "StatGraphs"
uninstallPackage "GraphicalModels"
uninstallPackage "EigenSolver"
uninstallPackage "GraphicalModelsMLE"

installPackage "Graphs"
check Graphs
installPackage "StatGraphs"
check StatGraphs
installPackage "GraphicalModels"
check GraphicalModels
installPackage "EigenSolver"
check EigenSolver
installPackage "GraphicalModelsMLE"
check GraphicalModelsMLE


-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------
--                                                                                       --
--                          EXAMPLES IN JSAG PAPER                                       --
--                                                                                       --
-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------


-------------------------------------------------------------------------------------------
-- Undirected 4-cycle                                                                    --
-------------------------------------------------------------------------------------------
-- Connection to PD matrix completion problem

restart

-- Example 3.2
loadPackage "GraphicalModelsMLE";    
G=graph{{1,2},{2,3},{3,4},{4,1}}; 
S=matrix{{.105409, -.0745495, -.0186132, .0621907},{-.0745495, .0783734,-.00844503,-.0872842},{-.0186132, -.00844503, .128307, .0230245},{.0621907, -.0872842, .0230245,.109849}};
solverMLE(G,S,SampleData=>false)

-- Example 4.2
U=matrix{{3,5,9,5},{1,6,1,5},{2,9,6,6},{2,5,0,4}};
J=scoreEquations(gaussianRing G,U);
dim J

-- Example 5.2
MLdegree(gaussianRing G)


-------------------------------------------------------------------------------------------
-- Multiples edges                                                                       --
-------------------------------------------------------------------------------------------
-- Example with a multiple edge and positive dimension of score eq ideal

restart

--Example 4.4
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph{{1,3},{1,2},{2,4},{3,4}},graph{{1,2}});
R = gaussianRing G;
U = random(RR^4,RR^4);
J=scoreEquations(R,U);
dim J

--Example 5.3
MLdegree(gaussianRing G)


-------------------------------------------------------------------------------------------
-- Mixed graph with all type of edges                                                    --
-------------------------------------------------------------------------------------------
-- Example with 2 local maxima in the log-likelihood function

restart

--Example 3.3
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph{{1,3},{2,4}},bigraph {{3,4}},graph{{1,2}});
S=matrix {{34183/50000, 716539/10000000, 204869/250000, 12213/25000}, 
        {716539/10000000, 112191/500000, 309413/1000000, 1803/4000}, 
        {204869/250000, 309413/1000000, 3849/3125,15172/15625}, 
        {12213/25000, 1803/4000, 15172/15625, 4487/4000}};
solverMLE(G,S,SampleData=>false)

--Example 4.3
R = gaussianRing G;
(J,Sigma)=scoreEquations(R,S,SampleData=>false,CovarianceMatrix=>true);
dim J, degree J
sols=zeroDimSolve(J);netList sols
checkPD(apply(sols,i->sub(Sigma,matrix{coordinates(i)})))
-- compute Jacobian matrix (i.e. score equations)
scoreEq=-1/det Sigma*jacobianMatrixOfRationalFunction(det Sigma)-jacobianMatrixOfRationalFunction(trace(S*(inverse Sigma)));
-- compute Hessian matrix
Hessian=matrix for f in flatten entries scoreEq list flatten entries jacobianMatrixOfRationalFunction(f);
--compute eigenvalues of the Hessian matrix evaluated at real points in sols
apply({sols_0,sols_1,sols_4},i->eigenvalues sub(Hessian,matrix{coordinates(i)}))

-- Extra info (not in the paper)
MLdegree R

-- Example 6.1
restart
loadPackage "StatGraphs";
G = mixedGraph(digraph {{1,3},{2,4}},bigraph{{3,4}},graph{{1,2}});
partitionLMG G

-- Example 6.2
restart
loadPackage "GraphicalModels";
G = mixedGraph(digraph {{1,3},{2,4}},bigraph{{3,4}},graph{{1,2}});
R=gaussianRing G;
undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R
-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------


-- OLD TESTS 
restart
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}},graph {{1,2}});
S =  matrix {{7/20, 13/50, -3/50, -19/100}, {13/50, 73/100, -7/100, -9/100},
      {-3/50, -7/100, 2/5, 3/50},{-19/100, -9/100, 3/50, 59/100}};
scoreEquations(gaussianRing G,S);
dim oo, degree oo   
scoreEquations(gaussianRing G,S,SampleData=>false);
dim oo, degree oo   
solverMLE(G,S,SampleData=>false,ConcentrationMatrix=>true)



restart
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}},graph {{1,2}});
R = gaussianRing G;
MLdegree(R)
R=gaussianRing G
loadPackage "DeterminantalRepresentations"
scoreEquations(R,roundMatrix(6,random(RR^4,RR^4)));
dim oo, degree oo   
scoreEquations(R,roundMatrix(53,random(RR^4,RR^4)));
dim oo, degree oo   
scoreEquations(R,random(RR^4,RR^4));
dim oo, degree oo   



S=matrix {{.105409, -.0745495, -.0186132, .0621907}, {-.0745495, .0783734, -.00844503,
     -.0872842}, {-.0186132, -.00844503, .128307, .0230245}, {.0621907, -.0872842, .0230245,
     .109849}};
J=scoreEquations(R,S,SampleData=>false);
dim J, degree J


S=matrix {{34183/50000, 716539/10000000, 204869/250000, 12213/25000}, 
        {716539/10000000, 112191/500000, 309413/1000000, 1803/4000}, 
        {204869/250000, 309413/1000000, 3849/3125,15172/15625}, 
        {12213/25000, 1803/4000, 15172/15625, 4487/4000}};
J=scoreEquations(R,S,SampleData=>false);
dim J, degree J



loadPackage "DeterminantalRepresentations"
scoreEquations(R,roundMatrix(6,random(RR^4,RR^4)));
dim oo, degree oo   
scoreEquations(R,roundMatrix(53,random(RR^4,RR^4)));
dim oo, degree oo   
scoreEquations(R,random(RR^4,RR^4));
dim oo, degree oo   

random(RR^4,RR^4)
roundMatrix(53,oo)



S=matrix {{34183/50000, 716539/10000000, 204869/250000, 12213/25000}, 
        {716539/10000000, 112191/500000, 309413/1000000, 1803/4000}, 
        {204869/250000, 309413/1000000, 3849/3125,15172/15625}, 
        {12213/25000, 1803/4000, 15172/15625, 4487/4000}};
--Compute critical points
J=scoreEquations(R,S,SampleData=>false);
dim J, degree J
sols=zeroDimSolve(J);netList sols
--Compute Sigma
K=undirectedEdgesMatrix R;L = directedEdgesMatrix R,P = bidirectedEdgesMatrix R;
FR=frac(coefficientRing(R)[gens R-set support covarianceMatrix R]);
IdL = inverse (id_(FR^4)-sub(L,FR));W=directSum(inverse sub(K,FR),sub(P,FR));
Sigma = (transpose IdL) * W * IdL
--check which solutions give PD matrices
PDmatrices=checkPD(apply(sols,i->sub(Sigma,matrix{coordinates(i)})))
--check what kind of critical point is each solution
scoreEq=-1/det Sigma*jacobianMatrixOfRationalFunction(det Sigma)- jacobianMatrixOfRationalFunction(trace(S*(inverse Sigma)));
Hessian=matrix for f in flatten entries scoreEq list flatten entries jacobianMatrixOfRationalFunction(f);
apply({sols_0,sols_3,sols_4},i->eigenvalues sub(Hessian,matrix{coordinates(i)}))

solverMLE(G,S,SampleData=>false)

restart
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}},graph {{1,2}});
R = gaussianRing G;
S=matrix {{34183/50000, 716539/10000000, 204869/250000, 12213/25000}, 
        {716539/10000000, 112191/500000, 309413/1000000, 1803/4000}, 
        {204869/250000, 309413/1000000, 3849/3125,15172/15625}, 
        {12213/25000, 1803/4000, 15172/15625, 4487/4000}};
J=scoreEquations(R,S,SampleData=>false);
dim J, degree J
sols=zeroDimSolve(J);netList sols

K=undirectedEdgesMatrix R;L = directedEdgesMatrix R,P = bidirectedEdgesMatrix R;
FR=frac(coefficientRing(R)[gens R-set support covarianceMatrix R]);
IdL = inverse (id_(FR^4)-sub(L,FR));
W=directSum(inverse sub(K,FR),sub(P,FR));
Sinv = inverse ((transpose IdL) * W * IdL)

Sigma = (transpose IdL) * W * IdL

L=checkPD(apply(sols,i->sub(Sinv,matrix{coordinates(i)})))
LSigma=checkPD(apply(sols,i->sub(Sigma,matrix{coordinates(i)})))

partial=jacobianMatrixOfRationalFunction(det Sinv)*matrix{{1/det(Sinv)}}-jacobianMatrixOfRationalFunction(trace(S*Sinv))
partialSigma=-jacobianMatrixOfRationalFunction(det Sigma)*matrix{{1/det Sigma}}-jacobianMatrixOfRationalFunction(trace(S*(inverse Sigma)))


H=matrix for f in flatten entries partial list flatten entries jacobianMatrixOfRationalFunction(f)

HSigma=matrix for f in flatten entries partialSigma list flatten entries jacobianMatrixOfRationalFunction(f)


realsols={sols_0,sols_1,sols_4}
netList realsols
apply(realsols,i->eigenvalues sub(H,matrix{coordinates(i)}))


apply(realsols,i->eigenvalues sub(HSigma,matrix{coordinates(i)}))


solverMLE(G,S,SampleData=>false)



restart
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}});

loadPackage "DeterminantalRepresentations"
R=gaussianRing G;
scoreEquations(R,roundMatrix(6,random(RR^4,RR^4)));
dim oo, degree oo   
U=roundMatrix(6,random(RR^4,RR^4))

restart
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}});
U=matrix {{98899/1000000, 703827/1000000, 79253/200000, 577943/1000000}, {85811/100000, 94273/125000, 154629/250000, 94847/125000}, {110481/125000, 1499/100000, 509767/1000000,
      36417/40000}, {23673/62500, 273067/500000, 181627/200000, 16943/40000}}
J=scoreEquations(gaussianRing G,U);
dim oo, degree oo   
solverMLE(G,U)
sols=zeroDimSolve J
netList sols

restart
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}});
U=random(QQ^4,QQ^4)
J=scoreEquations(gaussianRing G,U);
dim oo, degree oo   
solverMLE(G,U)
sols=zeroDimSolve J
netList sols

restart
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}});
loadPackage "DeterminantalRepresentations"
U=roundMatrix(53,random(RR^4,RR^4))
J=scoreEquations(gaussianRing G,U);
dim oo, degree oo   
netList zeroDimSolve J
solverMLE(G,U)

restart
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}});
loadPackage "DeterminantalRepresentations"
U=random(RR^4,RR^4);
V=roundMatrix(53,sampleCovarianceMatrix(U));
J=scoreEquations(gaussianRing G,U);
dim oo, degree oo   
netList zeroDimSolve J
solverMLE(G,U)

restart
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}});
loadPackage "DeterminantalRepresentations"
U=random(RR^4,RR^4)
V=roundMatrix(4,sampleCovarianceMatrix(U));
J=scoreEquations(gaussianRing G,V,SampleData=>false);
dim oo, degree oo   
netList zeroDimSolve J
solverMLE(G,U)




scoreEquations(R,random(RR^4,RR^4));
dim oo, degree oo   

R = gaussianRing G;L = directedEdgesMatrix R,P = bidirectedEdgesMatrix R;
FR=frac(coefficientRing(R)[gens R-set support covarianceMatrix R]);
IdL = inverse (id_(FR^4)-sub(L,FR));
Sinv = inverse ((transpose IdL) * sub(P,FR) * IdL)
S=matrix {{34183/50000, 716539/10000000, 204869/250000, 12213/25000}, 
        {716539/10000000, 112191/500000, 309413/1000000, 1803/4000}, 
        {204869/250000, 309413/1000000, 3849/3125,15172/15625}, 
        {12213/25000, 1803/4000, 15172/15625, 4487/4000}};
J=scoreEquations(R,S,SampleData=>false);
dim J, degree J
sols=zeroDimSolve(J);netList sols
L=checkPD(apply(sols,i->sub(Sinv,matrix{coordinates(i)})))
partial=jacobianMatrixOfRationalFunction(det Sinv)*(1/det(Sinv))-jacobianMatrixOfRationalFunction(trace(S*Sinv))
H=matrix for f in flatten entries partial list flatten entries jacobianMatrixOfRationalFunction(f);
apply({sols_2,sols_3},i->eigenvalues sub(H,matrix{coordinates(i)}))



restart
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}},graph {{1,2}});
S=matrix {{.68366, .0716539, .819476, .48852}, {.0716539, .224382, .309413, .45075}, {.819476, .309413, 1.23168, .971008}, {.48852, .45075, .971008, 1.12175}}
(J,Sinv)=scoreEquationsInternal(gaussianRing G,S,SampleData=>false);
dim J,degree J
sols=zeroDimSolve(J)
netList sols
--3 real solutions
M=genListMatrix(sols,Sinv)
L=checkPD(M)
--3 PD matrices
maxMLE({L_2},S)
maxMLE({L_0},S)
maxMLE({L_1},S)

--Build Hessian matrix to determine whether they are local minima, maxima or saddle points
V=roundMatrix(53,S)
JAC=jacobianMatrixOfRationalFunction(det Sinv)*matrix{{1/det(Sinv)}}-jacobianMatrixOfRationalFunction(trace(V*Sinv))
H=for f in flatten entries JAC list flatten entries jacobianMatrixOfRationalFunction(f)
H=matrix H
entries H==entries transpose H
netList sols
MM=genListMatrix({sols_0,sols_3,sols_4},H)
for m in MM list eigenvalues m
-- 2 local maxima (all eigenvalues are negative) and 1 saddle point (both negative and positive) 



-- Example with multiple local maxima I: SUR (seemingly unrelated regressions) Drton 2004
restart
debug loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}});
S=matrix {{.68366, .0716539, .819476, .48852}, {.0716539, .224382, .309413, .45075}, {.819476, .309413, 1.23168, .971008}, {.48852, .45075, .971008, 1.12175}}
(J,Sinv)=scoreEquationsInternal(gaussianRing G,S,SampleData=>false);
dim J,degree J
sols=zeroDimSolve(J)
netList sols
--3 real solutions
M=genListMatrix(sols,Sinv)
L=checkPD(M)
--3 PD matrices
maxMLE({L_2},S)
maxMLE({L_0},S)
maxMLE({L_1},S)

--Build Hessian matrix to determine whether they are local minima, maxima or saddle points
V=roundMatrix(53,S)
JAC=jacobianMatrixOfRationalFunction(det Sinv)*matrix{{1/det(Sinv)}}-jacobianMatrixOfRationalFunction(trace(V*Sinv))
H=for f in flatten entries JAC list flatten entries jacobianMatrixOfRationalFunction(f)
H=matrix H
entries H==entries transpose H
netList sols
MM=genListMatrix({sols_0,sols_3,sols_4},H)
for m in MM list eigenvalues m
-- 2 local maxima (all eigenvalues are negative) and 1 saddle point (both negative and positive) 

restart
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}});
R = gaussianRing G;L = directedEdgesMatrix R,P = bidirectedEdgesMatrix R;
FR=frac(coefficientRing(R)[gens R-set support covarianceMatrix R]);
IdL = inverse (id_(FR^4)-sub(L,FR));
Sinv = inverse ((transpose IdL) * sub(P,FR) * IdL)
U=random(RR^4,RR^4);
J=scoreEquations(R,S,SampleData=>false);
dim J, degree J
sols=zeroDimSolve(J);netList sols
L=checkPD(apply(sols,i->sub(Sinv,matrix{coordinates(i)})))
partial=jacobianMatrixOfRationalFunction(det Sinv)*matrix{{1/det(Sinv)}}-jacobianMatrixOfRationalFunction(trace(S*Sinv))
H=matrix for f in flatten entries partial list flatten entries jacobianMatrixOfRationalFunction(f)
realsols={sols_2,sols_3,sols_4}
apply(realsols,i->eigenvalues sub(H,matrix{coordinates(i)}))


restart
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}});
R = gaussianRing G;L = directedEdgesMatrix R,P = bidirectedEdgesMatrix R;
FR=frac(coefficientRing(R)[gens R-set support covarianceMatrix R]);
IdL = inverse (id_(FR^4)-sub(L,FR));
Sinv = inverse ((transpose IdL) * sub(P,FR) * IdL)
S=matrix {{34183/50000, 716539/10000000, 204869/250000, 12213/25000}, 
        {716539/10000000, 112191/500000, 309413/1000000, 1803/4000}, 
        {204869/250000, 309413/1000000, 3849/3125,15172/15625}, 
        {12213/25000, 1803/4000, 15172/15625, 4487/4000}};
J=scoreEquations(R,S,SampleData=>false);
dim J, degree J
sols=zeroDimSolve(J);netList sols
L=checkPD(apply(sols,i->sub(Sinv,matrix{coordinates(i)})));
partial=jacobianMatrixOfRationalFunction(det Sinv)*matrix{{1/det(Sinv)}}-jacobianMatrixOfRationalFunction(trace(S*Sinv));
H=matrix for f in flatten entries partial list flatten entries jacobianMatrixOfRationalFunction(f);
apply({sols_0,sols_1,sols_4},i->eigenvalues sub(H,matrix{coordinates(i)}))


-- Example without MLE
restart
loadPackage "GraphicalModelsMLE";    
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}});
R = gaussianRing G;L = directedEdgesMatrix R,P = bidirectedEdgesMatrix R;
FR=frac(coefficientRing(R)[gens R-set support covarianceMatrix R]);
IdL = inverse (id_(FR^4)-sub(L,FR));
Sinv = inverse ((transpose IdL) * sub(P,FR) * IdL)
S=matrix {{13067/1600, -1147/400, -229/320, 1221/400}, {-1147/400, 2779/200, 271/320, -707/1200}, {-229/320, 271/320, 37/384, -587/2880}, {1221/400, -707/1200, -587/2880, 863/600}}
J=scoreEquations(R,S,SampleData=>false);
dim J, degree J
sols=zeroDimSolve(J);netList sols
solverMLE(G,S,SampleData=>false)
MLdegree(R)



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




-- 4-cycle example
G=graph{{1,2},{2,3},{3,4},{4,1}};
--U=random(RR^4,RR^4);
--S=sampleCovarianceMatrix(U)
--toString S
S=matrix {{.105409, -.0745495, -.0186132, .0621907}, {-.0745495, .0783734, -.00844503,
     -.0872842}, {-.0186132, -.00844503, .128307, .0230245}, {.0621907, -.0872842, .0230245,
     .109849}};
solverMLE(G,S,SampleData=>false)

R=gaussianRing G
loadPackage "DeterminantalRepresentations"
scoreEquations(R,roundMatrix(6,random(RR^4,RR^4)));
dim oo, degree oo   
scoreEquations(R,roundMatrix(53,random(RR^4,RR^4)));
dim oo, degree oo   
scoreEquations(R,random(RR^4,RR^4));
dim oo, degree oo   

