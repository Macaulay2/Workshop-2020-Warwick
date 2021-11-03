----------------------------------
-- Correspondence variety
---------------------------------
restart
debug loadPackage "GraphicalModelsMLE";
G = mixedGraph(graph{{1,2}},digraph{{1,3},{2,4}},bigraph{{3,4}});
R = gaussianRing G;
L = directedEdgesMatrix R;
P = bidirectedEdgesMatrix R;
K = undirectedEdgesMatrix R;
sup=set support covarianceMatrix R
R=R**QQ[v_(11)..v_(14),v_(21)..v_(24),v_(31)..v_(34),v_(41)..v_(44),t]

lpR=coefficientRing(R)[gens R-sup];
FR = frac(lpR);
K=sub(K, FR);
Kinv=inverse K;
P=sub(P,FR);
W = directSum(Kinv,P)
L= sub(L,FR);
d=numcols L;
IdL := inverse (id_(FR^d)-L);
S = (transpose IdL) * W * IdL
Sinv = inverse S;

--FR1=(ring Sinv)[v_(11)..v_(14),v_(22)..v_(24),v_(33),v_(34),v_(44)]
Sinv=sub(Sinv,FR)
V=matrix{{34183/50000,716539/1000000,204869/250000,v_(14)},{716539/1000000, 112191/500000,v_(23), 1803/4000},
    {204869/250000,v_(23),3849/3125,15172/15625},{v_(14),1803/4000,15172/15625,4487/4000}}
--V=matrix {{34183/50000, 716539/10000000, 204869/250000, 12213/25000}, 
--              {716539/10000000, 112191/500000, 309413/1000000, 1803/4000}, 
--              {204869/250000, 309413/1000000, 3849/3125,15172/15625}, 
--              {12213/25000, 1803/4000, 15172/15625, 4487/4000}};
C1 = trace(Sinv * V);
C1derivative = jacobianMatrixOfRationalFunction(C1);
LL =jacobianMatrixOfRationalFunction (det Sinv)*matrix{{1/det(Sinv)}} - C1derivative;
LL=flatten entries(LL);
denoms = apply(#LL, i -> lift(denominator(LL_i), lpR));
J=(ideal apply(#LL, i -> lift(numerator(LL_i),lpR)))
for i from 0 to (#denoms-1) do (
if degree denoms_i =={0} then J=J else  
J=saturate(J,denoms_i))
Imax=ideal(det Sinv- trace (V*Sinv)-t)
JImax=J+Imax

--listVar= {k_(1,1),k_(2,2),k_(1,2),l_(1,3),l_(2,4),p_(3,3),p_(4,4),p_(3,4),s_(1,1),s_(1,2),s_(1,3),s_(1,4),s_(2,2),s_(2,3),s_(2,4),s_(3,3),s_(3,4),s_(4,4)}
--listVar= {k_(1,1),k_(2,2),k_(1,2),l_(1,3),l_(2,4),p_(3,3),p_(4,4),p_(3,4)}
listVar= {k_(1,1),k_(2,2),k_(1,2),l_(1,3),l_(2,4),p_(3,3),p_(4,4),p_(3,4)}
Iaux= sub(ideal listVar, lpR)
Jnew=eliminate( flatten entries gens Iaux,J)
JnewImax=Jnew+Imax

jiR=coefficientRing(FR)[ support JnewImax]
subJnewImax=sub(JnewImax,jiR)

Iaux2=sub(ideal(v_(14),v_(23)),lpR )
Jaux=eliminate(flatten entries gens Iaux2,J)

matList={matrix {{.68366, .0716539, 1.00282, .234375}, {.0716539, .224382, .105105, .733937}, {1.00282, .105105, 1.76955, -.0700599},
      {.234375, .733937, -.0700599, 2.97432}}, matrix {{.68366, .0716539, .467724, .070198}, {.0716539, .224382, .0490218, .219823},
      {.467724, .0490218, .75038, .429714}, {.070198, .219823, .429714, .66928}}}
mat1=matList_0
mat1=roundMatrix(10,mat1)
mat2=matList_1
mat2=roundMatrix(10,mat2)
eq1=-13574062/1000000- trace (V*mat1)-936624/100000
eq2=-4910951/1000000- trace (V*mat2)-936624/100000

Jopt=Jnew+ideal(eq1)+ideal(eq2)
--for Sinv in L list log det Sinv- trace (V*Sinv)
----------------------------------
-- Original example
----------------------------------
restart
debug loadPackage "GraphicalModelsMLE";
G = mixedGraph(graph{{1,2}},digraph{{1,3},{2,4}},bigraph{{3,4}});
S=matrix {{34183/50000, 716539/10000000, 204869/250000, 12213/25000}, 
              {716539/10000000, 112191/500000, 309413/1000000, 1803/4000}, 
              {204869/250000, 309413/1000000, 3849/3125,15172/15625}, 
              {12213/25000, 1803/4000, 15172/15625, 4487/4000}};
solverMLE(G,S,SampleData=>false)
R = gaussianRing G;
(J,Sigma)=scoreEquations(R,S,SampleData=>false,CovarianceMatrix=>true);
dim J, degree J
sols=zeroDimSolve(J);netList sols
checkPD(apply(sols,i->sub(Sigma,matrix{coordinates(i)})))
scoreEq=-1/det Sigma*jacobianMatrixOfRationalFunction(det Sigma)-
jacobianMatrixOfRationalFunction(trace(S*(inverse Sigma)));
Hessian=matrix for f in flatten entries scoreEq list 
flatten entries jacobianMatrixOfRationalFunction(f);
apply({sols_0,sols_3,sols_4},i->eigenvalues sub(Hessian,matrix{coordinates(i)}))	  
----------------------------------
-- Original example
----------------------------------
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
MM=genListMatrix(sols,H)
for m in MM list eigenvalues m
-- 2 local maxima (all eigenvalues are negative) and 1 saddle point (both negative and positive) 

-----------------------------------
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
--(J,Sinv)=scoreEquationsInternal(gaussianRing G,S,SampleData=>false);
(J,Sigma)=scoreEquations(gaussianRing G,S,SampleData=>false, CovarianceMatrix=>true);
sols=zeroDimSolve(J)
netList sols
--3 real solutions
M=genListMatrix(sols,inverse Sigma)
L=checkPD(M)
solverMLE(G,S,SampleData=>false)

 S1= ( L_0+L_1)*(1/2)
 S1= matrix {{.68366, 0, .571755, 0}, {0, .224382, 0, .294633}, {.571755, 0, .848995, .238487}, {0, .294633, .238487, .761644}}
 N=matrix{{-0.2208423, -0.05311186,  0.52218624, -1.5876211},{-0.2051235,  0.94012152 ,-0.07092418, -0.3880292},
     {-0.2756983, -0.28463992, -0.64418473, -0.5294548},{-0.4108377, 0.39453094, -0.15861016,  1.9259545}}
 S2=S1+sampleCovarianceMatrix N
 S2=matrix {{.68366, .0716539, .571755, .48852}, {.0716539, .224382, .309413, .294633}, {.571755, .309413, .848995, .238487}, {.48852, .294633, .238487, .761644}}
 solverMLE(G,S2,SampleData=>false)
 (J,Sinv)=scoreEquations(gaussianRing G,S2,SampleData=>false, CovarianceMatrix=>true);
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

