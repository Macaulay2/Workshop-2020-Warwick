-- -*- coding: utf-8-unix -*-

newPackage(
     "GraphicalModelsMLERoser",
     Version => "0.3",
     Date => "June 5, 2020",
     Authors => {
	  {Name => "Luis David Garcia Puente",
	   Email => "lgarcia@shsu.edu",
	   HomePage => "http://www.shsu.edu/~ldg005"},
          {Name=> "David Swinarski", 
	   Email=> "dswinarski@fordham.edu",
	   HomePage=>"http://faculty.fordham.edu/dswinarski"}, 
          {Name=> "Elina Robeva", 
	   Email=> "erobeva@gmail.com",
	   HomePage=>"http://math.berkeley.edu/~erobeva"},
   	  {Name=> "Harshit J Motwani", 
	   Email=> "harshitmotwani2015@gmail.com"},
   	  {Name=> "Olga Kuznetsova", 
	   Email=> "kuznetsova.olga@gmail.com",
	   HomePage=>"https://okuznetsova.com"},
   	  {Name=> "Roser Homs Pons", 
	   Email=> "roserhp@gmail.com",
	   HomePage=>"https://personal-homepages.mis.mpg.de/homspons/index.html"}
	  },
     Headline => "maximum likelihood estimates for structural equation models",
     DebuggingMode => true,
     PackageExports => {"GraphicalModels","Graphs","EigenSolver"}
     -- WARNING: PLACE EIGENSOLVER PACKAGE IN THE USUAL M2 PACKAGE FOLDER (check path command)
     -- EigenSolver is used for computing the MLE. May need to be replaced or complemented 
     -- by a more efficient solver in the future
     )
export {"sampleCovarianceMatrix",
    "JacobianMatrixOfRationalFunction",
    "scoreEquationsFromCovarianceMatrix",
    "scoreEquationsFromCovarianceMatrixUndir",
    "evalConcentrationMatrix",
    "PDcheck",
    "MLEsolver"
       	} 
     
     
--**************************--
--  INTERNAL ROUTINES       --
--**************************--

--*************************************--
--  Functions (local) used throughout  --
--*************************************--




--------------------------------------------
-- turn the entries of an integer matrix into rational numbers
--------------------------------------------

matZZtoQQ = (M) -> (
    E:=entries(M);    
    return matrix apply(#E, i -> apply(#(E_i), j -> (1/1)*E_i_j))    
);

--------------------------------------------
-- change the entries of matrix M from ring R to ring lpR 
-- via the map F:R-->lpR
--------------------------------------------

matRtolpR = (M,F) -> (
    E:=entries(M);    
    return matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)))    
);



--**************************--
--  METHODS 	      	   	  --
--**************************--
sampleCovarianceMatrix = method();
-- why List and not Matrix as input?
sampleCovarianceMatrix(List) := (U) -> (
    n := #U;
    U = apply(#U, i -> if ring(U_i)===ZZ then matZZtoQQ(U_i) else U_i);
    Ubar := matrix{{(1/n)}} * sum(U);
    return ((1/n)*(sum apply(n, i -> (transpose (U#i-Ubar))*(U#i-Ubar))));        
);


JacobianMatrixOfRationalFunction = method();
JacobianMatrixOfRationalFunction(RingElement) := (F) -> (
    f:=numerator(F);
    g:=denominator(F);
    R:=ring(f);
    answer:=diff(vars(R), f) * g - diff(vars(R), g)*f;
    answer=substitute(answer, ring(F));
    return transpose(matrix({{(1/g)^2}})*answer)
);

scoreEquationsFromCovarianceMatrix = method();
scoreEquationsFromCovarianceMatrix(Ring,List) := (R, U) -> (
    V := sampleCovarianceMatrix(U);   
 --   R := gaussianRing(G);  
    -- Lambda
    L := directedEdgesMatrix R;
    -- d is equal to the number of vertices in G
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
    L = matRtolpR(L,F);
    W = matRtolpR(W,F);
    FR := frac(lpR);
    K := inverse (id_(lpR^d)-L);
    S := (transpose K) * W * K;
    Sinv := inverse substitute(S, FR);    
    C1 := trace(Sinv * V)/2;
    C1derivative := JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
    LL := (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-1/(2*det(S)))}} - (C1derivative);
    LL=flatten entries(LL);
    denoms := apply(#LL, i -> lift(denominator(LL_i), lpR));
    prod := product(denoms);
    J:=ideal apply(#LL, i -> lift(numerator(LL_i),lpR));
    J = saturate(J, prod);
    return J;
);



scoreEquationsFromCovarianceMatrixUndir = method();
scoreEquationsFromCovarianceMatrixUndir(Ring,Matrix) := (R, U) -> (
    -- check which method of computing transpose is correct.
    S:=U*transpose(U);
    -- need to fix conversion to rationals
    --S=matZZtoQQ(S);
    --V := sampleCovarianceMatrix(U);
    -- Concentration matrix K
    K:=undirectedEdgesMatrix R;
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
    K=F(K);
    I:=ideal{jacobian ideal{determinant(K)}-determinant(K)*jacobian(ideal{trace(K*S)})};
    J:=saturate(I,ideal{determinant(K)});
    return J;
 );


evalConcentrationMatrix = method();
evalConcentrationMatrix(List,Ring):= (L,R) ->
(
    T := {};
    K:=undirectedEdgesMatrix R;
    --ring mapping begins
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
    K=F(K);
    --ring mapping ends 
    for l in L do
    (
    	T = T|{coordinates(l)};	
    );
    M := {};
    for t in T do
    (
    	m := substitute(K,matrix{t});	
    	M = M|{m}
    );    
    return M
);


PDcheck = method();
PDcheck(List) := (L) -> (
   mat := {};
    for l in L do
    (
    	flag := 0;
    	-- Compute eigenvalues for each matrix
	L1 := eigenvalues l;
    	--Check whether all of them are positive
	for t in L1 do 
    	(	 
	    if 0 >= t then flag = 1;
     	);
        if flag == 0 then mat = mat | {l} ;
    );
    if mat == {} then print("none of the matrices are pd");
    return mat;
);



MLEsolver = method();
MLEsolver(Ideal,Ring):= (J,R) -> (
    --solve system with eigensolver
    sols:=zeroDimSolve(J);
    --evaluate concentration matrix on solutions
    M:=evalConcentrationMatrix(sols,R);
    --consider only PD matrices
    
    L:=PDcheck M;
    --check there is only one PD matrix
    if #L>1 then error("It is not an undirected graph"); 
    --return MLE: inverse of the PD matrix
    E:=inverse L_0;
    return E;    
);


--******************************************--
-- DOCUMENTATION     	       	    	    -- 
--******************************************--

beginDocumentation()

doc ///
    Key
        GraphicalModelsMLERoser
    Headline
        a package for MLE estimates of parameters for statistical graphical models 
    Description        
        Text
            Add some text and an example.   
	    
	    In the example below, we create the score equations (defining the critical points of the log likelihood function written in terms of the covariance matrix) associated to the four data vectors $(1,2,1,-1)$, $(2,1,3,0)$, $(-1,0,1,1)$, $(-5,3,4,-6)$ for a graphical model with four vertices, five directed edges, and one bidirected edge.
	    
        Example
	    needsPackage("Graphs");
	    needsPackage("GraphicalModels");
	    needsPackage("GraphicalModelsMLE");
	    G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
	    R = gaussianRing(G)
	    U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
            scoreEquationsFromCovarianceMatrix(R,U)
--	    scoreEquationsCovariance1(G, U)
    Caveat
        GraphicalModelsMLE requires Graphs.m2 and GraphicalModels.m2. 
///

--------------------------------
-- Documentation
--------------------------------

doc /// 
    Key
        sampleCovarianceMatrix
        (sampleCovarianceMatrix,List) 
    Headline
        compute the sample covariance matrix of a list of data vectors
    Usage
        sampleCovarianceMatrix U
    Inputs
        U:List
    Outputs
         :Matrix
    Description 
        Text
	    The sample covariance matrix is $S = \frac{1}{n} \sum_{i=1}^{n} (X^{(i)}-\bar{X}) (X^{(i)}-\bar{X})^T$.  The entries here are not the unbiased estimators of the variance/covariance; that is, the entries here correspond to the outputs of the commands VAR.P and COVARIANCE.P in Excel, not VAR and COVARIANCE in Excel.
	    
	    We assume that the data vectors are entered as a list of row matrices, all with the same width.
        Example
            M = {matrix{{1, 2, 0}}, matrix{{-1, 0, 5/1}}, matrix{{3, 5, 2/1}}, matrix{{-1, -4, 1/1}}};
	    sampleCovarianceMatrix(M)
	    U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}};
	    sampleCovarianceMatrix(U)	    
     ///

doc /// 
    Key
        JacobianMatrixOfRationalFunction
        (JacobianMatrixOfRationalFunction,RingElement) 
    Headline
        compute the Jacobian matrix of a rational function
    Usage
        JacobianMatrixOfRationalFunction(F)
    Inputs
        F:RingElement
    Outputs
         :Matrix
    Description 
        Text
	    This function computes the Jacobian matrix of a rational function
        Example
	    R=QQ[x,y];
	    FR=frac R;
	    F=1/(x^2+y^2);
            JacobianMatrixOfRationalFunction(F)
	    R=QQ[t_1,t_2,t_3];
	    FR=frac R;
	    JacobianMatrixOfRationalFunction( (t_1^2*t_2)/(t_1+t_2^2+t_3^3) )
   ///


doc /// 
    Key
        scoreEquationsFromCovarianceMatrix
        (scoreEquationsFromCovarianceMatrix,Ring,List) 
    Headline
        computes the score equations that arise from the log likelihood formula in terms of the covariance matrix Sigma
    Usage
        scoreEquationsFromCovarianceMatrix(R,U)
    Inputs
        R:Ring
	U:List
    Outputs
         :Ideal
    Description 
        Text
	    This function computes the score equations that arise from the log likelihood formula in terms of the covariance matrix $\Sigma$.  
	    Suppose we are given a list of data vectors, each a row matrix.  
	    The log likelihood function we want to maximize is given in Prop. 2.1.12 of Sturmfels's lecture notes (to do: update this reference to the printed version).  
	    
        Example
	    needsPackage("Graphs");
	    needsPackage("GraphicalModels");
	    G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
	    R = gaussianRing(G)
	    U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
            scoreEquationsFromCovarianceMatrix(R,U)
    ///

doc ///
    Key
    	scoreEquationsFromCovarianceMatrixUndir
	(scoreEquationsFromCovarianceMatrixUndir,Ring,Matrix)
    Headline
    	computes the score equations for undirected graphs
    Usage
    	scoreEquationsFromCovarianceMatrixUndir(R,U)
    Inputs
    	R:Ring
	U:List
    Outputs
    	:Ideal
    Description
    	Text
	    This function computes the score equations for undirected graphs. 
	    See Example 2.1.13 of Sturmfels' lecture notes
	Example
	    G=graph{{1,2},{2,3},{3,4},{1,4}}
	    R=gaussianRing(G)
	    U=random(ZZ^4,ZZ^4)
	    scoreEquationsFromCovarianceMatrixUndir(R,U)				
     ///

doc   ///
    Key
    	PDcheck
	(PDcheck,List)
    Headline
    	checks which matrices from the list are positive definite
    Usage
    	PDcheck(L)
    Inputs
    	L: List  
    Outputs
    	 : List
    Description
    	Text
	   This function takes a list of matrices and returns another list with
	   only positive definite matrices
      	Example
	    L={matrix{{1,0},{0,1}},matrix{{-2,0},{0,1}}}				
    	    PDcheck(L)
     	 ///


doc   ///
    Key
    	evalConcentrationMatrix
	(evalConcentrationMatrix,List,Ring)
    Headline
    	evaluates a list of points on the concentration matrix of the Gaussian ring
    Usage
    	evalConcentrationMatrix(L)
    Inputs
    	L: List
	R: Ring  
    Outputs
    	 : List
    Description
    	Text
	   This function takes a list of points and returns a list with
	   the concentration matrix of R evaluated on the points
      	Example
	    G=graph{{1,2},{2,3},{3,4},{1,4}}
	    R=gaussianRing(G)
	    U=random(ZZ^4,ZZ^4)
	    J=scoreEquationsFromCovarianceMatrixUndir(R,U)
	    sols=zeroDimSolve(J)
	    M=evalConcentrationMatrix(sols,R)	
     	 ///


doc ///
    Key
    	MLEsolver
	(MLEsolver,Ideal,Ring)
    Headline
    	computes MLE from score equations
    Usage
    	MLEsolver(J,R)
    Inputs
    	J: Ideal
      	R: Ring 
    Outputs
    	:Matrix
    Description
    	Text
	    This function computes the MLE from score equations. At the moment it
	    only works fine for undirected graphs. 
	    See Example 2.1.13 of Sturmfels' lecture notes
	Example
	    G=graph{{1,2},{2,3},{3,4},{1,4}}
	    R=gaussianRing(G)
	    U=random(ZZ^4,ZZ^4)
	    J=scoreEquationsFromCovarianceMatrixUndir(R,U)
	    MLEsolver(J,R)				
///

--******************************************--
-- TESTS     	       	    	      	    --
--******************************************--

TEST /// 
R=QQ[x,y];
FR=frac R;
F=1/(x^2+y^2);
M=entries JacobianMatrixOfRationalFunction(F);
N=transpose {{-2*x/(x^2 + y^2)^2,-2*y/(x^2 + y^2)^2 }};
assert(M === N)
///

TEST ///
R=QQ[x_1,x_2,x_3];
FR=frac R;
M=entries JacobianMatrixOfRationalFunction( (x_1^2*x_2)/(x_1+x_2^2+x_3^3) );
N=transpose {{2*x_1*x_2/(x_2^2 + x_3^3 + x_1) - x_1^2*x_2/(x_2^2 + x_3^3 + x_1)^2, -2*x_1^2*x_2^2/(x_2^2 + x_3^3 + x_1)^2 + x_1^2/(x_2^2 + x_3^3 + x_1) , -3*x_1^2*x_2*x_3^2/(x_2^2 + x_3^3 + x_1)^2 }};
assert(M === N)
/// 

TEST ///
M = {matrix{{1, 2, 0}}, matrix{{-1, 0, 5/1}}, matrix{{3, 5, 2/1}}, matrix{{-1, -4, 1/1}}};
N = sampleCovarianceMatrix(M);
A = matrix {{11/4, 39/8, -1}, {39/8, 171/16, 0}, {-1, 0, 7/2}};
assert(N===A)	
///

TEST ///
X = {matrix {{36, -3, -25, -36}}, matrix {{-10, 11, -29, -20}}, matrix {{-18, 33, -15, -11}}, matrix {{-42, 0, 20, 43}}, matrix {{-30, -26, 32, 2}}, matrix {{2, -38, -24, -43}} };
Y = sampleCovarianceMatrix(X);
B = matrix matrix {{5621/9, -1037/18, -7835/18, -10565/18}, {-1037/18, 19505/36, -4897/36, 5147/36}, {-7835/18, -4897/36, 20465/36, 18941/36}, {-10565/18, 5147/36, 18941/36, 28889/36}};
assert(Y===B)	
///

TEST ///
G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
R=gaussianRing(G)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
J=scoreEquationsFromCovarianceMatrix(R,U);
I=ideal(20*p_(3,4)+39,50*p_(4,4)-271,440104*p_(3,3)-742363,230*p_(2,2)-203,16*p_(1,1)-115,5*l_(3,4)+2,110026*l_(2,3)-2575,55013*l_(1,3)-600,115*l_(1,2)+26);
assert(J===I)
///     

TEST ///
G=graph{{1,2},{2,3},{3,4},{1,4}}
R=gaussianRing(G)
U=random(ZZ^4,ZZ^4)
J=scoreEquationsFromCovarianceMatrixUndir(R,U)
assert(dim J===0)
assert(degree J===5)
MLEsolver(J,R)
///   

TEST ///
G=graph{{1,2},{2,4},{3,4},{2,3}}
R=gaussianRing(G)
U=random(ZZ^4,ZZ^4)
J=scoreEquationsFromCovarianceMatrixUndir(R,U)
assert(dim J===0)
MLEsolver(J,R)
///
--------------------------------------
--------------------------------------
end--
