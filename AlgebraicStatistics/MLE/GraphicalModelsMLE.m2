-- -*- coding: utf-8-unix -*-

newPackage(
     "GraphicalModelsMLE",
     Version => "0.3",
     Date => "June 5, 2020",
     Authors => {
	  {Name=> "Carlos Amendola", 
	   Email=> "carlos.amendola@tum.de",
	   HomePage=>"http://www.carlos-amendola.com/"},
       
	  {Name => "Luis David Garcia Puente",
	   Email => "lgarcia@shsu.edu",
	   HomePage => "http://www.shsu.edu/~ldg005"},
       
          {Name=> "Roser Homs Pons", 
	   Email=> "roserhp@gmail.com",
	   HomePage=>"https://personal-homepages.mis.mpg.de/homspons/index.html"},
       
          {Name=> "Olga Kuznetsova", 
	   Email=> "kuznetsova.olga@gmail.com",
	   HomePage=>"https://okuznetsova.com"},
       
          {Name=> "Harshit J Motwani", 
	   Email=> "harshitmotwani2015@gmail.com"},
       
          {Name=> "Elina Robeva", 
	   Email=> "erobeva@gmail.com",
	   HomePage=>"http://math.berkeley.edu/~erobeva"},
       
          {Name=> "David Swinarski", 
	   Email=> "dswinarski@fordham.edu",
	   HomePage=>"http://faculty.fordham.edu/dswinarski"}
	  },
     Headline => "maximum likelihood estimates for graphical statistical models",
     DebuggingMode => true,
     PackageExports => {"GraphicalModels","Graphs","EigenSolver","StatGraphs"}
     )
export {
    "checkPD",
    "checkPSD",
    "doSaturate",-- optional argument in scoreEquationsFromCovarianceMatrix and solverMLE
    "jacobianMatrixOfRationalFunction",
    "sampleCovarianceMatrix",
    "sampleData",-- optional argument in scoreEquationsFromCovarianceMatrix and solverMLE
    "saturateOptions", -- optional argument in scoreEquationsFromCovarianceMatrix and solverMLE
    "scoreEquationsFromCovarianceMatrix",
    "scoreEquationsFromCovarianceMatrixUndir",
    "solverMLE"
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


-----------------------------------------
-- move to a new ring, lpR, which does not have the s variables
-- When the function gaussianRing is applied to an object of 
-- class Graph, Digraph and MixedGraph, then it outputs
-- a ring with indeterminates s(i,j) for 1 ≤i ≤j ≤n, and
-- additionally l(i,j), p(i,j) for mixed graphs or k(i,j) for 
-- graphs with 1 ≤i ≤j ≤n where n is the number of vertices in G
-- (see its documentation for details) In our computations, we
-- do not need the indeterminates s(i,j) .  
-- Output is a map F:R-->lpR and the ring lpR.
-----------------------------------------
changeRing=(d,R)->(
    -- count the number of S variables. d is the number of
    -- rows in a relevant matrix (clarify!!)
    numSvars:=lift(d*(d+1)/2,ZZ);
    --lp ring is the ring without the s variables
    lpRvarlist:=apply(numgens(R)-numSvars,i->(gens(R))_i);
    KK:=coefficientRing(R);
    lpR:=KK[lpRvarlist];
    -- here i is taken to numgens(R)-numSvars-1 because
    -- indexing starts from 0. But for subscripts it
    -- starts from 1.
    lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
    F:=map(lpR,R,lpRTarget);
    return (F,lpR)
    );


-- alternative function to changeRing, removes s variables and returns a list instead of sequence, I was getting some errors with changeRing
removeSvar = (R) ->(
    numSvars := #set support covarianceMatrix R;
    lpRvarlist := gens R - set support covarianceMatrix R;
    KK := coefficientRing(R);
    lpR := KK[lpRvarlist];
    lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
    F:=map(lpR,R,lpRTarget);
    return {F,lpR};
);
------------------------------------------------------
-- Substitues a list of points on a list of matrices
    -- input -  list of points from solves
    -- output - list of matrices after substituting these values
------------------------------------------------------
genListMatrix = (L,R) ->
(
    T := {};
    K:=undirectedEdgesMatrix R;--checks for Graph, should be changed to allow for MixedGraphs
    --ring mapping begins
     -- d is equal to the number of vertices in G
    d := numRows K;
    -- move to a new ring, lpR, which does not have the s variables
   -- F:=changeRing(d,R);
   F := removeSvar(R);
   K = F_0(K);
   -- K = matRtolpR(K,F);
    
    
    --ring mapping ends 
    for l in L do
    (
    	T = T|{coordinates(l)};	
    );
    M := {};
    for t in T do
    (
    	m := substitute(K,matrix{t});	
    	M = M|{m};
    );    
    return M
);

----------------------------------------------
-- Selects all argmax for log det K- trace(S*K),
-- where K is an element of the list L.
-- We assume that L is the intersection of the 
-- variety of the ideal generated by the Jacobian
-- of critical equations and the cone of PD matrices.

-- input - list L of candidate Sinv matrices (Sinv is Sigma^{-1}) and
-- sample covariance matrix V. Notation in line with scoreEquationsFromCovarianceMatrix
-- output -list of argmax
----------------------------------------------
maxMLE=(L,V)->(
    if #L==0 then  error("No critical points to evaluate");
    if #L==1 then  E:=inverse L_0;
    if #L>=1 then 
    	eval:=for Sinv in L list log det Sinv- trace (V*Sinv);
	evalReal:={};
	for pt in eval do (if isReal pt then evalReal=evalReal  | {pt});
	if #evalReal==0 then  error("No critical point evaluates to a real solution");
	maxPt:=max evalReal;
	indexOptimal:=positions(eval, i ->i== maxPt);
	E={};
        for i in indexOptimal do E=E | {L_i};
    return (maxPt, E) 
    );


--**************************--
--  METHODS 	      	   	  --
--**************************--
sampleCovarianceMatrix = method(TypicalValue =>Matrix);
sampleCovarianceMatrix(List) := (U) -> (
    n := #U;
    U = apply(#U, i -> if ring(U_i)===ZZ then matZZtoQQ(U_i) else U_i);
    Ubar := matrix{{(1/n)}} * sum(U);
    return ((1/n)*(sum apply(n, i -> (transpose (U#i-Ubar))*(U#i-Ubar))));        
);

sampleCovarianceMatrix(Matrix) := (U) -> (
   X := {};
   n := numRows U;
   -- converting it to list of matrix; rows of matrix correponds to the elements of the list
   X = for i to n-1 list U^{i};
   return sampleCovarianceMatrix(X);
);



jacobianMatrixOfRationalFunction = method(TypicalValue =>Matrix);
jacobianMatrixOfRationalFunction(RingElement) := (F) -> (
    f:=numerator(F);
    g:=denominator(F);
    R:=ring(f);
    answer:=diff(vars(R), f) * g - diff(vars(R), g)*f;
    answer=substitute(answer, ring(F));
    return transpose(matrix({{(1/g)^2}})*answer)
);

scoreEquationsFromCovarianceMatrix = method(TypicalValue =>Ideal, Options =>{doSaturate => true, saturateOptions => options saturate, sampleData=>true});
scoreEquationsFromCovarianceMatrix(Ring,List) := opts ->(R, U) -> ( 
    ----------------------------------------------------
    -- Extract information about the graph
    ---------------------------------------------------- 
    -- Lambda
    L := directedEdgesMatrix R;
    -- K 
    K:= undirectedEdgesMatrix R;
    -- Psi
    P := bidirectedEdgesMatrix R;
    
    ----------------------------------------------------
    -- Create an auxiliary ring and its fraction field
    -- which do not have the s variables
    ----------------------------------------------------
    -- d is equal to the number of vertices in G
    d := numRows L;
    -- create a new ring, lpR, which does not have the s variables
    (F,lpR):=changeRing(d,R);
    -- create its fraction field
    FR := frac(lpR);
    
    -----------------------------------------------------
    -- Construct Omega
    -----------------------------------------------------
    -- Kinv
    K=substitute(K, FR);
    Kinv:=inverse K;
    P=substitute(P,FR);
       
     --Omega
    if K==0 then W:=P else (if P==0 then W=Kinv else W = directSum(Kinv,P));
    
    -- move to FR, the fraction field of lpR
    L= substitute (L,FR);
    
    -- Sigma
    if L==0 then S:=W else (
	IdL := inverse (id_(FR^d)-L);
    	S = (transpose IdL) * W * IdL
	);
    if S == Kinv then Sinv:= K else Sinv = inverse S; 
    
    -- Sample covariance matrix
    if opts.sampleData then V := sampleCovarianceMatrix(U) else 
    (V=U_0;
     for i from 1 to  #U-1 do V= V||U_i); -- if U is a list (even if it is inputted as a matrix, see the method scoreEquationsFromCovarianceMatrix(Ring,Matrix))
     
    -- Compute ideal J   
    C1 := trace(Sinv * V);
    C1derivative := jacobianMatrixOfRationalFunction(C1);
    LL :=jacobianMatrixOfRationalFunction (det Sinv)*matrix{{1/det(Sinv)}} - C1derivative;
    LL=flatten entries(LL);
    denoms := apply(#LL, i -> lift(denominator(LL_i), lpR));
    J:=ideal apply(#LL, i -> lift(numerator(LL_i),lpR));
    --Saturate
    if opts.doSaturate then 
    (   argSaturate:=opts.saturateOptions  >>newOpts-> args ->(args, newOpts);
    	for i from 0 to (#denoms-1) do (
	    if degree denoms_i =={0} then J=J else  
	    J=saturate(argSaturate(J,denoms_i))); 
	);
    return J;
);

scoreEquationsFromCovarianceMatrix(Ring,Matrix) := opts -> (R, U) -> (
   X := {};
   n := numRows U;
   -- converting U to list of matrices; rows of matrix correponds to the elements of the list
   X = for i to n-1 list U^{i};
   return scoreEquationsFromCovarianceMatrix(R,X,opts);
);


scoreEquationsFromCovarianceMatrixUndir = method();
scoreEquationsFromCovarianceMatrixUndir(Ring,Matrix) := (R, U) -> (
    --update to not assume zero mean variables
   n := numRows U;
   S:=sampleCovarianceMatrix U;
    -- Concentration matrix K
    K:=undirectedEdgesMatrix R;
    -- d is equal to the number of vertices in G
    d := numRows K;
    -- move to a new ring, lpR, which does not have the s variables
    (F,lpR):=changeRing(d,R);
    K=F(K);
    I:=ideal{jacobian ideal{determinant(K)}-determinant(K)*jacobian(ideal{trace(K*S)})};
    J:=saturate(I,ideal{determinant(K)});
    return J;
 );

scoreEquationsFromCovarianceMatrixUndir(Ring,List) := (R, U) -> (
    n := #U;
    L := for i from 0 to n-1 list {U_i};
    M := matrix(L);
    return scoreEquationsFromCovarianceMatrixUndir(R,M);
);

checkPD = method(TypicalValue =>List);
checkPD(List) := (L) -> (
   mat := {};
    for l in L do
    (
    	flag := 0;
    	-- Compute eigenvalues for each matrix
	L1 := eigenvalues l;
    	--Check whether all of them are positive and real
	for t in L1 do 
    	(	 
	    if 0 >= t then flag = 1;
	    if not isReal t then flag=1;
     	);
        if flag == 0 then mat = mat | {l} ;
    );
    if mat == {} then print("none of the matrices are pd");
    return mat;
);

checkPSD = method(TypicalValue =>List);
checkPSD(List) := (L) -> (
   mat := {};
    for l in L do
    (
    	flag := 0;
    	-- Compute eigenvalues for each matrix
	L1 := eigenvalues l;
    	--Check whether all of them are non-negative and real
	for t in L1 do 
    	(	 
	    if 0 > t then flag = 1;
	    if not isReal t then flag=1;
     	);
        if flag == 0 then mat = mat | {l} ;
    );
    if mat == {} then print("none of the matrices are psd");
    return mat;
);


solverMLE = method(TypicalValue =>List, Options =>{doSaturate => true, saturateOptions => options saturate, sampleData=>true});-- add a choice of and sigmaInverse solver e,sigmaInverse=>true
solverMLE(MixedGraph,List):=  opts ->(G,U) -> (
    -- generate the Gaussian ring of the MixedGraph
    R:= gaussianRing(G);
    -- compute the sample covariance matrix
    if opts.sampleData then V := sampleCovarianceMatrix(U) else 
    (V=U_0;
     for i from 1 to  #U-1 do V= V||U_i); -- if U is a list (even if it is inputted as a matrix, see the method scoreEquationsFromCovarianceMatrix(Ring,Matrix))
    -- generate the ideal of the score equations
    if opts.doSaturate then (
	 argSaturate1:=opts.saturateOptions  >>newOpts-> args ->(args, saturateOptions=>newOpts,sampleData=>false);
         J:=scoreEquationsFromCovarianceMatrix(argSaturate1(R,V));)
    else J= scoreEquationsFromCovarianceMatrix(R,V,doSaturate=>false, sampleData=>false);
    -- check that the system has finitely many solutions
    if dim J =!= 0 then (
	print ("the ideal is not zero-dimensional");
	return J;
	) else (
    -- solve system with eigensolver
    sols:=zeroDimSolve(J);
    --evaluate matrices on solutions
    M:=genListMatrix(sols,R);
    --consider only PD matrices    
    L:=checkPD M;
    --find the optimal points
    optSols:=maxMLE(L,V);
    return optSols);    
);

solverMLE(Graph,List) := opts -> (G, U) -> (
    return solverMLE(mixedGraph (G),U, opts);
    );

solverMLE(Digraph,List) := opts -> (G, U) -> (
    return solverMLE(mixedGraph (G),U, opts);
    );
solverMLE(Bigraph,List) := opts -> (G, U) -> (
    return solverMLE(mixedGraph (G),U, opts);
    );
solverMLE(Graph,Digraph,List) := opts -> (G,D,U) -> (
    return solverMLE(mixedGraph (G,D),U, opts);
    );
solverMLE(Digraph,Graph,List) := opts -> (D,G,U) -> (
    return solverMLE(mixedGraph (D,G),U, opts);
    );
solverMLE(Digraph,Bigraph,List) := opts -> (D,B,U) -> (
    return solverMLE(mixedGraph (D,B),U, opts);
    );
solverMLE(Bigraph,Digraph,List) := opts -> (B,D,U) -> (
    return solverMLE(mixedGraph (B,D),U, opts);
    );
solverMLE(Graph, Bigraph,List) := opts -> (G,B,U) -> (
    return solverMLE(mixedGraph (G,B),U, opts);
    );
solverMLE(Bigraph,Graph,List) := opts -> (B,G,U) -> (
    return solverMLE(mixedGraph (B,G),U, opts);
    );
solverMLE(Graph, Digraph, Bigraph, List) := opts -> (G,D,B,U) -> (
    return solverMLE(mixedGraph (G,D,B),U, opts);
    );
solverMLE(Digraph, Bigraph, Graph, List) := opts -> (D,B,G,U) -> (
    return solverMLE(mixedGraph (D,B,G),U, opts);
    );
solverMLE(Bigraph, Graph, Digraph, List) := opts -> (B,G,D,U) -> (
    return solverMLE(mixedGraph (B,G,D),U, opts);
    );
solverMLE(Graph,Bigraph, Digraph, List) := opts -> (G,B,D,U) -> (
    return solverMLE(mixedGraph (G,B,D),U, opts);
    );
solverMLE(Bigraph, Digraph,Graph, List) := opts -> (B,D,G,U) -> (
    return solverMLE(mixedGraph (B,D,G),U, opts);
    );
solverMLE(Digraph, Graph, Bigraph, List) := opts -> (D,G,B,U) -> (
    return solverMLE(mixedGraph (D,G,B),U, opts);
    );

--******************************************--
-- DOCUMENTATION     	       	    	    -- 
--******************************************--

beginDocumentation()

doc ///
    Key
        GraphicalModelsMLE
    Headline
        a package for MLE estimates of parameters for statistical graphical models 
    Description        
        Text
            Add some text and an example.   
	    
	    In the example below, we create the score equations (defining the critical points of the log likelihood function written in terms of the covariance matrix) associated to the four data vectors $(1,2,1,-1)$, $(2,1,3,0)$, $(-1,0,1,1)$, $(-5,3,4,-6)$ for a graphical model with four vertices, five directed edges, and one bidirected edge.
	    
        Example	   
	    G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
	    R = gaussianRing(G)
	    U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
            scoreEquationsFromCovarianceMatrix(R,U)

    Caveat
        GraphicalModelsMLE requires  @TO Graphs@,  @TO StatGraphs@ and  @TO GraphicalModels@.
///

--------------------------------
-- Documentation
--------------------------------

doc /// 
    Key
        sampleCovarianceMatrix
        (sampleCovarianceMatrix,List) 
    Headline
        compute the sample covariance matrix of a list of data vectors or a matrix
    Usage
        sampleCovarianceMatrix U
    Inputs
        U:List
	  each element corresponds to an observations vector  
	U:Matrix
	  each row corresponds to an observations vector
    Outputs
         :Matrix
    Description 
        Text
	    The sample covariance matrix is $S = \frac{1}{n} \sum_{i=1}^{n} (X^{(i)}-\bar{X}) (X^{(i)}-\bar{X})^T$.  
	    The entries here are not the unbiased estimators of the variance/covariance; that is, the entries here 
	    correspond to the outputs of the commands VAR.P and COVARIANCE.P in Excel, not VAR and COVARIANCE in Excel.
	    
	    We assume that the data vectors are entered as a list of row matrices, all with the same width, or as a matrix.
	    Each element of the list (or a row of the matrix) correspond to an observations vector
        Example
          M = {matrix{{1, 2, 0}}, matrix{{-1, 0, 5/1}}, matrix{{3, 5, 2/1}}, matrix{{-1, -4, 1/1}}};
	  sampleCovarianceMatrix(M)
	  U= matrix{{1,2,1,-1},{2,1,3,0},{-1, 0, 1, 1},{-5, 3, 4, -6}};
	  sampleCovarianceMatrix(U)	    
     ///

doc /// 
    Key
        jacobianMatrixOfRationalFunction
        (jacobianMatrixOfRationalFunction,RingElement) 
    Headline
        compute the Jacobian matrix of a rational function
    Usage
        jacobianMatrixOfRationalFunction(F)
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
            jacobianMatrixOfRationalFunction(F)
	    R=QQ[t_1,t_2,t_3];
	    FR=frac R;
	    jacobianMatrixOfRationalFunction( (t_1^2*t_2)/(t_1+t_2^2+t_3^3) )
   ///

-------------------------------------------------------
-- Documentation scoreEquationsFromCovarianceMatrix ---
-------------------------------------------------------

doc /// 
    Key
        scoreEquationsFromCovarianceMatrix
        (scoreEquationsFromCovarianceMatrix,Ring,List) 
	(scoreEquationsFromCovarianceMatrix,Ring,Matrix)
    Headline
        computes the score equations that arise from the log-likelihood function of $\Sigma^{-1}$ of a Gaussian model
    Usage
        scoreEquationsFromCovarianceMatrix(R,U)
    Inputs
        R:Ring
	  the Gaussian ring of an underlying loopless mixed graph;
	U:List 
	  If sampleData=>true (the default setting), then this is data given as a list. 
	  Can be a list of lists or a list of matrices.
	  Each element of the list corresponds to an observations vector. 
	  If sampleData=>false, then this is a sample covariance matrix given as a list, 
	  i.e., each element of the list corresponds to a row in the covariance matrix.
	U:Matrix
	  If sampleData=>true (the default setting), then data given as  a matrix.
	  Each row of the matrix corresponds to an observations vector. 
	  If sampleData=>false, then this is a sample covariance matrix.
	   
    Outputs
         :Ideal
	  ideal generated by the Jacobian of the log-likelihood function of $\Sigma^{-1}$
    Description 
        Text
	    This function computes the score equations that arise from the 
	    maximization of the log-likelihood function of a graphical Gaussian
	    statistical model given in Proposition 7.1.10 (Sullivant, 2018):
	    
	    $\ell(\Sigma)=-n/2 log det \Sigma - n/2 tr (S\Sigma^{-1})$,
	    
	    as a functio of $\Sigma^{-1}$.
	    
	    The underlying graph is assumed to be a loopless mixed graph G (Sadeghi 
	    and Lauritzen, 2014) given by type @TO MixedGraph@. The nodes of $G$ are partitioned as $V = U\cup W$, such that:
	
 	    - if $i-j$ in $G$ then $i,j\in U$
	    
  	    - if $i\leftarrow \rightarrow j$ in $G$ then $i,j\in W$ 
	    
	    -  there is no directed edge $i\to j$ in $G$ such that $i\in W$ and $j\in U$.
	    
	    Vertices that are not adjacent to an undirected or a bidirected edge are assumed to be in U.
	    
	    We require that the directed part is a DAG, i.e., there should not be any
	    directed cycles after the identification of the connected 
	    undirected and bidirected components. 
	    
	    The covariance matrix  $\Sigma$ is given by		    	    
	    
	    $\Sigma=(I-\Lambda)^{-T} diag(K^{-1}, \Psi) (I-\Lambda)^{-1}$, where
	    
	    - $\Lambda $ is a  $\mathbb{R}^{V x V }$ matrix such that $\lambda_{i,j}=0$
	    whenever $i \rightarrow j$ is not in  $E$;
	    
	    - $diag(K^{-1}, \Psi) $ is a block-diagonal $\mathbb{R}^{V x V }$ matrix;
	    
	    - $K$ is a  $\mathbb{R}^{U x U }$ matrix such that $k_{i,j}=0$
	    whenever $i - j$ is not in  $E$;
	    
	    - $\Psi$ is a  $\mathbb{R}^{W x W }$ matrix such that $ \psi_{i,j}=0$
	    whenever $i \leftarrow  \rightarrow  j$ is not in  $E$;
	    
	    References
	    
	    Sadeghi, K. and Lauritzen, S., 2014. Markov properties for mixed graphs. Bernoulli, 20(2), pp.676-696.
	    
	    Sullivant, S., 2018. Algebraic statistics (Vol. 194). American Mathematical Soc.
	
	Example
	    G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
	    R = gaussianRing(G)
	    U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
            scoreEquationsFromCovarianceMatrix(R,U)
    ///
    
doc ///
  Key
    doSaturate
  Headline
    optional input to cancel the saturation of the ideal generated of the Jacobian by its denominators
  SeeAlso
     scoreEquationsFromCovarianceMatrix
   ///
doc ///
  Key
    [scoreEquationsFromCovarianceMatrix, doSaturate]
  Headline
     optional input to cancel the saturation of the ideal generated of the Jacobian by its denominators
  Usage
    scoreEquationsFromCovarianceMatrix(R,U,doSaturate=>true)  
  Inputs 
    true: Boolean
    
  Description
    Example
     G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}});
     R=gaussianRing(G);
     U=matrix{{1,2,1,-1},{2,1,3,0},{-1, 0, 1, 1},{-5, 3, 4, -6}};
     JnoSat=scoreEquationsFromCovarianceMatrix(R,U,doSaturate=>false);
     dim JnoSat  
     degree JnoSat
  SeeAlso
     scoreEquationsFromCovarianceMatrix   	
///

doc ///
  Key
    saturateOptions
  Headline
    optional input to set up saturation, use any option from saturate
  SeeAlso
     scoreEquationsFromCovarianceMatrix
     doSaturate
     saturate
   ///
doc ///
  Key
    [scoreEquationsFromCovarianceMatrix,  saturateOptions]
  Headline
    optional input to set up saturation, use any option from saturate
  Usage
    scoreEquationsFromCovarianceMatrix(R,U,saturateOptions=>{options saturate})  
  Inputs 
    L: List
     list of options to set up the saturation. Accepts any option from the function
     @TO saturate@.
    
  Description
    Example
     G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}});
     R=gaussianRing(G);
     U=matrix{{1,2,1,-1},{2,1,3,0},{-1, 0, 1, 1},{-5, 3, 4, -6}};
     J=scoreEquationsFromCovarianceMatrix(R,U,saturateOptions => {DegreeLimit=>1, MinimalGenerators => false});

  SeeAlso
     scoreEquationsFromCovarianceMatrix   
     doSaturate
     saturate	
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
    	checkPD
	(checkPD,List)
    Headline
    	checks which matrices from a list are positive definite
    Usage
    	checkPD(L)
    Inputs
    	L: List  
	    list of matrices
    Outputs
    	 : List
	   list of positive definite matrices
    Description
    	Text
	   This function takes a list of matrices and returns another list with
	   only positive definite matrices
      	Example
	    L={matrix{{1,0},{0,1}},matrix{{-2,0},{0,1}},matrix{{sqrt(-1),0},{0,sqrt (-1)}}}				
    	    checkPD(L)
     	 ///

doc   ///
    Key
    	checkPSD
	(checkPSD,List)
    Headline
    	checks which matrices from a list are positive semidefinite
    Usage
    	checkPSD(L)
    Inputs
    	L: List  
	    list of matrices
    Outputs
    	 : List
	   list of positive semidefinite matrices
    Description
    	Text
	   This function takes a list of matrices and returns another list with
	   only positive semidefinite matrices
      	Example
	    L={matrix{{1,0},{0,1}},matrix{{-2,0},{0,1}},matrix{{sqrt(-1),0},{0,sqrt (-1)}},matrix{{0,0},{0,0}}}				
    	    checkPSD(L)
     	 ///

doc ///
    Key
        solverMLE
	(solverMLE,MixedGraph,List)
    Headline
    	computes MLE from a @TO MixedGraph@ and a data sample (or sample covariance matrix)
    Usage
    	solverMLE(G,U)
    Inputs
    	G: MixedGraph
      	U: List 
    Outputs
        : List
	
    Description
    	Text
	    This function computes the critical points from the score equations and 
	    selects those that lie in the cone of positive-definite matrices.
	    See Example 2.1.13 of Sturmfels' lecture notes
	Example
	    G=graph{{1,2},{2,3},{3,4},{1,4}}
	    U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
	    --U=random(ZZ^4,ZZ^4)
	    solverMLE(G,U)				
///

--******************************************--
-- TESTS     	       	    	      	    --
--******************************************--

TEST /// 
R=QQ[x,y];
FR=frac R;
F=1/(x^2+y^2);
M=entries jacobianMatrixOfRationalFunction(F);
N=transpose {{-2*x/(x^2 + y^2)^2,-2*y/(x^2 + y^2)^2 }};
assert(M === N)
///

TEST ///
R=QQ[x_1,x_2,x_3];
FR=frac R;
M=entries jacobianMatrixOfRationalFunction( (x_1^2*x_2)/(x_1+x_2^2+x_3^3) );
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
B = matrix {{5621/9, -1037/18, -7835/18, -10565/18}, {-1037/18, 19505/36, -4897/36, 5147/36}, {-7835/18, -4897/36, 20465/36, 18941/36}, {-10565/18, 5147/36, 18941/36, 28889/36}};
assert(Y===B)	
///

TEST ///
X = matrix{{48,89,27,28},{23,19,29,94},{135,23,44,71},{91,75,24,98}};
Y = sampleCovarianceMatrix(X);
B = matrix {{29147/16, -1313/8, 220, 1609/16}, {-1313/8, 3827/4, -155, -3451/8}, {220, -155, 119/2, -63/4}, {1609/16, -3451/8, -63/4, 12379/16}};
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
///   

TEST ///
L={matrix{{1,0},{0,1}},matrix{{-2,0},{0,1}},matrix{{sqrt(-1),0},{0,sqrt (-1)}},matrix{{0.0001*sqrt(-1),0},{0,0.0000001*sqrt (-1)}}};
Y = checkPD(L);
B = {matrix{{1, 0}, {0, 1}}};
assert(Y===B)	
///

TEST ///
L={matrix{{1,0},{0,1}},matrix{{-2,0},{0,1}},matrix{{sqrt(-1),0},{0,sqrt (-1)}},matrix{{0.0001*sqrt(-1),0},{0,0.0000001*sqrt (-1)}},matrix{{0,0},{0,0}}};
Y = checkPD(L);
B = {matrix{{1, 0}, {0, 1}},matrix{{0,0},{0,0}}};
assert(Y===B)	
///
    
    
--------------------------------------
--------------------------------------
end--
--------------------------------------
--------------------------------------


--blank documentation node:
doc /// 
    Key
       gaussianMatrix
       (gaussianMatrix,Digraph,Matrix,List) 
    Headline
    Usage
    Inputs
    Outputs
    Description 
        Text
        Example
    SeeAlso
///

restart
uninstallPackage "GraphicalModelsMLE"
restart
--installPackage("Graphs", UserMode=>true)
installPackage ("GraphicalModelsMLE", RemakeAllDocumentation => true, UserMode=>true)
installPackage("GraphicalModelsMLE",UserMode=>true,DebuggingMode => true)
installPackage("GraphicalModelsMLE",UserMode=>true,DebuggingMode => true, FileName =>"/Users/lgp/Software/Macaulay2/Workshop-2020-Warwick/AlgebraicStatistics/GraphicalModelsMLE.m2")
check GraphicalModelsMLE

viewHelp "GraphicalModelsMLE"
help GraphicalModelsMLE

----------------------
-- Parameterization -- ????????????????????????????????????????????????????????????????????????
---------------------- 
---- We need this for both directed and undirected graphs:

----  parameterizations and for toric varieties the corresponding matrix. 
----  In the case of toric varieties the matrix is easy.  Here is the code, 
----  commented out to be used later when we are ready. 
---- 
----  toAMatrix = method()
----  toAMatrix List := Matrix => (M) -> (
----      if any(M,isMonomial)
----         then error "this parameterization does not correspond to a toric ideal." 
----         else (
----              Mexp := apply(M, exponents);
----              transpose matrix apply(Mexp, flatten)))
----
---- isMonomial = method()
---- isMonomial RingElement := Boolean => (m) -> (
----      termList := terms m;
----      if #termList == 1 then true else false)

---- isMonomial works well as long as m is actually a polynomial or monomial and not 
---- an element of ZZ, QQ, RR, etc.


end;
restart
installPackage("GraphicalModelsMLE")
check "GraphicalModelsMLE"


