-------------
--FUNCTIONS--
-------------
needsPackage "GraphicalModelsMLE"

-----------------
-- DUAL VARIETY
-----------------
-- Source: Algorithm 5.1 from Blekherman, Parrilo, and Thomas, eds. 
-- Semidefinite optimization and convex algebraic geometry. 
-- SIAM, 2012.

-- Input:
-- IX - ideal of variety X in C^{n}
-- n - number of coordinates of the ambient space of X 
--     (i.e., all indices start from 1)
-- x - variables for X
-- u - variables for the dual of X
--restart
dualVariety=(IX,n,x,u)->(
    c:=codim IX;
    --JacX:=submatrix(transpose jacobian IX,toList(0..n-1));
    JacX:=diff(matrix{toList(x_1..x_n)}, transpose gens IX);
    AugJacX:=matrix{toList(u_1..u_n)}||JacX;
    SingX:=IX+minors(c,jacobian IX);
    conormalX:=saturate(IX+minors(c+1,AugJacX),SingX);
    dualX:=eliminate(toList(x_1..x_n),conormalX);
    dualX
    )

-----------------------------
-- ALGEBRAIC BOUNDARY and boundary components
-----------------------------
--algebraic boundary H_G
-- See p. 7 from Sturmfels and Uhler paragraph after Proposition 2.4
-- Algorithm for K a s x s matrix
-- 1. For all 1 <=p<=s compute I:= the ideal of the p-minors of K
-- 2. Compute the dual variety of each minimal prime of I
-- 3. Keep only those varieties, whose ideal is principal
-- 4. H_L is the product of these principal generators

boundaryComponents=(K,p,n)->(
    I:=minors(p,K);
    minPrimes:=minimalPrimes I;
    m:= length minPrimes;
    allComponents:=for i to  m-1 list (dualComponent:=dualVariety(minPrimes_i,n,l,t),numgens trim dualComponent);
    boundaryComponents:= for i in allComponents list (if i_1==1 then i_0)
    )

algBoundary=(V,n,K)->(
    delete (null, flatten (for p from 1 to V list boundaryComponents(K,p,n))) 
    )

--NOTE: use as input for algBoundary embededdK(K)

--QUESTION: shouldn't we consider p from 2 to V? the ideal of 1-minors doesn't make much sense...

-----------------------------------------
-- auxiliary function for algebraic boundary and related functions
-- returns K in the right ring
-----------------------------------------
  
embeddedK=(K)->(
    (V,n,Rtotal,S):=coloredData(K);
    R:=ring K;
    R2:=QQ[gens R,gens Rtotal];
    K=sub(K,R2);
    return(V,n,K);
    )


-------------------------------------------------
-- PROJECTIONS OF MATRICES OF GIVEN RANK ON GRAPH
-------------------------------------------------
--INPUT:
-- stats - ideal of sufficient statistics
-- rk - rank of matrices
-- S - variable matrix for samle covariance
rankProjection=(stats,rk,S)->(
    I:=minors(rk+1,S);
    varList:= support I;
    rankProjection:=eliminate(varList,I+stats)
    )

-------------------------------------------------
-- EXTENSION OF PARTIAL SOLUTIONS
-------------------------------------------------
--INPUT:
-- stats - ideal of sufficient statistics
-- rk - rank of matrices
-- S - variable matrix for samle covariance
-- numDrop - number of variables to drop for s_{i,j} (in the monomial order)
eliminationIdeal=(stats,rk,S,numDrop)->(
    I:=minors(rk+1,S);
    varList:= support I;
    return eliminate(drop(varList,-numDrop),I+stats);
    )


------------------------------------------------
-- IDEAL OF SUFFICIENT STATISTICS
------------------------------------------------
-- auxiliary function
naiveMatrixProduct=(M1,M2)->(
    m=numcols M1-1;
    matrix toList apply(0..m, i -> toList apply(0..m, j -> M1_(i,j)*M2_(i,j))))

-- ring and other basic data
-- INPUT: concentration matrix K
-- OUTPUT: sequence (Rtotal,S,stats) with
------- 1. p: number of vertices of the graph, i.e. size of matrix K
--------2. n: dimension of the space of concentration matrices
------  3. Rtotal: ring with variables s_ij and t_i with the suitable ordering for extension theorem  
------  4. S: Sample covariance matrix with variables s_ij
coloredData=(K)->(
    p:=numcols K;
    n:=#support K;
    w:=flatten toList apply(1..p, i -> toList apply(i..p, j -> (i,j)));
    v:=apply (w, ij -> s_ij);
    Rtotal:=QQ[v,reverse toList apply(1..n, i->t_i)];
    S:=genericSymmetricMatrix(Rtotal,s_(1,1),p);
    return(p,n,Rtotal,S);
    )

-- automatically generate the ideal of sufficient stats from K:
-- INPUT: concentration matrix K
-- OUTPUT: stats: ideal of sufficient statistics
suffStat=(K)->(
    (p,n,Rtotal,S):=coloredData(K);
    M:=mutableMatrix{{apply(1..n,i->0)}};
    stats:=ideal{};
    for i from 0 to n-1 do (
    	M_(0,i)=1;
    	aux=matrix M;
    	stats=stats+(t_(i+1)-sum (flatten entries naiveMatrixProduct(sub(K,aux),S)));
    	M_(0,i)=0;
      	);
    return stats;
    )

------------------------------------------------
-- VANISHING POLYNOMIALS AND INTERIOR POINTS
------------------------------------------------
-- find points that evaluate to different signs of a polynomial
--INPUT: sequence (f,m,k,p,n,stats)
-- f: polynomial we wish to evaluate
-- m: rank of empirical covariance matrix
-- k: maximum number of points we want to try
-- p: number of vertices in the graph (size of covariance matrix)
-- n: dimension of the linear space of concentration matrices
-- stats: ideal of sufficient statistics
differentSign=(f,m,k,p,n,stats)->(
     X:=random(QQ^p,QQ^m);    
     SE:=(1/m)*X*transpose(X);
     aux:=matrix {join(flatten toList apply(0..(p-1),i->toList apply(i..(p-1),j->-SE_(i,j))),
	toList apply(1..n,i->0))};
     suffE:=sub(stats,aux);
     subs:=toList apply(1..n,i->t_i=>suffE_(i-1));
     eval:=sub(f,subs);
     L1:={suffE,eval};
     --initialize vars
     sign:=1;
     suffEnew:=suffE;
     evalnew:=eval;
     i:=0;
     L2:={};
     while (i<k and sign==1) list (
          i=i+1;
          X=random(QQ^p,QQ^m);
          SE=(1/m)*X*transpose(X);
          aux=matrix {join(flatten toList apply(0..(p-1),i->toList apply(i..(p-1),j->-SE_(i,j))),
	   toList apply(1..n,i->0))};
          suffEnew=sub(stats,aux);
          subs=toList apply(1..n,i->t_i=>suffEnew_(i-1));
          evalnew=sub(f,subs);
          if sub(evalnew*eval,QQ)<0 then sign=-1;
          if sign==1 then continue; (suffEnew,evalnew,i));
     if sign==-1 then L2={suffEnew,evalnew,i};
     return(L1,L2); 
     )

--check empirically whether the generators of an ideal change signs for
--covariance matrices of a given rank
--INPUT: sequence (I,m,k,p,n,stats)
-- I: rank projection ideal
-- m: rank of empirical covariance matrix
-- k: maximum number of points we want to try
-- p: number of vertices in the graph (size of covariance matrix)
-- n: dimension of the linear space of concentration matrices
-- stats: ideal of sufficient statistics
empiricalVanishingPolynomials=(I,m,k,p,n,stats)->(
    L:={};
    for f in flatten entries gens I do (
	(L1,L2)=differentSign(f,m,k,p,n,stats);
	if L2!={} then L=append(L,f););
    return L;
    )
	
-- generate points in the interior from two points with different sign
-- auxiliary functions
genListMatrix = (L,A) ->
(
		T:= for l in L list coordinates(l);
		M:= for t in T list substitute(A,matrix{t});
    return M
);

maxMLE=(L,V)->(
    if #L==0 then  error("No critical points to evaluate");
    if #L==1 then  (E:=L_0; maxPt:=log det L_0- trace (V*L_0))
    else
    	(eval:=for Sinv in L list log det Sinv- trace (V*Sinv);
	evalReal:=for pt in eval when isReal pt list pt;
	if #evalReal==0 then  error("No critical point evaluates to a real solution");
	maxPt=max evalReal;
	indexOptimal:=positions(eval, i ->i== maxPt);
		E= for i in indexOptimal list L_i;);
    return (maxPt, E)
    );

adjacencyMat=(A)->(
    matrix toList apply(0..(p-1),i->toList apply(0..(p-1),j-> if A_(i,j)!=0 then 1_QQ else 0_QQ))
)

realPartMatrix = A -> matrix apply(entries A, r -> r/realPart)

checkRealReg = method(TypicalValue =>List);
checkRealReg(List):=(L) -> (
		for l in L
		list (if not length (select(eigenvalues l, i->abs(realPart i)<0.000000001 
			             or abs(imaginaryPart i )>0.0000000001))==0
		      then continue; 
		      realPartMatrix l)
);
checkRealReg(Matrix):=(L)->{
    return checkRealReg(L);
};


--m: rk of empirical covariance matrix
--k: number of points we want to try	
empiricalMLEexistence=(m,k,K)->(
    count:=0;
    for i from 1 to k do (    
    	p=numcols K;
    	--X=random(QQ^p,QQ^m);              
    	X=random(RR^p,RR^m);  
	X=matrix toList apply(0..p-1,i->toList apply(0..m-1,j->promote(X_(i,j),QQ)));  
	S=(1/m)*X*transpose(X);
	--S=matrix toList apply(0..2,i->toList apply(0..2,j->lift(S_(i,j),QQ)));
    	I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
    	J=saturate(I,det K);
	if dim J!=0 then return(count,"The ideal does not have the right dimension ");
    	criticalPoints=zeroDimSolve(J);
        criticalMatrices=genListMatrix(criticalPoints,K);
	if checkPD(criticalMatrices,ZeroTolerance=>1e-7)!={} then count=count+1;
        );
        return(count);
     )	

empiricalMLENoPD=(m,k,K)->(
    count:=0;
    for i from 1 to k do (    
    	p=numcols K;
    	X=random(QQ^p,QQ^m);              
    	S=(1/m)*X*transpose(X);
    	I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
    	J=saturate(I,det K);
	if dim J!=0 then return("The ideal does not have the right dimension ");
    	criticalPoints=zeroDimSolve(J);
        criticalMatrices=genListMatrix(criticalPoints,K);
	if length(criticalMatrices)!=0 then count=count+1;
        );
        return(count);
     )	

empiricalMLERealReg=(m,k,K)->(
    count:=0;
    for i from 1 to k do (    
    	p=numcols K;
    	X=random(QQ^p,QQ^m);              
    	S=(1/m)*X*transpose(X);
    	I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
    	J=saturate(I,det K);
	if dim J!=0 then return("The ideal does not have the right dimension ");
    	criticalPoints=zeroDimSolve(J);
        criticalMatrices=genListMatrix(criticalPoints,K);
	if checkRealReg(criticalMatrices)!={} then count=count+1;
        );
        return(count);
     )	

MLE=(K,coord)->(
        score:=ideal{jacobian(matrix{{det K}})-det(K)*sub(transpose(matrix{coord}),ring K)};
	scoreEq=saturate(score,det K);
	if dim scoreEq!=0 then return("The ideal does not have the right dimension ");
	criticalPoints=zeroDimSolve(scoreEq);
	criticalMatrices=genListMatrix(criticalPoints,K);
	existence=checkPD(criticalMatrices);
	if existence!={} then return(degree scoreEq)
	else return("MLE doesn't exist");
     )

	
interiorPoint=(v1,v2,n,f,K)->(
    Rt:=QQ[t_1..t_n,lambda];
    tol:=0.000000000000001; --only applies to discard complex sols
    v:=sub(vector delete(lambda,gens Rt),Rt);
    P:=sub(v1,Rt);
    Q:=sub(v2,Rt);
    paramr:=ideal flatten entries (v-P-lambda*(Q-P));
    r=eliminate(lambda,paramr);
    Rtt:=QQ[t_1..t_n];
    r=sub(r,Rtt);
    f=sub(f,Rtt);
    I:=ideal{f,r};
    sols:=zeroDimSolve I;
    basicMatSuf:={};
    M:=mutableMatrix {{apply(1..n,i->0)}};
    for i from 1 to n do (
        M_(0,i-1)=1;
        aux=matrix M;
        A=adjacencyMat naiveMatrixProduct(sub(K,aux),S);
        basicMatSuf=append(basicMatSuf,{A,sum(flatten entries A)});
        M_(0,i-1)=0;
        );
    Rl=ring(K);
    intPoints:={};
    for k from 0 to (#sols -1) do(
	--next two lines allow us to omit complex solutions
        listIm=apply(coordinates sols_k,i->abs(imaginaryPart i));
        if sum apply(listIm,i->if i>tol then 1 else 0)>0 then continue;
	--now we should only work with the real part to avoid problems
	coord=apply(coordinates sols_k,j->lift(j,QQ));
	MatSuf=sum toList apply(0..(n-1),j->(coord_j/(basicMatSuf_j)_1)*(basicMatSuf_j)_0);
	scoreEq=saturate(ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(MatSuf*K)}})},det K);
        criticalPoints=zeroDimSolve(scoreEq);
        criticalMatrices=genListMatrix(criticalPoints,K);
	if checkPD(criticalMatrices)!={} then intPoints=append(intPoints,coord);
        );	
    return intPoints;
)


----------------------------------------------------------
-- RANK m COMPLETION OF A GIVEN SUFF STATS
----------------------------------------------------------

rankCompletion=(intP,p,rk,stats,S,Rtotal,tol)->(
    pp:=lift(p*(p+1)/2,ZZ);
    dropList:=toList (1..pp);
    varList:=reverse flatten toList apply(1..p,i->toList apply(i..p,j->s_(i,j)));
    subs:=toList apply(1..n,i->sub(t_i,Rtotal)=>intP_(i-1));
    for i from 0 to pp-1 do(
         partialEq=eliminationIdeal(stats,rk,S,dropList_i);
         use Rtotal;
         partialSol=sub(partialEq,subs);      	 
	 degList=apply(flatten entries gens partialSol,i->degree i);
         pos=positions(flatten degList,i->i>0);
         if pos=={} then return("No partial solution");   
         pos1=positions(flatten degList,i->i==1);
         if pos1=={} 
	     then( 
		 sols=zeroDimSolve(sub(ideal partialSol_(pos_0),QQ[support partialSol_(pos_0)]));
		 sol={};
		 k=0;
		 while (k<#sols and sol=={}) do(
	         --next two lines allow us to omit complex solutions
                 listIm=apply(coordinates sols_k,i->abs(imaginaryPart i));
                 if sum apply(listIm,i->if i>tol then 1 else 0)>0 then continue else sol=coordinates sols_k;
		 k=k+1;);
                 if sol=={} then return("No real partial solutions") else subs=join({sub(varList_i,Rtotal)=>lift(sol_0,QQ)},subs);
	         )
	     else(  
         	 f=sub(partialSol_(pos1_0),QQ[support partialSol_(pos1_0)]);
         	 A=-matrix{{sub(((coefficients f)_1)_(0,0),QQ)}};
         	 b=matrix{{sub(((coefficients f)_1)_(1,0),QQ)}};
         	 a=(solve(A,b))_(0,0);
          	 subs=join({sub(varList_i,Rtotal)=>a},subs);
		 );
	 --check this is really a solution with certain tolerance;
	 if (sum apply(flatten entries gens sub(partialSol,subs),i->if abs(sub(i,RR))<tol then 0 else 1))>0 then return ("No partial solution");    
         );    
         return sub(S,drop(subs,-n));
      )
  
checkRankCompletion=(compl,rk,tol)->(  
    L=apply(flatten entries gens minors(rk+1,compl),i->sub(i,RR));
    if sum apply(L,i->if i>tol then 1 else 0)>0 then return("This is not a rank completion")
    else return(eigenvalues sub(compl,QQ));	
    )

------------------------------------------------
-- TBD
------------------------------------------------
-- generate point on the intersection of complement  
