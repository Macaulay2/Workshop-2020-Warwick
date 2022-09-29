restart
needsPackage "GraphicalModelsMLE"

--FUNCTIONS

--m: rk of empirical covariance matrix
--k: number of points we want to try	
empiricalMLEexistence=(m,k,K)->(
    count:=0;
    for i from 1 to k do (    
    	p=numcols K;
    	--X=random(QQ^p,QQ^m);              
    	X=random(RR^p,RR^m);  
	X=transpose matrix {toList apply(0..p-1,i->promote(X_(i,0),QQ))};            
	S=(1/m)*X*transpose(X);
    	I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
    	J=saturate(I,det K);
	if dim J!=0 then return(count,"The ideal does not have the right dimension ");
    	criticalPoints=zeroDimSolve(J);
        criticalMatrices=genListMatrix(criticalPoints,K);
	if checkPD(criticalMatrices,ZeroTolerance=>1e-7)!={} then count=count+1;
        );
        return(count);
     )	

adjugateMLEexistence=(m,k,K)->(
    count:=0;
    for i from 1 to k do (    
    	p=numcols K;         
    	X=random(QQ^m,QQ^p);            
	S=exteriorPower(m,transpose(X)*X);
    	I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
    	J=saturate(I,det K);
	if dim J!=0 then return(count,"The ideal does not have the right dimension ");
    	criticalPoints=zeroDimSolve(J);
        criticalMatrices=genListMatrix(criticalPoints,K);
	if checkPD(criticalMatrices,ZeroTolerance=>1e-7)!={} then count=count+1;
        );
        return(count);
     )	
 
--auxiliary functions (internal in package "GraphicalModelsMLE")
genListMatrix = (L,A) ->
(
		T:= for l in L list coordinates(l);
		M:= for t in T list substitute(A,matrix{t});
    return M
);
 
-------------------------------------------------
--EXAMPLE
-------------------------------------------------

-- APPROACH 1 (DOESN'T WORK): 
-- 1) Generate k  observations (reals promoted to rationals)
-- 2) Compute the sample covariance matrix S=X*transpose(X)
-- 3) Compute the MLE of S

-- this approach works when all vertices are colored (30%) but never when 
-- all edges are colored

R=QQ[l_1..l_4]
k=100
K=matrix{{l_1,l_4,l_4},{l_4,l_2,l_4},{l_4,l_4,l_3}} -- all edges
--K=matrix{{l_1,l_2,l_3},{l_2,l_1,l_4},{l_3,l_4,l_1}} -- all vertices
empiricalMLEexistence(1,k,K) 

-- APPROACH 1 (WORKS approx. 50% of the time): 
-- 1) Generate k rank 1 PSD matrices that are adjugates of rank 2 PSD matrices
-- 2) Compute the MLE of S

adjugateMLEexistence(2,k,K)

-- Specific example of 2x3 matrix whose adjugate is a sample covariance matrix
-- for which MLE exists for all edges:
X=matrix {{8/5, 2/3, 8/9}, {7/6, 1/3, 5}}
S=exteriorPower(2,transpose(X)*X)
rank S

--compute the MLE
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
J=saturate(I,det K);
if dim J!=0 then return(count,"The ideal does not have the right dimension ");
criticalPoints=zeroDimSolve(J);
criticalMatrices=genListMatrix(criticalPoints,K);
sol=checkPD(criticalMatrices,ZeroTolerance=>1e-7)
SigmaHat=inverse sol_0

-- confirm correctness of solution
S_(0,0)-SigmaHat_(0,0)<1e-7
S_(1,1)-SigmaHat_(1,1)<1e-7
S_(2,2)-SigmaHat_(2,2)<1e-7
(S_(0,1)+S_(0,2)+S_(1,2))-(SigmaHat_(0,1)+SigmaHat_(0,2)+SigmaHat_(1,2))<1e-7
(sol_0)_(0,1)==(sol_0)_(0,2)
(sol_0)_(0,1)==(sol_0)_(1,2)

--for reference
-- sol_0= matrix {{18.4865, -.143707, -.143707}, {-.143707, .0754344, -.143707},
--      {-.143707, -.143707, .391758}}
--matrix {{.0597531, .516614, .211426}, {.516614, 48.4829, 17.9743},
--      {.211426, 17.9743, 9.22359}}