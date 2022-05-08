-------------
--FUNCTIONS--
-------------

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
    JacX:=submatrix(transpose jacobian IX,toList(0..n-1));
    AugJacX:=matrix{toList(u_1..u_n)}||JacX;
    SingX:=minors(c,JacX);
    conormalX:=saturate(IX+minors(c+1,AugJacX),SingX);
    dualX:=eliminate(toList(x_1..x_n),conormalX);
    dualX
    )

-----------------------------
-- ALGEBRAIC BOUNDARY
-----------------------------
--algebraic boundary H_G
-- See p. 7 from Sturmfels and Uhler paragraph after Proposition 2.4
-- Algorithm for K a s x s matrix
-- 1. For all 1 <=p<=s compute I:= the ideal of the p-minors of K
-- 2. Compute the dual variety of each minimal prime of I
-- 3. Keep only those varieties, whose ideal is principal
-- 4. H_L is the product of these principal generators

boundaryComponents=(K,p)->(
    I:=minors(p,K);
    minPrimes:=minimalPrimes I;
    m:= length minPrimes;
    allComponents:=for i to  m-1 list (dualComponent:=dualVariety(minPrimes_i,n,l,t),numgens trim dualComponent);
    boundaryComponents:= for i in allComponents list (if i_1==1 then i_0)
    )


algBoundary=(K)->(
    s=numgens target K;       
    delete (null, flatten (for p from 1 to s list boundaryComponents(K,p))) 
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