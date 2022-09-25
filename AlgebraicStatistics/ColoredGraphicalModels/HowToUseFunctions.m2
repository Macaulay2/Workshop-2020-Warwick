-- THIS FILE IS MEANT TO BE USED AS A TEMPLATE TO USE THE PACKAGE functions.m2
-- FOR ANY COLORED GRAPHICAL MODEL THAT WE WANT TO STUDY
-- It has the following content:

-- PART I: GENERAL SETTING
-- (elimination ideals, boundary components)
-- PART II: FIND RK 1 POINTS IN THE INTERSECTION OF THE ALGEBRAIC BOUNDARY AND THE INTERIOR
-- (change of signs, points in the intersection, rank completion)
-- PART III: EMPIRICAL CHECK OF MLE EXISTENCE OF RK 1 POINTS 
-- (random rk 1 points and those coming from adjugates of random rk 2 points)
-- PART IV: DEEPER INTO CRITICAL POINTS 
-- (ML-degree depending on rk)


restart
load "functions.m2"
load "setupReciprocalVarieties.m2"

--PART I: GENERAL SETTING 
R=QQ[l_1..l_4]
K=matrix{{l_1,l_4,l_4},{l_4,l_2,l_4},{l_4,l_4,l_3}}
(p,n,Rtotal,S)=coloredData(K)

--Sufficient statistics
stats=sub(suffStat(K),Rtotal)

--Computation of elimination ideals IG_n
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)));
netList IG1_*

--Computation of boundary components and algebraic boundary
(V,n,K2)=embeddedK(K);
boundaryComponents(K2,3,n)
boundaryComponents(K2,2,n)
boundaryComponents(K2,1,n)
algBoundary(V,n,K2)



--PART II: FIND RK 1 POINTS IN THE INTERSECTION OF THE ALGEBRAIC BOUNDARY AND THE INTERIOR

--Step 1: Generate rank 1 matrices whose sufficient statistics lie 
--in the intersection between the algebraic boundary and the interior
--of the cone of sufficient statistics (i.e. rk 1 points for which MLE exists)
use ring IG1
(L1,L2)=differentSign(IG1_0,3,100,p,n,stats);
netList L1 --initial point (suff statistics, evaluation)
netList L2 --if exists, point with different sign (suff statistics, evaluation, number of points we had to try)
v1=vector(flatten entries gens L1_0);
v2=vector(flatten entries gens L2_0);
P=interiorPoint(v1,v2,n,IG1_0,K)
intP=P_0 
use ring IG1
sub(sub(IG1,{t_1=>intP_0,t_2=>intP_1,t_3=>intP_2,t_4=>intP_3}),RR) 
--Careful with the tolerance on the rank 1 matrix completion: 
--must be higher than the zero tolerance of intP in IG1
rk1=rankCompletion(intP,p,1,stats,S,Rtotal,0.000001)
--check whether it is rank 1 and positive definite
rank rk1
eigenvalues sub(rk1,QQ)

--Step 2: Check whether the previous rank 1 matrix lies in the
--subvariety of rk 1 matrices of K_G^{-1}
aux=flatten toList apply(1..3,i->(toList apply(i..3,j->sub(s_(i,j),Rtotal)=>rk1_(i-1,j-1))));
netList (sub(sub(I1,aux),RR))_*

--Step 3: Check whether the MLE exists
MLE(K,intP)
--the output is either the ML-degree, if it exists, or "MLE doesn't exist", if it doesn't

--NOTE: we need to be careful here because it can be that we get the MLE because
--those matrices are actually rk3 (approximately rk1 but not exactly, 
-- check substitution of coordinates in IG1 and/or eigenvalues)
--This seems to be supported by the fact that the ML-degree here is usually 4, 
--instead of 3, as will happen with adjugates of rk2 matrices


--PART III: EMPIRICAL CHECK OF MLE EXISTENCE OF RK 1 POINTS

--Check MLE existence of random rk 1 points
empiricalMLEexistence(1,500,K)  --0
--Only checking if there are critical points
empiricalMLENoPD(1,500,K)  --500
--Enforce real and non-singular to the previous case
empiricalMLERealReg(1,500,K)  --498,499

--Check MLE existence of rk 1 matrices coming from adjugates of rk 2 matrices
k=20
for i to k do(
    U=random(QQ^2,QQ^3);
    A=transpose(U)*U;
    adjA=exteriorPower(2,A);
    adjA=adjugate(A);
    aux=flatten toList apply(1..3,i->(toList apply(i..3,j->sub(s_(i,j),Rtotal)=>adjA_(i-1,j-1))));
    suff=sub(stats,aux);
    use Rtotal;
    coord=-flatten entries gens sub(suff,{t_1=>0,t_2=>0,t_3=>0,t_4=>0});
    print MLE(K,coord);
)

--one iteration line-by-line
    U=random(QQ^2,QQ^3)
    A=transpose(U)*U
    adjA=adjugate(A)
    aux=flatten toList apply(1..3,i->(toList apply(i..3,j->sub(s_(i,j),Rtotal)=>adjA_(i-1,j-1))))
    suff=sub(stats,aux)
    use Rtotal
    coord=-flatten entries gens sub(suff,{t_1=>0,t_2=>0,t_3=>0,t_4=>0})
    print MLE(K,coord)
--the output is either the ML-degree, if it exists, or "MLE doesn't exist", if it doesn't

--NOTE: both in empiricalMLEexistence and adjugates, we are considering truly rk 1 matrices,
--as opposed to PART II, and the ML-degree seems to be always 3.

--example for MLE to exist and MLdegree=3
{441/256, 729/400, 289/100, -4407/800}


--PART IV: DEEPER INTO CRITICAL POINTS 

--Compute all details for given sufficient statistics
use R
coord --using last one, define it yourself if desired
score=ideal{jacobian(matrix{{det K}})-det(K)*sub(transpose(matrix{coord}),R)};
dim score, degree score
scoreEq=saturate(score,det K);
dim scoreEq, degree scoreEq
criticalPoints=zeroDimSolve(scoreEq);
criticalMatrices=genListMatrix(criticalPoints,K)
--note that all previous matrices indeed satisfy edge equality and have the same sufficient statistics
--for different graphs we need to change the expression to compute sufficient statistics
netList apply(criticalMatrices,i->{M_(0,0),M_(1,1),M_(2,2),2*M_(0,1)+2*M_(0,2)+2*M_(1,2)})
--however, at most one is positive definite
apply(criticalMatrices,i->eigenvalues i)
checkPD(criticalMatrices)
MLEK=criticalMatrices_2 --pick the only PD matrix, if it exists
M=inverse MLEK --MLE for the covariance matrix

--Test MLE for different ranks
use R
p=3
n=1
X=random(QQ^p,QQ^n);              
S=(1/n)*X*transpose(X);
rank S
eigenvalues S
minors(2,S) --that's what matters for M2 in terms of computing the rank

I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
netList I_*
dim I, degree I
--(1,4)
J=saturate(I,det K);
dim J, degree J
--n=3 MLdeg=4
--n=2 MLdeg=4
--n=1 MLdeg=3
