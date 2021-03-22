restart

-- STEP 1: Setup concentration matrix associated to colored graph (LSSM)
m=4
L=flatten for i from 1 to binomial(m+1,2) list l_i
Raux=QQ[L]

--MULTIVARIATE GAUSSIANS paper (Sturmfels-Uhler), TABLE 2:

-- Graph 1 
K=matrix{{l_1,l_2,0,l_2},{l_2,l_1,l_3,0},{0,l_3,l_1,l_2},{l_2,0,l_2,l_1}}
R=QQ[support K]
K=sub(K,R)
n=dim R   --dim d in table

-- Graph 2
K=matrix{{l_1,l_3,0,l_3},{l_3,l_2,l_4,0},{0,l_4,l_1,l_3},{l_3,0,l_3,l_2}}
R=QQ[support K]
K=sub(K,R)
n=dim R

-- STEP 2: ML-degree computation
X=random(QQ^m,QQ^m);              
S=(1/m)*X*transpose(X);
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(S*K)}})};
J=saturate(I,det K);
degree J -- ML-degree in table

-- STEP 3: degree of the "covariance" variety (reciprocal variety of the LSSM)
Lcov=flatten flatten for i from 1 to m list for j from i to m list s_{i,j}
FullRing=QQ[gens R,Lcov]
Lcon=take(gens FullRing, n)
Cov=genericSymmetricMatrix(FullRing,s_{1,1},m)
K=sub(K,FullRing)
ICov=ideal flatten flatten (K*Cov-id_(FullRing^m))
JCov=eliminate(Lcon,ICov)
CovRing=QQ[Lcov]
JCov=sub(JCov,CovRing)
dim JCov==n
degree(JCov)  -- degree in table
betti trim JCov  -- mingens P_G in table
