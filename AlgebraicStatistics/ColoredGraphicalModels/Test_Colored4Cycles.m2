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

--STEP 4: elimination criterion

Cov=sub(Cov,CovRing)
I3=minors(4,Cov)
eliminate({s_{1,3},s_{2,4}},I3) 
-- Graph 1: 0 ideal => with a sample of size n=3 the MLE exists with probability 1

I2=minors(3,Cov)

eliminate({s_{1,3},s_{2,4}},I2) 
isPrime I2
netList I2_*
-- Graph 1: non-zero ideal, elimination criterion says nothing about sample of size n=2

--STEP 5: cone of sufficient statistics


--Compute the dual variety of det(K)=0 in the sense of projective algebraic geometry
T=QQ[l_1,l_2,l_3]
K=sub(K,T)
boundaryK=minors(4,K)
primaryDecomposition(boundaryK)
minimalPrimes boundaryK
I1=(minimalPrimes boundaryK)_0
c1=codim I1
T2=QQ[gens T,t_1,t_2,t_3]
jacI1=sub(jacobian(I1),T2)
matI1=jacI1|matrix{{t_1},{t_2},{t_3}}
I1=sub(I1,T2)
J1=I1+minors(1,matI1)
J1=saturate(J1,minors(1,jacI1))

locus1=minors(2,K)
primaryDecomposition(locus1)
minimalPrimes locus1
locus2=minors(3,K)
primaryDecomposition(locus2)
minimalPrimes locus2

-----Tests
R2=QQ[Lcov,t_1,t_2,t_3,l_1,l_2,l_3]
gens R2

K=sub(K,R2)
K1=sub(K,{l_1=>1,l_2=>0,l_3=>0})
K2=sub(K,{l_1=>0,l_2=>1,l_3=>0})
K3=sub(K,{l_1=>0,l_2=>0,l_3=>1})
Cov=sub(Cov,R2)

suf=ideal{t_1-trace(Cov*K1),t_2-trace(Cov*K2),t_3-trace(Cov*K3)}

R3=QQ[Lcov,l_1,l_2,l_3,p_1,p_2,p_3,p_4]
K=sub(K,R3)
Cov=sub(Cov,R3)
--K PD matrix: ppal leading minors are positive
-- l_1>0
l_1*p_1-1
-- l_1^2-l_2^2>0
(l_1^2-l_2^2)*p_2-1
-- 3rd leading minor
det(submatrix'(K,{3},{3}))*p_3-1
-- det>0
det(K)*p_4-1
coneSuf=ideal{trace(Cov*K),l_1*p_1-1,(l_1^2-l_2^2)*p_2-1,det(submatrix'(K,{3},{3}))*p_3-1,det(K)*p_4-1}
coneSufStat=eliminate({p_1,p_2,p_3,p_4},coneSuf)
eliminate({l_1,l_2,l_3},coneSuf)





--Graph 1
netList primaryDecomposition ideal{det K}
netList JCov_*
isPrime JCov

R2=QQ[Lcov,t_1,t_2,t_3]
suf=ideal{t_1-s_{1,1},
t_2-(2*(s_{1,2}+s_{1,4}+s_{3,4})),
t_3-2*s_{2,3}}+sub(JCov,R2)
isPrime suf

A=sub(suf,{t_1=>3,t_2=>5,t_3=>-2})
ACov=sub(A,CovRing)
dim ACov,degree ACov

sufel=eliminate({t_1,t_2,t_3},suf)
sufel==sub(JCov,R2)

toString Lcov
sufel=eliminate({s_{1, 1}, s_{1, 2}, s_{1, 3}, s_{1, 4}, s_{2, 2}, s_{2, 3}, s_{2, 4}, s_{3, 3},
      s_{4, 4}},suf)

sufel==sub(JCov,R2)
dim suf, degree suf

S2=QQ[Lcov,{t_1,t_2,t_3},MonomialOrder=>Lex]
describe S2
suf=ideal{t_1-s_{1,1},
t_2-(2*(s_{1,2}+s_{1,4}+s_{3,4})),
t_3-2*s_{2,3}}+sub(JCov,S2)
aux=flatten entries gens gb suf
netList aux

R3=(QQ[t_1..t_3])[Lcov]
gens R3
suf=ideal{t_1-s_{1,1},
t_2-(2*(s_{1,2}+s_{1,4}+s_{3,4})),
t_3-2*s_{2,3}}+sub(JCov,R3)
dim suf, degree suf
