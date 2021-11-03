restart
R=QQ[l_1..l_3,s_11..s_13,s_22,s_23,s_33,t_1..t_3,x_11..x_13,x_22,x_23,x_33]
R=QQ[l_1..l_3,t_1..t_3,x_11..x_13,x_22,x_23,x_33]
gens R
K=matrix{{l_1+l_2+l_3,l_3,l_2},{l_3,l_1+l_2+l_3,l_1},{l_2,l_1,l_1+l_2+l_3}}
--F=frac R
--K=sub(K,F)
--inverse K
S=matrix{{s_11,s_12,s_13},{s_12,s_22,s_23},{s_13,s_23,s_33}}-- sample covariance
Sigma=matrix{{x_11,x_12,x_13},{x_12,x_22,x_23},{x_13,x_23,x_33}}--estimate
-- compute P_L
I=ideal(K*Sigma-id_(R^3))
eliminate({l_1,l_2,l_3},I)

-- compute Sigma hat
poly1=t_1-trace Sigma-2*x_23
poly2=t_2-trace Sigma- 2*x_13
poly3=t_3-trace Sigma - 2*x_12
--poly4=trace(Sigma*K)-trace(S*K)
--I=ideal(K*Sigma-id_(R^3),poly1,poly2,poly3,poly4)
I=ideal(K*Sigma-id_(R^3),poly1,poly2,poly3)
J=eliminate({l_1,l_2,l_3},I)
--J=eliminate({l_1,l_2,l_3,s_11,s_12,s_13,s_22,s_23,s_33},I)
eliminate({x_11,x_12,x_13,x_22,x_23},J)
--I=ideal(M_(0,0)-1,M_(0,1),M_(0,2),M_(1,0),M_(1,1)-1,M_(1,2),M_(2,0),M_(2,1),M_(2,2)-1)
--I=ideal(poly1,poly2,poly3,poly4)

--algebraic boundary
R1=QQ[l_1..l_3,t_1..t_3]
K=matrix{{l_1+l_2+l_3,l_3,l_2},{l_3,l_1+l_2+l_3,l_1},{l_2,l_1,l_1+l_2+l_3}}
--I=minors(2,K)
I=minors(3,K)
J=(minimalPrimes I)_0
codim J
jacobian J
M=matrix{{t_1,t_2,t_3}}|| submatrix(transpose jacobian J,{0,1,2})
conormal=I+minors(codim J+1,M)
conormal=saturate(conormal, minors(codim I,submatrix(transpose jacobian J,{0,1,2})))
eliminate({l_1,l_2,l_3},conormal)
