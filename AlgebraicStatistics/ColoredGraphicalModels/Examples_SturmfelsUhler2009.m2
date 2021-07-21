restart
R=QQ[l_1..l_3,s_11..s_13,s_22,s_23,s_33,t_1..t_3,x_11..x_13,x_22,x_23,x_33]
gens R
K=matrix{{l_1+l_2+l_3,l_3,l_2},{l_3,l_1+l_2+l_3,l_1},{l_2,l_1,l_1+l_2+l_3}}
--F=frac R
S=matrix{{s_11,s_12,s_13},{s_12,s_22,s_23},{s_13,s_23,s_33}}
Sigma=matrix{{x_11,x_12,x_13},{x_12,x_22,x_23},{x_13,x_23,x_33}}
poly1=t_1-trace S+2*s_23
poly2=t_2-trace S+ 2*s_13
poly3=t_2-trace S + 2*s_12
poly4=trace(Sigma*K)-trace(S*K)
M=Sigma*K
--I=ideal(M_(0,0)-1,M_(0,1),M_(0,2),M_(1,0),M_(1,1)-1,M_(1,2),M_(2,0),M_(2,1),M_(2,2)-1)
--I=ideal(poly1,poly2,poly3,poly4)
