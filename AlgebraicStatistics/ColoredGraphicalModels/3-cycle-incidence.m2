restart
R=QQ[sigma_11..sigma_33,s_11..s_33]
I=ideal(sigma_11+sigma_12-s_11-s_22,sigma_13-s_13,sigma_23-s_23,sigma_33-s_33)
J=ideal(s_11*s_22-s_12^2)
K=eliminate(toList(s_11..s_33),I+J)

A=matrix{{s_11,s_12,s_13},{s_12,s_22,s_23},{s_13,s_23,s_33}}
J=ideal(det A)
