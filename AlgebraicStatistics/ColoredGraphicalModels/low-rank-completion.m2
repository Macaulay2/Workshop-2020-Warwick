restart

--compute conormal variety
conormalVariety=(IX,n,x,u)->(
    c:=codim IX;
    JacX:=diff(matrix{toList(x_1..x_n)}, transpose gens IX);
    AugJacX:=matrix{toList(u_1..u_n)}||JacX;
    SingX:=IX+minors(c,jacobian IX);
    conormalX:=saturate(IX+minors(c+1,AugJacX),SingX);
    conormalX
    )

--Uhler 2010 Graph 7
n=5
R=QQ[l_1..l_n,t_1..t_n]
I=ideal(4*t_2^2*t_3^2-4*t_1*t_2*t_4^2+t_4^4-8*t_1*t_2*t_4*t_5+2*t_4^3*t_5-4*t_1*t_2*t_5^2+2*t_4^2*t_5^2+2*t_4*t_5^3+t_5^4)
CN=conormalVariety(I,n,t,l)  

diff(matrix{toList(t_1..t_n)}, transpose gens I)
submatrix(transpose jacobian I,toList(0..n-1))

R=QQ[l_1..l_n,t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
K=matrix{{l_1,l_3,0,l_3},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_3,0,l_5,l_1}}
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

boundaryComponents(K,1)
boundaryComponents(K,2)
boundaryComponents(K,3)
I=(boundaryComponents(K,4))_0
minimalPrimes (minors(4,K))
dualVariety(I,n,t,l)
dim I
codim I

I=dualVariety(ideal(l_1^2+l_2^2),n,l,t)
J=dualVariety(I,n,t,l)
I==dualVariety(J,n,l,t)
