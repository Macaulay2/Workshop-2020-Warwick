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
restart
dualVariety=(IX,n,x,u)->(
    c:=codim IX;
    JacX:=submatrix(transpose jacobian IX,toList(0..n-1));
    AugJacX:=matrix{toList(u_1..u_n)}||JacX;
    SingX:=minors(c,JacX);
    conormalX:=saturate(IX+minors(c+1,AugJacX),SingX);
    dualX:=eliminate(toList(x_1..x_n),conormalX);
    dualX
    )

--Example 6.2 in Uhler (Example 5.1 in Sturmfels and Uhler)
--restart
n=5
R=QQ[l_1..l_n,t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
K=matrix{{l_1,l_3,0,l_4},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_4,0,l_5,l_2}}
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--sufficient statistics
stat1=t_1-s_11-s_22
stat2=t_2-s_33-s_44
stat3=t_3-2*s_12
stat4=t_4-2*(s_23+s_14)
stat5=t_5-2*s_34

-- compute P_G
P_G=eliminate({l_1,l_2,l_3,l_4,l_5},ideal(K*S-id_(R^4)))

--algebraic boundary
p=2
I2=minors(p,K)
minPrimesI2=minimalPrimes I2
m=length minPrimesI2
boundaryP2=for i to m-1 list dualVariety(minPrimesI2_i,n,l,t)

p=3
I3=minors(p,K)
--jacI3 = diff(matrix{{l_1..l_5}}, transpose gens I3)
--singI3 = I3 + minors(codim I3,jacI3)
--saturate(I3,singI3)
trueI3=saturate(saturate(I3,I2),minors(1,K))
--trueI3=saturate(I3,I2)
--decompI3=primaryDecomposition trueI3
--length decompI3
--minPrimesI3=minimalPrimes trueI3
--n=length minPrimesI3
--boundaryP3=for i to n-1 list{
--J=(minimalPrimes I)_i;
--codim J;
--jacobian J;
--M=matrix{{t_1,t_2,t_3,t_4,t_5}}|| submatrix(transpose jacobian J,{0,1,2,3,4});
--conormal=I+minors(codim J+1,M);
--conormal=saturate(conormal, minors(codim I,submatrix(transpose jacobian J,{0,1,2,3,4})));
--eliminate({l_1,l_2,l_3,l_4,l_5},conormal)}

J=(minimalPrimes trueI3)_0
codim J
jacobian J
M=matrix{{t_1,t_2,t_3,t_4,t_5}}|| submatrix(transpose jacobian J,{0,1,2,3,4})
conormal=trueI3+minors(codim J+1,M)
conormal=saturate(conormal, minors(codim trueI3,submatrix(transpose jacobian J,{0,1,2,3,4})))
eliminate({l_1,l_2,l_3,l_4,l_5},conormal)

-- ideal(t_1-t_3)*(t_1+t_3)*(t_2-t_5)*(t_2+t_5)
--I_G2=eliminate({s_11,s_12,s_13,s_24},minors(3,S)+ideal(stat1,stat2,stat3,stat4,stat5))
Iminors=minors(3,S)
Iinverse=ideal(K*S-id_(R^4))
I_G2=eliminate({s_13,s_24},Iminors+Iinverse)
I_G2=eliminate({s_11,s_12,s_13,s_14,s_22,s_23,s_24,s_33,s_34,s_44},I_G2)
--leadTerm minors(3,S)

gens gb (minors(3,S)+ideal(stat1,stat2,stat3,stat4,stat5))

Itest=ideal(K*S-id_(R^4),stat1,stat2,stat3,stat4,stat5)+minors(3,S)
Jtest=eliminate({l_1,l_2,l_3,l_4,l_5},Itest)
eliminate({s_13,s_24},Jtest)
