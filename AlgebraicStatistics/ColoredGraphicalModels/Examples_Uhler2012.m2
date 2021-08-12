--Example 6.2 in Uhler (Example 5.1 in Sturmfels and Uhler)
restart
R=QQ[l_1..l_5,t_1..t_5,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
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
decompI2=primaryDecomposition I2
length decompI2
minPrimesI2=minimalPrimes I2
n=length minPrimesI2
boundaryP2=for i to n-1 list{
J=(minimalPrimes I2)_i;
codim J;
jacobian J;
M=matrix{{t_1,t_2,t_3,t_4,t_5}}|| submatrix(transpose jacobian J,{0,1,2,3,4});
conormal=I2+minors(codim J+1,M);
conormal=saturate(conormal, minors(codim I2,submatrix(transpose jacobian J,{0,1,2,3,4})));
eliminate({l_1,l_2,l_3,l_4,l_5},conormal)}

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
