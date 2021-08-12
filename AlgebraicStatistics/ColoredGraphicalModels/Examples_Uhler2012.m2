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
I=minors(p,K)
decompI=primaryDecomposition I
length decompI
minPrimesI=minimalPrimes I
n=length minPrimesI
boundaryP2=for i to n-1 list{
J=(minimalPrimes I)_i;
codim J;
jacobian J;
M=matrix{{t_1,t_2,t_3,t_4,t_5}}|| submatrix(transpose jacobian J,{0,1,2,3,4});
conormal=I+minors(codim J+1,M);
conormal=saturate(conormal, minors(codim I,submatrix(transpose jacobian J,{0,1,2,3,4})));
eliminate({l_1,l_2,l_3,l_4,l_5},conormal)}

p=3
I=minors(p,K)
decompI=primaryDecomposition I
length decompI
minPrimesI=minimalPrimes I
n=length minPrimesI
boundaryP3=for i to n-1 list{
J=(minimalPrimes I)_i;
codim J;
jacobian J;
M=matrix{{t_1,t_2,t_3,t_4,t_5}}|| submatrix(transpose jacobian J,{0,1,2,3,4});
conormal=I+minors(codim J+1,M);
conormal=saturate(conormal, minors(codim I,submatrix(transpose jacobian J,{0,1,2,3,4})));
eliminate({l_1,l_2,l_3,l_4,l_5},conormal)}