-----------------------------------
-----------------------------------
--G={(1,2,3),(12,13,23)}
-----------------------------------
-----------------------------------
restart
n=6
R=QQ[t_1..t_n,s_11..s_13,s_22..s_23,s_33]
S=matrix{{s_11,s_12,s_13},{s_12,s_22,s_23},{s_13,s_23,s_33}}--estimate

--sufficient statistics
--uncolored {12,13,23}
varList=flatten {toList(s_11..s_13), toList(s_22..s_23),s_33}
Istat=ideal(t_1-s_11,t_2-s_22,t_3-s_33,t_4-2*s_12,t_5-2*s_13,t_6-2*s_23)

I_G3=eliminate(varList,minors(4,S)+Istat)
I_G2=eliminate(varList,minors(3,S)+Istat)
I_G1=eliminate(varList,minors(2,S)+Istat)

restart
n=5
R=QQ[t_1..t_n,s_11..s_13,s_22..s_23,s_33]
S=matrix{{s_11,s_12,s_13},{s_12,s_22,s_23},{s_13,s_23,s_33}}--estimate

--sufficient statistics
--color vertices 2 and 3 
varList=flatten {toList(s_11..s_13), toList(s_22..s_23),s_33}
Istat=ideal(t_1-s_11,t_2-s_22-s_33,t_3-2*s_12,t_4-2*s_13,t_5-2*s_23)

I_G3=eliminate(varList,minors(4,S)+Istat)
I_G2=eliminate(varList,minors(3,S)+Istat)
I_G1=eliminate(varList,minors(2,S)+Istat)

restart
n=5
R=QQ[t_1..t_n,s_11..s_13,s_22..s_23,s_33]
S=matrix{{s_11,s_12,s_13},{s_12,s_22,s_23},{s_13,s_23,s_33}}--estimate

--sufficient statistics
--color edges  12 and 13 
varList=flatten {toList(s_11..s_13), toList(s_22..s_23),s_33}
Istat=ideal(t_1-s_11,t_2-s_22,t_3-s_33,t_4-2*s_12-2*s_13,t_5-2*s_23)

I_G3=eliminate(varList,minors(4,S)+Istat)
I_G2=eliminate(varList,minors(3,S)+Istat)
I_G1=eliminate(varList,minors(2,S)+Istat)

restart
n=4
R=QQ[t_1..t_n,s_11..s_13,s_22..s_23,s_33]
S=matrix{{s_11,s_12,s_13},{s_12,s_22,s_23},{s_13,s_23,s_33}}--estimate

--sufficient statistics
--color all vertices 
varList=flatten {toList(s_11..s_13), toList(s_22..s_23),s_33}
Istat=ideal(t_1-s_11-s_22-s_33,t_2-2*s_12,t_3-2*s_13,t_4-2*s_23)

I_G3=eliminate(varList,minors(4,S)+Istat)
I_G2=eliminate(varList,minors(3,S)+Istat)
I_G1=eliminate(varList,minors(2,S)+Istat)

restart
n=4
R=QQ[t_1..t_n,s_11..s_13,s_22..s_23,s_33]
S=matrix{{s_11,s_12,s_13},{s_12,s_22,s_23},{s_13,s_23,s_33}}--estimate

--sufficient statistics
--color all edges
varList=flatten {toList(s_11..s_13), toList(s_22..s_23),s_33}
Istat=ideal(t_1-s_11,t_2-s_22,t_3-s_33,t_4-2*s_12-2*s_13-2*s_23)

I_G3=eliminate(varList,minors(4,S)+Istat)
I_G2=eliminate(varList,minors(3,S)+Istat)
I_G1=eliminate(varList,minors(2,S)+Istat)

restart
n=4
R=QQ[t_1..t_n,s_11..s_13,s_22..s_23,s_33]
S=matrix{{s_11,s_12,s_13},{s_12,s_22,s_23},{s_13,s_23,s_33}}--estimate

--sufficient statistics
--color vertices 2 and 3 and edges 13 and 23
varList=flatten {toList(s_11..s_13), toList(s_22..s_23),s_33}
Istat=ideal(t_1-s_11,t_2-s_22-s_33,t_3-2*s_12,t_4-2*s_13-2*s_23)

I_G3=eliminate(varList,minors(4,S)+Istat)
I_G2=eliminate(varList,minors(3,S)+Istat)
I_G1=eliminate(varList,minors(2,S)+Istat)

restart
n=4
R=QQ[t_1..t_n,s_11..s_13,s_22..s_23,s_33]
S=matrix{{s_11,s_12,s_13},{s_12,s_22,s_23},{s_13,s_23,s_33}}--estimate

--sufficient statistics
--color vertices 2 and 3 and edges 12 and 13
varList=flatten {toList(s_11..s_13), toList(s_22..s_23),s_33}
Istat=ideal(t_1-s_11,t_2-s_22-s_33,t_3-2*s_12-2*s_13,t_4-2*s_23)

I_G3=eliminate(varList,minors(4,S)+Istat)
I_G2=eliminate(varList,minors(3,S)+Istat)
I_G1=eliminate(varList,minors(2,S)+Istat)
------------------------------
------------------------------
--G={(1,2,3,4),(12,13,23,34)}
------------------------------
------------------------------
restart
n=8
R=QQ[t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--sufficient statistics
--uncolored {12,13,23,34}
varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
Istat=ideal(t_1-s_11,t_2-s_22,t_3-s_33,t_4-s_44,t_5-2*s_12,t_6-2*s_13,t_7-2*s_23,t_8-2*s_34)

I_G3=eliminate(varList,minors(4,S)+Istat)
I_G2=eliminate(varList,minors(3,S)+Istat)
I_G1=eliminate(varList,minors(2,S)+Istat)

restart
n=7
R=QQ[t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--colored {12,13,23,34} with vertices 2 and 4 with the same color
varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
Istat=ideal(t_1-s_11,t_2-s_22-s_44,t_3-s_33,t_4-2*s_12,t_5-2*s_13,t_6-2*s_23,t_7-2*s_34)

I_G3=eliminate(varList,minors(4,S)+Istat)
I_G2=eliminate(varList,minors(3,S)+Istat)
I_G1=eliminate(varList,minors(2,S)+Istat)


restart
n=7
R=QQ[t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--colored {12,13,23,34} with edges  12  and 34 with the same color
varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
Istat=ideal(t_1-s_11,t_2-s_22,t_3-s_33,t_4-s_44,t_5-2*s_12-2*s_34,t_6-2*s_13,t_7-2*s_23)

I_G3=eliminate(varList,minors(4,S)+Istat)
I_G2=eliminate(varList,minors(3,S)+Istat)
I_G1=eliminate(varList,minors(2,S)+Istat)

restart
n=7
R=QQ[t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--colored {12,13,23,34} with edges  13  and 34 with the same color
varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
Istat=ideal(t_1-s_11,t_2-s_22,t_3-s_33,t_4-s_44,t_5-2*s_12,t_6-2*s_13-2*s_34,t_7-2*s_23)

I_G3=eliminate(varList,minors(4,S)+Istat)
I_G2=eliminate(varList,minors(3,S)+Istat)
I_G1=eliminate(varList,minors(2,S)+Istat)


restart
R=QQ[x_1..x_4]
X=matrix{{x_1,x_2},{x_2,x_1},{x_3,x_4},{x_4,x_3}}
S=X*transpose X

restart
n=4
R=QQ[l_1..l_n,t_1..t_n,s_11..s_13,s_22..s_23,s_33]
K=matrix{{l_1,l_4,l_4},{l_4,l_2,l_4},{l_4,l_4,l_3}}
S=matrix{{s_11,s_12,s_13},{s_12,s_22,s_23},{s_13,s_23,s_33}}--estimate

BC=boundaryComponents(K,3)
algBoundary(K)

 I2=minors(3,K)
 minPrimes2=minimalPrimes I2
dualVariety(minPrimes2_0,n,l,t)

decompose I2
