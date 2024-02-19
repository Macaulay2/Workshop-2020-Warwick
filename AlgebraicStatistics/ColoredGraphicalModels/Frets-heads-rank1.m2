--Example 6.2 in Uhler (Example 5.1 in Sturmfels and Uhler)
-- adding rank 1 constraints
restart

n=7
R=QQ[t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
--K=matrix{{l_1,l_3,0,l_4},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_4,0,l_5,l_2}}
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--sufficient statistics
stat1=t_1-s_11-s_22
stat2=t_2-s_33-s_44
stat3=t_3-2*s_12
stat4=t_4-2*(s_23+s_14)
stat5=t_5-2*s_34
stat6=t_6-2*s_13
stat7=t_7-2*s_24

varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
Istat=ideal(stat1,stat2,stat3,stat4,stat5,stat6,stat7)
Istatorig=ideal(stat1,stat2,stat3,stat4,stat5)

I_G2=eliminate(varList,minors(3,S)+Istat)
I_G1=eliminate(varList,minors(2,S)+Istat)
I_G1orig=eliminate(varList,minors(2,S)+Istatorig)

eliminate({t_6,t_7},I_G1)==I_G1orig

--uncolored

restart
n=10
R=QQ[t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
--K=matrix{{l_1,l_3,0,l_4},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_4,0,l_5,l_2}}
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--sufficient statistics
stat1=t_1-s_11
stat2=t_2-s_22
stat3=t_3-s_33
stat4=t_4-s_44
stat5=t_5-2*s_12
stat6=t_6-2*s_23
stat7=t_7-2*s_14
stat8=t_8-2*s_34
stat9=t_9-2*s_13
stat10=t_10-2*s_24

varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
Istat=ideal(stat1,stat2,stat3,stat4,stat5,stat6,stat7,stat8,stat9,stat10)
Istatorig=ideal(stat1,stat2,stat3,stat4,stat5,stat6,stat7,stat8)

I_G1=eliminate(varList,minors(2,S)+Istat)
I_G1orig=eliminate(varList,minors(2,S)+Istatorig)