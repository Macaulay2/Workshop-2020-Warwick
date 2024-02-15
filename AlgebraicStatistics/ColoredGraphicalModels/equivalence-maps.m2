--Fret's heads
--standard computation of the cone of sufficient statistics
restart
n=5
R=QQ[t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--sufficient statistics
stat1=t_1-s_11-s_22
stat2=t_2-s_33-s_44
stat3=t_3-2*s_12
stat4=t_4-2*(s_23+s_14)
stat5=t_5-2*s_34

varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
Istat=ideal(stat1,stat2,stat3,stat4,stat5)

I_G2=eliminate(varList,minors(3,S)+Istat)--0
I_G1=eliminate(varList,minors(2,S)+Istat)

--map composition

restart
n=7
R=QQ[t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--sufficient statistics
stat1=t_1-s_11-s_22
stat2=t_2-s_33-s_44
stat3=t_3-2*s_12
stat4=t_4-2*(s_23+s_14)
stat5=t_5-2*s_34

varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
Istat=ideal(stat1,stat2,stat3,stat4,stat5)

--eliminate missing edges
I_G2proj=eliminate({s_13,s_24},minors(3,S))
I_G2stats=I_G2proj+Istat
I_G2composition=eliminate({s_11,s_12,s_14,s_22,s_23,s_33,s_34,s_44},I_G2stats) --0


I_G1proj=eliminate({s_13,s_24},minors(2,S))
I_G1stats=I_G2proj+Istat
I_G1composition=eliminate({s_11,s_12,s_14,s_22,s_23,s_33,s_34,s_44},I_G2stats) 

I_G1composition=ideal(4*t_2^2*t_3^2-4*t_1*t_2*t_4^2+t_4^4+8*t_1*t_2*t_3*t_5-4*t_3*t_4^2*t_5+4*t_1^2*t_5^2)

I_G1== I_G1composition
--auxiliary stats
stat6=t_6-2*s_13
stat7=t_7-2*s_24


IstatAll=ideal(stat1,stat2,stat3,stat4,stat5,stat6,stat7)
I1=minors(2,S)
I1stats=I_1+IstatAll
varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
I1statsOnly=eliminate(varList,I1stats)
eliminate({t_6,t_7},I1statsOnly)
