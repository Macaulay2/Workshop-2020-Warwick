-- Example 4.10 in Uhler
restart
load "./functions.m2"
n=12
R=QQ[l_1..l_n,t_1..t_n,s_11..s_15,s_22..s_25,s_33..s_35,s_44..s_45,s_55]
K=matrix{{l_1,l_6,0,l_9,l_10},{l_6,l_2,l_7, 0,l_11},{0,l_7,l_3,l_8,l_12},{l_9,0,l_8,l_4,0},{l_10,l_11,l_12,0,l_5}}
S=matrix{{s_11,s_12,s_13,s_14,s_15},{s_12,s_22,s_23,s_24,s_25},{s_13,s_23,s_33,s_34,s_35},{s_14,s_24,s_34,s_44,s_45},{s_15,s_25,s_35,s_45,s_55}}--estimate

p=2
J1=minors(p,K)
minPrimes1=minimalPrimes J1
dualVariety(minPrimes1_0,n,l,t)

p=3
J2=minors(p,K)
minPrimes2=minimalPrimes J2
minPrimes2={ideal(l_7*l_10-l_6*l_12,-l_7*l_11+l_2*l_12,l_2*l_10-l_6*l_11,l_5*l_6-l_
      10*l_11,l_5*l_7-l_11*l_12,l_2*l_5-l_11^2,l_4*l_6*l_7+l_2*l_8*l_9,l_8*l_9
      *l_11+l_4*l_6*l_12,l_5*l_8*l_9+l_4*l_10*l_12,l_3*l_9*l_11+l_6*l_8*l_12-l
      _7*l_9*l_12,l_6*l_7*l_8+l_2*l_3*l_9-l_7^2*l_9,l_3*l_4*l_11-l_8^2*l_11-l_
      4*l_7*l_12,l_3*l_4*l_6-l_6*l_8^2+l_7*l_8*l_9,l_3*l_4*l_10-l_8^2*l_10+l_8
      *l_9*l_12,l_2*l_3*l_4-l_4*l_7^2-l_2*l_8^2,l_3*l_5*l_9+l_8*l_10*l_12-l_9*
      l_12^2,l_3*l_4*l_5-l_5*l_8^2-l_4*l_12^2,-l_6*l_8*l_10+l_1*l_8*l_11+l_6*l
      _9*l_12,l_1*l_7*l_8+l_3*l_6*l_9,l_3*l_9*l_10+l_1*l_8*l_12,l_1*l_2*l_8-l_
      6^2*l_8+l_6*l_7*l_9,-l_4*l_6*l_10+l_1*l_4*l_11-l_9^2*l_11,l_1*l_4*l_7+l_
      6*l_8*l_9-l_7*l_9^2,l_8*l_9*l_10+l_1*l_4*l_12-l_9^2*l_12,l_1*l_2*l_4-l_4
      *l_6^2-l_2*l_9^2,l_1*l_5*l_8-l_8*l_10^2+l_9*l_10*l_12,l_1*l_4*l_5-l_5*l_
      9^2-l_4*l_10^2,-l_3*l_6*l_10+l_1*l_3*l_11-l_1*l_7*l_12,l_1*l_2*l_3-l_3*l
      _6^2-l_1*l_7^2,l_1*l_3*l_4-l_1*l_8^2-l_3*l_9^2,l_1*l_3*l_5-l_3*l_10^2-l_
      1*l_12^2), ideal(l_8,l_9,l_4,l_5*l_6*l_7-l_7*l_10*l_11+l_2*l_10*l_12-l_6
      *l_11*l_12,l_2*l_3*l_10-l_7^2*l_10-l_3*l_6*l_11+l_6*l_7*l_12,l_3*l_5*l_6
      -l_3*l_10*l_11+l_7*l_10*l_12-l_6*l_12^2,l_2*l_3*l_5-l_5*l_7^2-l_3*l_11^2
      +2*l_7*l_11*l_12-l_2*l_12^2,l_6*l_7*l_10-l_1*l_7*l_11+l_1*l_2*l_12-l_6^2
      *l_12,l_1*l_5*l_7-l_7*l_10^2+l_6*l_10*l_12-l_1*l_11*l_12,l_1*l_2*l_5-l_5
      *l_6^2-l_2*l_10^2+2*l_6*l_10*l_11-l_1*l_11^2,-l_3*l_6*l_10+l_1*l_3*l_11-
      l_1*l_7*l_12,l_1*l_2*l_3-l_3*l_6^2-l_1*l_7^2,l_1*l_3*l_5-l_3*l_10^2-l_1*
      l_12^2)}

p=4
J3=minors(p,K)
minPrimes2=minimalPrimes J3

--sufficient statistics
varList=flatten {toList(s_11..s_15), toList(s_22..s_25),toList(s_33..s_35),s_44,s_45,s_55}
Istat=ideal(t_1-s_11, t_2-s_22, t_3-s_33,t_4-s_44,t_5-s_55,t_6-2*s_12,t_7-2*s_23,t_8-2*s_34,t_9-2*s_14,t_10-2*s_15,t_11-2*s_25,t_12-2*s_35)

I_G3=eliminate(varList,minors(4,S)+Istat)
I_G2=eliminate(varList,minors(3,S)+Istat)
I_G1=eliminate(varList,minors(2,S)+Istat)
