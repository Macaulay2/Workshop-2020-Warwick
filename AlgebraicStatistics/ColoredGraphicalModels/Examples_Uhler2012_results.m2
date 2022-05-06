----------------------------------------------------------------------
----------------------------------------------------------------------
-- Uhler 2012, Table 2
----------------------------------------------------------------------
----------------------------------------------------------------------

-- GRAPH 1
n=3
R=QQ[l_1..l_n,t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
K=matrix{{l_1,l_2,0,l_2},{l_2,l_1,l_3,0},{0,l_3,l_1,l_2},{l_2,0,l_2,l_1}}
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--ideal of variety of L^{-1}
P_G= ideal(s_22-s_33,s_13-s_24,s_12-s_34,s_11-s_44,s_24^2+s_24*s_33-s_23*s_
      34-s_34^2,s_23*s_24+s_14*s_33-s_24*s_34-s_33*s_34,s_14*s_24+s_14*s_33-
      s_24*s_34-s_34*s_44,s_23^2-s_33^2-s_34^2+s_33*s_44,s_14*s_23+s_24*s_33
      +s_14*s_34-s_23*s_34-s_34^2-s_24*s_44)
--algebraic boundary  
H_G={ideal(t_2+t_3), ideal t_3, ideal(2*t_1^2-t_2^2+2*t_2*t_3-t_3^2),
      ideal(t_1^2-2*t_1*t_2+t_2^2-6*t_1*t_3+2*t_2*t_3+5*t_3^2),
      ideal(t_1^2+2*t_1*t_2+t_2^2+6*t_1*t_3+2*t_2*t_3+5*t_3^2)}
  
--sufficient statistics
stat1=t_1-s_11-s_22-s_33-s_44
stat2=t_2-2*(s_12+s_14+s_34)
stat3=t_3-2*s_23

--I_Gn
I_G3=ideal()
I_G2=ideal()
I_G1=ideal()  

-- GRAPH 7
restart
n=5
R=QQ[l_1..l_n,t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
K=matrix{{l_1,l_3,0,l_3},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_3,0,l_5,l_1}}
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--ideal of variety of L^{-1}
P_G=ideal(s_23^2-s_22*s_33-s_34^2+s_33*s_44,s_14*s_23-s_13*s_24+s_14*s_34-
      s_13*s_44,s_13*s_23-s_12*s_33+s_14*s_33-s_13*s_34,s_14*s_22-s_12*s_24+
      s_14*s_24-s_12*s_44,s_13*s_22-s_12*s_23+s_13*s_24-s_12*s_34,s_13*s_14-
      s_23*s_24-s_11*s_34+s_22*s_34,s_13^2-s_11*s_33+s_24*s_33-s_23*s_34-s_
      34^2+s_33*s_44,s_12*s_13-s_11*s_23-s_24*s_34+s_23*s_44,s_14^2*s_33-s_
      24^2*s_33-s_11*s_34^2+s_22*s_34^2-s_24*s_33*s_44+s_23*s_34*s_44,s_11*s
      _12*s_14+s_11*s_14^2-s_11^2*s_24-s_12^2*s_24-s_14^2*s_24+s_11*s_22*s_
      24+s_24^3-s_11^2*s_44-s_12^2*s_44-s_14^2*s_44+s_11*s_22*s_44+s_11*s_24
      *s_44-s_22*s_24*s_44+s_24^2*s_44+s_11*s_44^2-s_22*s_44^2,s_11*s_12^2-s
      _11*s_14^2-s_11^2*s_22-s_12^2*s_22+s_11*s_22^2-s_12*s_14*s_24+s_14^2*s
      _24+s_22*s_24^2+s_11^2*s_44+s_12^2*s_44-s_12*s_14*s_44+s_14^2*s_44-s_
      22^2*s_44-s_24^2*s_44-s_11*s_44^2+s_22*s_44^2)
  
-- algebraic boundary 
H_G=  {ideal t_2, ideal(64*t_1^2*t_2^2-32*t_2^2*t_3^2-16*t_1*t_2*t_4^2+t_4^4
      +32*t_1*t_2*t_4*t_5-4*t_4^3*t_5-16*t_1*t_2*t_5^2+6*t_4^2*t_5^2-4*t_4*t
      _5^3+t_5^4), ideal(4*t_2^2*t_3^2-4*t_1*t_2*t_4^2+t_4^4-8*t_1*t_2*t_4*t
      _5+2*t_4^3*t_5-4*t_1*t_2*t_5^2+2*t_4^2*t_5^2+2*t_4*t_5^3+t_5^4)}

--I_Gn
I_G3=ideal()
I_G2=ideal()
I_G1=ideal(4*t_2^2*t_3^2-4*t_1*t_2*t_4^2+t_4^4-8*t_1*t_2*t_4*t_5+2*t_4^3*t_
      5-4*t_1*t_2*t_5^2+2*t_4^2*t_5^2+2*t_4*t_5^3+t_5^4)
  
-- GRAPH 9
restart
n=6
R=QQ[l_1..l_n,t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
K=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

-- compute P_G (ideal of variety of L^{-1})
P_G= ideal(s_23^2-s_22*s_33-s_34^2+s_33*s_44,s_14*s_23-s_13*s_24+s_14*s_34-
      s_13*s_44,s_13*s_23-s_12*s_33+s_14*s_33-s_13*s_34,s_14*s_22-s_12*s_24+
      s_14*s_24-s_12*s_44,s_13*s_22-s_12*s_23+s_13*s_24-s_12*s_34,s_12*s_14*
      s_33-s_11*s_24*s_33-s_12*s_13*s_34+s_13*s_14*s_34+s_11*s_23*s_34-s_13^
      2*s_44)

-- algebraic boundary 
H_G= {ideal t_1, ideal t_3, ideal(4*t_2*t_3-t_5^2-t_6^2),
      ideal(t_3*t_4^2-t_1*t_5^2-2*t_1*t_5*t_6-t_1*t_6^2),
      ideal(8*t_1*t_2*t_3-t_3*t_4^2-t_1*t_5^2+2*t_1*t_5*t_6-t_1*t_6^2)}
  
--I_Gn
I_G3=ideal()
I_G2=ideal()
I_G1=ideal(4*t_2*t_3-t_5^2-t_6^2,t_3*t_4^2-t_1*t_5^2-2*t_1*t_5*t_6-t_1*t_6^2
     )


-- GRAPH 8
restart
n=5
R=QQ[t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--sufficient statistics
varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
Istat=ideal(t_1-s_11-s_33,t_2-s_22-s_44,t_3-2*s_12-2*s_14,t_4-2*s_23,t_5-2*s_34)

--I_Gn
I_G3=ideal()
I_G2=ideal()
I_G1=ideal(4*t_1*t_2*t_4^2-t_3^2*t_4^2-t_4^4+8*t_1*t_2*t_4*t_5-2*t_4^3*t_5+4*t_1*t_2*t_5^2-t_3^2*t_5^2-2*t_4^2*t_5^2-2*t_4*t_5^3-t_5^4)

-- GRAPH 11
restart
n=5
R=QQ[l_1..l_n,t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
K=matrix{{l_1,l_3,0,l_4},{l_3,l_2,l_4,0},{0,l_4,l_1,l_5},{l_4,0,l_5,l_2}}
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate
--sufficient statistics
varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
Istat=ideal(t_1-s_11-s_33,t_2-s_22-s_44,t_3-2*s_12,t_4-2*s_23-2*s_14,t_5-2*s_34)


--algebraic boundary
H_G= {ideal(4*t_1*t_2-t_3^2-2*t_3*t_5-t_5^2),
      ideal(4*t_1*t_2-t_3^2-t_4^2+2*t_3*t_5-t_5^2),
      ideal(t_4^2-4*t_3*t_5)}

--I_Gn
I_G3=()
I_G2=()
I_G1= ideal(4*t_1*t_2-t_3^2-t_4^2+2*t_3*t_5-t_5^2)

-- GRAPH 13
restart
n=6
R=QQ[l_1..l_n,t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
K=matrix{{l_1,l_3,0,l_6},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_6,0,l_5,l_2}}
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--sufficient statistics
varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
Istat=ideal(t_1-s_11-s_22,t_2-s_33-s_44,t_3-2*s_12,t_4-2*s_23,t_5-2*s_34,t_6-2*s_14)

--algebraic boundary
H_G={ideal(t_2+t_5), ideal(t_2-t_5), ideal(t_1+t_3),
     ideal(t_1-t_3), ideal(t_3*t_5-t_4*t_6),
     ideal(4*t_2^2*t_3^2-4*t_1*t_2*t_4^2+t_4^4+8*t_1*t_2*t_3*t_5-4*t
     _3*t_4^2*t_5+4*t_1^2*t_5^2-8*t_1*t_2*t_4*t_6+4*t_4^3*t_6-8*t_3*
     t_4*t_5*t_6-4*t_1*t_2*t_6^2+6*t_4^2*t_6^2-4*t_3*t_5*t_6^2+4*t_4
     *t_6^3+t_6^4), ideal(4*t_2^2*t_3^2-4*t_1*t_2*t_4^2+t_4^4-8*t_1*
     t_2*t_3*t_5+4*t_3*t_4^2*t_5+4*t_1^2*t_5^2+8*t_1*t_2*t_4*t_6-4*t
     _4^3*t_6-8*t_3*t_4*t_5*t_6-4*t_1*t_2*t_6^2+6*t_4^2*t_6^2+4*t_3*
     t_5*t_6^2-4*t_4*t_6^3+t_6^4)}

--I_Gn
I_G3=ideal()
I_G2=ideal()
I_G1=ideal(t_3*t_5-t_4*t_6,4*t_2^2*t_3^2-4*t_1*t_2*t_4^2+t_4^4+4*t_1
     ^2*t_5^2-4*t_1*t_2*t_6^2-2*t_4^2*t_6^2+t_6^4)
 
  
-- GRAPH 16
restart
n=7
R=QQ[l_1..l_n,t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
K=matrix{{l_1,l_4,0,l_7},{l_4,l_1,l_5,0},{0,l_5,l_2,l_6},{l_7,0,l_6,l_3}}
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

-- compute P_G (ideal of variety of L^{-1})
P_G=ideal(s_13*s_14-s_23*s_24-s_11*s_34+s_22*s_34,s_14^2*s_33-s_24^2*s_33-s
     _11*s_34^2+s_22*s_34^2+s_13^2*s_44-s_23^2*s_44-s_11*s_33*s_44+s_22*s_33
     *s_44,s_14*s_23*s_24-s_13*s_24^2-s_14*s_22*s_34+s_12*s_24*s_34+s_13*s_
     22*s_44-s_12*s_23*s_44,s_13^2*s_24-s_23^2*s_24+s_12*s_14*s_33-s_11*s_24
     *s_33-s_12*s_13*s_34+s_22*s_23*s_34)

-- algebraic boundary 
H_G=  

--I_Gn
I_G3= ideal()
I_G2= ideal()
I_G1= ideal(t_4*t_6-t_5*t_7,4*t_2*t_3-t_6^2,t_3*t_5^2-t_1*t_6^2+t_2*t
     _7^2,4*t_1*t_2*t_5*t_6-t_5^3*t_6-4*t_2^2*t_4*t_7,4*t_3^2*t_4*t_
     5-4*t_1*t_3*t_6*t_7+t_6*t_7^3,4*t_3^2*t_4^2-4*t_1*t_3*t_7^2+t_7
     ^4,4*t_2^2*t_4^2-4*t_1*t_2*t_5^2+t_5^4)
 
-- GRAPH 17
restart
n=7
R=QQ[l_1..l_n,t_1..t_n]
K=matrix{{l_1,l_4,0,l_7},{l_4,l_2,l_5,0},{0,l_5,l_1,l_6},{l_7,0,l_6,l_3}}  

restart
n=7
R=QQ[t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--sufficient statistics
varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
Istat=ideal(t_1-s_11-s_33,t_2-s_22,t_3-s_44,t_4-2*s_12,t_5-2*s_23,t_6-2*s_34,t_7-2*s_14)

--I_Gn
I_G3= ideal()
I_G2= ideal()
I_G1=ideal(t_4*t_6-t_5*t_7,4*t_1*t_3-t_6^2-t_7^2,4*t_1*t_2-t_4^2-t_5
     ^2,t_3*t_5^2-t_2*t_6^2,t_3*t_4*t_5-t_2*t_6*t_7)
 
-- GRAPH 18 
restart
n=8
R=QQ[l_1..l_n,t_1..t_n]
K=matrix{{l_1,l_5,0,l_8},{l_5,l_2,l_6,0},{0,l_6,l_3,l_7},{l_8,0,l_7,l_4}}

restart
n=8
R=QQ[t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--sufficient statistics
varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
Istat=ideal(t_1-s_11,t_2-s_22,t_3-s_33,t_4-s_44,t_5-2*s_12,t_6-2*s_23,t_7-2*s_34,t_8-2*s_14)

--I_Gn
I_G3= ideal()
I_G2=  ideal(t_3^2*t_4^2*t_5^4-2*t_1*t_3*t_4^2*t_5^2*t_6^2+t_1^2*t_4^2
     *t_6^4-2*t_1*t_2*t_3*t_4*t_5^2*t_7^2-2*t_1^2*t_2*t_4*t_6^2*t_7^
     2+t_1*t_4*t_5^2*t_6^2*t_7^2+t_1^2*t_2^2*t_7^4+8*t_1*t_2*t_3*t_4
     *t_5*t_6*t_7*t_8-t_3*t_4*t_5^3*t_6*t_7*t_8-t_1*t_4*t_5*t_6^3*t_
     7*t_8-t_1*t_2*t_5*t_6*t_7^3*t_8-2*t_2*t_3^2*t_4*t_5^2*t_8^2-2*t
     _1*t_2*t_3*t_4*t_6^2*t_8^2+t_3*t_4*t_5^2*t_6^2*t_8^2-2*t_1*t_2^
     2*t_3*t_7^2*t_8^2+t_2*t_3*t_5^2*t_7^2*t_8^2+t_1*t_2*t_6^2*t_7^2
     *t_8^2-t_2*t_3*t_5*t_6*t_7*t_8^3+t_2^2*t_3^2*t_8^4)
 
I_G1=ideal(t_5*t_7-t_6*t_8,4*t_3*t_4-t_7^2,4*t_1*t_4-t_8^2,4*t_2*t_
      3-t_6^2,4*t_1*t_2-t_5^2,t_1*t_6*t_7-t_3*t_5*t_8,t_4*t_5*t_6-t_
      2*t_7*t_8)