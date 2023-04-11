restart
R=QQ[s_11..s_33,d_11..d_33]
--1. 2 colored vertices in a 3-cycle
I=ideal(d_11+d_22-s_11-s_22,
d_33-s_33,d_12-s_12,d_13-s_13,d_23-s_23,d_22*d_33-d_23^2-d_11*d_33-d_13^2)
Irank=ideal(s_11*s_22-s_12^2, s_11*s_33-s_13^2, s_22*s_33-s_23^2)
eliminate(toList(s_11..s_33),I+Irank)
-- contains d11d33, so the matrix cannot be PD

--2. Graph 11
restart
R=QQ[s_11..s_44,d_11..d_44]
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}
adjS=exteriorPower_3 S
D=matrix{{d_11,d_12,d_13,d_14},{d_12,d_22,d_23,d_24},{d_13,d_23,d_33,d_34},{d_14,d_24,d_34,d_44}}
adjD=exteriorPower_3 D

I=ideal(d_11-s_11 + d_33-s_33,d_22+d_44-s_22-s_44,d_12-s_12,
d_14+d_23-s_14-s_23, d_34-s_34, adjD_(0,2), adjD_(1,3),
adjD_(3,3)-adjD_(1,1), adjD_(2,2)-adjD_(0,0), adjD_(1,2)-adjD_(0,3))
Irank=minors(2,S)
J=eliminate(toList(s_11..s_44),I+Irank)
Jdecomp=decompose J
isSubset(minors(4,D),Jdecomp_0) --true
isSubset(minors(4,D),Jdecomp_1) --false
isSubset(minors(3,D),Jdecomp_1) --false
isSubset(minors(2,D),Jdecomp_1) --false
eliminate({d_12,d_13,d_14,d_23,d_24,d_34},Jdecomp_1)
eliminate({d_11,d_22,d_33,d_44},Jdecomp_1)
eliminate({d_14,d_24,d_34,d_44},Jdecomp_1)

J=eliminate(toList(d_11..d_44),I+Irank)
--results for elimination of s variables
ideal(d_12^2+d_14^2-d_11*d_22+2*d_14*d_23+d_23^2-d_22*d_33-2*d_12*d_34+d_
     34^2-d_11*d_44-d_33*d_44,d_14*d_23*d_24-d_13*d_24^2-d_14*d_22*d_34+d_12*d_
     24*d_34+d_13*d_22*d_44-d_12*d_23*d_44,d_14^2*d_23-d_14*d_23^2-d_13*d_14*d_
     24+d_13*d_23*d_24+d_14*d_22*d_33-d_12*d_24*d_33-d_12*d_14*d_34-d_13*d_22*d
     _34+d_12*d_23*d_34+d_11*d_24*d_34+d_12*d_13*d_44-d_11*d_23*d_44,d_13*d_14*
     d_23-d_13^2*d_24-d_12*d_14*d_33+d_11*d_24*d_33+d_12*d_13*d_34-d_11*d_23*d_
     34,d_14^2*d_22-2*d_12*d_14*d_24+d_11*d_24^2-d_24^2*d_33+2*d_23*d_24*d_34-d
     _22*d_34^2-d_14^2*d_44-2*d_14*d_23*d_44-2*d_23^2*d_44+2*d_22*d_33*d_44+2*d
     _12*d_34*d_44-d_34^2*d_44+d_11*d_44^2+d_33*d_44^2,d_13^2*d_22-2*d_12*d_13*
     d_23+d_11*d_23^2-2*d_14^2*d_33-2*d_14*d_23*d_33-d_23^2*d_33+d_22*d_33^2+2*
     d_13*d_14*d_34+2*d_12*d_33*d_34-d_11*d_34^2-d_33*d_34^2-d_13^2*d_44+2*d_11
     *d_33*d_44+d_33^2*d_44)

--primal decomposition
Jdecomp={ideal(d_12^2+d_14^2-d_11*d_22+2*d_14*d_23+d_23^2-d_22*d_33-2*d_12*d_34+d
      _34^2-d_11*d_44-d_33*d_44,d_24^2*d_33-2*d_23*d_24*d_34+d_22*d_34^2+d_23^2
      *d_44-d_22*d_33*d_44,d_14*d_24*d_33-d_14*d_23*d_34-d_13*d_24*d_34+d_12*d_
      34^2+d_13*d_23*d_44-d_12*d_33*d_44,d_14^2*d_33-2*d_13*d_14*d_34+d_11*d_34
      ^2+d_13^2*d_44-d_11*d_33*d_44,d_14*d_23*d_24-d_13*d_24^2-d_14*d_22*d_34+d
      _12*d_24*d_34+d_13*d_22*d_44-d_12*d_23*d_44,d_14*d_23^2-d_13*d_23*d_24-d_
      14*d_22*d_33+d_12*d_24*d_33+d_13*d_22*d_34-d_12*d_23*d_34,d_14^2*d_23-d_
      13*d_14*d_24-d_12*d_14*d_34+d_11*d_24*d_34+d_12*d_13*d_44-d_11*d_23*d_44,
      d_13*d_14*d_23-d_13^2*d_24-d_12*d_14*d_33+d_11*d_24*d_33+d_12*d_13*d_34-d
      _11*d_23*d_34,d_14^2*d_22-2*d_12*d_14*d_24+d_11*d_24^2-d_14^2*d_44-2*d_14
      *d_23*d_44-d_23^2*d_44+d_22*d_33*d_44+2*d_12*d_34*d_44-d_34^2*d_44+d_11*d
      _44^2+d_33*d_44^2,d_13*d_14*d_22-d_12*d_14*d_23-d_12*d_13*d_24+d_11*d_23*
      d_24-d_14^2*d_34-2*d_14*d_23*d_34-d_23^2*d_34+d_22*d_33*d_34+2*d_12*d_34^
      2-d_34^3+d_11*d_34*d_44+d_33*d_34*d_44,d_13^2*d_22-2*d_12*d_13*d_23+d_11*
      d_23^2-2*d_14*d_23*d_33-d_23^2*d_33+d_22*d_33^2-2*d_13*d_14*d_34+2*d_12*d
      _33*d_34+d_11*d_34^2-d_33*d_34^2+d_13^2*d_44+d_33^2*d_44),
      ideal(d_14-d_23,d_24*d_33-d_13*d_44,d_22*d_33-d_11*d_44,d_22*d_23-d_12*d_
      24+d_24*d_34-d_23*d_44,d_13*d_22-d_11*d_24,d_12*d_13-d_11*d_23+d_23*d_33-
      d_13*d_34,d_12^2-d_11*d_22+4*d_23^2-2*d_12*d_34+d_34^2-2*d_11*d_44-d_33*d
      _44,2*d_23*d_24*d_34+d_22*d_34^2+2*d_13*d_24*d_44-2*d_12*d_34*d_44+d_34^2
      *d_44-3*d_11*d_44^2-d_33*d_44^2,2*d_13*d_23*d_34-2*d_12*d_33*d_34+d_11*d_
      34^2+d_33*d_34^2+2*d_13^2*d_44-3*d_11*d_33*d_44-d_33^2*d_44,2*d_23*d_24^2
      +d_22*d_24*d_34-3*d_12*d_24*d_44+4*d_24*d_34*d_44-4*d_23*d_44^2,2*d_13*d_
      24^2+d_12*d_24*d_34-3*d_24*d_34^2-3*d_11*d_24*d_44+4*d_23*d_34*d_44-d_13*
      d_44^2,2*d_12*d_24^2+d_22^2*d_34-2*d_24^2*d_34-3*d_12*d_22*d_44+3*d_22*d_
      34*d_44-d_12*d_44^2,2*d_11*d_24^2+d_12*d_22*d_34-3*d_22*d_34^2-3*d_11*d_
      22*d_44-2*d_13*d_24*d_44+3*d_12*d_34*d_44-d_34^2*d_44+2*d_11*d_44^2+d_33*
      d_44^2,2*d_23^2*d_24+d_12*d_24*d_34-d_24*d_34^2-2*d_12*d_23*d_44-d_11*d_
      24*d_44+2*d_23*d_34*d_44-d_13*d_44^2,2*d_13*d_23*d_24+d_11*d_24*d_34-3*d_
      11*d_23*d_44-d_23*d_33*d_44+d_13*d_34*d_44,2*d_12*d_23*d_24+d_12*d_22*d_
      34-3*d_11*d_22*d_44+8*d_23^2*d_44+2*d_13*d_24*d_44-3*d_12*d_34*d_44+2*d_
      34^2*d_44-7*d_11*d_44^2-2*d_33*d_44^2,2*d_11*d_23*d_24+d_11*d_22*d_34-3*d
      _11*d_12*d_44-2*d_13*d_23*d_44-d_12*d_33*d_44+3*d_11*d_34*d_44,2*d_13^2*d
      _24+d_11*d_23*d_34+3*d_23*d_33*d_34-2*d_13*d_34^2-3*d_11*d_13*d_44-d_13*d
      _33*d_44,2*d_11*d_13*d_24+d_11*d_12*d_34+3*d_12*d_33*d_34-3*d_11*d_34^2-d
      _33*d_34^2-3*d_11^2*d_44-2*d_13^2*d_44+2*d_11*d_33*d_44+d_33^2*d_44,2*d_
      13*d_23^2-2*d_12*d_23*d_33+d_11*d_23*d_34+d_23*d_33*d_34-d_11*d_13*d_44-d
      _13*d_33*d_44,2*d_11*d_23^2+6*d_23^2*d_33+d_11*d_12*d_34-d_12*d_33*d_34-d
      _11*d_34^2+d_33*d_34^2-3*d_11^2*d_44-4*d_11*d_33*d_44-d_33^2*d_44,2*d_13^
      2*d_23-3*d_11*d_23*d_33-d_23*d_33^2+d_11*d_13*d_34+d_13*d_33*d_34,2*d_11*
      d_13*d_23-3*d_11*d_12*d_33-2*d_13*d_23*d_33-d_12*d_33^2+d_11^2*d_34+3*d_
      11*d_33*d_34)}

J11=ideal(d_14-d_23,d_24*d_33-d_13*d_44,d_22*d_33-d_11*d_44,d_22*d_23-d_12*d_
      24+d_24*d_34-d_23*d_44,d_13*d_22-d_11*d_24,d_12*d_13-d_11*d_23+d_23*d_33-
      d_13*d_34,d_12^2-d_11*d_22+4*d_23^2-2*d_12*d_34+d_34^2-2*d_11*d_44-d_33*d
      _44,2*d_23*d_24*d_34+d_22*d_34^2+2*d_13*d_24*d_44-2*d_12*d_34*d_44+d_34^2
      *d_44-3*d_11*d_44^2-d_33*d_44^2,2*d_13*d_23*d_34-2*d_12*d_33*d_34+d_11*d_
      34^2+d_33*d_34^2+2*d_13^2*d_44-3*d_11*d_33*d_44-d_33^2*d_44,2*d_23*d_24^2
      +d_22*d_24*d_34-3*d_12*d_24*d_44+4*d_24*d_34*d_44-4*d_23*d_44^2,2*d_13*d_
      24^2+d_12*d_24*d_34-3*d_24*d_34^2-3*d_11*d_24*d_44+4*d_23*d_34*d_44-d_13*
      d_44^2,2*d_12*d_24^2+d_22^2*d_34-2*d_24^2*d_34-3*d_12*d_22*d_44+3*d_22*d_
      34*d_44-d_12*d_44^2,2*d_11*d_24^2+d_12*d_22*d_34-3*d_22*d_34^2-3*d_11*d_
      22*d_44-2*d_13*d_24*d_44+3*d_12*d_34*d_44-d_34^2*d_44+2*d_11*d_44^2+d_33*
      d_44^2,2*d_23^2*d_24+d_12*d_24*d_34-d_24*d_34^2-2*d_12*d_23*d_44-d_11*d_
      24*d_44+2*d_23*d_34*d_44-d_13*d_44^2,2*d_13*d_23*d_24+d_11*d_24*d_34-3*d_
      11*d_23*d_44-d_23*d_33*d_44+d_13*d_34*d_44,2*d_12*d_23*d_24+d_12*d_22*d_
      34-3*d_11*d_22*d_44+8*d_23^2*d_44+2*d_13*d_24*d_44-3*d_12*d_34*d_44+2*d_
      34^2*d_44-7*d_11*d_44^2-2*d_33*d_44^2,2*d_11*d_23*d_24+d_11*d_22*d_34-3*d
      _11*d_12*d_44-2*d_13*d_23*d_44-d_12*d_33*d_44+3*d_11*d_34*d_44,2*d_13^2*d
      _24+d_11*d_23*d_34+3*d_23*d_33*d_34-2*d_13*d_34^2-3*d_11*d_13*d_44-d_13*d
      _33*d_44,2*d_11*d_13*d_24+d_11*d_12*d_34+3*d_12*d_33*d_34-3*d_11*d_34^2-d
      _33*d_34^2-3*d_11^2*d_44-2*d_13^2*d_44+2*d_11*d_33*d_44+d_33^2*d_44,2*d_
      13*d_23^2-2*d_12*d_23*d_33+d_11*d_23*d_34+d_23*d_33*d_34-d_11*d_13*d_44-d
      _13*d_33*d_44,2*d_11*d_23^2+6*d_23^2*d_33+d_11*d_12*d_34-d_12*d_33*d_34-d
      _11*d_34^2+d_33*d_34^2-3*d_11^2*d_44-4*d_11*d_33*d_44-d_33^2*d_44,2*d_13^
      2*d_23-3*d_11*d_23*d_33-d_23*d_33^2+d_11*d_13*d_34+d_13*d_33*d_34,2*d_11*
      d_13*d_23-3*d_11*d_12*d_33-2*d_13*d_23*d_33-d_12*d_33^2+d_11^2*d_34+3*d_
      11*d_33*d_34)
isSubset(ideal(d_33*d_44-d_34^2),J11)
isSubset(ideal(d_22*d_44-d_24^2),J11)
isSubset(ideal(d_22*d_33-d_23^2),J11)
isSubset(ideal(d_11*d_44-d_14^2),J11)
isSubset(ideal(d_11*d_22-d_12^2),J11)
isSubset(ideal(d_11*d_33-d_13^2),J11)
isSubset(ideal det submatrix'(D,{0},{0}),J11)
isSubset(ideal det submatrix'(D,{1},{1}),J11)
isSubset(ideal det submatrix'(D,{2},{2}),J11)
isSubset(ideal det submatrix'(D,{3},{3}),J11)

restart
R=QQ[s_11..s_44,d_11..d_44,t_1..t_5]
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}
adjS=exteriorPower_3 S
D=matrix{{d_11,d_12,d_13,d_14},{d_12,d_22,d_23,d_24},{d_13,d_23,d_33,d_34},{d_14,d_24,d_34,d_44}}
adjD=exteriorPower_3 D

I=ideal(d_11-s_11 + d_33-s_33,d_22+d_44-s_22-s_44,d_12-s_12,
d_14+d_23-s_14-s_23, d_34-s_34, adjD_(0,2), adjD_(1,3),
adjD_(3,3)-adjD_(1,1), adjD_(2,2)-adjD_(0,0), adjD_(1,2)-adjD_(0,3))
J11=ideal(d_14-d_23,d_24*d_33-d_13*d_44,d_22*d_33-d_11*d_44,d_22*d_23-d_12*d_
      24+d_24*d_34-d_23*d_44,d_13*d_22-d_11*d_24,d_12*d_13-d_11*d_23+d_23*d_33-
      d_13*d_34,d_12^2-d_11*d_22+4*d_23^2-2*d_12*d_34+d_34^2-2*d_11*d_44-d_33*d
      _44,2*d_23*d_24*d_34+d_22*d_34^2+2*d_13*d_24*d_44-2*d_12*d_34*d_44+d_34^2
      *d_44-3*d_11*d_44^2-d_33*d_44^2,2*d_13*d_23*d_34-2*d_12*d_33*d_34+d_11*d_
      34^2+d_33*d_34^2+2*d_13^2*d_44-3*d_11*d_33*d_44-d_33^2*d_44,2*d_23*d_24^2
      +d_22*d_24*d_34-3*d_12*d_24*d_44+4*d_24*d_34*d_44-4*d_23*d_44^2,2*d_13*d_
      24^2+d_12*d_24*d_34-3*d_24*d_34^2-3*d_11*d_24*d_44+4*d_23*d_34*d_44-d_13*
      d_44^2,2*d_12*d_24^2+d_22^2*d_34-2*d_24^2*d_34-3*d_12*d_22*d_44+3*d_22*d_
      34*d_44-d_12*d_44^2,2*d_11*d_24^2+d_12*d_22*d_34-3*d_22*d_34^2-3*d_11*d_
      22*d_44-2*d_13*d_24*d_44+3*d_12*d_34*d_44-d_34^2*d_44+2*d_11*d_44^2+d_33*
      d_44^2,2*d_23^2*d_24+d_12*d_24*d_34-d_24*d_34^2-2*d_12*d_23*d_44-d_11*d_
      24*d_44+2*d_23*d_34*d_44-d_13*d_44^2,2*d_13*d_23*d_24+d_11*d_24*d_34-3*d_
      11*d_23*d_44-d_23*d_33*d_44+d_13*d_34*d_44,2*d_12*d_23*d_24+d_12*d_22*d_
      34-3*d_11*d_22*d_44+8*d_23^2*d_44+2*d_13*d_24*d_44-3*d_12*d_34*d_44+2*d_
      34^2*d_44-7*d_11*d_44^2-2*d_33*d_44^2,2*d_11*d_23*d_24+d_11*d_22*d_34-3*d
      _11*d_12*d_44-2*d_13*d_23*d_44-d_12*d_33*d_44+3*d_11*d_34*d_44,2*d_13^2*d
      _24+d_11*d_23*d_34+3*d_23*d_33*d_34-2*d_13*d_34^2-3*d_11*d_13*d_44-d_13*d
      _33*d_44,2*d_11*d_13*d_24+d_11*d_12*d_34+3*d_12*d_33*d_34-3*d_11*d_34^2-d
      _33*d_34^2-3*d_11^2*d_44-2*d_13^2*d_44+2*d_11*d_33*d_44+d_33^2*d_44,2*d_
      13*d_23^2-2*d_12*d_23*d_33+d_11*d_23*d_34+d_23*d_33*d_34-d_11*d_13*d_44-d
      _13*d_33*d_44,2*d_11*d_23^2+6*d_23^2*d_33+d_11*d_12*d_34-d_12*d_33*d_34-d
      _11*d_34^2+d_33*d_34^2-3*d_11^2*d_44-4*d_11*d_33*d_44-d_33^2*d_44,2*d_13^
      2*d_23-3*d_11*d_23*d_33-d_23*d_33^2+d_11*d_13*d_34+d_13*d_33*d_34,2*d_11*
      d_13*d_23-3*d_11*d_12*d_33-2*d_13*d_23*d_33-d_12*d_33^2+d_11^2*d_34+3*d_
      11*d_33*d_34)
  Istat=ideal(t_1-s_11-s_33,t_2-s_22-s_44,t_3-2*s_12,t_4-2*s_23-2*s_14,t_5-2*s_34)
  IG1= ideal(4*t_1*t_2-t_3^2-t_4^2+2*t_3*t_5-t_5^2)
  JointVariety=I+J11+Istat+IG1
  IS=eliminate({d_11,d_12,d_13,d_14,d_22,d_23,d_24,d_33,d_34,d_44,t_1,t_2,t_3,t_4,t_5},JointVariety)
Irank=minors(2,S)
IS==Irank
isSubset(Irank,IS)
isSubset(IS,Irank)

Ifiber=ideal(d_11-s_11 + d_33-s_33,d_22+d_44-s_22-s_44,d_12-s_12,
d_14+d_23-s_14-s_23, d_34-s_34)
decompose Irank
isPrime Irank
isPrime Ifiber
decompose (Irank+Ifiber)
isPrime (Irank+Ifiber)
bdry=eliminate(toList(s_11..s_44),Irank+Ifiber)
polyBdry=(flatten entries gens bdry)_0
loadPackage "SemidefiniteProgramming"
loadPackage "SumsOfSquares"
sol=solveSOS (-polyBdry)
--NM=lcm drop (flatten entries jacobian sub(bdry,QQ[d_11,d_12,d_13,d_14,d_22,d_23,d_24,d_33,d_34,d_44]),{2,6})
--minimalPrimes (sub(bdry,ring NM)+NM)
isPrime bdry 
delete (0,flatten entries jacobian bdry)
Iminors=ideal(d_11*d_22-d_12^2+d_22*d_33-d_23^2+d_33*d_44-d_34^2+d_11*d_44-d_14^2)
eliminate(toList(s_11..s_44),Irank+Ifiber+Iminors)

-- Graph 9
restart
R=QQ[s_11..s_44,d_11..d_44]
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}
adjS=exteriorPower_3 S
D=matrix{{d_11,d_12,d_13,d_14},{d_12,d_22,d_23,d_24},{d_13,d_23,d_33,d_34},{d_14,d_24,d_34,d_44}}
adjD=exteriorPower_3 D
I=ideal(d_11-s_11,d_22+d_44-s_22-s_44,d_33-s_33,d_12+d_14-s_12-s_14,d_23-s_23,d_34-s_34,
    adjD_(0,2), adjD_(1,3),adjD_(2,2)-adjD_(0,0),adjD_(0,3)-adjD_(2,3))
Irank=minors(2,S)
J=eliminate(toList(s_11..s_44),I+Irank)
Jdecomp=decompose J
isSubset(minors(4,D),Jdecomp_0) --true
isSubset(minors(4,D),Jdecomp_1) --true
isSubset(minors(4,D),Jdecomp_2) --true
isSubset(minors(4,D),Jdecomp_3) --true
isSubset(minors(4,D),Jdecomp_4) --false
eliminate({d_12,d_13,d_14,d_23,d_24,d_34},Jdecomp_4)
eliminate({d_11,d_22,d_33,d_44},Jdecomp_1)
eliminate({d_14,d_24,d_34,d_44},Jdecomp_1)
Jdecomp_4

--SOS decomposition when diagonal entries are squared
restart
R=QQ[s_11..s_44,d_11..d_44]
S=matrix{{s_11^2,s_12,s_13,s_14},{s_12,s_22^2,s_23,s_24},{s_13,s_23,s_33^2,s_34},{s_14,s_24,s_34,s_44^2}}
D=matrix{{d_11^2,d_12,d_13,d_14},{d_12,d_22^2,d_23,d_24},{d_13,d_23,d_33^2,d_34},{d_14,d_24,d_34,d_44^2}}
Irank=minors(2,S)
Ifiber=ideal(d_11^2-s_11^2 + d_33^2-s_33^2,d_22^2+d_44^2-s_22^2-s_44^2,d_12-s_12,
d_14+d_23-s_14-s_23, d_34-s_34)
bdry=eliminate(toList(s_11..s_44),Irank+Ifiber)
polyBdry=(flatten entries gens bdry)_0
loadPackage "SemidefiniteProgramming"
loadPackage "SumsOfSquares"
sol=solveSOS (-polyBdry) --no solution
peek sol
s=sosPoly sol
peek s
coefs=s#coefficients
gene=s#generators
netList gene
prod=apply(coefs,gene,(i,j)->i*j^2)
sos=sum(prod)
sos===-polyBdry
