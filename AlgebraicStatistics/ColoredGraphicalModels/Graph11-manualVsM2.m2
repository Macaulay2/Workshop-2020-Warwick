restart
R=QQ[s_11..s_44,d_11..d_44,t_1..t_5]
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}
adjS=exteriorPower_3 S
D=matrix{{d_11,d_12,d_13,d_14},{d_12,d_22,d_23,d_24},{d_13,d_23,d_33,d_34},{d_14,d_24,d_34,d_44}}
adjD=exteriorPower_3 D

Irank=minors(2,S)
Ifiber=ideal(d_11-s_11 + d_33-s_33,d_22+d_44-s_22-s_44,d_12-s_12,
d_14+d_23-s_14-s_23, d_34-s_34)
bdry=eliminate(toList(s_11..s_44),Irank+Ifiber)
--bdry=ideal(d_12^2+d_14^2-d_11*d_22+2*d_14*d_23+d_23^2-d_22*d_33-2*d_12*d_34+d_34^2-d_11*d_44-d_33*d_44)

IPrincipalMinors=ideal(s_11*s_22-s_12^2,s_22*s_33-s_23^2,s_33*s_44-s_34^2,s_11*s_44-s_14^2)
bdryPrincipalMinors=eliminate(toList(s_11..s_44),IPrincipalMinors+Ifiber)
IAdditionalMinors=ideal(s_12*s_34-s_23*s_14)
bdryPrincipalMinors=eliminate(toList(s_11..s_44),IPrincipalMinors+IAdditionalMinors+Ifiber)
J=ideal(d_12*d_34-d_14*d_23)
isSubset(J,IPrincipalMinors+IAdditionalMinors+Ifiber)
gens gb (IPrincipalMinors+IAdditionalMinors+Ifiber) --14
gens gb IPrincipalMinors --4
gens gb IAdditionalMinors --1
gens gb Ifiber --5
gens gb (IPrincipalMinors+IAdditionalMinors)  --9 (Irank has 20 gens)

