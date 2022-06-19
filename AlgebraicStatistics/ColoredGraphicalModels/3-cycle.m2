restart
R=QQ[s_11,s_12,s_13,s_22,s_23,s_33]
S=genericSymmetricMatrix(R,3)
I=minors(2,S)
(dim I,degree I) --(3,4)


R1=QQ[join(gens R,x_1..x_6)]
use R1
X=matrix{{x_1..x_3},{x_4..x_6}}
SHat=transpose(X)*X
Adj=exteriorPower (2,SHat)
rank Adj
S1=genericSymmetricMatrix(R1,3)
J=eliminate(toList(x_1..x_6),ideal(Adj-S1))
J==minors(2,S1)

for i to 100 list random RR
