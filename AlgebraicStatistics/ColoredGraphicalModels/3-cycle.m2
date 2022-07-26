restart
R=QQ[s_11,s_12,s_13,s_22,s_23,s_33]
S=genericSymmetricMatrix(R,3)
I=minors(2,S)
(dim I,degree I) --(3,4)


R1=QQ[join(x_1..x_6,gens R)]
use R1
X=matrix{{x_1..x_3},{x_4..x_6}}
SHat=transpose(X)*X
Adj=exteriorPower (2,SHat)
rank Adj
S1=genericSymmetricMatrix(R,3)
S1=sub(S1,R1)
AdjS1=ideal(Adj-S1) 
Elim1=eliminate(toList(x_2..x_6), AdjS1)
J=eliminate(toList(x_1..x_6),ideal(Adj-S1))
J==minors(2,S1)

X1=matrix{{x_1..x_3}}
SHat1=transpose(X1)*X1
Adj1=exteriorPower (2,SHat1)
rank Adj1

AdjAdj=exteriorPower (2,Adj)
rank AdjAdj

rank jacobian ideal flatten entries Adj

Y=random(QQ^3,QQ^1)
S=Y*transpose(Y)

I=ideal(Adj-S)
dim I

R2=QQ[x_1..x_6]
use R2
I2=sub (I,R2)
(dim I2, codim I2)
degree I2

GenH=ideal(x_1^2+3*x_3+7*x_5,11*x_1+3*x_2+5*x_6,13*x_2^2+5*x_3^3+x_4+5*x_6) -- not generic!!!

dim (I2+GenH)
degree (I2+GenH)
-- intersect with generic 9 dimensional variety
--for i to 100 list random R
restart
R=QQ[l1,l2,l3,l4]
Kv=matrix{{l1,l2,l3},{l2,l1,l4},{l3,l4,l1}} -- all vertices
Ke=matrix{{l1,l4,l4},{l4,l2,l4},{l4,l4,l3}} -- all edges
K2e=matrix{{l1,l4,0},{l4,l2,l4},{0,l4,l3}} -- 2 edges, one removed
fv=det Kv
fe=det Ke  
f2e=det K2e  
Hv=jacobian ideal flatten entries jacobian ideal fv
He=jacobian ideal flatten entries jacobian ideal fe
H2e=jacobian ideal flatten entries jacobian ideal f2e
rank Hv
rank He
rank H2e
