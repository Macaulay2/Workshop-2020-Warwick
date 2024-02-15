
-----------------------------------------------------------
--4-cycle with equal vertices
-----------------------------------------------------------
restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_2,0,l_3},{l_2,l_1,l_4,0},{0,l_4,l_1,l_5},{l_3,0,l_5,l_1}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

--Elimination ideal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
netList IG1_*

IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2



restart
R=QQ[a_1..a_4,l_1..l_5]

L2=transpose(matrix{{a_1,a_2,a_3,a_4}})*matrix{{a_1,a_2,a_3,a_4}}
K1=sub(matrix{{1,0,0,0},{0,-1,0,0},{0,0,0,0},{0,0,0,0}},R)
K2=sub(matrix{{1,0,0,0},{0,0,0,0},{0,0,-1,0},{0,0,0,0}},R)
K3=sub(matrix{{1,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,-1}},R)
K4=sub(matrix{{0,0,1,0},{0,0,0,0},{1,0,0,0},{0,0,0,0}},R)
K5=sub(matrix{{0,0,0,0},{0,0,0,1},{0,0,0,0},{0,1,0,0}},R)
K={0,K1,K2,K3,K4,K5}
KT=sum toList apply(1..5,i->l_i*K_i)

P=L2+KT

det P

toList apply(0..(numcols P),i->toList(0..i))

P_{0,1}^{0,1}

leadingPrincipalMinors=M->(
mm:=numcols M;
ss:=toList apply(0..(numcols M-1),i->toList(0..i));
for s in ss list M_s^s
);

LKT=leadingPrincipalMinors(KT)
netList apply(LKT,i->det i)
rank KT

LP=leadingPrincipalMinors(P);
netList apply(LP,i->det i)

toExternalString apply(LKT,i->det i)
FindInstance[l_1+l_2+l_3>0 && -l_1^2-l_1*l_2-l_1*l_3>0 && l_1^2*l_2+l_1*l_2^2+l_1*l_2*l_3+l_1*l_4^2>0 && -l_1^2*l_2*l_3-l_1*l_2^2*l_3-l_1*l_2*l_3^2-l_1*l_3*l_4^2+l_1*l_2*l_5^2+l_2^2*l_5^2+l_2*l_3*l_5^2+l_4^2*l_5^2>0]


l

-----------------------------------------------------------
--5-cycle with equal vertices
-----------------------------------------------------------
restart
load "functions.m2"
R=QQ[l_1..l_6]
K=matrix{{l_1,l_2,0,0,l_3},{l_2,l_1,l_4,0,0},{0,l_4,l_1,l_5,0},{0,0,l_5,l_1,l_6},{l_3,0,0,l_6,l_1}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

--Elimination ideal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
netList IG1_*

IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2
