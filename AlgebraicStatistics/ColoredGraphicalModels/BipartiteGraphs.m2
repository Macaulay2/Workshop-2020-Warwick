restart
load "functions.m2"

G=graph{{1,6},{1,7},{1,8},{1,9},{1,10},{2,6},{2,7},{2,8},{2,9},{2,10},{3,6},{3,7},{3,8},{3,9},{3,10},{4,6},{4,7},{4,8},{4,9},{4,10},{5,6},{5,7},{5,8},{5,9},{5,10}}
RG=gaussianRing G
KG=undirectedEdgesMatrix RG
gens RG
support KG

R=QQ[support KG]
K=sub(KG,R);

I5=time minors(6,K);
time rank K  --10
--eliminating after computing the minors is not the same as
--setting to 0 and computing the minors

varList={s_(1,2),s_(1,3),s_(1,4),s_(1,5),s_(2,3),s_(2,4),s_(2,5),s_(3,4),s_(3,5),s_(4,5),
         s_(6,7),s_(6,8),s_(6,9),s_(6,10),s_(7,8),s_(7,9),s_(7,10),s_(8,9),s_(8,10),s_(9,10)};

M6=time minors(6,S);
I5aux=time eliminate(varList,M6);


(p,n,Rtotal,S)=coloredData(K)

suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*


-- Computation of the elimination ideals
P5=time rankProjection(stats,5,S);

IG5=time sub(rankProjection(stats,5,S),ring(suffStat(K)));
IG4=time sub(rankProjection(stats,4,S),ring(suffStat(K)));
--according to corollary 2.11 in Bleckerman-Sinn: IG5=0 and IG4 not 0
--because the generic completion rank is 5 but the mle threshold is 4

