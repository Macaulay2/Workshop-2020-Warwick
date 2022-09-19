-- 3-cycle: 2 vertices and 2 edges including edge between the two vertices
restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_3,l_3},{l_3,l_1,l_4},{l_3,l_4,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

Raux=QQ[gens Rtotal,d_11,d_12,d_13,d_22,d_23,d_33]
gens Raux
D=genericSymmetricMatrix(Raux,d_11,3)
netList 

load "setupReciprocalVarieties.m2"
Dadj=adjugate D
f1=Dadj_(0,0)-Dadj_(1,1)
f2=Dadj_(0,1)-Dadj_(0,2)

stats=sub(stats,Raux)
eq=trim sub(stats,{t_1=>d_11+d_22,t_2=>d_33,t_3=>2*d_12+2*d_13,t_4=>2*d_23})
netList eq_*

MLEeq=ideal{f1}+ideal{f2}+eq
netList MLEeq_*    
netList (trim MLEeq)_*

mingens MLEeq
netList minimalPrimes MLEeq
isPrime MLEeq
codim MLEeq

S=sub(S,Raux)
MLEeq1=MLEeq+minors(2,S)
netList (trim MLEeq1)_*
netList minimalPrimes MLEeq1
isPrime MLEeq1

toString gens Raux
MLEeq1d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(2,2), s_(2,3), s_(3,3)},MLEeq1)
netList MLEeq1d_*
codim MLEeq1d, degree MLEeq1d
L=minimalPrimes MLEeq1d
netList (L_2)_*
codim L_2,degree L_2
det D % gb L_2


MLEeq2=MLEeq+minors(3,S)
netList (trim MLEeq2)_*
netList minimalPrimes MLEeq2
isPrime MLEeq2
MLEeq2d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(2,2), s_(2,3), s_(3,3)},MLEeq2)
det D % gb MLEeq2d

--nice grobner basis

-- 3-cycle with 3 edges equal
restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_4,l_4},{l_4,l_2,l_4},{l_4,l_4,l_3}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

Raux=QQ[gens Rtotal,d_11,d_12,d_13,d_22,d_23,d_33]
gens Raux
D=genericSymmetricMatrix(Raux,d_11,3)
netList 

load "setupReciprocalVarieties.m2"
Dadj=adjugate D
f1=Dadj_(0,1)-Dadj_(1,2)
f2=Dadj_(0,1)-Dadj_(0,2)

stats=sub(stats,Raux)
eq=trim sub(stats,{t_1=>d_11,t_2=>d_22,t_3=>d_33,t_4=>2*d_12+2*d_13+2*d_23})
netList eq_*

MLEeq=ideal{f1}+ideal{f2}+eq
netList MLEeq_*    
netList (trim MLEeq)_*

mingens MLEeq
netList minimalPrimes MLEeq
isPrime MLEeq
codim MLEeq

S=sub(S,Raux)
MLEeq1=MLEeq+minors(2,S)
netList (trim MLEeq1)_*
netList minimalPrimes MLEeq1
isPrime MLEeq1

toString gens Raux
MLEeq1d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(2,2), s_(2,3), s_(3,3)},MLEeq1)
netList MLEeq1d_*
codim MLEeq1d, degree MLEeq1d
L=minimalPrimes MLEeq1d
length L
netList (L_1)_*
netList (L_0)_*
codim L_1,degree L_1
det D % gb L_1


MLEeq1s=trim eliminate({d_11, d_12, d_13, d_22, d_23, d_33},MLEeq1)
netList MLEeq1s_* 
MLEeq1s==minors(2,S)

MLEeq2=MLEeq+minors(3,S)
netList (trim MLEeq2)_*
netList minimalPrimes MLEeq2
isPrime MLEeq2
MLEeq2d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(2,2), s_(2,3), s_(3,3)},MLEeq2)
det D % gb MLEeq2d

--nice grobner basis

-- 3-cycle with 2 vertices equal
restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_3,l_4},{l_3,l_1,l_5},{l_4,l_5,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

Raux=QQ[gens Rtotal,d_11,d_12,d_13,d_22,d_23,d_33]
gens Raux
D=genericSymmetricMatrix(Raux,d_11,3)

load "setupReciprocalVarieties.m2"
Dadj=adjugate D
f=Dadj_(0,0)-Dadj_(1,1)

stats=sub(stats,Raux)
eq=trim sub(stats,{t_1=>d_11+d_22,t_2=>d_33,t_3=>2*d_12,t_4=>2*d_13,t_5=>2*d_23})
netList eq_*

MLEeq=ideal{f}+eq
netList MLEeq_*    
netList (trim MLEeq)_*

mingens MLEeq
netList minimalPrimes MLEeq
isPrime MLEeq
codim MLEeq

S=sub(S,Raux)
MLEeq1=MLEeq+minors(2,S)
netList (trim MLEeq1)_*
netList minimalPrimes MLEeq1

toString gens Raux
MLEeq1d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(2,2), s_(2,3), s_(3,3)},MLEeq1)
netList MLEeq1d_*
codim MLEeq1d, degree MLEeq1d
L=minimalPrimes MLEeq1d
length L
netList (L_1)_*
netList (L_0)_*
codim L_1,degree L_1
det D % gb L_1


MLEeq1s=trim eliminate({d_11, d_12, d_13, d_22, d_23, d_33},MLEeq1)
netList MLEeq1s_* 
MLEeq1s==minors(2,S)

MLEeq2=MLEeq+minors(3,S)
netList (trim MLEeq2)_*
netList minimalPrimes MLEeq2
isPrime MLEeq2
MLEeq2d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(2,2), s_(2,3), s_(3,3)},MLEeq2)
det D % gb MLEeq2d



--4-cycle: GRAPH 9
restart
load "functions.m2"
R=QQ[l_1..l_6]
K=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

Raux=QQ[gens Rtotal,d_11,d_12,d_13,d_14,d_22,d_23,d_24,d_33,d_34,d_44]
gens Raux
D=genericSymmetricMatrix(Raux,d_11,4)

load "setupReciprocalVarieties.m2"
Dadj=adjugate D
f1=Dadj_(1,1)-Dadj_(3,3)
f2=Dadj_(0,1)-Dadj_(0,3)

stats=sub(stats,Raux)
eq=trim sub(stats,{t_1=>d_11,t_2=>d_22+d_44,t_3=>d_33,t_4=>2*d_12+2*d_14,t_5=>2*d_23,t_6=>2*d_34})
netList eq_*

MLEeq=ideal{f1}+ideal{f2}+eq
netList MLEeq_*    
netList (trim MLEeq)_*

mingens MLEeq
netList minimalPrimes MLEeq
codim MLEeq

S=sub(S,Raux)
MLEeq1=MLEeq+minors(2,S)
netList (trim MLEeq1)_*
netList minimalPrimes MLEeq1

toString gens Raux
MLEeq1d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(1,4),s_(2,2), s_(2,3),s_(2,4), s_(3,3),s_(3,4),s_(4,4)},MLEeq1)
netList MLEeq1d_*
codim MLEeq1d, degree MLEeq1d
L=minimalPrimes MLEeq1d
length L
netList (L_0)_*
netList (L_1)_*
netList (L_2)_*
netList (L_3)_*
netList (L_4)_*
netList (L_5)_*

codim L_1,degree L_1
det D % gb L_5


MLEeq1s=trim eliminate({d_11, d_12, d_13, d_22, d_23, d_33},MLEeq1)
netList MLEeq1s_* 
MLEeq1s==minors(2,S)

MLEeq2=MLEeq+minors(3,S)
netList (trim MLEeq2)_*
netList minimalPrimes MLEeq2
isPrime MLEeq2
MLEeq2d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(1,4),s_(2,2), s_(2,3),s_(2,4), s_(3,3),s_(3,4),s_(4,4)},MLEeq2)
det D % gb MLEeq2d
