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
L=minimalPrimes MLEeq
netList L
codim MLEeq

S=sub(S,Raux)
MLEeq1=MLEeq+minors(2,S)
netList (trim MLEeq1)_*
L1=time minimalPrimes MLEeq1
netList L1
length L1
netList (L1_0)_*
netList (L1_1)_*
netList (L1_2)_*
netList (L1_3)_*
netList (L1_4)_*
netList (L1_5)_*
det D % gb L1_5

adjeq=ideal{f1,f2}
adjeq1=adjeq+minors(2,S)
adjeq1==MLEeq1
L1adj=time minimalPrimes adjeq1
length L1adj
netList (L1adj_0)_*
netList (L1adj_1)_*
det D % gb L1adj_0
det D % gb L1adj_1
det D % gb adjeq
det D % gb adjeq1

toString gens Raux
MLEeq1d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(1,4),s_(2,2), s_(2,3),s_(2,4), s_(3,3),s_(3,4),s_(4,4)},MLEeq1)
netList MLEeq1d_*
codim MLEeq1d, degree MLEeq1d
L1d=minimalPrimes MLEeq1d
length L1d
netList (L1d_0)_*
netList (L1d_1)_*
netList (L1d_2)_*
netList (L1d_3)_*
netList (L1d_4)_*
netList (L1d_5)_*

codim L1d_1,degree L1d_1
codim L1d_5,degree L1d_5
det D % gb L1d_5
det D % gb L1d_4
det D % gb L1d_3
det D % gb L1d_2
d_33*d_44-d_34^2 % gb L1d_2
d_33 % gb L1d_1
d_33 % gb L1d_0

MLEeq1s=trim eliminate({d_11, d_12, d_13,d_14, d_22, d_23,d_24,d_33,d_34,d_44},MLEeq1)
netList MLEeq1s_* 
MLEeq1s==minors(2,S)

MLEeq2=MLEeq+minors(3,S)
netList (trim MLEeq2)_*
netList minimalPrimes MLEeq2
isPrime MLEeq2
MLEeq2d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(1,4),s_(2,2), s_(2,3),s_(2,4), s_(3,3),s_(3,4),s_(4,4)},MLEeq2)
det D % gb MLEeq2d


-- Frets' head
restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_3,0,l_4},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_4,0,l_5,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

Raux=QQ[gens Rtotal,d_11,d_12,d_13,d_14,d_22,d_23,d_24,d_33,d_34,d_44]
gens Raux
D=genericSymmetricMatrix(Raux,d_11,4)

load "setupReciprocalVarieties.m2"
Dadj=adjugate D
f1=Dadj_(0,0)-Dadj_(1,1); --lambda1 from Uhler's paper
f2=Dadj_(2,2)-Dadj_(3,3); --lambda2
f3=Dadj_(0,3)-Dadj_(1,2); --lambda4

stats=sub(stats,Raux)
eq=trim sub(stats,{t_1=>d_11+d_22,t_2=>d_33+d_44,t_3=>2*d_12,t_4=>2*d_23+2*d_14,t_5=>2*d_34})
netList eq_*

MLEeq=ideal{f1,f2,f3}+eq
netList (trim MLEeq)_*

S=sub(S,Raux)
MLEeq1=MLEeq+minors(2,S);
MLEeq1d=time trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(1,4),s_(2,2), s_(2,3),s_(2,4), s_(3,3),s_(3,4),s_(4,4)},MLEeq1);
-- used 822.141 seconds
netList MLEeq1d_*
codim MLEeq1d, degree MLEeq1d
--4,108
betti MLEeq1d
--3 cubics, 1 quartic

ideal(d_14^2*d_33-d_24^2*d_33-2*d_13*d_14*d_34+2*d_23*d_24*d_34+d_11*d_34^2-d_22*d_34^2+d_13^2*d_44-d_23^2*d_44-d_11*d_33*d_44+d_22*d_33*d_44,d_14^2*d_23-d_14*d_23^2-d_13*d_14*d_24+d_
      13*d_23*d_24+d_14*d_22*d_33-d_12*d_24*d_33-d_12*d_14*d_34-d_13*d_22*d_34+d_12*d_23*d_34+d_11*d_24*d_34+d_12*d_13*d_44-d_11*d_23*d_44,d_13^2*d_22-d_14^2*d_22-2*d_12*d_13*d_23+d_11*d_23
      ^2+2*d_12*d_14*d_24-d_11*d_24^2+d_12^2*d_33-d_11*d_22*d_33-d_12^2*d_44+d_11*d_22*d_44,d_14^4+14*d_14*d_23^3+d_23^4+4*d_13*d_14^2*d_24+6*d_13*d_14*d_23*d_24-10*d_13*d_23^2*d_24-2*d_11*
      d_14*d_23*d_33-12*d_14*d_22*d_23*d_33-d_11*d_23^2*d_33-d_22*d_23^2*d_33+4*d_12*d_14*d_24*d_33+10*d_12*d_23*d_24*d_33-d_11*d_24^2*d_33-5*d_22*d_24^2*d_33+d_12^2*d_33^2-2*d_11*d_13*d_14
      *d_34-6*d_13*d_14*d_22*d_34-2*d_12*d_14*d_23*d_34+10*d_13*d_22*d_23*d_34-14*d_12*d_23^2*d_34-4*d_11*d_14*d_24*d_34-8*d_11*d_23*d_24*d_34+10*d_22*d_23*d_24*d_34+2*d_11*d_12*d_33*d_34+2
      *d_12*d_22*d_33*d_34+2*d_11^2*d_34^2+6*d_11*d_22*d_34^2-4*d_22^2*d_34^2+d_11*d_13^2*d_44-4*d_12*d_13*d_14*d_44-d_11*d_14^2*d_44+4*d_14^2*d_22*d_44+2*d_11*d_14*d_23*d_44-2*d_14*d_22*d_
      23*d_44+3*d_11*d_23^2*d_44-6*d_22*d_23^2*d_44-10*d_12*d_14*d_24*d_44+5*d_11*d_24^2*d_44-d_11^2*d_33*d_44-3*d_12^2*d_33*d_44+d_11*d_22*d_33*d_44+5*d_22^2*d_33*d_44+2*d_11*d_12*d_34*d_
      44+2*d_12*d_22*d_34*d_44+6*d_12^2*d_44^2-5*d_11*d_22*d_44^2)

L1d=time minimalPrimes MLEeq1d;
length L1d
netList (L1d_0)_*
netList (L1d_1)_*
netList (L1d_2)_*
netList (L1d_3)_*
netList (L1d_4)_*
netList (L1d_5)_*

codim L1d_1,degree L1d_1
codim L1d_5,degree L1d_5
det D % gb L1d_5
det D % gb L1d_4
det D % gb L1d_3
det D % gb L1d_2
d_33*d_44-d_34^2 % gb L1d_2
d_33 % gb L1d_1
d_33 % gb L1d_0

MLEeq1s=trim eliminate({d_11, d_12, d_13,d_14, d_22, d_23,d_24,d_33,d_34,d_44},MLEeq1)
netList MLEeq1s_* 
MLEeq1s==minors(2,S)

MLEeq2=MLEeq+minors(3,S)
netList (trim MLEeq2)_*
netList minimalPrimes MLEeq2
isPrime MLEeq2
MLEeq2d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(1,4),s_(2,2), s_(2,3),s_(2,4), s_(3,3),s_(3,4),s_(4,4)},MLEeq2)
det D % gb MLEeq2d
